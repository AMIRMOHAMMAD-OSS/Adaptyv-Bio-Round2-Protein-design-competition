import argparse
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time
import zipfile


RF_MODELS = [
    ("Base_ckpt.pt", 483616107, "6f5902ac237024bdd0c176cb93063dc4"),
    ("Complex_base_ckpt.pt", 483619179, "e29311f6f1bf1af907f9ef9f44b8328b"),
]
AF_ARCHIVE = "alphafold_params_colab_2022-12-06.tar"
AF_SIZE = 4099624960
AF_NAMES = {f"params_model_{i}_multimer_v3.npz" for i in range(1, 6)}
AF_MARKER = "download_complexes_multimer_v3_finished.txt"


@contextmanager
def lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f"Another download is using {path}")
        yield


def check_space(directory, needed):
    directory.mkdir(parents=True, exist_ok=True)
    locations = [directory]
    host_drive = os.environ.get("BINDER_HOST_DRIVE")
    if host_drive:
        locations.append(Path(host_drive))
    elif "microsoft" in os.uname().release.lower() and Path("/mnt/c").is_dir():
        locations.append(Path("/mnt/c"))
    for location in locations:
        free = shutil.disk_usage(location).free
        print(f"Free at {location}: {free / 1024**3:.1f} GiB; needed: {needed / 1024**3:.1f} GiB", flush=True)
        if free < needed:
            raise RuntimeError(f"Free more space on {location}. Download progress is retained.")


def missing_bytes(destination, expected):
    partial = Path(str(destination) + ".part")
    existing = destination if destination.exists() else partial
    return max(0, expected - (existing.stat().st_size if existing.exists() else 0))


def download(url, destination, expected, attempts=30, retry_delay=5):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(str(destination) + ".part")
    if destination.exists():
        if destination.stat().st_size != expected:
            raise RuntimeError(f"Unexpected size; file preserved: {destination}")
        print(f"Already downloaded: {destination}", flush=True)
        return
    for attempt in range(1, attempts + 1):
        size = partial.stat().st_size if partial.exists() else 0
        if size == expected:
            os.replace(partial, destination)
            return
        if size > expected:
            raise RuntimeError(f"Oversized partial file preserved: {partial}")
        print(f"Attempt {attempt}/{attempts}: {destination.name}, {size}/{expected} bytes", flush=True)
        result = subprocess.run([
            "curl", "--fail", "--location", "--http1.1", "--retry", "0",
            "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "90",
            "--continue-at", "-", "--output", str(partial), url,
        ])
        size = partial.stat().st_size if partial.exists() else 0
        if result.returncode == 0 and size == expected:
            os.replace(partial, destination)
            return
        if result.returncode in {23, 26, 60, 77}:
            raise RuntimeError(f"curl failed with code {result.returncode}; partial file retained.")
        if attempt < attempts:
            time.sleep(retry_delay)
    raise RuntimeError(f"Download incomplete. Rerun to resume {partial}")


def valid_npz(path):
    try:
        with zipfile.ZipFile(path) as archive:
            return bool(archive.namelist()) and archive.testzip() is None
    except (OSError, zipfile.BadZipFile, EOFError):
        return False


def extract_params(archive, destination):
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:") as bundle:
        members = {}
        for member in bundle.getmembers():
            name = Path(member.name).name
            if name in AF_NAMES:
                if not member.isfile() or name in members:
                    raise RuntimeError(f"Unexpected archive member: {member.name}")
                members[name] = member
        if set(members) != AF_NAMES:
            raise RuntimeError(f"Missing parameter files: {sorted(AF_NAMES - set(members))}")
        for name in sorted(AF_NAMES):
            member = members[name]
            target = destination / name
            if target.is_file() and target.stat().st_size == member.size and valid_npz(target):
                print(f"Verified existing parameter file: {name}", flush=True)
                continue
            temporary = destination / (name + ".extracting")
            print(f"Extracting and checking: {name}", flush=True)
            with bundle.extractfile(member) as source, temporary.open("wb") as output:
                shutil.copyfileobj(source, output, length=2**20)
            if temporary.stat().st_size != member.size or not valid_npz(temporary):
                raise RuntimeError(f"Parameter file failed validation: {temporary}")
            os.replace(temporary, target)
    (destination / AF_MARKER).touch()


def fetch_rf(root):
    root = Path(root) / "models"
    with lock(root / ".download.lock"):
        needed = sum(missing_bytes(root / name, size) for name, size, _ in RF_MODELS)
        if needed:
            check_space(root, needed + 1024**3)
        for name, size, key in RF_MODELS:
            download(f"https://files.ipd.uw.edu/pub/RFdiffusion/{key}/{name}", root / name, size)


def fetch_af(root):
    root = Path(root)
    destination = root / "params"
    archive = root / "downloads" / AF_ARCHIVE
    with lock(root / ".binderflow-download.lock"):
        complete = all(valid_npz(destination / name) for name in sorted(AF_NAMES))
        if complete:
            (destination / AF_MARKER).touch()
            print("All five AlphaFold Multimer-v3 parameter files are verified.", flush=True)
            return
        needed = missing_bytes(archive, AF_SIZE) + AF_SIZE + 2 * 1024**3
        check_space(root, needed)
        download(f"https://storage.googleapis.com/alphafold/{AF_ARCHIVE}", archive, AF_SIZE)
        extract_params(archive, destination)
        archive.unlink()
        print("All five AlphaFold Multimer-v3 parameter files are verified.", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rfd-root", required=True, type=Path)
    parser.add_argument("--colabfold-params", required=True, type=Path)
    args = parser.parse_args()
    fetch_rf(args.rfd_root)
    fetch_af(args.colabfold_params)


if __name__ == "__main__":
    main()
