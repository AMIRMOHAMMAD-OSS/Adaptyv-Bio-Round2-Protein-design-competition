from functools import partial
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
from pathlib import Path
import socket
import tarfile
import tempfile
import threading
import unittest
import zipfile

import numpy as np

from scripts.download_models import AF_MARKER, AF_NAMES, download, extract_params


class DownloadTests(unittest.TestCase):
    def test_resume_after_interrupted_http_download(self):
        payload = bytes(range(256)) * 32
        requests = []

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                request_range = self.headers.get("Range")
                requests.append(request_range)
                if request_range is None:
                    self.send_response(200)
                    self.send_header("Content-Length", str(len(payload)))
                    self.end_headers()
                    self.wfile.write(payload[:1024])
                    self.wfile.flush()
                    self.connection.shutdown(socket.SHUT_RDWR)
                    self.connection.close()
                else:
                    offset = int(request_range.removeprefix("bytes=").removesuffix("-"))
                    self.send_response(206)
                    self.send_header("Content-Length", str(len(payload) - offset))
                    self.send_header("Content-Range", f"bytes {offset}-{len(payload)-1}/{len(payload)}")
                    self.end_headers()
                    self.wfile.write(payload[offset:])

        server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        thread = threading.Thread(target=partial(server.serve_forever, poll_interval=0.05), daemon=True)
        thread.start()
        try:
            with tempfile.TemporaryDirectory() as temporary:
                destination = Path(temporary) / "model.bin"
                url = f"http://127.0.0.1:{server.server_port}/model.bin"
                download(url, destination, len(payload), attempts=3, retry_delay=0)
                self.assertEqual(destination.read_bytes(), payload)
                self.assertEqual(requests, [None, "bytes=1024-"])
                download(url, destination, len(payload), attempts=1, retry_delay=0)
                self.assertEqual(len(requests), 2)
        finally:
            server.shutdown()
            server.server_close()
            thread.join()

    def make_archive(self, path, invalid=False, duplicate=False):
        buffer = io.BytesIO()
        np.savez(buffer, weights=np.array([1.0, 2.0]))
        data = buffer.getvalue()
        with tarfile.open(path, "w") as archive:
            names = sorted(AF_NAMES)
            if duplicate:
                names.append(names[0])
            for index, name in enumerate(names):
                content = b"invalid" if invalid and index == 0 else data
                member = tarfile.TarInfo(name)
                member.size = len(content)
                archive.addfile(member, io.BytesIO(content))
            member = tarfile.TarInfo("../../outside.txt")
            member.size = 6
            archive.addfile(member, io.BytesIO(b"unused"))

    def test_only_verified_parameters_are_extracted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, destination = root / "models.tar", root / "params"
            self.make_archive(archive)
            extract_params(archive, destination)
            self.assertEqual({p.name for p in destination.iterdir()}, AF_NAMES | {AF_MARKER})
            self.assertFalse((root / "outside.txt").exists())
            for name in AF_NAMES:
                with zipfile.ZipFile(destination / name) as weights:
                    self.assertIsNone(weights.testzip())

    def test_invalid_or_duplicate_parameters_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for case in ("invalid", "duplicate"):
                archive, destination = root / f"{case}.tar", root / case
                self.make_archive(archive, **{case: True})
                with self.assertRaises(RuntimeError):
                    extract_params(archive, destination)
                self.assertFalse((destination / AF_MARKER).exists())


if __name__ == "__main__":
    unittest.main()
