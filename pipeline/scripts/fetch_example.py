from pathlib import Path
from urllib.request import urlopen
import hashlib
import json
root = Path(__file__).resolve().parents[1]
out = root / 'inputs/1MOX.pdb'
url = 'https://files.rcsb.org/download/1MOX.pdb'
expected = 'afb29731066559e9d4c437875b58a87be3a297d19623641474276cbe42102a46'
if not out.exists():
    data = urlopen(url, timeout=120).read()
    if hashlib.sha256(data).hexdigest() != expected:
        raise SystemExit('RCSB file changed from the checked example; inspect the new structure before adopting it.')
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
if hashlib.sha256(out.read_bytes()).hexdigest() != expected:
    raise SystemExit('Existing 1MOX file differs from the checked example')
print(out)
