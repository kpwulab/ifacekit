"""Pack the Top8000 Ramachandran grids into an npz file (development comparison only; not shipped).
Usage: python tools/build_rama_data.py /path/to/rotarama_data OUT.npz
Source tables: Top8000 Ramachandran distributions (Richardson lab), as distributed with MolProbity."""
import sys
import numpy as np

FILES = {"general": "rama8000-general-noGPIVpreP.data", "glycine": "rama8000-gly-sym.data",
         "cispro": "rama8000-cispro.data", "transpro": "rama8000-transpro.data",
         "prepro": "rama8000-prepro-noGP.data", "ileval": "rama8000-ileval-nopreP.data"}
out = {}
for k, fn in FILES.items():
    g = np.zeros((180, 180), dtype=np.float32)
    for line in open(f"{sys.argv[1]}/{fn}"):
        if line.startswith("#") or not line.strip():
            continue
        x, y, v = line.split()
        g[int((float(x) + 180) // 2), int((float(y) + 180) // 2)] = float(v)
    out[k] = g
np.savez_compressed(sys.argv[2], **out)
print({k: (v.shape, float(v.max())) for k, v in out.items()})
