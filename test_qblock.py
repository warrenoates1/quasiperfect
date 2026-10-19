"""BlockSearch must find exactly what Search finds, on many (c, omega)."""
import sys, time
from qp import Search
import qblock

calls = [0]
_orig = qblock.BlockSearch.block_root
def _br(self, *a, **k):
    calls[0] += 1
    return _orig(self, *a, **k)
qblock.BlockSearch.block_root = _br

CMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
W = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [3, 4]
bad = sols = 0
t0 = time.time()
for w in W:
    for c in range(1, CMAX + 1):
        a = Search(c=c); a.run(w)
        b = qblock.BlockSearch(c=c); b.run(w)
        sols += len(a.solutions)
        if sorted(a.solutions) != sorted(b.solutions) or bool(a.unresolved) != bool(b.unresolved):
            bad += 1
            print("MISMATCH", c, w, a.solutions, b.solutions, len(a.unresolved), len(b.unresolved))
print(f"{bad} mismatches, {sols} solutions compared, {calls[0]} block roots, {time.time()-t0:.0f}s")
