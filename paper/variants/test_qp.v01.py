"""Cross-check the tree search against brute force on sigma(n) = 2n + c.

For every odd k <= K, n = k^2 is checked directly. For each c > 0 that
brute force hits, the solver must find every brute-force solution, and every
solver solution <= K^2 must appear in the brute-force list.
"""
import sys
import time
from collections import defaultdict

from qp import Search, factorize, sig

K = int(sys.argv[1]) if len(sys.argv) > 1 else 20001
CMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 5000

brute = defaultdict(list)          # (c, omega) -> [n]
for k in range(3, K + 1, 2):
    f = factorize(k)
    n = k * k
    s = 1
    for p, e in f.items():
        s *= sig(p, 2 * e)
    c = s - 2 * n
    if 0 < c <= CMAX:
        brute[(c, len(f))].append(n)

print(f"brute force: {sum(map(len, brute.values()))} hits over {len(brute)} (c, omega) pairs")
bad = 0
t0 = time.time()
for (c, w), ns in sorted(brute.items()):
    s = Search(c=c)
    s.run(w)
    found = sorted(s.solutions)
    small = [n for n in found if n <= K * K]
    ok = sorted(ns) == small and not s.unresolved
    if not ok:
        bad += 1
    print(f"c={c:5d} omega={w}: brute={sorted(ns)} solver={found} "
          f"unresolved={len(s.unresolved)} nodes={s.nodes} {'OK' if ok else 'MISMATCH'}")
print(f"{bad} mismatches, {time.time() - t0:.1f}s")
