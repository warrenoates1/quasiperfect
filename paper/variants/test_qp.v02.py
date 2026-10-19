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

# Phase 2: planted solutions with large exponents and large last primes,
# beyond brute-force reach. Each must be rediscovered from its c alone.
import random
from qp import is_probable_prime

rng = random.Random(1)
odd_primes = [p for p in range(3, 400) if is_probable_prime(p)]
tried = missed = 0
while tried < int(sys.argv[3]) if len(sys.argv) > 3 else tried < 150:
    w = rng.randint(2, 4)
    ps = sorted(rng.sample(odd_primes[:12], w - 1))
    ps.append(rng.choice([p for p in odd_primes if p > ps[-1]] or [ps[-1] + 2]))
    if not is_probable_prime(ps[-1]) or len(set(ps)) != w:
        continue
    exps = [rng.choice([1, 1, 2, 3, 5, 9]) for _ in ps]
    n, s = 1, 1
    for p, a in zip(ps, exps):
        n *= p ** (2 * a)
        s *= sig(p, 2 * a)
    c = s - 2 * n
    if c <= 0 or c * 50 > n:          # keep c small relative to n, like the real problem
        continue
    tried += 1
    S = Search(c=c)
    S.run(w)
    if n not in S.solutions or S.unresolved:
        missed += 1
        print("MISS", ps, exps, c, S.solutions, S.unresolved[:2])
print(f"phase 2: {tried} planted, {missed} missed")
