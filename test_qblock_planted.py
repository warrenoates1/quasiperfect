"""Planted solutions: n = (m q q')^2, c = sigma(n) - 2n; both searches must find n."""
import time
from itertools import combinations
from qp import Search, sig, factorize
import qblock

def sigma_sq(k):
    s = 1
    for p, e in factorize(k).items():
        s *= sig(p, 2 * e)
    return s

t0 = time.time(); bad = tried = roots = 0
for base in (3 * 5 * 11, 3 * 5 * 13, 3 * 5 * 17):
    for q, q2 in combinations([19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103], 2):
        if base % q == 0 or base % q2 == 0:
            continue
        k = base * q * q2; n = k * k
        c = sigma_sq(k) - 2 * n
        if c <= 0 or c > 10**18:
            continue
        w = len(factorize(k))
        a = Search(c=c); a.run(w)
        b = qblock.BlockSearch(c=c); b.run(w); roots += b.splits
        tried += 1
        okb = n in b.solutions and sorted(a.solutions) == sorted(b.solutions)
        if not okb:
            bad += 1
            print("BAD", k, c, n in a.solutions, n in b.solutions, a.solutions, b.solutions)
print(f"{tried} planted, {bad} bad, {roots} splits, {time.time()-t0:.0f}s")
