"""Unit tests for the two leaf solvers on planted solutions of sigma(n) = 2n + c.

Brute force cannot reach large exponents or large last primes, so plant n from
a chosen factorisation, set c = sigma(n) - 2n, and check that each leaf solver
recovers the missing piece from the rest:

  last_prime:  every other component fixed, the largest prime unknown.
  single_open: every prime known, one exponent open from 2.
"""
import random
import sys

from qp import Search, is_probable_prime, sig

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
odd_primes = [p for p in range(3, 2000) if is_probable_prime(p)]

done = {"last_prime": 0, "single_open": 0}
fails = 0
while min(done.values()) < N:
    w = rng.randint(2, 5)
    ps = sorted(rng.sample(odd_primes[:15], w - 1))
    big = rng.choice([p for p in odd_primes if p > ps[-1]])
    ps.append(big)
    exps = [rng.choice([1, 1, 1, 2, 3, 4, 7, 12]) for _ in ps]
    n, s = 1, 1
    for p, a in zip(ps, exps):
        n *= p ** (2 * a)
        s *= sig(p, 2 * a)
    c = s - 2 * n
    if c <= 0:
        continue
    fixed = tuple((p, a, True) for p, a in zip(ps, exps))

    # last prime: the rest must be deficient (h < 2) for the leaf to apply
    m, sm = 1, 1
    for p, a in zip(ps[:-1], exps[:-1]):
        m *= p ** (2 * a)
        sm *= sig(p, 2 * a)
    if sm < 2 * m and done["last_prime"] < N:
        S = Search(c=c)
        S.last_prime(fixed[:-1], __import__("fractions").Fraction(sm, m), m, ps[-2])
        done["last_prime"] += 1
        if n not in S.solutions or S.unresolved:
            fails += 1
            print("last_prime MISS", ps, exps, c, S.solutions, S.unresolved[:1])

    # single open exponent, any position
    if done["single_open"] < N:
        i = rng.randrange(w)
        comps = fixed[:i] + ((ps[i], 1, False),) + fixed[i + 1:]
        S = Search(c=c)
        S.single_open(comps, i)
        done["single_open"] += 1
        if n not in S.solutions:
            fails += 1
            print("single_open MISS", ps, exps, i, c, S.solutions)

print(f"{done} planted, {fails} failures")
