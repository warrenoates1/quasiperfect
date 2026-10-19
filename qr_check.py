"""Verify the quadratic-residue constraint on sigma(n) = 2n + c, and measure its bite.

Claim. Let n = m^2 be an odd square with sigma(n) = 2n + c. If a prime r divides
2n + c then r does not divide 2m (else r | c and r | 2n, handled separately), and
(2m)^2 = 4m^2 = 2(2n) = 2(-c) = -2c  (mod r),
so -2c must be a quadratic residue mod r. For c = 1 that says every prime factor
of sigma(n) is congruent to 1 or 3 mod 8.

Since sigma(p^(2a)) divides sigma(n) whenever p^(2a) is the exact power of p in n,
the test applies to each component on its own, which is what makes it useful as a
prune: it rejects a (prime, exponent) pair outright, with no reference to the
rest of the state.

This script (1) tests the claim by brute force, including for c != 1, and
(2) lists which exponents survive for small primes.
"""
import sys
from math import gcd

from qp import factorize, is_probable_prime, sig


def jacobi(a, n):
    """Jacobi symbol (a/n) for odd n > 0."""
    assert n > 0 and n % 2 == 1
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def _rho(n):
    if n % 2 == 0:
        return 2
    for c in range(1, 40):
        x = y = 2
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d
    raise RuntimeError(f"rho failed on {n}")


def factor(n, out=None):
    """Full prime factorisation via Pollard rho (list of primes with multiplicity)."""
    out = [] if out is None else out
    if n == 1:
        return out
    if is_probable_prime(n):
        out.append(n)
        return out
    d = _rho(n)
    factor(d, out)
    factor(n // d, out)
    return out


def allowed(s, c):
    """Does every prime factor r of s satisfy (-2c/r) = 1?"""
    for r in set(factor(s)):
        if c % r == 0:
            continue          # r | c: the argument does not apply
        if jacobi(-2 * c % r, r) != 1:
            return False, r
    return True, None


def test_claim(limit=30001, cmax=10**6):
    """Brute force: every prime factor of sigma(m^2) - 2m^2 ... check the real
    identity on actual solutions, and the general claim on 2m^2 + c."""
    print("(1) The claim on 2m^2 + c, for every odd m and the c it actually produces")
    bad = checked = sols = 0
    for m in range(3, limit, 2):
        n = m * m
        s = 1
        for p, e in factorize(m).items():
            s *= sig(p, 2 * e)
        c = s - 2 * n
        if not (0 < c <= cmax):
            continue
        sols += 1
        ok, r = allowed(s, c)          # s = sigma(n) = 2n + c, so the claim applies
        checked += 1
        if not ok:
            bad += 1
            print(f"    COUNTEREXAMPLE m={m} n={n} c={c}: prime {r} of sigma(n) has "
                  f"(-2c/r) != 1")
    print(f"    {checked} genuine solutions of sigma(n) = 2n + c checked, "
          f"{bad} counterexamples")

    print("\n(2) The same claim stated per component, on those solutions")
    bad = checked = 0
    for m in range(3, limit, 2):
        n = m * m
        f = factorize(m)
        s = 1
        for p, e in f.items():
            s *= sig(p, 2 * e)
        c = s - 2 * n
        if not (0 < c <= cmax):
            continue
        for p, e in f.items():
            ok, r = allowed(sig(p, 2 * e), c)
            checked += 1
            if not ok:
                bad += 1
                print(f"    COUNTEREXAMPLE m={m} c={c} component {p}^{2*e}: "
                      f"prime {r} has (-2c/r) != 1")
    print(f"    {checked} components of genuine solutions checked, {bad} rejected "
          f"(any rejection would make the prune UNSOUND)")


def survivors():
    print("\n(3) For c = 1: which exponents survive, per prime")
    print("    (a prime with no surviving exponent cannot divide a quasiperfect number)")
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 257, 65537]:
        ok = []
        for a in range(1, 11):
            s = sig(p, 2 * a)
            if s.bit_length() > 320:
                ok.append(f">{2*a-2}?")
                break
            good, _ = allowed(s, 1)
            if good:
                ok.append(2 * a)
        print(f"    p={p:<7} surviving exponents 2a<=20: {ok if ok else 'NONE'}")


if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else 30001
    test_claim(lim)
    survivors()
