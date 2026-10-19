"""Close a two-free-prime leaf by the discriminant criterion, with a residue sieve.

For a leaf n = m p^(2a) q^2 with m fixed, writing A = sigma(m),
S = A sigma(p^(2a)) and T = 2 m p^(2a), the equation sigma(n) = 2n + c is
quadratic in q:

    (S - T) q^2 + S q + (S - c) = 0

so Delta = S^2 - 4(S-T)(S-c) must be a perfect square and q follows as
(S +/- sqrt(Delta)) / (2(T-S)). Derived and tested in the paper; S != T always,
since S is odd and T even.

Two points about where the saving comes from, because they are easy to get wrong.
The criterion applies for ANY exponent on p -- S and T simply absorb p^(2a) --
but only for exponent 2 on q, since sigma(q^(2b)) has degree 2b. And the existing
last-prime lemma already closes every exponent of q, so this is a faster path
rather than a new capability: it replaces a polynomial bisection per candidate p
with one perfect-square test.

The real speed comes from sieving p first. A perfect square is a quadratic
residue to every modulus, so Delta(p) failing to be a residue modulo any small r
rejects p outright, with no large arithmetic. Each modulus rejects about half,
and the test is a table lookup.

Soundness: the sieve only ever rejects p for which Delta(p) cannot be square, so
it cannot discard a solution. Exponents of q above 2 are not handled here and are
left to the caller.
"""
from math import isqrt

# Moduli for the sieve. Small, pairwise coprime-ish, and cheap to index.
SIEVE_MODULI = (64, 63, 65, 11, 13, 17, 19, 23, 29, 31)


def _qr_table(r):
    """Bit set of quadratic residues modulo r."""
    s = 0
    for x in range(r):
        s |= 1 << (x * x % r)
    return s


_QR = {r: _qr_table(r) for r in SIEVE_MODULI}


def sieve_rejects(delta):
    """True if delta is definitely not a perfect square, by residue alone."""
    for r in SIEVE_MODULI:
        if not (_QR[r] >> (delta % r)) & 1:
            return True
    return False


def close_leaf(m, A, sig_pa, c, want_prime=None):
    """Solve for the last prime q at exponent 2, given the rest of the leaf fixed.

    m       the fixed part times p^(2a), i.e. everything except q
    A       sigma of the fixed part times sigma(p^(2a)), i.e. S
    sig_pa  unused placeholder kept for call-site clarity; S is passed as A
    Returns a list of integer q solving the quadratic, unfiltered by primality.
    """
    S, T = A, 2 * m
    if S == T:                      # impossible: S odd, T even. Guard anyway.
        return []
    delta = S * S - 4 * (S - T) * (S - c)
    if delta < 0 or sieve_rejects(delta):
        return []
    root = isqrt(delta)
    if root * root != delta:
        return []
    out = []
    den = 2 * (T - S)
    for num in (S + root, S - root):
        if den and num % den == 0:
            q = num // den
            if q > 1:
                out.append(q)
    return out


def scan(m0, A0, c, p_lo, p_hi, exps, primes, verify):
    """Scan candidate intermediate primes p over [p_lo, p_hi].

    For each p in `primes` within range and each exponent 2a with a in `exps`,
    close the leaf for q at exponent 2 and hand any candidate to `verify`.
    Returns (candidates_tested, sieve_rejections, solutions).
    """
    from qp import sig
    tested = rejected = 0
    sols = []
    for p in primes:
        if p < p_lo:
            continue
        if p > p_hi:
            break
        for a in exps:
            pa = p ** (2 * a)
            m = m0 * pa
            S = A0 * sig(p, 2 * a)
            T = 2 * m
            if S == T:
                continue
            delta = S * S - 4 * (S - T) * (S - c)
            tested += 1
            if delta < 0 or sieve_rejects(delta):
                rejected += 1
                continue
            root = isqrt(delta)
            if root * root != delta:
                continue
            den = 2 * (T - S)
            for num in (S + root, S - root):
                if den and num % den == 0:
                    q = num // den
                    if q > p and verify(p, a, q):
                        sols.append((p, a, q))
    return tested, rejected, sols
