"""Search the last two primes with the smaller one held as a block of primes.

Why: at a node with r = 2 and several open exponents, Search.enumerate_next
tries the next prime q one at a time, and for each q re-refines every open
exponent. When hi is within ~1e-9 of 2 the q-range holds tens of millions of
primes and each costs thousands of states (omega = 17, 3 and 5 barred, left one
such node unresolved: ~2 s per q, ~15,000 core-hours). Neighbouring q give
almost the same tree, so here q is an interval [qa, qb] of candidate primes,
carried through the search with bounds valid for every prime in it:

    f(q, a) = sigma(q^2a)/q^2a  >=  f(qb, 1)      (decreasing in q, rising in a)
    f(q, a)                      <   qa/(qa - 1)
    q^2a                         >=  qa^2

Every prune applied to the block is therefore implied by the same prune for
each q in it, so the block search visits a superset of what the per-q search
would have kept. The block is split in half whenever its own spread is the
larger source of slack (or an exact q is needed), and a block holding a single
prime is handed back to Search.search with that prime as an ordinary open
component, so every leaf is decided by the original, tested code.

Two prunes of the per-q path are dropped as unnecessary for soundness (they
are re-applied once q is exact): the coprimality of sigma(p^2a) with q when a
component is fixed, and the per-q lower check on lo f(q).
"""
from fractions import Fraction as F

from qp import Search, is_probable_prime, next_prime


def prime_at_least(n):
    return n if n == 2 or (n > 2 and n % 2 and is_probable_prime(n)) else next_prime(n)


def prime_at_most(n, floor):
    if n % 2 == 0:
        n -= 1
    while n >= floor and not is_probable_prime(n):
        n -= 2
    return n


class BlockSearch(Search):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.splits = 0

    def enumerate_next(self, comps, r, lo, hi, nmin, last):
        if r != 2:
            return super().enumerate_next(comps, r, lo, hi, nmin, last)
        c = self.c
        bound = lo / (2 - lo + F(c, nmin))
        qa = next_prime(max(last, bound.numerator // bound.denominator - 1))
        self.block_root(comps, qa, lo, hi)

    def block_root(self, comps, qa, lo=None, hi=None):
        if hi is None:
            lo, hi, _ = self.bounds(comps)
        # (3) at r = 2: q < q2 so hi (q/(q-1))^2 > hi q/(q-1) q2/(q2-1) > 2 is
        # necessary; qend is the largest integer satisfying it.
        n, d = hi.numerator, hi.denominator
        s = (2 * d / n) ** 0.5
        qend = int(s / (s - 1)) + 2
        while not n * qend * qend > 2 * d * (qend - 1) ** 2:
            qend -= 1
        while n * (qend + 1) ** 2 > 2 * d * qend ** 2:
            qend += 1
        if qend >= qa:
            self.bsearch(comps, (qa, qend), 1)

    def bbounds(self, comps, blk):
        lo, hi, n = self.bounds(comps)
        qa, qb = blk
        return (lo * F(qb * qb + qb + 1, qb * qb), hi * F(qa, qa - 1), n * qa * qa,
                lo, hi)

    def bsearch(self, comps, blk, r):
        self.nodes += 1
        self.current = comps
        if self.log and self.nodes % 200_000 == 0:
            self.log(self)
        c = self.c
        qa, qb = blk
        if r == 0:
            qb = min(qb, comps[-1][0] - 1)      # q < q', q' was added last
        qa = prime_at_least(qa)
        qb = prime_at_most(qb, qa)
        if qa > qb:
            return
        if qa == qb:
            return self.finalize(comps, qa, r)
        blk = (qa, qb)
        lo, hi, nmin, clo, chi = self.bbounds(comps, blk)
        opens = [i for i, t in enumerate(comps) if not t[2]]

        if r == 0:
            if hi <= 2 or lo > 2 + F(c, nmin):
                return
            return self.choose(comps, blk, r, clo, chi, opens)

        # r == 1: q' > q >= qa remains
        if lo >= 2:
            if c == 1:                         # lo(1 + 1/q') < 2 + 1/(n q'^2) impossible
                return
            return self.split(comps, blk, r)
        q2 = next_prime(qa)
        if hi * F(q2, q2 - 1) <= 2:
            return
        if hi >= 2:
            return self.choose(comps, blk, r, clo, chi, opens)
        width = F(2) / (2 - hi) - lo / (2 - lo + F(c, nmin))
        if width > self.refine_width:
            return self.choose(comps, blk, r, clo, chi, opens)
        bound = lo / (2 - lo + F(c, nmin))
        q = next_prime(max(qa, bound.numerator // bound.denominator - 1))
        count = 0
        while hi * F(q, q - 1) > 2:
            if (q not in self.exclude
                    and lo * F(q * q + q + 1, q * q) <= 2 + F(c, nmin * q * q)
                    and self.coprime_ok_add(comps, q)):
                self.bsearch(comps + ((q, 1, False),), (qa, min(qb, q - 1)), 0)
            count += 1
            if count > self.q_limit:
                self.unresolved.append(("block-q2", comps, blk, q))
                return
            q = next_prime(q)

    def choose(self, comps, blk, r, clo, chi, opens):
        qa, qb = blk
        block_gap = F(qa, qa - 1) / F(qb * qb + qb + 1, qb * qb)
        if opens and chi / clo > block_gap:
            i = min(opens, key=lambda j: (comps[j][1], comps[j][0]))
            p, a, _ = comps[i]
            if self.coprime_ok_fix(comps, p, a):
                self.bsearch(comps[:i] + ((p, a, True),) + comps[i + 1:], blk, r)
            self.bsearch(comps[:i] + ((p, a + 1, False),) + comps[i + 1:], blk, r)
        else:
            self.split(comps, blk, r)

    def split(self, comps, blk, r):
        self.splits += 1
        qa, qb = blk
        mid = (qa + qb) // 2
        self.bsearch(comps, (qa, mid), r)
        self.bsearch(comps, (mid + 1, qb), r)

    def finalize(self, comps, q, r):
        if q in self.exclude or not self.coprime_ok_add(comps, q):
            return
        if r == 1:
            self.search(comps + ((q, 1, False),), 1)
        else:                                   # keep q' last
            self.search(comps[:-1] + ((q, 1, False), comps[-1]), 0)
