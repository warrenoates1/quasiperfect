"""Exhaustive search for odd squares n with sigma(n) = 2n + c and omega(n) = k.

c = 1 is the quasiperfect case. Any quasiperfect number is an odd square
(Cattaneo 1951), so n = prod p_i^(2 a_i) over odd primes p_i, and the search
branches on the primes in increasing order and on their exponents.

Soundness rests on these facts, each used as a prune or a finite bound:

  h(x) = sigma(x)/x is multiplicative, and h(p^(2a)) rises to p/(p-1) as a grows.
  A component whose exponent is not yet fixed ("open", exponent >= 2A) is carried
  as the interval [h(p^(2A)), p/(p-1)).

  (1) Target: h(n) = 2 + c/n, with 2 < 2 + c/n <= 2 + c/nmin.
  (2) Remaining primes only raise h: h(n) > lo(partial). The first remaining prime q
      raises it by at least 1 + 1/q + 1/q^2, and n >= nmin * q^2.
  (3) Upper bound: the r remaining primes, all > last, contribute at most the
      product of q/(q-1) over the r consecutive primes after last.
  (4) gcd(n, sigma(n)) = gcd(n, c): with c = 1 no prime of n divides sigma(p^(2a))
      for any fixed component (applied only when c = 1).
  (5) One open exponent, everything else fixed: sigma(n) = 2n + c = c (mod p^(2a))
      gives p^(2a) | sigma(n') + c(p - 1), where n' = n / p^(2a). Finite.
  (6) One unknown last prime q, everything else fixed (m, S = sigma(m),
      D = 2m - S > 0): q^(2a) (2m - D q) = c q + S - c, solved exactly per a.
  (7) All primes known: prod over open components of p/(p-1) times fixed parts has
      odd numerator and even denominator, so it never equals 2. Splitting open
      exponents therefore always ends: the interval eventually excludes the target.
"""
import sys
import time
from fractions import Fraction as F

sys.setrecursionlimit(1_000_000)

_SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]


def is_probable_prime(n):
    # Never rejects a prime, so a false "prime" only adds a branch: sound.
    if n < 2:
        return False
    for p in _SMALL:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime(n):
    q = n + 1 if n % 2 == 0 else n + 2
    if n < 2:
        return 2
    while not is_probable_prime(q):
        q += 2
    return q


def sig(p, e):
    return (p ** (e + 1) - 1) // (p - 1)


def isqrt_ceil_div(a, b):
    return -(-a // b)


class Search:
    def __init__(self, c=1, q_limit=10**7, log=None):
        self.c = c
        self.q_limit = q_limit   # max primes enumerated at one node before flagging
        self.nodes = 0
        self.leaf_checks = 0
        self.solutions = []
        self.unresolved = []
        self.log = log

    # comps: tuple of (p, a, fixed); exponent 2a, or >= 2a when not fixed.
    def bounds(self, comps):
        lo = F(1)
        hi = F(1)
        nmin = 1
        for p, a, fx in comps:
            pe = p ** (2 * a)
            h = F(sig(p, 2 * a), pe)
            lo *= h
            nmin *= pe
            hi *= h if fx else F(p, p - 1)
        return lo, hi, nmin

    def coprime_ok_fix(self, comps, p, a):
        if self.c != 1:
            return True
        s = sig(p, 2 * a)
        return all(s % q for q, _, _ in comps if q != p)

    def coprime_ok_add(self, comps, q):
        if self.c != 1:
            return True
        return all(sig(p, 2 * a) % q for p, a, fx in comps if fx)

    def record(self, n):
        # callers have already checked sigma(n) = 2n + c exactly from the factorisation
        self.solutions.append(n)

    def run(self, k):
        self.search((), k)

    def search(self, comps, r):
        self.nodes += 1
        if self.log and self.nodes % 200_000 == 0:
            self.log(self)
        c = self.c
        lo, hi, nmin = self.bounds(comps)
        last = comps[-1][0] if comps else 2
        opens = [i for i, t in enumerate(comps) if not t[2]]

        if r == 0:
            if not opens:
                self.leaf_checks += 1
                if lo * nmin == 2 * nmin + c:
                    self.record(nmin)
                return
            if hi <= 2:                       # (1): open parts stay strictly below hi
                return
            if lo > 2 + F(c, nmin):           # (1)
                return
            if len(opens) == 1:
                self.single_open(comps, opens[0])
                return
            self.refine(comps, opens, r)
            return

        # r >= 1
        if lo >= 2:
            # (1)+(2): lo(1 + 1/q) < 2 + c/(nmin q^2)  =>  q < c/(nmin lo)
            qb = F(c) / (nmin * lo)
            q = next_prime(last)
            while q < qb:                     # never entered when c = 1
                if self.coprime_ok_add(comps, q):
                    self.search(comps + ((q, 1, False),), r - 1)
                q = next_prime(q)
            return
        # (3)
        ub = hi
        q = last
        for _ in range(r):
            q = next_prime(q)
            ub *= F(q, q - 1)
        if ub <= 2:
            return
        if hi >= 2:
            self.refine(comps, opens, r)
            return
        if r == 1 and not opens:
            self.last_prime(comps, lo, nmin, last)
            return
        self.enumerate_next(comps, r, lo, hi, nmin, last)

    def refine(self, comps, opens, r):
        # split the open component with the smallest exponent floor (round robin)
        i = min(opens, key=lambda j: (comps[j][1], comps[j][0]))
        p, a, _ = comps[i]
        if self.coprime_ok_fix(comps, p, a):
            self.search(comps[:i] + ((p, a, True),) + comps[i + 1:], r)
        self.search(comps[:i] + ((p, a + 1, False),) + comps[i + 1:], r)

    def single_open(self, comps, i):
        # (5)
        c = self.c
        p, A, _ = comps[i]
        s_rest = 1
        n_rest = 1
        for j, (q, a, _) in enumerate(comps):
            if j != i:
                s_rest *= sig(q, 2 * a)
                n_rest *= q ** (2 * a)
        T = s_rest + c * (p - 1)
        v = 0
        while T % p == 0:
            T //= p
            v += 1
        for a in range(A, v // 2 + 1):
            self.leaf_checks += 1
            pe = p ** (2 * a)
            if s_rest * sig(p, 2 * a) == 2 * n_rest * pe + c:
                self.record(n_rest * pe)

    def last_prime(self, comps, lo, m, last):
        # (6): q^(2a) (2m - Dq) = c q + S - c, with q > last and q < 2m/D.
        c = self.c
        S = lo.numerator * (m // lo.denominator)
        assert F(S, m) == lo
        D = 2 * m - S
        assert D > 0
        qmax = (2 * m - 1) // D               # 2m - Dq >= 1
        qmin = last + 1
        if qmax < qmin:
            return
        a = 1
        while True:
            e = 2 * a
            if qmin ** e > c * qmax + S:      # g(q) >= q^e > RHS for every q
                break

            def P(q):
                return q ** e * (2 * m - D * q) - c * q - (S - c)

            peak = (e * 2 * m) // ((e + 1) * D)   # g rises below ~peak, falls above
            # rising branch: prove no root, else scan it (flagged)
            if qmin <= peak:
                g_lo = qmin ** e * (2 * m - D * qmin)
                if not g_lo > c * (peak + 1) + S:
                    self.unresolved.append(("rising-branch", comps, a))
            # falling branch: P strictly decreasing on [max(qmin, peak+1), qmax]
            L = max(qmin, peak + 1)
            R = qmax
            if L <= R and P(L) >= 0 >= P(R):
                while L < R:
                    mid = (L + R) // 2
                    if P(mid) > 0:
                        L = mid + 1
                    else:
                        R = mid
                if P(L) == 0:
                    self.leaf_checks += 1
                    self.record(m * L ** e)
            a += 1

    def enumerate_next(self, comps, r, lo, hi, nmin, last):
        c = self.c
        # (2): lo (1 + 1/q) < 2 + c/nmin  =>  q > lo / (2 - lo + c/nmin)
        bound = lo / (2 - lo + F(c, nmin))
        start = max(last, bound.numerator // bound.denominator - 1)
        q = next_prime(start)
        count = 0
        while True:
            # (3) with q as the smallest remaining prime
            ub = hi
            t = q
            ub *= F(t, t - 1)
            for _ in range(r - 1):
                t = next_prime(t)
                ub *= F(t, t - 1)
            if ub <= 2:
                break
            if (lo * F(q * q + q + 1, q * q) <= 2 + F(c, nmin * q * q)
                    and self.coprime_ok_add(comps, q)):
                self.search(comps + ((q, 1, False),), r - 1)
            count += 1
            if count > self.q_limit:
                self.unresolved.append(("q-range", comps, r, q))
                return
            q = next_prime(q)


def factorize(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def sigma_of(n):
    s = 1
    for p, e in factorize(n).items():
        s *= sig(p, e)
    return s


if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    c = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    for k in range(1, kmax + 1):
        t0 = time.time()
        s = Search(c=c, log=lambda s: print(f"  ... {s.nodes:,} nodes", flush=True))
        s.run(k)
        print(f"c={c} omega={k}: nodes={s.nodes:,} leaf checks={s.leaf_checks:,} "
              f"solutions={s.solutions} unresolved={len(s.unresolved)} "
              f"({time.time() - t0:.1f}s)", flush=True)
        for u in s.unresolved[:5]:
            print("   unresolved:", u)
