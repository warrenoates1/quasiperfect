"""Compiled search of one r = 2 node whose next-prime range is too long for Search.

Built for the single node omega = 17 (3 and 5 barred) left unresolved
("q-range", stopped by q_limit at q = 571,287,683): fifteen components, nine
with open exponents, two primes q < q' still to choose. Search spends ~8,000
states per q there, and the remaining range holds ~27 million primes.

Method. Write the abundancy in logarithms. With hi the node's upper bound (open
exponents at infinity) and K = log(2/hi) > 0, an open component p with exponent
2a contributes log(p/(p-1)) - d_p(a), where d_p(a) = -log(1 - p^-(2a+1)), so

    log f(q) + log f(q') = K + sum_p d_p + log(1 + 1/(2N)),

with f(x) = sigma(x^2b)/x^2b in [log1p((x+1)/x^2), -log1p(-1/x)) for every b >= 1.
For each q the search keeps a box over the open exponents (each either fixed at
a or open at >= a, contributing [d_p(a), d_p(a)] or [0, d_p(a)]), turns the box
into the range of q' it allows, splits the open component with the widest
contribution until that range is short, then tries each prime q' in it and
splits on until the box is narrower than the target the pair (q, q') leaves.

Sums are in double-double arithmetic (~32 significant digits; the big terms
are ~1.8e-9, so absolute error ~1e-40), every table value is rounded from an
exact rational, and every comparison carries a slack S = 1e-35 far above the
rounding error, so the search can only keep too much, never prune a real
solution. (A first plain-double version needed S = 2e-23 and left ~0.56 near
misses per q for the exact code, ~75 ms each.) Whatever survives to the bottom is written
out as a candidate and settled by the original exact code (Search.search from
that box with q and q' as open components). The term log(1 + 1/(2N)) is below
1e-140 and inside S.

Prunes of Search not used here (coprimality with q and q', the residue test on
q and q') only shrink its tree, so dropping them is sound. The residue/mod-8
test on fixing an open exponent is kept, precomputed per (p, a) with
coprime_ok_fix against the node's own primes.
"""
import math
import numpy as np
from numba import njit

S_SLACK = 1e-35
W_MAX = 32
D_FLOOR = 1e-50          # exponents whose d_p(a) is below this are left open


@njit(cache=True)
def mulmod(a, b, n):                      # n < 2^38
    return (((a * (b >> 19)) % n << 19) + a * (b & 524287)) % n


@njit(cache=True)
def powmod(a, e, n):
    r = 1
    a %= n
    while e:
        if e & 1:
            r = mulmod(r, a, n)
        a = mulmod(a, a, n)
        e >>= 1
    return r


@njit(cache=True)
def isprime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):         # deterministic below 3.4e14
        x = powmod(a, d, n)
        if x == 1 or x == n - 1:
            continue
        comp = True
        for _ in range(s - 1):
            x = mulmod(x, x, n)
            if x == n - 1:
                comp = False
                break
        if comp:
            return False
    return True


@njit(cache=True)
def dd_add(ah, al, bh, bl):
    s = ah + bh
    v = s - ah
    e = (ah - (s - v)) + (bh - v) + al + bl
    h = s + e
    return h, e - (h - s)


@njit(cache=True)
def two_prod(a, b):
    p = a * b
    # Dekker split
    c = 134217729.0 * a
    ah = c - (c - a)
    al = a - ah
    c = 134217729.0 * b
    bh = c - (c - b)
    bl = b - bh
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


@njit(cache=True)
def dd_recip(x):
    # 1/x for an integer x < 2^53, as a double-double
    h = 1.0 / x
    p, pe = two_prod(h, x)          # h x = 1 + r
    r = (p - 1.0) + pe
    return h, -h * r


@njit(cache=True)
def log_f_bounds(x):
    # (log(1 + 1/x + 1/x^2), -log(1 - 1/x)) as double-doubles, x ~ 1e8..1e11
    yh, yl = dd_recip(float(x))
    # y^2 and u^2 are ~3e-18, so they need their own low parts; cubes and up
    # are below 1e-26 and a double is ample.
    s2, s2l = two_prod(yh, yh)
    s2l += 2.0 * yh * yl
    y3 = s2 * yh
    y4 = s2 * s2
    # u = y + y^2; log1p(u) = u - u^2/2 + u^3/3 - u^4/4
    uh, ul = dd_add(yh, yl, s2, s2l)
    u2, u2l = two_prod(uh, uh)
    u2l += 2.0 * uh * ul
    u3 = u2 * uh
    u4 = u2 * u2
    loh, lol = dd_add(uh, ul, -0.5 * u2, -0.5 * u2l)
    loh, lol = dd_add(loh, lol, u3 / 3.0 - 0.25 * u4, 0.0)
    # -log1p(-y) = y + y^2/2 + y^3/3 + y^4/4
    hih, hil = dd_add(yh, yl, 0.5 * s2, 0.5 * s2l)
    hih, hil = dd_add(hih, hil, y3 / 3.0 + 0.25 * y4, 0.0)
    return loh, lol, hih, hil


@njit(cache=True)
def run_q(q, K, KL, dtab, dtl, allowed, floor, amax, out, nout, stats):
    """Search one q. Candidates appended to out[nout:] as (q, q', a_0..a_k, fx_0..fx_k)."""
    nop = floor.shape[0]
    maxd = 4096
    st_a = np.empty((maxd, nop), np.int64)
    st_f = np.empty((maxd, nop), np.int64)
    st_p = np.empty(maxd, np.int64)
    top = 0
    for i in range(nop):
        st_a[0, i] = floor[i]
        st_f[0, i] = 0
    st_p[0] = 0
    top = 1
    qlo, qlol, qhi, qhil = log_f_bounds(q)
    a = np.empty(nop, np.int64)
    fx = np.empty(nop, np.int64)
    while top > 0:
        top -= 1
        stats[0] += 1
        for i in range(nop):
            a[i] = st_a[top, i]
            fx[i] = st_f[top, i]
        qp = st_p[top]
        Dlo = 0.0
        Dll = 0.0
        Dhi = 0.0
        Dhl = 0.0
        best = -1
        bestd = 0.0
        for i in range(nop):
            if a[i] > amax[i]:
                Dhi, Dhl = dd_add(Dhi, Dhl, D_FLOOR, 0.0)
                continue
            d = dtab[i, a[i]]
            dl = dtl[i, a[i]]
            if fx[i]:
                Dlo, Dll = dd_add(Dlo, Dll, d, dl)
                Dhi, Dhl = dd_add(Dhi, Dhl, d, dl)
            else:
                Dhi, Dhl = dd_add(Dhi, Dhl, d, dl)
                if d > bestd:
                    bestd = d
                    best = i
        if qp == 0:
            h, l = dd_add(K, KL, Dhi, Dhl)
            h, l = dd_add(h, l, -qlo, -qlol)
            Rhi = h + l + S_SLACK
            h, l = dd_add(K, KL, Dlo, Dll)
            h, l = dd_add(h, l, -qhi, -qhil)
            Rlo = h + l - S_SLACK
            if Rhi <= 0.0:
                continue
            y = math.expm1(Rhi)
            Lb = (1.0 + math.sqrt(1.0 + 4.0 * y)) / (2.0 * y)
            lo_i = max(int(math.floor(Lb)) - 1, q + 1)
            if Rlo <= 0.0:
                hi_i = 1 << 40
            else:
                hi_i = int(math.ceil(1.0 / (-math.expm1(-Rlo)))) + 1
            if hi_i < lo_i:
                continue
            if hi_i - lo_i > W_MAX and best >= 0:
                pass                                   # split below
            else:
                if hi_i - lo_i > 10_000_000 or hi_i >= (1 << 38):
                    stats[2] += 1                      # unresolved
                    continue
                for qq in range(lo_i, hi_i + 1):
                    if qq % 2 == 1 and isprime(qq):
                        if top >= maxd:
                            stats[3] += 1
                            return nout
                        for i in range(nop):
                            st_a[top, i] = a[i]
                            st_f[top, i] = fx[i]
                        st_p[top] = qq
                        top += 1
                continue
        else:
            plo, plol, phi, phil = log_f_bounds(qp)
            # need Dhi >= qlo + plo - K and Dlo <= qhi + phi - K (within S)
            h, l = dd_add(Dhi, Dhl, K, KL)
            h, l = dd_add(h, l, -qlo, -qlol)
            h, l = dd_add(h, l, -plo, -plol)
            if h + l < -S_SLACK:
                continue
            h, l = dd_add(Dlo, Dll, K, KL)
            h, l = dd_add(h, l, -qhi, -qhil)
            h, l = dd_add(h, l, -phi, -phil)
            if h + l > S_SLACK:
                continue
            tw = (qhi - qlo) + (phi - plo) + (qhil - qlol) + (phil - plol) + 2 * S_SLACK
            if best < 0 or (Dhi - Dlo) + (Dhl - Dll) <= tw:
                if nout < out.shape[0]:
                    out[nout, 0] = q
                    out[nout, 1] = qp
                    for i in range(nop):
                        out[nout, 2 + i] = a[i]
                        out[nout, 2 + nop + i] = fx[i]
                nout += 1
                continue
        # split the open component with the widest contribution
        if top + 2 > maxd:
            stats[3] += 1
            return nout
        for i in range(nop):
            st_a[top, i] = a[i]
            st_f[top, i] = fx[i]
        st_a[top, best] = a[best] + 1
        st_p[top] = qp
        top += 1
        if allowed[best, a[best]]:
            for i in range(nop):
                st_a[top, i] = a[i]
                st_f[top, i] = fx[i]
            st_f[top, best] = 1
            st_p[top] = qp
            top += 1
        stats[1] += 1
    return nout


@njit(cache=True)
def run_block(qs, K, KL, dtab, dtl, allowed, floor, amax, out, stats):
    nout = 0
    for j in range(qs.shape[0]):
        nout = run_q(qs[j], K, KL, dtab, dtl, allowed, floor, amax, out, nout, stats)
    return nout


def prepare(comps, search):
    """Tables for a node's open components: (opens, K, KL, dtab, dtl, allowed, floor, amax).

    K = log(2/hi) and d_p(a) = -log(1 - p^-(2a+1)) are summed exactly as
    rationals from their series (truncation below 1e-43) and split into
    double-double (high, low) pairs."""
    from fractions import Fraction as F
    lo, hi, _ = search.bounds(comps)

    def dd(x):
        h = float(x)
        return h, float(x - F(h))

    u = (2 - hi) / hi                       # K = log1p(u), u ~ 1.8e-9
    K, KL = dd(u - u**2 / 2 + u**3 / 3 - u**4 / 4 + u**5 / 5)
    opens = [(p, a) for p, a, fx in comps if not fx]
    amax = []
    for p, a in opens:
        m = a
        while p ** -(2 * (m + 1) + 1) > D_FLOOR:
            m += 1
        amax.append(m)
    A = max(amax) + 2
    dtab = np.zeros((len(opens), A))
    dtl = np.zeros((len(opens), A))
    allowed = np.zeros((len(opens), A), np.bool_)
    for i, (p, a0) in enumerate(opens):
        for a in range(a0, amax[i] + 1):
            t = F(1, p ** (2 * a + 1))      # t <= 7^-11 ~ 5e-10
            dtab[i, a], dtl[i, a] = dd(t + t**2 / 2 + t**3 / 3 + t**4 / 4 + t**5 / 5)
            allowed[i, a] = search.coprime_ok_fix(comps, p, a)
    floor = np.array([a for _, a in opens], np.int64)
    return opens, K, KL, dtab, dtl, allowed, floor, np.array(amax, np.int64)


def primes_between(a, b):
    """Primes in [a, b] by a segmented sieve."""
    import numpy as np
    r = int(math.isqrt(b)) + 1
    small = np.ones(r + 1, bool)
    small[:2] = False
    for i in range(2, int(r ** 0.5) + 1):
        if small[i]:
            small[i * i::i] = False
    ps = np.nonzero(small)[0]
    seg = np.ones(b - a + 1, bool)
    for p in ps:
        st = max(p * p, (a + p - 1) // p * p)
        seg[st - a::p] = False
    if a <= 1:
        seg[:2 - a] = False
    return (np.nonzero(seg)[0] + a).astype(np.int64)
