"""Compiled search for any node with two primes left to choose (r = 2).

Generalises qfast.py, which closed one such node for omega = 17. There the
last two primes were huge, so their own exponents could be ignored; here they
may be small (omega = 7 has nodes like 3^a 5^b 17^c 257^d q5^2 . q . q'), so
the exponents of q and q' are components of the box like any other.

Log form. Every component x with exponent 2b contributes
log f(x, b) = L(x) - d_x(b), with L(x) = -log(1 - 1/x) and
d_x(b) = -log(1 - x^-(2b+1)), so sigma(N) = 2N + c becomes

    L(q) + L(q') = K + sum_x d_x(b_x) + eps,   K = log 2 - sum_{node} L(p),
    0 < eps = log(1 + c/(2N)) <= c / (2 N_min).

A box gives each component an exponent either fixed (d exact) or open from a
floor (d in (0, d(floor)]). The search goes r = 2 (choose q), r = 1 (choose q'),
r = 0 (check), splitting the open component with the widest contribution when a
range or a box is too wide. All sums are double-double; node constants come from
exact rationals or 60-digit Decimal logs; every comparison carries slack S far
above the rounding error, so the search keeps too much and never too little.
Survivors are written out and settled by Search.search, the original exact code.

Prunes kept: the residue test when a node exponent is fixed (precomputed with
Search.coprime_ok_fix) and the mod-8 test when q's or q''s exponent is fixed.
Prunes dropped (coprimality with q and q', the residue test on q and q') only
shrink the tree, so dropping them is sound.
"""
import math
from decimal import Decimal, getcontext
from fractions import Fraction as F

import numpy as np
from numba import njit

from qp import is_probable_prime

S_SLACK = 1e-35
D_FLOOR = 1e-58          # bound on d for exponents past the table
W_Q2 = 32                # enumerate q' when the range holds at most this many integers
Q2_MAX = 4e18            # beyond this q' goes to the exact code (int64 limit)
W_LAST = 4096           # window allowed when q' is the only unknown left
SIEVE_SEG = 8192        # q is walked in sieved segments of this many integers

getcontext().prec = 80


# ---------------------------------------------------------------- double-double
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
    c = 134217729.0 * a
    ah = c - (c - a)
    al = a - ah
    c = 134217729.0 * b
    bh = c - (c - b)
    bl = b - bh
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


@njit(cache=True)
def dd_mul(ah, al, bh, bl):
    p, e = two_prod(ah, bh)
    e += ah * bl + al * bh
    h = p + e
    return h, e - (h - p)


@njit(cache=True)
def dd_recip(x):
    h = 1.0 / x
    p, pe = two_prod(h, x)
    r = (p - 1.0) + pe
    return h, -h * r


@njit(cache=True)
def neglog1m(xh, xl):
    """-log(1 - x) = sum x^k / k for 0 < x <= 1/7, double-double."""
    sh, sl = xh, xl
    ph, pl = xh, xl
    k = 1
    while True:
        k += 1
        ph, pl = dd_mul(ph, pl, xh, xl)
        t = ph / k
        if t < 1e-48:
            break
        # t as double-double: ph/k with its remainder
        th = t
        rh, rl = two_prod(th, float(k))
        tl = ((ph - rh) - rl + pl) / k
        sh, sl = dd_add(sh, sl, th, tl)
    return sh, sl


@njit(cache=True)
def L_of(x):
    yh, yl = dd_recip(float(x))
    return neglog1m(yh, yl)


@njit(cache=True)
def d_of(x, b):
    """d_x(b) = -log(1 - x^-(2b+1)) as double-double."""
    yh, yl = dd_recip(float(x))
    ph, pl = yh, yl
    for _ in range(2 * b):
        ph, pl = dd_mul(ph, pl, yh, yl)
    if ph < 1e-300:
        return 0.0, 0.0
    return neglog1m(ph, pl)


# ---------------------------------------------------------------- primes
@njit(cache=True)
def mulmod(a, b, n):
    return (((a * (b >> 19)) % n << 19) + a * (b & 524287)) % n


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
    for a in (2, 3, 5, 7, 11, 13, 17):
        x = 1
        base = a % n
        e = d
        while e:
            if e & 1:
                x = mulmod(x, base, n)
            base = mulmod(base, base, n)
            e >>= 1
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
def mod8_ok(p, b):
    # sigma(p^2b) = 1 + b(p+1) mod 8 must be 1 or 3 (residue lemma, c = 1 only)
    r = (1 + b * (p + 1)) % 8
    return r == 1 or r == 3


# ---------------------------------------------------------------- kernel
@njit(cache=True)
def box_sum(n0, a, fx, dtab, dtl, amax):
    """Node part of the box: (Dlo, Dlo_l, Dhi, Dhi_l, widest index, width, lnN floor)."""
    Dlo = 0.0
    Dll = 0.0
    Dhi = 0.0
    Dhl = 0.0
    best = -1
    bw = 0.0
    for i in range(n0):
        if a[i] > amax[i]:
            Dhi, Dhl = dd_add(Dhi, Dhl, D_FLOOR, 0.0)
            continue
        d = dtab[i, a[i]]
        dl = dtl[i, a[i]]
        Dhi, Dhl = dd_add(Dhi, Dhl, d, dl)
        if fx[i]:
            Dlo, Dll = dd_add(Dlo, Dll, d, dl)
        elif d > bw:
            bw = d
            best = i
    return Dlo, Dll, Dhi, Dhl, best, bw


@njit(cache=True)
def run_node(n0, primes, floor0, fixed0, dtab, dtl, allowed, amax, lnp,
             K, KL, c, use_mod8, last, sp, out, stats, qfrom=0, qto=1 << 62):
    """Search one r = 2 node, for the smaller prime q in [qfrom, qto] only (so a
    dear node can be split across processes). Candidates go to out as rows:
    q, q', bq, fq, bq', fq', a_0..a_{n0-1}, fx_0..fx_{n0-1}. Returns count."""
    maxd = 1 << 14
    # stack: level (2,1,0), node a/fx, q, bq, fq, q', bq', fq'
    st_lv = np.empty(maxd, np.int64)
    st_a = np.empty((maxd, n0), np.int64)
    st_f = np.empty((maxd, n0), np.int64)
    st_q = np.empty((maxd, 6), np.int64)
    top = 0
    st_lv[0] = 2
    for i in range(n0):
        st_a[0, i] = floor0[i]
        st_f[0, i] = fixed0[i]
    for j in range(6):
        st_q[0, j] = 0
    top = 1
    nout = 0
    a = np.empty(n0, np.int64)
    fx = np.empty(n0, np.int64)
    seg = np.empty(SIEVE_SEG, np.bool_)
    while top > 0:
        top -= 1
        stats[0] += 1
        lv = st_lv[top]
        for i in range(n0):
            a[i] = st_a[top, i]
            fx[i] = st_f[top, i]
        if lv == 3:
            # cursor over the q-range: sieve one segment, leave a cursor for the rest
            stats[0] -= 1
            lo_s = st_q[top, 0]
            qmax = st_q[top, 1]
            hi_s = min(lo_s + SIEVE_SEG - 1, qmax)
            if hi_s < qmax:
                st_q[top, 0] = hi_s + 1
                top += 1                          # cursor stays, primes go above it
            for j in range(hi_s - lo_s + 1):
                seg[j] = True
            for k in range(sp.shape[0]):
                p = sp[k]
                if p * p > hi_s:
                    break
                st = max(p * p, ((lo_s + p - 1) // p) * p)
                for m in range(st, hi_s + 1, p):
                    seg[m - lo_s] = False
            for j in range(hi_s - lo_s, -1, -1):  # push high first, so low q runs first
                qq = lo_s + j
                if seg[j] and qq > 2:
                    if top >= maxd:
                        stats[3] += 1
                        return nout
                    st_lv[top] = 1
                    for i in range(n0):
                        st_a[top, i] = a[i]
                        st_f[top, i] = fx[i]
                    st_q[top, 0] = qq
                    st_q[top, 1] = 1
                    st_q[top, 2] = 0
                    st_q[top, 3] = 0
                    st_q[top, 4] = 0
                    st_q[top, 5] = 0
                    top += 1
            continue
        q = st_q[top, 0]
        bq = st_q[top, 1]
        fq = st_q[top, 2]
        q2 = st_q[top, 3]
        bq2 = st_q[top, 4]
        fq2 = st_q[top, 5]
        Dlo, Dll, Dhi, Dhl, best, bw = box_sum(n0, a, fx, dtab, dtl, amax)
        lnN = 0.0
        for i in range(n0):
            lnN += 2.0 * a[i] * lnp[i]
        lnM = lnN                              # node part only
        # widest open among q, q' exponents (index n0 = q, n0 + 1 = q')
        if lv <= 1:
            dh, dl = d_of(q, bq)
            Dhi, Dhl = dd_add(Dhi, Dhl, dh, dl)
            if fq:
                Dlo, Dll = dd_add(Dlo, Dll, dh, dl)
            elif dh > bw:
                bw = dh
                best = n0
            lnN += 2.0 * bq * math.log(q)
        if lv == 0:
            dh, dl = d_of(q2, bq2)
            Dhi, Dhl = dd_add(Dhi, Dhl, dh, dl)
            if fq2:
                Dlo, Dll = dd_add(Dlo, Dll, dh, dl)
            elif dh > bw:
                bw = dh
                best = n0 + 1
            lnN += 2.0 * bq2 * math.log(q2)
        # eps bound: c / (2 N_min), with q, q' at their least possible size
        if lv == 2:
            lnN += 4.0 * math.log(last + 1.0)
        elif lv == 1:
            lnN += 2.0 * math.log(q + 1.0)
        epsmax = 0.5 * c * math.exp(-lnN) * 1.000001
        split = False
        if lv == 2:
            # q <= q' gives 2 L(q) >= K + Dlo, so q <= qmax
            h, l = dd_add(K, KL, Dlo, Dll)
            X = 0.5 * (h + l) - S_SLACK
            if X <= 0.0 and best >= 0:
                split = True
            else:
                if X <= 0.0:
                    # all fixed, node abundancy >= 2: only eps can balance, and
                    # 1/q < eps <= c/(2 M q^4) gives q^3 < c/(2M)
                    qmax = int((0.5 * c * math.exp(-lnM)) ** (1.0 / 3.0) * 1.000001) + 2
                else:
                    qmax = int(1.0 / (-math.expm1(-X))) + 2
                # q' > q gives log f(q,1) < K + Dhi + eps, so q >= qlow
                h, l = dd_add(K, KL, Dhi, Dhl)
                Y = h + l + epsmax + S_SLACK
                if Y <= 0.0:
                    continue
                if Y > 30.0:                      # huge eps (c > 1 only): no lower bound
                    qlow = last + 1
                else:
                    yy = math.expm1(Y)
                    qlow = max(int((1.0 + math.sqrt(1.0 + 4.0 * yy)) / (2.0 * yy)) - 2, last + 1)
                qlow = max(qlow, qfrom)
                qmax = min(qmax, qto)
                if qmax < qlow:
                    continue
                if best >= 0 and qmax > 5000 and qmax > 8 * qlow:
                    split = True
                elif qmax >= (1 << 38):
                    stats[2] += 1
                    continue
                else:
                    # walk primes q in [qlow, qmax] by a cursor state
                    if top >= maxd:
                        stats[3] += 1
                        return nout
                    st_lv[top] = 3
                    for i in range(n0):
                        st_a[top, i] = a[i]
                        st_f[top, i] = fx[i]
                    st_q[top, 0] = qlow
                    st_q[top, 1] = qmax
                    top += 1
                    continue
        elif lv == 1:
            Lh, Ll = L_of(q)
            # R = L(q') - d_q' = K + D + eps - L(q); q' needs log f(q',1) <= Rhi, L(q') > Rlo
            h, l = dd_add(K, KL, Dhi, Dhl)
            h, l = dd_add(h, l, -Lh, -Ll)
            Rhi = h + l + epsmax + S_SLACK
            h, l = dd_add(K, KL, Dlo, Dll)
            h, l = dd_add(h, l, -Lh, -Ll)
            Rlo = h + l - S_SLACK
            if Rhi <= 0.0:
                continue
            # bounds as floats first: q' can be far beyond int64 when R is tiny
            if Rhi > 30.0:                        # huge eps (c > 1 only): no lower bound
                flo = 0.0
            else:
                yy = math.expm1(Rhi)
                flo = (1.0 + math.sqrt(1.0 + 4.0 * yy)) / (2.0 * yy)
            if Rlo <= 0.0 and best < 0:
                # all fixed and K + D - L(q) <= 0: 1/q' < R <= eps <= c/(2 N1 q'^2)
                fhi = 0.5 * c * math.exp(-(lnM + 2.0 * bq * math.log(q))) * 1.000001
            elif Rlo <= 0.0:
                fhi = 1e300
            else:
                fhi = 1.0 / (-math.expm1(-Rlo))
            if fhi < flo - 4.0:
                continue
            if fhi > Q2_MAX and best >= 0:
                split = True
                lo_i = 0
                hi_i = 0
            elif fhi > Q2_MAX:
                lo_i = 0
                hi_i = 1 << 62                    # forces the hand-off below
            else:
                # doubles carry ~16 digits: widen by a relative margin as well
                lo_i = max(int(flo * (1.0 - 4e-15)) - 2, q + 1)
                hi_i = int(fhi * (1.0 + 4e-15)) + 2
                if hi_i < lo_i:
                    continue
            if split:
                pass
            elif hi_i - lo_i > W_Q2 and best >= 0:
                split = True
            elif hi_i - lo_i > W_LAST:
                # nothing left open but q': hand the node to the exact code,
                # whose last-prime routine solves for q' directly (q' = 0 row)
                if nout < out.shape[0]:
                    out[nout, 0] = q
                    out[nout, 1] = 0
                    out[nout, 2] = bq
                    out[nout, 3] = fq
                    out[nout, 4] = 0
                    out[nout, 5] = 0
                    for i in range(n0):
                        out[nout, 6 + i] = a[i]
                        out[nout, 6 + n0 + i] = fx[i]
                nout += 1
                stats[4] += 1
                continue
            else:
                for qq in range(lo_i, hi_i + 1):
                    # below 2^38 test primality here; above it the window is a
                    # few integers wide and the exact code tests any survivor
                    if qq % 2 == 1 and (qq >= (1 << 38) or isprime(qq)):
                        if top >= maxd:
                            stats[3] += 1
                            return nout
                        st_lv[top] = 0
                        for i in range(n0):
                            st_a[top, i] = a[i]
                            st_f[top, i] = fx[i]
                        st_q[top, 0] = q
                        st_q[top, 1] = bq
                        st_q[top, 2] = fq
                        st_q[top, 3] = qq
                        st_q[top, 4] = 1
                        st_q[top, 5] = 0
                        top += 1
                continue
        else:
            Lh, Ll = L_of(q)
            L2h, L2l = L_of(q2)
            h, l = dd_add(Lh, Ll, L2h, L2l)
            # need K + Dlo <= L + L' <= K + Dhi + eps
            u, ul = dd_add(h, l, -K, -KL)
            v, vl = dd_add(u, ul, -Dlo, -Dll)
            if v + vl < -S_SLACK:
                continue
            v, vl = dd_add(u, ul, -Dhi, -Dhl)
            if v + vl > epsmax + S_SLACK:
                continue
            if best < 0 or (Dhi - Dlo) + (Dhl - Dll) <= epsmax + 2 * S_SLACK:
                if nout < out.shape[0]:
                    out[nout, 0] = q
                    out[nout, 1] = q2
                    out[nout, 2] = bq
                    out[nout, 3] = fq
                    out[nout, 4] = bq2
                    out[nout, 5] = fq2
                    for i in range(n0):
                        out[nout, 6 + i] = a[i]
                        out[nout, 6 + n0 + i] = fx[i]
                nout += 1
                continue
            split = True
        if not split:
            continue
        # split the widest open component: raised child, then fixed child
        if top + 2 > maxd:
            stats[3] += 1
            return nout
        stats[1] += 1
        for rep in range(2):
            st_lv[top] = lv
            for i in range(n0):
                st_a[top, i] = a[i]
                st_f[top, i] = fx[i]
            st_q[top, 0] = q
            st_q[top, 1] = bq
            st_q[top, 2] = fq
            st_q[top, 3] = q2
            st_q[top, 4] = bq2
            st_q[top, 5] = fq2
            if best < n0:
                if rep == 0:
                    st_a[top, best] = a[best] + 1
                else:
                    if not allowed[best, a[best]]:
                        break
                    st_f[top, best] = 1
            elif best == n0:
                if rep == 0:
                    st_q[top, 1] = bq + 1
                else:
                    if use_mod8 and not mod8_ok(q, bq):
                        break
                    st_q[top, 2] = 1
            else:
                if rep == 0:
                    st_q[top, 4] = bq2 + 1
                else:
                    if use_mod8 and not mod8_ok(q2, bq2):
                        break
                    st_q[top, 5] = 1
            top += 1
    return nout


# ---------------------------------------------------------------- Python side
def _dec_neglog1m(x):
    """-log(1 - x) for a Fraction 0 < x < 1, to 80 digits."""
    return -(Decimal(1) - Decimal(x.numerator) / Decimal(x.denominator)).ln()


def _dd(dec):
    h = float(dec)
    return h, float(dec - Decimal(h))


_SP = None


def small_primes():
    global _SP
    if _SP is None:
        n = 1 << 20
        s = np.ones(n, bool)
        s[:2] = False
        for i in range(2, int(n ** 0.5) + 1):
            if s[i]:
                s[i * i::i] = False
        _SP = np.nonzero(s)[0].astype(np.int64)
    return _SP


def prepare(comps, search):
    n0 = len(comps)
    amax = []
    for p, a, fx in comps:
        m = a
        if not fx:
            while F(1, p ** (2 * (m + 1) + 1)) > F(1, 10 ** 60):
                m += 1
        amax.append(m)
    A = max(amax) + 2
    dtab = np.zeros((n0, A))
    dtl = np.zeros((n0, A))
    allowed = np.zeros((n0, A), np.bool_)
    for i, (p, a0, fx) in enumerate(comps):
        for a in range(a0, amax[i] + 1):
            dtab[i, a], dtl[i, a] = _dd(_dec_neglog1m(F(1, p ** (2 * a + 1))))
            allowed[i, a] = search.coprime_ok_fix(comps, p, a)
    Ksum = Decimal(2).ln()
    for p, _, _ in comps:
        Ksum -= _dec_neglog1m(F(1, p))
    K, KL = _dd(Ksum)
    return dict(n0=n0, primes=np.array([p for p, _, _ in comps], np.int64),
                floor0=np.array([a for _, a, _ in comps], np.int64),
                fixed0=np.array([1 if fx else 0 for _, _, fx in comps], np.int64),
                dtab=dtab, dtl=dtl, allowed=allowed, amax=np.array(amax, np.int64),
                lnp=np.array([math.log(p) for p, _, _ in comps]), K=K, KL=KL)


def search_node(comps, search, cap=1_000_000, qfrom=0, qto=1 << 62):
    """Run one r = 2 node; settle survivors exactly. Returns (kernel_stats, rows)."""
    P = prepare(comps, search)
    last = max(p for p, _, _ in comps)
    n0 = P["n0"]
    out = np.zeros((cap, 6 + 2 * n0), np.int64)
    stats = np.zeros(5, np.int64)
    use_mod8 = 1 if search.c == 1 else 0
    m = run_node(n0, P["primes"], P["floor0"], P["fixed0"], P["dtab"], P["dtl"],
                 P["allowed"], P["amax"], P["lnp"], P["K"], P["KL"], float(search.c),
                 use_mod8, last, small_primes(), out, stats, qfrom, qto)
    rows = out[:min(m, cap)]
    if m > cap or stats[2] or stats[3]:
        search.unresolved.append(("qnode", comps, (qfrom, qto), int(m), stats.tolist()))
    for row in rows:
        q, q2, bq, fq, bq2, fq2 = (int(v) for v in row[:6])
        if q in search.exclude or q2 in search.exclude:
            continue
        new = tuple((p, int(row[6 + i]), bool(row[6 + n0 + i]))
                    for i, (p, _, _) in enumerate(comps))
        if q2 >= (1 << 38) and not is_probable_prime(q2):
            continue                              # kernel only tests primality below 2^38
        if q2 == 0:                               # q' left to Search's last_prime
            search.search(new + ((q, bq, bool(fq)),), 1)
        else:
            search.search(new + ((q, bq, bool(fq)), (q2, bq2, bool(fq2))), 0)
    return stats, len(rows)


class NodeSearch:
    """Mixin: route every r = 2 enumeration through the compiled kernel."""

    def enumerate_next(self, comps, r, lo, hi, nmin, last):
        if r != 2:
            return super().enumerate_next(comps, r, lo, hi, nmin, last)
        st, n = search_node(comps, self)
        self.kstates = getattr(self, "kstates", 0) + int(st[0])
        self.kcands = getattr(self, "kcands", 0) + n
