"""Measure the double-double rounding error of qnode's building blocks against
80-digit Decimal, as a check on the proved bound: |error| <= REL * magnitude,
REL = 2048 u^2 (u = 2^-53). Reports the worst observed error in units of u^2 * |value|."""
import random
from decimal import Decimal, getcontext
from fractions import Fraction as F
import qnode

getcontext().prec = 90
U2 = Decimal(2) ** -106


def dec_L(x):
    return -(Decimal(1) - Decimal(1) / Decimal(x)).ln()


def dec_d(x, b):
    # series in exact rationals: Decimal's 1 - t loses t entirely once t < 1e-90
    t = F(1, x ** (2 * b + 1))
    s, p, k = F(0), t, 1
    while True:
        s += p / k
        if p < F(1, 10 ** 100) * t:
            break
        k += 1
        p *= t
    return Decimal(s.numerator) / Decimal(s.denominator)


random.seed(7)
worst = {"L": 0, "d": 0, "chain": 0}
absw = {"L": 0, "d": 0}
TRUNC = Decimal("2e-48")      # series truncation, an absolute term in the bound
for _ in range(4000):
    x = random.choice([random.randrange(7, 1000), random.randrange(1000, 10**6),
                       random.randrange(10**6, 10**12)])
    h, l = qnode.L_of(x)
    ex = dec_L(x)
    err = abs(Decimal(h) + Decimal(l) - ex)
    absw["L"] = max(absw["L"], err)
    worst["L"] = max(worst["L"], max(err - TRUNC, 0) / (U2 * ex))
    b = random.randrange(1, 6)
    h2, l2 = qnode.d_of(x, b)
    ed = dec_d(x, b)
    if ed > Decimal("1e-280"):
        err = abs(Decimal(h2) + Decimal(l2) - ed)
        absw["d"] = max(absw["d"], err)
        worst["d"] = max(worst["d"], max(err - TRUNC, 0) / (U2 * ed))
    # a chain like the r = 0 check: L(q) + L(q') - K - D, D a sum of ~10 terms
    y = random.randrange(x + 2, x * 50 + 10)
    terms_pos = [(qnode.L_of(x), dec_L(x)), (qnode.L_of(y), dec_L(y))]
    terms_neg = [(qnode.d_of(random.randrange(3, 300), random.randrange(1, 8)), None) for _ in range(10)]
    hs, ls = 0.0, 0.0
    exact = Decimal(0)
    mag = Decimal(0)
    for (th, tl), _ in terms_pos:
        hs, ls = qnode.dd_add(hs, ls, th, tl)
        exact += Decimal(th) + Decimal(tl)
        mag += abs(Decimal(th) + Decimal(tl))
    for (th, tl), _ in terms_neg:
        hs, ls = qnode.dd_add(hs, ls, -th, -tl)
        exact -= Decimal(th) + Decimal(tl)
        mag += abs(Decimal(th) + Decimal(tl))
    worst["chain"] = max(worst["chain"], abs(Decimal(hs) + Decimal(ls) - exact) / (U2 * mag))
for k, v in worst.items():
    print(f"{k:6s} worst error beyond 2e-48 truncation = {float(v):8.3f} u^2 x magnitude   (bound used: 2048)")
for k, v in absw.items():
    print(f"{k:6s} worst absolute error = {float(v):.2e}")
