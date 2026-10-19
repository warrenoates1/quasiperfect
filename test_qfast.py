"""Planted-solution test for qfast.run_q, then timing on the real node."""
import ast, json, math, random, time
import numpy as np
from qp import Search
import qfast

rec = [json.loads(l) for l in open('results/cond_3-5_17.jsonl')]
_, comps, r, qstop = ast.literal_eval(rec[-1]['unresolved'][0])
s = Search(c=1, exclude={3, 5})
opens, K, KL, dtab, dtl, allowed, floor, amax = qfast.prepare(comps, s)
from fractions import Fraction as F
print("opens", opens, "K", K, "amax", list(amax))
nop = len(opens)
out = np.zeros((200000, 2 + 2 * nop), np.int64)
everything = np.ones_like(allowed)
random.seed(1)
ok = 0
T = 300
for t in range(T):
    q = qfast.primes_between(random.randrange(571_300_000, 1_130_000_000), 0) if False else None
    # pick q, q' primes and a fully fixed exponent vector, then choose K to fit exactly
    while True:
        q = random.randrange(571_300_000, 1_130_000_000) | 1
        if qfast.isprime(q): break
    while True:
        qq = random.randrange(q + 2, q * (3 if t % 3 else 50)) | 1
        if qfast.isprime(qq): break
    a = [f + random.choice([0, 0, 1, 2, 5]) for f in floor]
    D = sum(F(dtab[i, a[i]]) + F(dtl[i, a[i]]) for i in range(nop))
    b1 = qfast.log_f_bounds(q); b2 = qfast.log_f_bounds(qq)
    Ke = F(b1[0]) + F(b1[1]) + F(b2[0]) + F(b2[1]) - D
    Kp = float(Ke); KpL = float(Ke - F(Kp))
    stats = np.zeros(4, np.int64)
    n = qfast.run_q(q, Kp, KpL, dtab, dtl, everything, floor, amax, out, 0, stats)
    hit = any(out[j, 1] == qq and all((out[j, 2 + i] == a[i]) if out[j, 2 + nop + i] else (out[j, 2 + i] <= a[i]) for i in range(nop)) for j in range(min(n, len(out))))
    ok += hit
    if not hit:
        print("MISS", q, qq, a, n, list(stats))
print(f"planted: {ok}/{T} found")

qs = qfast.primes_between(571_287_684, 571_287_684 + 20000)
stats = np.zeros(4, np.int64)
t0 = time.time()
n = qfast.run_block(qs[:5], K, KL, dtab, dtl, allowed, floor, amax, out, stats)   # compile
t0 = time.time(); stats[:] = 0
n = qfast.run_block(qs, K, KL, dtab, dtl, allowed, floor, amax, out, stats)
dt = time.time() - t0
nf = sum(all(out[j, 2 + nop:]) for j in range(min(n, len(out))))
print("all-fixed candidates", nf)
print(f"{len(qs)} q in {dt:.2f}s ({dt/len(qs)*1e3:.2f} ms/q), states {stats[0]}, splits {stats[1]}, unresolved {stats[2]}, overflow {stats[3]}, candidates {n}")
