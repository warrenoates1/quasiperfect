"""NodeSearch (compiled r = 2 kernel) against Search, on real solutions of
sigma(n) = 2n + c and on the omega <= 6 reproduction."""
import sys, time
from qp import Search
import qnode

class NS(qnode.NodeSearch, Search):
    pass

mode = sys.argv[1] if len(sys.argv) > 1 else "c"
if mode == "c":
    CMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    W = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [3, 4]
    bad = sols = 0; t0 = time.time()
    for w in W:
        for c in range(1, CMAX + 1):
            a = Search(c=c); a.run(w)
            b = NS(c=c); b.run(w)
            sols += len(a.solutions)
            if sorted(a.solutions) != sorted(b.solutions) or bool(a.unresolved) != bool(b.unresolved):
                bad += 1
                print("MISMATCH", c, w, a.solutions, b.solutions, a.unresolved[:1], b.unresolved[:1], flush=True)
    print(f"{bad} mismatches, {sols} solutions compared, {time.time()-t0:.0f}s")
elif mode == "w":
    for w in range(1, int(sys.argv[2]) + 1):
        t0 = time.time()
        b = NS(c=1); b.run(w)
        print(w, b.nodes, getattr(b, 'kstates', 0), getattr(b, 'kcands', 0), b.solutions, len(b.unresolved), f"{time.time()-t0:.1f}s", flush=True)

if mode == "b":
    # every (c, omega) with a brute-force solution n = k^2, k <= K, c <= CMAX
    from collections import defaultdict
    from qp import factorize, sig
    K = int(sys.argv[2]); CMAX = int(sys.argv[3])
    brute = defaultdict(list)
    for k in range(3, K + 1, 2):
        f = factorize(k); n = k * k; s = 1
        for p, e in f.items():
            s *= sig(p, 2 * e)
        cc = s - 2 * n
        if 0 < cc <= CMAX:
            brute[(cc, len(f))].append(n)
    bad = 0; t0 = time.time(); found = 0
    print(f"{len(brute)} cases to run", flush=True)
    STEP = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    for idx, ((cc, w), ns) in enumerate(sorted(brute.items())[::STEP]):
        t1 = time.time()
        b = NS(c=cc); b.run(w)
        if time.time() - t1 > 5:
            print(f"slow: c={cc} omega={w} {time.time()-t1:.0f}s", flush=True)
        if idx % 200 == 0:
            print(f"  {idx} cases, c up to {cc}, {time.time()-t0:.0f}s", flush=True)
        small = sorted(x for x in b.solutions if x <= K * K)
        ok = small == sorted(ns) and not b.unresolved
        found += len(ns)
        if not ok:
            bad += 1
            print("BAD", cc, w, ns, b.solutions, b.unresolved[:1], flush=True)
    print(f"{len(brute)} (c, omega) cases, {found} brute-force solutions, {bad} bad, {time.time()-t0:.0f}s")
