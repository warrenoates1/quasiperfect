"""Is subtree cost predictable from the subtree itself, before searching it?

Every sampled subtree's identity and node count is recorded, so this is an
empirical question rather than a guess. Candidate predictors, all computable from
the subtree in microseconds:

  gap_hi = 2 - hi    how much room the abundancy ceiling leaves. Lemma range(3)
                     bounds the next prime by r*hi/gap_hi, so a small gap means a
                     huge window of primes to rule out.
  gap_lo = 2 - lo    the same for the floor. When lo >= 2 the subtree is dead at
                     once by Lemma range(4), so these are the cheapest of all.

Reports which predictor correlates best with measured cost, the structural
signature of the dearest subtrees, and what it would cost to clear the cheap
ones.

Usage: python analyse_cost.py [results/strat_omega7.jsonl]
"""
import json
import math
import sys
from fractions import Fraction as F

from qp import sig

TOTAL = 80458          # depth-5 subtrees at omega = 7


def features(comps, r):
    lo = hi = F(1)
    nmin = 1
    for p, a, fixed in comps:
        h = F(sig(p, 2 * a), p ** (2 * a))
        lo *= h
        nmin *= p ** (2 * a)
        hi *= h if fixed else F(p, p - 1)
    gh, gl = 2 - hi, 2 - lo
    return {
        "lo": float(lo), "hi": float(hi), "r": r,
        "gap_hi": float(gh), "gap_lo": float(gl),
        "dead": gl <= 0,
        # window, finite and meaningful only when hi < 2
        "window": float(r * hi / gh) if gh > 0 else float("inf"),
        # The admissible interval for the next prime, from Lemma range(2) below
        # and range(3) above. When the lower end exceeds the upper end nothing
        # fits and the subtree dies at once, however close hi is to 2 -- which is
        # why "hi near 2" alone mispredicts. Cost should track the interval's
        # width, not either end.
        "q_lo": float(lo / (gl + F(1, nmin))) if gl > 0 else float("inf"),
        "q_hi": float(r * hi / gh) if gh > 0 else float("inf"),
        "score": max(0.0, (float(r * hi / gh) if gh > 0 else float("inf"))
                     - (float(lo / (gl + F(1, nmin))) if gl > 0 else float("inf"))),
        "primes": [p for p, _, _ in comps],
        "exps": [2 * a for _, a, _ in comps],
    }


def corr(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return (sxy / math.sqrt(sxx * syy), sxy / sxx) if sxx and syy else (float("nan"),) * 2


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "results/strat_omega7.jsonl"
    best = {}
    for line in open(path):
        d = json.loads(line)
        comps, r = json.loads(d["task"])
        f = features([tuple(c) for c in comps], r)
        f.update(nodes=d["nodes"], censored=d["censored"], secs=d["secs"])
        k = (tuple(f["primes"]), tuple(f["exps"]))
        if k not in best or f["nodes"] > best[k]["nodes"]:
            best[k] = f
    rows = list(best.values())
    done = [f for f in rows if not f["censored"]]
    dead = [f for f in rows if f["dead"]]
    print(f"{len(rows)} distinct subtrees; {len(done)} completed; "
          f"{len(dead)} already dead by Lemma range(4)\n")

    print("predictor vs measured states (completed subtrees, log-log):")
    for name, key in (("window  = r*hi/(2-hi)", "window"),
                      ("interval width q_hi-q_lo", "score")):
        pts = [(f[key], f["nodes"]) for f in done
               if 0 < f[key] < float("inf") and f["nodes"] > 0]
        rho, slope = corr([math.log10(x) for x, _ in pts],
                          [math.log10(y) for _, y in pts])
        print(f"  {name:<26} n={len(pts):<4} rho = {rho:+.3f}   states ~ x^{slope:.2f}")

    empty = [f for f in rows if f["score"] == 0.0]
    if empty:
        ns = sorted(f["nodes"] for f in empty)
        print(f"\n  subtrees with an empty interval (no prime fits): {len(empty)}, "
              f"median {ns[len(ns)//2]} states, max {ns[-1]:,}")

    print("\nbinned by predicted score:")
    print(f"  {'score':>20}  {'subtrees':>8} {'censored':>8} {'median':>12} {'max':>13}"
          f" {'est. total states':>18}")
    edges = [0, 1, 1e3, 1e4, 1e5, 1e6, 1e7, float("inf")]
    plan = []
    for a, b in zip(edges, edges[1:]):
        grp = [f for f in rows if a <= f["score"] < b]
        if not grp:
            continue
        ns = sorted(f["nodes"] for f in grp)
        cz = sum(f["censored"] for f in grp)
        frac = len(grp) / len(rows)
        mean = sum(ns) / len(ns)
        est = TOTAL * frac * mean
        plan.append((a, b, len(grp), cz, mean, est))
        print(f"  {a:>9.3g} - {b:<8.3g}  {len(grp):>8} {cz:>8} {ns[len(ns)//2]:>12,} "
              f"{ns[-1]:>13,} {est:>18,.3g}")

    print("\nclearing the cheap strata first (cumulative, scaled to "
          f"{TOTAL:,} subtrees):")
    cs = cn = 0
    for a, b, k, cz, mean, est in plan:
        cs += est
        cn += TOTAL * k / len(rows)
        hrs = cs / 10021 / 3600
        flag = "" if cz == 0 else f"  <-- {cz}/{k} censored, cost is a lower bound"
        print(f"  score < {b:<9.3g} {cn:>9,.0f} subtrees ({cn/TOTAL:5.1%})  "
              f"{cs:>10.3g} states  {hrs:>8,.0f} core-hours  "
              f"{hrs/14/24:>6.1f} days on 14 cores{flag}")

    print("\nstructural signature of the ten dearest:")
    for f in sorted(rows, key=lambda f: -f["nodes"])[:10]:
        prs = " ".join(f"{p}^{e}" for p, e in zip(f["primes"], f["exps"]))
        print(f"  {f['nodes']:>13,}{'*' if f['censored'] else ' '} "
              f"2-lo {f['gap_lo']:.2e}  score {f['score']:>12,.0f}  {prs}")
    print("  (* = censored)")


if __name__ == "__main__":
    main()
