"""Characterise the subtrees the cost predictor cannot score: those with hi > 2.

The predictor of Section "Cost is predictable before searching" uses the width of
the admissible interval for the next prime. Its upper end comes from Lemma
range(3), which only bites when hi < 2; where hi > 2 there is no upper bound and
the width is formally infinite. That is about 18% of the depth-5 subtrees at
omega = 7, and their real cost is unknown: the few seen in early sampling were
trivial (3-6 states), which would mean the predictor has them exactly backwards,
but a serial probe of 300 was killed before finishing so nothing was established.

This samples them properly, in parallel and under a cap, and records every
result. Two outcomes matter: the censoring rate (are they cheap?) and whether any
solution or unresolved state appears (they are real coverage either way).

Usage: python probe_hi2.py K DEPTH SAMPLE [CAP] [WORKERS]
"""
import json
import os
import random
import sys
import time
from multiprocessing import Pool

from qp import Search
from run_ordered import run_one, score
from run_parallel import Splitter, key


def main():
    k, depth, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    cap = int(sys.argv[4]) if len(sys.argv) > 4 else 2_000_000
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else 14

    sp = Splitter(depth)
    sp.run(k)
    pool_tasks = [(c, r) for c, r in sp.tasks if score(c, r) == float("inf")]
    total = len(sp.tasks)
    print(f"omega={k} depth={depth}: {total:,} subtrees, "
          f"{len(pool_tasks):,} with hi > 2 ({len(pool_tasks)/total:.1%})", flush=True)
    rng = random.Random(4242)
    sample = rng.sample(pool_tasks, min(n, len(pool_tasks)))
    print(f"sampling {len(sample):,} at a cap of {cap:,} states on {workers} workers",
          flush=True)

    os.makedirs("results", exist_ok=True)
    recs = []
    t0 = time.time()
    with open(f"results/hi2_omega{k}.jsonl", "a") as out, Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(
                run_one, [(c, r, cap, float("inf")) for c, r in sample],
                chunksize=4), 1):
            out.write(json.dumps(rec) + "\n")
            recs.append(rec)
            if rec["solutions"] or rec["unresolved"]:
                print(f"  *** FLAG {rec}", flush=True)
            if i % 50 == 0 or i == len(sample):
                out.flush()
                ns = sorted(x["nodes"] for x in recs)
                print(f"  {i}/{len(sample)}  median {ns[len(ns)//2]:,}  "
                      f"max {ns[-1]:,}  censored {sum(x['censored'] for x in recs)}  "
                      f"{time.time()-t0:,.0f}s", flush=True)

    ns = sorted(x["nodes"] for x in recs)
    cz = sum(x["censored"] for x in recs)
    mean = sum(ns) / len(ns)
    print(f"\n{len(recs)} sampled, {cz} censored at {cap:,}")
    print(f"  states: min {ns[0]:,}  median {ns[len(ns)//2]:,}  mean {mean:,.0f}  "
          f"max {ns[-1]:,}")
    print(f"  deciles: {[f'{ns[int(len(ns)*f)]:,}' for f in (0.1,0.25,0.5,0.75,0.9)]}")
    print(f"  extrapolated cost of all {len(pool_tasks):,}: "
          f"{mean*len(pool_tasks):,.3g} states"
          f"{'  (LOWER BOUND, some censored)' if cz else ''}")
    print(f"  solutions {[s for x in recs for s in x['solutions']]}, "
          f"unresolved {sum(len(x['unresolved']) for x in recs)}")
    cheap = sum(1 for v in ns if v < 1000)
    print(f"\n  verdict: {cheap}/{len(ns)} ({cheap/len(ns):.0%}) finish under 1,000 "
          f"states; {cz/len(ns):.0%} hit the cap")


if __name__ == "__main__":
    main()
