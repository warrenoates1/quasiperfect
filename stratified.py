"""Stratified cost estimate for a full omega = K search.

Uniform sampling fails here: the subtree sizes are so heavy-tailed that a few
hundred random subtrees miss the ones that dominate, and capping them hides
exactly the quantity we need. At omega = 6 one subtree held 70% of the work; in
the uniform omega = 7 sample, 122 of 400 subtrees hit a 2e6 cap and supplied 83%
of the sampled states while truncated, so the scaled figure was a lower bound of
unknown quality.

This does it in two passes instead:

  pass 1 (probe)  a random sample, each capped at PROBE nodes, classifies
                  subtrees as light (finished) or heavy (hit the cap), and
                  measures the light mean exactly.
  pass 2 (heavy)  a random subsample of the heavy ones is run with a much
                  larger cap, giving the heavy mean. Any that still hit the
                  large cap are reported as censored, and the heavy mean is then
                  itself a lower bound -- but a far tighter one.

Estimate = total_subtrees * (frac_light * mean_light + frac_heavy * mean_heavy).

Every subtree run to completion is also a genuine partial verification, so
solutions and unresolved states are accumulated and reported.

Usage: python stratified.py K DEPTH NPROBE NHEAVY [PROBE] [BIG] [WORKERS]
"""
import json
import os
import random
import sys
import time
from multiprocessing import Pool

from qp import Search
from run_parallel import Splitter, key


def run_one(args):
    comps, r, cap = args

    class Capped(Search):
        def search(self, c, rr):
            if self.nodes > cap:
                raise StopIteration
            return super().search(c, rr)

    s = Capped(c=1)
    t = time.time()
    hit = False
    try:
        s.search(tuple(tuple(x) for x in comps), r)
    except (StopIteration, RecursionError):
        hit = True
    return {"task": key((tuple(tuple(x) for x in comps), r)),
            "nodes": s.nodes, "secs": round(time.time() - t, 2), "censored": hit,
            "unresolved": len(s.unresolved),
            "solutions": [str(x) for x in s.solutions]}


def main():
    k, depth = int(sys.argv[1]), int(sys.argv[2])
    nprobe, nheavy = int(sys.argv[3]), int(sys.argv[4])
    probe = int(sys.argv[5]) if len(sys.argv) > 5 else 2_000_000
    big = int(sys.argv[6]) if len(sys.argv) > 6 else 100_000_000
    workers = int(sys.argv[7]) if len(sys.argv) > 7 else 14

    sp = Splitter(depth)
    sp.run(k)
    total = len(sp.tasks)
    rng = random.Random(777)
    sample = rng.sample(sp.tasks, min(nprobe, total))
    print(f"omega={k} depth={depth}: {total:,} subtrees", flush=True)
    print(f"pass 1: probing {len(sample):,} at a cap of {probe:,} nodes", flush=True)

    os.makedirs("results", exist_ok=True)
    out = open(f"results/strat_omega{k}.jsonl", "a")
    t0 = time.time()
    with Pool(workers) as pool:
        recs = []
        for i, rec in enumerate(pool.imap_unordered(
                run_one, [(c, r, probe) for c, r in sample], chunksize=1), 1):
            rec["pass"] = 1
            out.write(json.dumps(rec) + "\n")
            recs.append(rec)
            if i % 40 == 0:
                print(f"  probed {i}/{len(sample)}, {time.time()-t0:,.0f}s", flush=True)
    out.flush()

    light = [r for r in recs if not r["censored"]]
    heavy = [r for r in recs if r["censored"]]
    fl, fh = len(light) / len(recs), len(heavy) / len(recs)
    mean_l = sum(r["nodes"] for r in light) / max(1, len(light))
    rate = sum(r["nodes"] for r in recs) / max(1e-9, sum(r["secs"] for r in recs))
    print(f"\npass 1 done in {time.time()-t0:,.0f}s: {len(light)} light "
          f"({fl:.1%}), {len(heavy)} heavy ({fh:.1%})")
    print(f"  mean light subtree: {mean_l:,.0f} nodes")
    print(f"  measured rate: {rate:,.0f} nodes per core-second")

    if not heavy:
        print("no heavy subtrees: the uniform estimate is already reliable.")
        return
    tasks = [(json.loads(r["task"])[0], json.loads(r["task"])[1]) for r in heavy]
    pick = rng.sample(tasks, min(nheavy, len(tasks)))
    print(f"\npass 2: running {len(pick)} heavy subtrees at a cap of {big:,} nodes",
          flush=True)
    t1 = time.time()
    hrecs = []
    with Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(
                run_one, [(c, r, big) for c, r in pick], chunksize=1), 1):
            rec["pass"] = 2
            out.write(json.dumps(rec) + "\n")
            out.flush()
            hrecs.append(rec)
            tag = "CENSORED" if rec["censored"] else "complete"
            print(f"  [{i}/{len(pick)}] {rec['nodes']:>13,} nodes {rec['secs']:>9,.0f}s "
                  f"{tag}", flush=True)
    out.close()

    cens = [r for r in hrecs if r["censored"]]
    mean_h = sum(r["nodes"] for r in hrecs) / len(hrecs)
    est = total * (fl * mean_l + fh * mean_h)
    cores = est / rate / 3600
    print(f"\npass 2 done in {time.time()-t1:,.0f}s: {len(cens)} of {len(hrecs)} "
          f"still censored at {big:,}")
    print(f"  mean heavy subtree: {mean_h:,.0f} nodes"
          f"{'  (LOWER BOUND, some censored)' if cens else '  (exact)'}")
    print(f"\nstratified estimate for omega={k}:")
    print(f"  nodes       {'>=' if cens else '~'} {est:,.3g}")
    print(f"  core-hours  {'>=' if cens else '~'} {cores:,.3g}")
    print(f"  on {workers} cores {'>=' if cens else '~'} {cores/workers/24:,.3g} days")
    allr = recs + hrecs
    print(f"\npartial verification from {len(allr)} subtrees actually searched: "
          f"{sum(r['nodes'] for r in allr):,} states, "
          f"solutions {[s for r in allr for s in r['solutions']]}, "
          f"unresolved {sum(r['unresolved'] for r in allr)}")


if __name__ == "__main__":
    main()
