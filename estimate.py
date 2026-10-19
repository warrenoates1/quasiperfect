"""Estimate the cost of a full omega = K search from a random sample of subtrees.

Runs a random sample of the depth-DEPTH subtrees to completion, under a node
cap, and scales up. The distribution of subtree sizes is very heavy-tailed -- at
omega = 6 a single subtree held 70% of the whole computation -- so the scaled
figure is a LOWER BOUND: a random sample of a few hundred will usually miss the
few subtrees that dominate. Capped subtrees are reported separately, since each
contributes an unknown amount above the cap.

Usage: python estimate.py K DEPTH SAMPLE [CAP] [WORKERS]
"""
import json
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
    return {"nodes": s.nodes, "secs": time.time() - t, "capped": hit,
            "unresolved": len(s.unresolved), "solutions": [str(x) for x in s.solutions]}


def main():
    k, depth, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    cap = int(sys.argv[4]) if len(sys.argv) > 4 else 2_000_000
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else 14

    sp = Splitter(depth)
    sp.run(k)
    total = len(sp.tasks)
    rng = random.Random(12345)
    sample = rng.sample(sp.tasks, min(n, total))
    print(f"omega={k} depth={depth}: {total:,} subtrees; sampling {len(sample):,}, "
          f"cap {cap:,} nodes/subtree", flush=True)

    t0 = time.time()
    with Pool(workers) as pool:
        recs = pool.map(run_one, [(c, r, cap) for c, r in sample], chunksize=1)

    nodes = sum(x["nodes"] for x in recs)
    cpu = sum(x["secs"] for x in recs)
    capped = [x for x in recs if x["capped"]]
    sols = [s for x in recs for s in x["solutions"]]
    unres = sum(x["unresolved"] for x in recs)
    scale = total / len(sample)

    print(f"sample: {nodes:,} nodes, {cpu:,.0f} core-seconds, "
          f"{len(capped)} of {len(recs)} hit the cap")
    print(f"sample rate: {nodes / cpu:,.0f} nodes per core-second")
    print(f"solutions in sample: {sols}   unresolved in sample: {unres}")
    print(f"\nscaling by {scale:,.1f}x -- LOWER BOUND, the tail is not sampled:")
    print(f"  nodes        >= {nodes * scale:,.3g}")
    print(f"  core-hours   >= {cpu * scale / 3600:,.3g}")
    print(f"  on 14 cores  >= {cpu * scale / 3600 / 14 / 24:,.3g} days")
    print(f"\nwall clock for this estimate: {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    main()
