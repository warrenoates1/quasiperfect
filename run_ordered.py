"""Search subtrees in predicted-cost order, cheapest first.

Measured on the omega = 7 sample, the width of the admissible interval for the
next prime predicts a subtree's cost well (log-log rho = 0.85, states growing as
width^1.37). The width is a product of fractions, so it costs microseconds to
compute and needs no searching. That makes it possible to do the cheap part of a
search that cannot be finished whole: sort every subtree by predicted cost and
work upwards, so coverage accumulates fast and the intractable states are left
clearly identified at the end rather than blocking progress.

Each finished subtree is appended to results/ordered_omega<K>.jsonl with both its
predicted and actual cost, so the predictor can be audited afterwards against
what really happened. Re-running resumes from that file, and the run can be
stopped at any time without losing completed work.

Usage: python run_ordered.py K DEPTH [MAX_SCORE] [CAP] [WORKERS]
  MAX_SCORE  only run subtrees predicted below this (default 1e4, about 42% of
             the omega = 7 tree and a few hours on 14 cores)
  CAP        abandon a subtree above this many states, recording it as censored
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from multiprocessing import Pool

from qp import Search, sig
from run_parallel import Splitter, key


def score(comps, r):
    """Predicted cost: the width of the admissible interval for the next prime."""
    lo = hi = F(1)
    nmin = 1
    for p, a, fixed in comps:
        h = F(sig(p, 2 * a), p ** (2 * a))
        lo *= h
        nmin *= p ** (2 * a)
        hi *= h if fixed else F(p, p - 1)
    gl, gh = 2 - lo, 2 - hi
    if gl <= 0:
        return 0.0                      # already dead by Lemma range(4)
    if gh <= 0:
        # hi > 2: UNSCORED, not scored as expensive. Lemma range(3) gives no upper
        # bound here, so the interval has no width to measure. Such a state must
        # have an open component (were all fixed, hi = lo > 2 and it would be dead
        # above), and raising every open exponent together drives lo past 2, which
        # kills that path within a step or two. But the search also fixes
        # components, and a fixed child has a smaller hi, so it can become an
        # ordinary state of any cost. Measured on 400 sampled subtrees: strongly
        # bimodal, most finishing in tens of states while a minority run for a long
        # time. Returning infinity keeps them out of any band selected as cheap,
        # which is the conservative choice; it does not assert they are dear.
        return float("inf")
    q_lo = float(lo / (gl + F(1, nmin)))
    return max(0.0, float(r * hi / gh) - q_lo)


def run_one(args):
    comps, r, cap, sc = args

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
    return {"task": key((tuple(tuple(x) for x in comps), r)), "predicted": sc,
            "nodes": s.nodes, "secs": round(time.time() - t, 2), "censored": hit,
            "unresolved": [repr(u) for u in s.unresolved],
            "solutions": [str(x) for x in s.solutions]}


def main():
    k, depth = int(sys.argv[1]), int(sys.argv[2])
    cutoff = float(sys.argv[3]) if len(sys.argv) > 3 else 1e4
    cap = int(sys.argv[4]) if len(sys.argv) > 4 else 50_000_000
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else 14

    sp = Splitter(depth)
    sp.run(k)
    total = len(sp.tasks)
    scored = [(score(c, r), c, r) for c, r in sp.tasks]
    scored.sort(key=lambda t: t[0])
    chosen = [t for t in scored if t[0] <= cutoff]
    print(f"omega={k} depth={depth}: {total:,} subtrees; {len(chosen):,} "
          f"({len(chosen)/total:.1%}) predicted at or below {cutoff:,.3g}", flush=True)
    print(f"head tree above the split: {sp.nodes:,} states, "
          f"{len(sp.unresolved)} unresolved, solutions {sp.solutions}", flush=True)

    os.makedirs("results", exist_ok=True)
    path = f"results/ordered_omega{k}.jsonl"
    done = set()
    bad = 0
    if os.path.exists(path):
        for line in open(path):
            # A machine losing power mid-write leaves a truncated final line.
            # Skipping it costs one subtree, which is redone; letting json.loads
            # raise here would make the whole run unresumable, which is worse.
            try:
                done.add(json.loads(line)["task"])
            except (ValueError, KeyError):
                bad += 1
    if bad:
        print(f"skipped {bad} unparseable record(s) - truncated by an "
              f"interrupted write; those subtrees will be redone", flush=True)
    todo = [(c, r, cap, sc) for sc, c, r in chosen if key((tuple(tuple(x) for x in c), r)) not in done]
    print(f"{len(done):,} already recorded; {len(todo):,} to run on {workers} workers",
          flush=True)

    t0 = time.time()
    nodes = cz = 0
    flagged = []
    last_sync = time.time()
    with open(path, "a") as out, Pool(workers) as pool:
        for i, rec in enumerate(pool.imap(run_one, todo, chunksize=8), 1):
            out.write(json.dumps(rec) + "\n")
            # Flush every record so that killing the process -- a shutdown, a
            # forced restart for updates -- loses nothing already computed, and
            # fsync periodically so that losing power loses only seconds of it.
            out.flush()
            if time.time() - last_sync > 5:
                os.fsync(out.fileno())
                last_sync = time.time()
            nodes += rec["nodes"]
            cz += rec["censored"]
            if rec["solutions"] or rec["unresolved"]:
                flagged.append(rec)
                print(f"  *** FLAG {rec}", flush=True)
            if i % 500 == 0 or i == len(todo):
                os.fsync(out.fileno())
                el = time.time() - t0
                print(f"  {i:,}/{len(todo):,} subtrees  {nodes:,} states  "
                      f"{cz} censored  {el/60:,.0f} min  "
                      f"(eta {el/i*(len(todo)-i)/60:,.0f} min)", flush=True)

    print(f"\nDONE: {len(todo):,} subtrees, {nodes:,} states, {cz} censored, "
          f"{(time.time()-t0)/3600:.2f} h")
    print(f"coverage of omega={k} at depth {depth}: "
          f"{(len(done)+len(todo))/total:.1%} of subtrees")
    print(f"solutions {[s for r in flagged for s in r['solutions']]}, "
          f"unresolved {sum(len(r['unresolved']) for r in flagged)}")


if __name__ == "__main__":
    main()
