"""Audit a completed run: does the recorded evidence actually cover the whole tree?

The claim "no solution with this many prime factors" rests on the union of the
recorded subtrees being the entire search tree. That is checkable without
re-running anything: re-derive the task list from the splitter, then require
every task to have a record, no record to be for an unknown task, no task to
be left unresolved, and no solution to have been found. Duplicate records for
the same task are reported and reconciled rather than silently summed.

Usage: python audit_results.py K DEPTH
"""
import json
import os
import sys
from collections import defaultdict

from run_parallel import Splitter, key


def main():
    k, depth = int(sys.argv[1]), int(sys.argv[2])
    sp = Splitter(depth)
    sp.run(k)
    expected = {key(t) for t in sp.tasks}

    path = f"results/omega{k}_d{depth}.jsonl"
    recs = defaultdict(list)
    for line in open(path):
        r = json.loads(line)
        recs[r["task"]].append(r)

    missing = expected - recs.keys()
    unknown = recs.keys() - expected
    dupes = {t: v for t, v in recs.items() if len(v) > 1}

    print(f"tasks expected from the splitter : {len(expected)}")
    print(f"tasks with at least one record   : {len(recs)}")
    print(f"records in file                  : {sum(len(v) for v in recs.values())}")
    print(f"missing tasks                    : {len(missing)}")
    print(f"records for unexpected tasks     : {len(unknown)}")
    print(f"tasks with duplicate records     : {len(dupes)}")
    for t, v in dupes.items():
        print(f"  duplicate: {t}")
        for r in v:
            tag = f"split into {r['split']['subtasks']}" if "split" in r else "single process"
            print(f"    nodes={r['nodes']:>13,}  secs={r['seconds']:>10,.0f}  {tag}")

    # One record per task for the totals: prefer the first, and check the
    # duplicates agree on the only thing that matters (nothing found, nothing left).
    chosen = {t: v[0] for t, v in recs.items()}
    nodes = sp.nodes + sum(r["nodes"] for r in chosen.values())
    leaves = sp.leaf_checks + sum(r["leaf_checks"] for r in chosen.values())
    sols = [s for r in recs.values() for x in r for s in x["solutions"]] + \
           [str(n) for n in sp.solutions]
    unres = [u for r in recs.values() for x in r for u in x["unresolved"]] + \
            [repr(u) for u in sp.unresolved]

    print(f"\nhead tree (above the split)      : {sp.nodes:,} nodes")
    print(f"total nodes (one record per task): {nodes:,}")
    print(f"exact leaf checks                : {leaves:,}")
    print(f"solutions, across ALL records    : {sols}")
    print(f"unresolved, across ALL records   : {len(unres)}")

    ok = not missing and not unknown and not sols and not unres
    print(f"\n{'COMPLETE: the recorded subtrees cover the whole tree, nothing found, nothing unresolved.' if ok else 'INCOMPLETE - see above.'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
