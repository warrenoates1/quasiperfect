"""Split the omega = k search into independent subtrees and run them in parallel.

The tree is walked normally until a node has `depth` primes; that node becomes
a task instead of being expanded. Every task is an ordinary Search.search call,
so the union of tasks is exactly the full tree. Finished tasks are appended to
results/omega<k>_d<depth>.jsonl, so an interrupted run resumes where it left off.

Usage: python run_parallel.py K [DEPTH] [WORKERS]
"""
import json
import os
import sys
import time
from multiprocessing import Pool

from qp import Search


class Splitter(Search):
    def __init__(self, depth):
        super().__init__(c=1)
        self.depth = depth
        self.tasks = []

    def search(self, comps, r):
        if len(comps) == self.depth and r > 0:
            self.tasks.append((comps, r))
            return
        super().search(comps, r)


def key(task):
    return json.dumps([[list(c) for c in task[0]], task[1]])


def run_task(task):
    comps, r = task
    t0 = time.time()
    s = Search(c=1)
    s.search(tuple(tuple(c) for c in comps), r)
    return {"task": key(task), "nodes": s.nodes, "leaf_checks": s.leaf_checks,
            "solutions": [str(n) for n in s.solutions],
            "unresolved": [repr(u) for u in s.unresolved],
            "seconds": round(time.time() - t0, 2)}


def main():
    k = int(sys.argv[1])
    depth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else max(1, os.cpu_count() - 2)

    sp = Splitter(depth)
    sp.run(k)
    head = {"nodes": sp.nodes, "solutions": sp.solutions, "unresolved": sp.unresolved}
    tasks = sp.tasks

    os.makedirs("results", exist_ok=True)
    path = f"results/omega{k}_d{depth}.jsonl"
    done = {}
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                rec = json.loads(line)
                done[rec["task"]] = rec
    todo = [t for t in tasks if key(t) not in done]
    print(f"omega={k} depth={depth}: {len(tasks)} tasks, {len(done)} already done, "
          f"{len(todo)} to run on {workers} workers; head tree {head}", flush=True)

    t0 = time.time()
    with open(path, "a") as out, Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_task, todo, chunksize=1), 1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            done[rec["task"]] = rec
            if rec["solutions"] or rec["unresolved"]:
                print("FLAG", rec, flush=True)
            if i % 50 == 0 or i == len(todo):
                print(f"  {len(done)}/{len(tasks)} tasks, {time.time() - t0:.0f}s", flush=True)

    recs = [done[key(t)] for t in tasks]
    nodes = head["nodes"] + sum(r["nodes"] for r in recs)
    sols = head["solutions"] + [s for r in recs for s in r["solutions"]]
    unres = len(head["unresolved"]) + sum(len(r["unresolved"]) for r in recs)
    print(f"DONE omega={k}: nodes={nodes:,} solutions={sols} unresolved={unres}", flush=True)


if __name__ == "__main__":
    main()
