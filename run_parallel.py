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

from priority import raise_priority
from qp import Search


class Splitter(Search):
    def __init__(self, depth, exclude=()):
        super().__init__(c=1, exclude=exclude)
        self.depth = depth
        self.tasks = []

    def search(self, comps, r):
        if len(comps) == self.depth and r > 0:
            self.tasks.append((comps, r))
            return
        super().search(comps, r)


def key(task):
    return json.dumps([[list(c) for c in task[0]], task[1]])


HERE = os.path.dirname(os.path.abspath(__file__))
PAUSE_FILE = os.path.join(HERE, "results", "PAUSE")
PROGRESS_DIR = os.path.join(HERE, "results", "progress")


def _write_progress(path, rec):
    """Replace a worker's progress file atomically, so a reader never sees half."""
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(rec, f)
    os.replace(tmp, path)


def _pause_hook(paused, task, t0, c0):
    """Search.log callback, called every 200,000 states.

    Does two jobs through the filesystem, the one channel shared with runs the
    S4U scheduled task launches in session 0 - where an ordinary session cannot
    open the processes at all, not to query, suspend or kill:

    - waits while results/PAUSE exists;
    - writes results/progress/<pid>.json with the subtree and the states searched
      so far. Without this a running subtree is a black box until it finishes:
      on 4 Oct 2026 nine omega = 8 workers had run 13.5 hours each with nothing
      to say whether they were 10% done or 0.1%.

    Costs one small file write per 200,000 states, which is noise.
    """
    os.makedirs(PROGRESS_DIR, exist_ok=True)
    path = os.path.join(PROGRESS_DIR, f"{os.getpid()}.json")

    def hook(search):
        wall = time.time() - t0 - paused[0]
        _write_progress(path, {
            "task": key(task), "pid": os.getpid(), "nodes": search.nodes,
            "seconds": round(wall, 1),
            "cpu_seconds": round(time.process_time() - c0, 1),
            "paused_seconds": round(paused[0], 1),
            "paused": os.path.exists(PAUSE_FILE), "updated": time.time()})
        if os.path.exists(PAUSE_FILE):
            p0 = time.time()
            while os.path.exists(PAUSE_FILE):
                time.sleep(5)
            paused[0] += time.time() - p0
    hook.path = path
    return hook


def run_task(task):
    # Every worker opts out of Windows' background power throttling. Without
    # this the pool runs at about a quarter speed; see priority.py.
    prio = raise_priority()
    # (comps, r) or (comps, r, exclude): the barred primes must reach the worker,
    # or a parallel conditional search would silently become unconditional.
    comps, r = task[0], task[1]
    excl = task[2] if len(task) > 2 else ()
    paused = [0.0]
    t0, c0 = time.time(), time.process_time()
    hook = _pause_hook(paused, task, t0, c0)
    s = Search(c=1, exclude=excl, log=hook)
    s.search(tuple(tuple(c) for c in comps), r)
    wall = time.time() - t0 - paused[0]
    cpu = time.process_time() - c0
    try:                                   # finished: its result goes in the jsonl
        os.remove(hook.path)
    except OSError:
        pass
    # cpu/wall is the throttling gauge: ~1.0 means this worker had a whole core,
    # ~0.3 means Windows was holding it back. Recorded per subtree, so it is
    # visible even for workers this session cannot open.
    return {"task": key(task), "nodes": s.nodes, "leaf_checks": s.leaf_checks,
            "solutions": [str(n) for n in s.solutions],
            "unresolved": [repr(u) for u in s.unresolved],
            "seconds": round(wall, 2), "cpu_seconds": round(cpu, 2),
            "paused_seconds": round(paused[0], 1),
            "priority_ok": prio["priority"], "ecoqos_off": prio["ecoqos_off"]}


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
