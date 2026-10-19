"""Split one task from a run_parallel results file one level deeper and run it in parallel.

Usage: python split_task.py K DEPTH TASK_JSON [WORKERS]

The sub-tasks come from walking TASK's subtree with Splitter(DEPTH), so together
they cover exactly TASK. Finished sub-tasks go to results/split_<hash>.jsonl
(resumable). When all are done, one record for TASK, marked "split", is appended
to results/omega<K>_d<depth of TASK>.jsonl so the main run treats it as done.
"""
import hashlib
import json
import os
import sys
import time
from multiprocessing import Pool

from run_parallel import Splitter, key, run_task


def main():
    k, depth, task_json = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 13
    comps, r = json.loads(task_json)
    comps = tuple(tuple(c) for c in comps)
    task = (comps, r)

    sp = Splitter(depth)
    sp.search(comps, r)            # walks the subtree, collecting sub-tasks
    head_nodes = sp.nodes
    subs = sp.tasks

    os.makedirs("results", exist_ok=True)
    tag = hashlib.sha1(key(task).encode()).hexdigest()[:10]
    path = f"results/split_{tag}.jsonl"
    done = {}
    if os.path.exists(path):
        for line in open(path):
            rec = json.loads(line)
            done[rec["task"]] = rec
    todo = [t for t in subs if key(t) not in done]
    print(f"{key(task)} -> {len(subs)} sub-tasks at depth {depth}, "
          f"{len(done)} done, head nodes {head_nodes}, unresolved {sp.unresolved}", flush=True)

    t0 = time.time()
    with open(path, "a") as out, Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_task, todo, chunksize=1), 1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            done[rec["task"]] = rec
            if rec["solutions"] or rec["unresolved"]:
                print("FLAG", rec, flush=True)
            if i % 25 == 0 or i == len(todo):
                print(f"  {len(done)}/{len(subs)} sub-tasks, {time.time() - t0:.0f}s", flush=True)

    recs = [done[key(t)] for t in subs]
    merged = {
        "task": key(task),
        "nodes": head_nodes + sum(x["nodes"] for x in recs),
        "leaf_checks": sp.leaf_checks + sum(x["leaf_checks"] for x in recs),
        "solutions": [str(n) for n in sp.solutions] + [s for x in recs for s in x["solutions"]],
        "unresolved": [repr(u) for u in sp.unresolved] + [u for x in recs for u in x["unresolved"]],
        "seconds": round(sum(x["seconds"] for x in recs), 2),
        "split": {"depth": depth, "subtasks": len(subs), "file": path},
    }
    main_path = f"results/omega{k}_d{len(comps)}.jsonl"
    existing = {json.loads(l)["task"] for l in open(main_path)} if os.path.exists(main_path) else set()
    if merged["task"] not in existing:
        with open(main_path, "a") as f:
            f.write(json.dumps(merged) + "\n")
    print("MERGED", {k2: v for k2, v in merged.items() if k2 != "task"}, flush=True)


if __name__ == "__main__":
    main()
