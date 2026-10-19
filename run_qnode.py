"""Full omega = K search with every r = 2 node done by the compiled kernel (qnode).

usage: python run_qnode.py K DEPTH WORKERS
Splits the tree at DEPTH components (run_parallel.Splitter), runs each subtree
with NodeSearch, appends one record per subtree to results/qnode_omega{K}_d{DEPTH}.jsonl
(resumable), and writes a .done summary when every subtree is in.
"""
import json, os, sys, time
from multiprocessing import Pool
from priority import raise_priority
from qp import Search
from run_parallel import Splitter, key
import qnode


class NS(qnode.NodeSearch, Search):
    pass


def run_task(task):
    prio = raise_priority()
    comps, r = task
    t0, c0 = time.time(), time.process_time()
    s = NS(c=1)
    s.search(tuple(tuple(c) for c in comps), r)
    return {"task": key(task), "nodes": s.nodes, "kernel_states": getattr(s, "kstates", 0),
            "kernel_cands": getattr(s, "kcands", 0), "leaf_checks": s.leaf_checks,
            "solutions": [str(n) for n in s.solutions],
            "unresolved": [repr(u) for u in s.unresolved],
            "seconds": round(time.time() - t0, 2), "cpu_seconds": round(time.process_time() - c0, 2),
            "priority_ok": prio["priority"]}


def main():
    k, depth, workers = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    sp = Splitter(depth)
    sp.run(k)
    head = {"nodes": sp.nodes, "solutions": sp.solutions, "unresolved": sp.unresolved}
    tasks = sp.tasks
    path = f"results/qnode_omega{k}_d{depth}.jsonl"
    done = {}
    if os.path.exists(path):
        for line in open(path):
            if line.strip():
                rec = json.loads(line)
                done[rec["task"]] = rec
    todo = [t for t in tasks if key(t) not in done]
    print(f"omega={k} depth={depth}: {len(tasks)} subtrees, {len(done)} done, {len(todo)} to run "
          f"on {workers} workers; head {head}", flush=True)
    t0 = time.time()
    with open(path, "a") as out, Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_task, todo, chunksize=4), 1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            done[rec["task"]] = rec
            if rec["solutions"] or rec["unresolved"]:
                print("FLAG", rec, flush=True)
            if i % 2000 == 0 or i == len(todo):
                print(f"  {len(done)}/{len(tasks)} subtrees, {time.time() - t0:.0f}s", flush=True)
    recs = [done[key(t)] for t in tasks]
    summary = {"omega": k, "depth": depth, "subtrees": len(tasks),
               "python_states": head["nodes"] + sum(r["nodes"] for r in recs),
               "kernel_states": sum(r["kernel_states"] for r in recs),
               "kernel_candidates": sum(r["kernel_cands"] for r in recs),
               "solutions": head["solutions"] + [x for r in recs for x in r["solutions"]],
               "unresolved": len(head["unresolved"]) + sum(len(r["unresolved"]) for r in recs),
               "core_seconds": round(sum(r["cpu_seconds"] for r in recs))}
    summary["verdict"] = "ELIMINATED" if not summary["solutions"] and not summary["unresolved"] else "NOT ELIMINATED"
    open(path.replace(".jsonl", ".done"), "w").write(json.dumps(summary) + "\n")
    print("DONE", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
