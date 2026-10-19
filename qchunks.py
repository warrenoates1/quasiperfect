"""Split dear r = 2 nodes by the range of their smaller prime q and run the pieces in parallel.

usage: python qchunks.py OUT.jsonl CHUNKS WORKERS TASK_JSON [TASK_JSON ...]
Each TASK_JSON is a subtree key from the omega run, e.g.
'[[[3, 5, false], ..., [65539, 1, false]], 2]'. The subtree is expanded by
Search down to its r = 2 nodes, and each r = 2 node's q-range [last+1, 2^38)
is cut at CHUNKS geometric points; every piece runs as its own job. Coverage:
the pieces of each node tile [last+1, 2^38) exactly, and q >= 2^38 is flagged
by the kernel itself as unresolved. Resumable by (node, piece).
"""
import json, math, os, sys, time
from multiprocessing import Pool
from priority import raise_priority
from qp import Search
import qnode


class Collect(Search):
    """Expand a subtree down to its r = 2 nodes without searching them."""
    def __init__(self):
        super().__init__(c=1)
        self.nodes2 = []

    def enumerate_next(self, comps, r, lo, hi, nmin, last):
        if r == 2:
            self.nodes2.append(comps)
            return
        return super().enumerate_next(comps, r, lo, hi, nmin, last)


class NS(qnode.NodeSearch, Search):
    pass


def piece(job):
    raise_priority()
    comps, qa, qb = job
    comps = tuple(tuple(c) for c in comps)
    s = NS(c=1)
    t0, c0 = time.time(), time.process_time()
    st, n = qnode.search_node(comps, s, qfrom=qa, qto=qb)
    return {"node": json.dumps(comps), "qa": qa, "qb": qb, "kernel_states": int(st[0]),
            "cands": n, "exact_nodes": s.nodes, "solutions": [str(x) for x in s.solutions],
            "unresolved": [repr(u) for u in s.unresolved],
            "seconds": round(time.time() - t0, 1), "cpu_seconds": round(time.process_time() - c0, 1)}


def jobs_for(task, chunks):
    comps, r = json.loads(task)
    c = Collect()
    c.search(tuple(tuple(x) for x in comps), r)
    head = {"nodes": c.nodes, "solutions": c.solutions, "unresolved": c.unresolved}
    out = []
    for node in c.nodes2:
        lo = max(p for p, _, _ in node) + 1
        hi = 1 << 38
        cuts = sorted({int(round(lo * (hi / lo) ** (i / chunks))) for i in range(chunks + 1)})
        cuts[0], cuts[-1] = lo, hi
        for a, b in zip(cuts, cuts[1:]):
            out.append((node, a, b - 1))
    return head, out


if __name__ == "__main__":
    path, chunks, workers = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    tasks = sys.argv[4:]
    if len(tasks) == 1 and tasks[0].startswith("@"):        # @file: one task per line
        tasks = [l.strip() for l in open(tasks[0][1:]) if l.strip()]
    jobs, heads = [], []
    for t in tasks:
        h, j = jobs_for(t, chunks)
        heads.append(h)
        jobs += j
    done = set()
    if os.path.exists(path):
        for l in open(path):
            if l.strip():
                r = json.loads(l)
                done.add((r["node"], r["qa"], r["qb"]))
    todo = [j for j in jobs if (json.dumps(tuple(tuple(c) for c in j[0])), j[1], j[2]) not in done]
    # dearest-looking first: small q pieces of each node
    print(f"{len(tasks)} subtrees -> {len(jobs)} pieces, {len(done)} done; heads {heads}", flush=True)
    t0 = time.time()
    with open(path, "a") as f, Pool(workers, initializer=raise_priority) as pool:
        for i, rec in enumerate(pool.imap_unordered(piece, todo, chunksize=1), 1):
            f.write(json.dumps(rec) + "\n")
            f.flush()
            if rec["solutions"] or rec["unresolved"]:
                print("FLAG", rec, flush=True)
            if i % 50 == 0 or i == len(todo):
                print(f"  {i}/{len(todo)} pieces, {time.time()-t0:.0f}s", flush=True)
    print("DONE", flush=True)
