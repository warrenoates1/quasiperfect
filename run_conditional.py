"""Eliminate one omega for one set of barred primes, in parallel.

The conditional searches of the paper are single searches and so single-threaded,
which wasted thirteen of fourteen cores. This splits one at a fixed prime depth
and runs the subtrees across cores, as run_parallel.py does for the unconditional
case, with the barred primes passed through to every worker -- omitting them there
would silently make a conditional search unconditional, which would be sound but
would answer a different and much harder question.

A value of omega counts as eliminated only if no solution is found AND no state is
left unresolved. Results append to results/cond_<barred>_<omega>.jsonl, so a run
can be resumed.

Usage: python run_conditional.py OMEGA DEPTH BARRED [WORKERS]
  BARRED  comma separated, e.g. 3,5
"""
import atexit
import json
import os
import subprocess
import sys
import time
from multiprocessing import Pool

from run_parallel import Splitter, key, run_task


def take_lock(lock):
    """Create lock holding our pid; refuse if a live process already holds it."""
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            pid = int(open(lock).read().strip() or 0)
        except (ValueError, OSError):
            pid = 0
        if pid and pid_alive(pid):
            return False
        os.remove(lock)                      # stale: holder died
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            return False
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    atexit.register(lambda: os.path.exists(lock) and os.remove(lock))
    return True


def pid_alive(pid):
    if os.name == "nt":
        out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"],
                             capture_output=True, text=True).stdout
        return str(pid) in out
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def main():
    k, depth = int(sys.argv[1]), int(sys.argv[2])
    barred = tuple(sorted(int(x) for x in sys.argv[3].split(",") if x))
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 14

    sp = Splitter(depth, exclude=barred)
    sp.search((), k)
    tasks = sp.tasks
    tag = "-".join(map(str, barred))
    print(f"omega={k}, barred {barred}, split at depth {depth}: {len(tasks):,} subtrees",
          flush=True)
    print(f"head tree: {sp.nodes:,} states, {len(sp.unresolved)} unresolved, "
          f"solutions {sp.solutions}", flush=True)
    if not tasks:
        print(f"DONE omega={k}: nothing below the split; head alone settles it", flush=True)
        return

    os.makedirs("results", exist_ok=True)
    path = f"results/cond_{tag}_{k}.jsonl"
    # One writer per results file. Two processes appending to the same jsonl would
    # duplicate work and interleave records, which is how a resume gets corrupted.
    # This matters because runs are queued back to back: a run started by hand while
    # a queued one is still pending would otherwise collide when the queue advances.
    lock = path + ".lock"
    if not take_lock(lock):
        sys.exit(f"another process is already writing {path}; refusing to start a second writer")
    done = {}
    if os.path.exists(path):
        for line in open(path):
            try:
                r = json.loads(line)
                done[r["task"]] = r
            except (ValueError, KeyError):
                pass
    todo = [(c, r, barred) for c, r in tasks if key((c, r)) not in done]
    print(f"{len(done):,} recorded, {len(todo):,} to run on {workers} workers", flush=True)

    t0 = time.time()
    with open(path, "a") as out, Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_task, todo, chunksize=1), 1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            done[rec["task"]] = rec
            if rec["solutions"] or rec["unresolved"]:
                print(f"  *** FLAG {rec}", flush=True)
            if i % 200 == 0 or i == len(todo):
                os.fsync(out.fileno())
                el = time.time() - t0
                print(f"  {i:,}/{len(todo):,}  {sum(x['nodes'] for x in done.values()):,}"
                      f" states  {el/60:,.1f} min  (eta {el/i*(len(todo)-i)/60:,.0f} min)",
                      flush=True)

    recs = [done[key((c, r))] for c, r in tasks]
    nodes = sp.nodes + sum(x["nodes"] for x in recs)
    sols = [s for x in recs for s in x["solutions"]] + [str(v) for v in sp.solutions]
    unres = sum(len(x["unresolved"]) for x in recs) + len(sp.unresolved)
    verdict = "ELIMINATED" if not sols and not unres else "NOT ELIMINATED"
    # Completion marker, so the resume task can tell "finished" from "interrupted"
    # without re-deriving the split. Written only once every subtree is recorded.
    with open(path + ".done", "w") as f:
        f.write(json.dumps({"omega": k, "barred": barred, "nodes": nodes,
                            "solutions": sols, "unresolved": unres,
                            "verdict": verdict}) + chr(10))
    print(f"\nDONE omega={k}, barred {barred}: {nodes:,} states, solutions {sols}, "
          f"unresolved {unres} -> {verdict}  ({(time.time()-t0)/60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
