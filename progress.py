"""Show what each running worker is doing, from the files it writes.

Workers write results/progress/<pid>.json every 200,000 states (see
run_parallel._pause_hook) and delete it when their subtree finishes. This reads
them, so progress is visible even for runs in session 0 that this session cannot
open. A file whose pid is no longer alive is reported as stale rather than
trusted: a worker that died leaves its last file behind.

Usage: python progress.py
"""
import glob
import json
import os
import subprocess
import time

from run_parallel import PAUSE_FILE, PROGRESS_DIR


def alive(pid):
    out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"],
                         capture_output=True, text=True).stdout
    return str(pid) in out


def main():
    files = sorted(glob.glob(os.path.join(PROGRESS_DIR, "*.json")))
    print("PAUSED" if os.path.exists(PAUSE_FILE) else "running", "-",
          len(files), "worker progress file(s)")
    now = time.time()
    for f in files:
        try:
            r = json.load(open(f))
        except (ValueError, OSError):
            continue
        rate = r["nodes"] / r["seconds"] if r["seconds"] else 0
        cw = r["cpu_seconds"] / r["seconds"] if r["seconds"] else 0
        primes = [c[0] for c in json.loads(r["task"])[0]]
        state = "stale (pid gone)" if not alive(r["pid"]) else (
            "paused" if r.get("paused") else "working")
        print(f"  pid {r['pid']:>6}  {state:16}  {r['nodes']:>14,} states  "
              f"{r['seconds']/3600:6.2f} h  {rate:8,.0f}/s  cpu/wall {cw:.2f}  "
              f"updated {now - r['updated']:.0f}s ago  primes {primes}")


if __name__ == "__main__":
    main()
