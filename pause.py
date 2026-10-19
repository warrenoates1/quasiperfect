"""Pause or resume the running quasiperfect searches, without losing work.

Killing a run costs whatever its workers hold in flight, since a subtree is only
recorded once it completes. Pausing costs nothing.

This used to suspend processes through NtSuspendProcess. That cannot reach runs
the S4U scheduled task relaunches: they live in session 0, and an ordinary
signed-in session gets access denied on every OpenProcess against them - query,
suspend and terminate alike (found 3 Oct 2026, after the original runs died with
the session they were started in and the task brought them back). So the pause
is now a file: every worker checks for results/PAUSE each 200,000 states (about
every ten seconds) and sleeps while it exists. Works from any session, for runs
started anywhere, and the time spent paused is recorded and excluded from each
subtree's timing.

Usage: python pause.py pause | resume | status
"""
import os
import sys

from run_parallel import PAUSE_FILE


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "status"
    if what in ("pause", "suspend"):
        open(PAUSE_FILE, "w").close()
        print("paused: workers stop at their next checkpoint, within ~10 seconds")
    elif what == "resume":
        if os.path.exists(PAUSE_FILE):
            os.remove(PAUSE_FILE)
        print("resumed")
    else:
        print("PAUSED" if os.path.exists(PAUSE_FILE) else "running (no PAUSE file)")


if __name__ == "__main__":
    main()
