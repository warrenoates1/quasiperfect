"""Print where the search is every N nodes, to see what dominates the tree."""
import sys

from qp import Search

k = int(sys.argv[1])
every = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
stop = int(sys.argv[3]) if len(sys.argv) > 3 else 400000


class Stop(Exception):
    pass


def log(s):
    print(s.nodes, s.current, flush=True)
    if s.nodes >= stop:
        raise Stop


s = Search(c=1, log=log)
import qp
qp_log_every = every
orig = Search.search


def search(self, comps, r):
    if self.nodes % every == 0 and self.nodes:
        print(self.nodes, r, comps, flush=True)
        if self.nodes >= stop:
            raise Stop
    return orig(self, comps, r)


Search.search = search
try:
    s.run(k)
except Stop:
    pass
