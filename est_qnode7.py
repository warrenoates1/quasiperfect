"""Whole-run estimate for omega = 7 with the node kernel, from the uniform probe sample."""
import json, time
from multiprocessing import Pool
from qp import Search
import qnode

class NS(qnode.NodeSearch, Search):
    pass

def one(t):
    comps, r = t
    s = NS(c=1)
    t0 = time.time()
    s.search(tuple(tuple(c) for c in comps), r)
    return time.time() - t0, getattr(s, 'kstates', 0), s.nodes, len(s.solutions), len(s.unresolved)

if __name__ == '__main__':
    from priority import raise_priority
    L = [json.loads(l) for l in open('results/strat_omega7.jsonl')]
    P = [l for l in L if l['pass'] in (1, 'probe')]
    tasks = [json.loads(l['task']) if isinstance(l['task'], str) else l['task'] for l in P]
    with Pool(8, initializer=raise_priority) as p:
        R = p.map(one, tasks, chunksize=1)
    secs = sum(r[0] for r in R)
    print(f"{len(R)} sampled subtrees: {secs:.0f} core-s, kernel states {sum(r[1] for r in R)}, "
          f"solutions {sum(r[3] for r in R)}, unresolved {sum(r[4] for r in R)}")
    print(f"max {max(r[0] for r in R):.1f}s; scaled to 80,458 subtrees: {secs/len(R)*80458/3600:.1f} core-hours")
    print(sorted(round(r[0], 2) for r in R)[-10:])
