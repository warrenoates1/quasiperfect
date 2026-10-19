"""Time the node kernel on the dearest sampled omega = 7 subtrees."""
import json, sys, time
from qp import Search
import qnode

class NS(qnode.NodeSearch, Search):
    pass

L = [json.loads(l) for l in open('results/strat_omega7.jsonl')]
H = sorted([l for l in L if l['nodes'] > 5e6], key=lambda l: -l['nodes'])
n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for l in H[:n]:
    comps, r = json.loads(l['task']) if isinstance(l['task'], str) else l['task']
    comps = tuple(tuple(c) for c in comps)
    s = NS(c=1)
    t0 = time.time()
    s.search(comps, r)
    print(f"python {l['nodes']:>11} {'(capped)' if l['censored'] else '        '} {l['secs']:>8.0f}s | "
          f"kernel {getattr(s,'kstates',0):>12} states, {getattr(s,'kcands',0)} cands, py {s.nodes}, "
          f"{time.time()-t0:8.1f}s sols={s.solutions} unres={len(s.unresolved)}", flush=True)
