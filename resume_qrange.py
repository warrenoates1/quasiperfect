"""Resume a 'q-range' unresolved node from where q_limit stopped it.

usage: python resume_qrange.py RESULTS.jsonl [--probe N]

Re-enters Search.enumerate_next at the recorded node with the start prime set
to the q where enumeration stopped (already searched; next_prime moves past it)
and q_limit lifted, so the remaining range runs to its natural end, where the
(3) bound ub <= 2 stops it. Same code path as the original run, so the result
is equivalent to having run that node with no limit.
"""
import ast, json, sys, time
from fractions import Fraction as F
from qp import Search

path = sys.argv[1]
probe = int(sys.argv[sys.argv.index("--probe") + 1]) if "--probe" in sys.argv else None
rec = [json.loads(l) for l in open(path) if l.strip()]
nodes = [ast.literal_eval(u) for r in rec for u in r["unresolved"]]
assert len(nodes) == 1 and nodes[0][0] == "q-range", nodes
_, comps, r, qstop = nodes[0]

s = Search(c=1, q_limit=probe if probe else 10**18, exclude={3, 5})
lo, hi, nmin = s.bounds(comps)
print("r", r, "qstop", qstop, "hi", float(hi), "lo", float(lo))
print("approx q end", float(r * hi / (2 - hi)))
t0 = time.time()
s.enumerate_next(comps, r, lo, hi, nmin, qstop)
out = {"comps": repr(comps), "r": r, "start_after": qstop, "probe": probe,
       "nodes": s.nodes, "solutions": s.solutions,
       "unresolved": [repr(u) for u in s.unresolved], "seconds": round(time.time() - t0, 1)}
print(json.dumps(out))
if not probe:
    open(path.replace(".jsonl", ".qrange.json"), "w").write(json.dumps(out) + "\n")
