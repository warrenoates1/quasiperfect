"""Combine the omega = 7 qnode run and the straggler pieces into one verdict, checking coverage."""
import json
from run_parallel import Splitter, key
from qchunks import jobs_for

sp = Splitter(5)
sp.run(7)
main = {}
for l in open('results/qnode_omega7_d5.jsonl'):
    if l.strip():
        r = json.loads(l)
        main[r['task']] = r
strag = [l.strip() for l in open('results/stragglers.txt') if l.strip()]
pieces = {}
for l in open('results/qchunks_omega7.jsonl'):
    if l.strip():
        r = json.loads(l)
        pieces[(r['node'], r['qa'], r['qb'])] = r

missing, sols, unres, core = [], [], 0, 0.0
pstates = ksstates = 0
heads_nodes = 0
for t in sp.tasks:
    k = key(t)
    if k in strag:
        head, jobs = jobs_for(k, 600)
        heads_nodes += head['nodes']
        sols += head['solutions']
        unres += len(head['unresolved'])
        for node, a, b in jobs:
            r = pieces.get((json.dumps(tuple(tuple(c) for c in node)), a, b))
            if r is None:
                missing.append((k, a, b))
                continue
            sols += r['solutions']
            unres += len(r['unresolved'])
            core += r['cpu_seconds']
            ksstates += r['kernel_states']
            pstates += r['exact_nodes']
    elif k in main:
        r = main[k]
        sols += r['solutions']
        unres += len(r['unresolved'])
        core += r['cpu_seconds']
        ksstates += r['kernel_states']
        pstates += r['nodes']
    else:
        missing.append((k,))
summary = {"omega": 7, "subtrees": len(sp.tasks), "head_states": sp.nodes,
           "straggler_head_states": heads_nodes, "exact_states": pstates,
           "kernel_states": ksstates, "core_hours": round(core / 3600, 1),
           "solutions": sols, "unresolved": unres, "missing": len(missing)}
summary["verdict"] = ("ELIMINATED" if not sols and not unres and not missing and not sp.solutions
                      and not sp.unresolved else "NOT ELIMINATED")
open('results/omega7_final.json', 'w').write(json.dumps(summary) + "\n")
print(json.dumps(summary))
if missing:
    print(missing[:5])
