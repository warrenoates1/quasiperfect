"""Run qfast over the rest of the unresolved q-range node of a results file.

usage: python run_qfast.py RESULTS.jsonl [CHUNK]
Writes RESULTS.qfast.jsonl (one line per chunk, resumable) and, when every
chunk is in, settles any candidates with the exact code and writes
RESULTS.qfast.done.
"""
import ast, json, os, sys, time
import numpy as np
from qp import Search
import qfast, qfast_verify
from priority import raise_priority
print("priority", raise_priority(), flush=True)

path = sys.argv[1]
CHUNK = int(sys.argv[2]) if len(sys.argv) > 2 else 10_000_000
rec = [json.loads(l) for l in open(path) if l.strip()]
nodes = [ast.literal_eval(u) for r in rec for u in r["unresolved"]]
assert len(nodes) == 1 and nodes[0][0] == "q-range"
_, comps, r, qstop = nodes[0]
assert r == 2
s = Search(c=1, exclude={3, 5})
opens, K, KL, dtab, dtl, allowed, floor, amax = qfast.prepare(comps, s)
lo, hi, _ = s.bounds(comps)
# (3) at r = 2: hi (q/(q-1))^2 > 2 is necessary for the smaller prime q
n, d = hi.numerator, hi.denominator
sq = (2 * d / n) ** 0.5
qend = int(sq / (sq - 1)) + 10
while not n * qend * qend > 2 * d * (qend - 1) ** 2:
    qend -= 1
assert not n * (qend + 1) ** 2 > 2 * d * qend ** 2
qa = qstop + 1            # qstop itself was searched by the original run
out_path = path.replace(".jsonl", ".qfast.jsonl")
done = set()
if os.path.exists(out_path):
    for l in open(out_path):
        if l.strip():
            done.add(json.loads(l)["start"])
print(f"node r={r}, q from {qa} to {qend}, K={K:.6e}", flush=True)
out = np.zeros((1_000_000, 2 + 2 * len(opens)), np.int64)
for start in range(qa, qend + 1, CHUNK):
    if start in done:
        continue
    end = min(start + CHUNK - 1, qend)
    qs = qfast.primes_between(start, end)
    st = np.zeros(4, np.int64)
    t0 = time.time()
    m = qfast.run_block(qs, K, KL, dtab, dtl, allowed, floor, amax, out, st)
    cands = out[:min(m, len(out))].tolist()
    line = {"start": start, "end": end, "primes": len(qs), "states": int(st[0]),
            "unresolved": int(st[2]), "overflow": int(st[3]), "candidates": m,
            "cand_rows": cands, "seconds": round(time.time() - t0, 1)}
    with open(out_path, "a") as f:
        f.write(json.dumps(line) + "\n")
    print(start, end, len(qs), int(st[0]), m, int(st[2]), int(st[3]), line["seconds"], flush=True)

lines = [json.loads(l) for l in open(out_path) if l.strip()]
starts = sorted(l["start"] for l in lines)
assert starts == list(range(qa, qend + 1, CHUNK)), "chunks missing"
rows = [row for l in lines for row in l["cand_rows"]]
bad = sum(l["unresolved"] + l["overflow"] for l in lines) + sum(l["candidates"] > len(l["cand_rows"]) for l in lines)
v = qfast_verify.verify(comps, opens, rows)
summary = {"q_from": qa, "q_to": qend, "primes": sum(l["primes"] for l in lines),
           "states": sum(l["states"] for l in lines), "candidates": len(rows),
           "exact_nodes": v.nodes, "solutions": v.solutions,
           "unresolved": len(v.unresolved) + bad,
           "verdict": "NODE CLOSED" if not v.solutions and not v.unresolved and not bad else "NOT CLOSED"}
open(path.replace(".jsonl", ".qfast.done"), "w").write(json.dumps(summary) + "\n")
print(json.dumps(summary))
