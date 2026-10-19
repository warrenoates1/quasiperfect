"""Spot check of the qfast-closed node with the original exact code.

Runs the body of Search.enumerate_next on N randomly sampled primes q from the range qfast
closed (one q at a time), in parallel.
Independent of qfast: shares only qp.py. Agreement is necessary, not sufficient.
usage: python spot_qrange.py N WORKERS
"""
import ast, json, random, sys, time
from multiprocessing import Pool
from qp import Search, next_prime

def one(q):
    rec = [json.loads(l) for l in open('results/cond_3-5_17.jsonl')]
    _, comps, r, _ = ast.literal_eval(rec[-1]['unresolved'][0])
    from fractions import Fraction as F
    s = Search(c=1, exclude={3, 5})
    lo, hi, n = s.bounds(comps)
    # the body of enumerate_next's loop for this one q
    if (lo * F(q * q + q + 1, q * q) <= 2 + F(1, n * q * q)
            and s.coprime_ok_add(comps, q)):
        s.search(comps + ((q, 1, False),), r - 1)
    return q, s.nodes, s.solutions, len(s.unresolved)

if __name__ == '__main__':
    from priority import raise_priority
    N, W = int(sys.argv[1]), int(sys.argv[2])
    random.seed(20261016)
    qs = sorted({next_prime(random.randrange(571_287_683, 1_137_878_224)) for _ in range(N)})
    t0 = time.time()
    with Pool(W, initializer=raise_priority) as p:
        res = p.map(one, qs, chunksize=4)
    out = {"sampled": len(res), "states": sum(r[1] for r in res),
           "solutions": [r[2] for r in res if r[2]], "unresolved": sum(r[3] for r in res),
           "q_min": qs[0], "q_max": qs[-1], "seconds": round(time.time() - t0)}
    open('results/cond_3-5_17.spot.json', 'w').write(json.dumps(out) + "\n")
    print(json.dumps(out))
