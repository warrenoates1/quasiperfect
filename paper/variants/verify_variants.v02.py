"""Check that the faulty variants really do fail, and the fixed ones really pass.

Section~\\ref{sec:history} of the paper claims specific defects in specific
versions. A claim like that should be tested, not asserted, so this runs each
named variant against the defect it is supposed to exhibit:

  Issue 1  qp.v01 (before) must leave unresolved states on the planted
           last-prime cases; qp.v07 (after) must leave none.
  Issue 8  qp.v15 (before) must fail to terminate on omega = 5 within a
           generous node budget; qp.v17 (after) must finish in 193,559 nodes.

Each variant runs in its own temporary directory so nothing shadows the
working copy. Usage: python verify_variants.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
VAR = os.path.join(HERE, "variants")
PROJ = os.path.abspath(os.path.join(HERE, ".."))

LEAF_PROBE = r"""
import random, sys
from fractions import Fraction
from qp import Search, sig, is_probable_prime
# The same planted last-prime construction as test_leaves.py, fixed seed, so the
# two variants are asked exactly the same questions.
rng = random.Random(1)
odd = [p for p in range(3, 2000) if is_probable_prime(p)]
flags = recovered = tried = 0
while tried < 200:
    w = rng.randint(2, 5)
    ps = sorted(rng.sample(odd[:15], w - 1)); ps.append(rng.choice([p for p in odd if p > ps[-1]]))
    exps = [rng.choice([1, 1, 1, 2, 3, 4, 7, 12]) for _ in ps]
    n = s = 1
    for p, a in zip(ps, exps):
        n *= p ** (2 * a); s *= sig(p, 2 * a)
    c = s - 2 * n
    if c <= 0: continue
    m = sm = 1
    for p, a in zip(ps[:-1], exps[:-1]):
        m *= p ** (2 * a); sm *= sig(p, 2 * a)
    if sm >= 2 * m: continue
    tried += 1
    S = Search(c=c)
    S.last_prime(tuple((p, a, True) for p, a in zip(ps[:-1], exps[:-1])), Fraction(sm, m), m, ps[-2])
    if S.unresolved: flags += 1
    if n in S.solutions: recovered += 1
print(f"tried={tried} recovered={recovered} unresolved_cases={flags}")
"""

TERM_PROBE = r"""
import sys
from qp import Search
BUDGET = 4_000_000
class Capped(Search):
    def search(self, comps, r):
        if self.nodes > BUDGET:
            raise SystemExit(f"DID NOT TERMINATE: over {BUDGET:,} nodes; "
                             f"deepest exponent floor seen = {self.deepest}")
        self.deepest = max(getattr(self, 'deepest', 0), max((t[1] for t in comps), default=0))
        return super().search(comps, r)
s = Capped(c=1)
s.deepest = 0
s.run(5)
print(f"TERMINATED nodes={s.nodes} solutions={s.solutions} unresolved={len(s.unresolved)}")
"""


def run(variant, probe, label):
    d = tempfile.mkdtemp(prefix="qpvar_")
    try:
        shutil.copy(os.path.join(VAR, variant), os.path.join(d, "qp.py"))
        open(os.path.join(d, "probe.py"), "w").write(probe)
        p = subprocess.run([sys.executable, "probe.py"], cwd=d,
                           capture_output=True, text=True, timeout=1800)
        out = (p.stdout + p.stderr).strip().splitlines()
        print(f"  {variant:16s} {label:34s} -> {out[-1] if out else '(no output)'}")
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    print("Issue 1 - last-prime root finding (planted cases, same seed):")
    run("qp.v01.py", LEAF_PROBE, "before: expect unresolved > 0")
    run("qp.v07.py", LEAF_PROBE, "after:  expect unresolved = 0")
    print("Issue 8 - termination of refinement (omega = 5):")
    # v15 introduces the faulty rule but references a parameter that only appears
    # in v16, so v16 is the first runnable version carrying the defect.
    run("qp.v16.py", TERM_PROBE, "before: expect no termination")
    run("qp.v17.py", TERM_PROBE, "after:  expect 193,559 nodes")
