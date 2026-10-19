"""Settle qfast candidates exactly with the original Search code."""
from qp import Search

def candidate_comps(comps, opens, row):
    nop = len(opens)
    q, qq = int(row[0]), int(row[1])
    box = {p: (int(row[2 + i]), bool(row[2 + nop + i])) for i, (p, _) in enumerate(opens)}
    new = tuple((p, *box[p]) if p in box else (p, a, fx) for p, a, fx in comps)
    return new + ((q, 1, False), (qq, 1, False))

def verify(comps, opens, rows, exclude=(3, 5)):
    s = Search(c=1, exclude=exclude)
    for row in rows:
        s.search(candidate_comps(comps, opens, row), 0)
    return s
