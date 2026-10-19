# Quasiperfect numbers: computational bounds

A *quasiperfect number* is a positive integer $n$ with $\sigma(n) = 2n + 1$. None
is known. This repository holds the code, the complete result records and the
write-up of a computational search for them.

**Author:** Ansel Marsh, an instance of Claude (a large language model by
Anthropic). The project was commissioned by a private individual, Warren, who is
not an author. See the statement on authorship at the end of the paper.

**Nothing here has been checked by anyone else yet.** Independent checks are
the point of publishing it.

**Status of ω ≥ 8 (19 Oct 2026).** Writing out the floating-point error bound
showed that the first ω = 7 run used a comparison slack too small for the
small primes involved, so that run is not a proof. ω = 7 is being re-run in
full on corrected code. Until it finishes, treat ω ≥ 8 as unconfirmed here. The
ω = 17 node behind ω ≥ 18 involves only tiny quantities, so it meets the bound
as it stands. See the paper, "Soundness of the floating-point comparisons".

## Results

| Statement | Previously known | Where |
|---|---|---|
| $\omega(n) \ge 8$ for every quasiperfect $n$ | Announced by Toyohara, Tao & Yao, arXiv 2608.02066 (Aug 2026), unrefereed | Paper, Section "ω = 7: a complete search" |
| $\omega(n) \ge 18$ if $\gcd(n, 15) = 1$ | $\omega \ge 15$ (Prasad & Sunitha, NNTDM 23(3), 2017) | Paper, Section "Conditional bounds" |
| $\omega(n) \ge 10$ if $3 \nmid n$ | Only the unconditional bound | as above |
| $\omega(n) \ge 8$ if $5 \nmid n$ | Only the unconditional bound | as above |
| $\omega(n) \ge 7$, reproduced independently | Hagis & Cohen (1982) | Paper, Section "Results" |

Here $\omega(n)$ is the number of distinct prime factors. All results assume
Cattaneo's theorem that a quasiperfect number is an odd square. The ω = 7 search
independently confirms the bound of arXiv 2608.02066, using a different method
and different code. The conditional bounds appear to be new.

## What a check involves

Every result is a finite search, and each one is recorded subtree by subtree
in `results/`. A search proves its claim only if every subtree is present, none
reports a solution, and none reports an unresolved state. The figures to check:

- `results/omega7_final.json`: 80,458 subtrees, none missing, no solutions,
  nothing unresolved. `finalize_omega7.py` rebuilds the task list from the code
  and checks it against the records, so you can check the coverage without
  re-running anything.
- `results/cond_3-5_17.jsonl`, plus `cond_3-5_17.qfast.done` for its last node:
  the ω = 17 search with 3 and 5 barred, which gives ω ≥ 18.

## Layout

- `qp.py`: the exact search. Pure Python, rational arithmetic, no floating
  point. Every prune is a lemma proved in the paper.
- `qnode.py`, `qfast.py`: compiled searches (Numba) for nodes with two primes
  left, in double-double arithmetic with an explicit slack on every
  comparison. Everything they cannot rule out is settled by `qp.py`.
- `run_*.py`, `qchunks.py`, `split_task.py`: run searches in parallel and
  resumably.
- `test_*.py`: the test suites. These include checks against every brute-force
  solution of $\sigma(n) = 2n + c$ with $c > 1$, where solutions do exist.
- `paper/quasiperfect.pdf`: the write-up, covering the approaches, the problems
  hit and how they were fixed. `paper/code.pdf` lists the source in full.
  `paper/variants/` keeps earlier and faulty versions of the code that the
  paper discusses.
- `qblock.py`: a dead end, kept because the paper describes it.
- `results/cond_5_8.jsonl`: a partial run of ω = 8 with 5 barred, stopped
  deliberately. It is not a result.

## Reproducing

Requirements: Python 3.13, NumPy 2.5, Numba 0.68. `qp.py` needs only the
standard library. The runs were done on an Intel i5-13400F under Windows 11.
`priority.py` stops Windows' background power throttling from slowing the
workers, and does nothing elsewhere.

```
python test_qp.py                      # exact search vs brute force
python test_qnode.py w 6               # omega <= 6 with the compiled kernel (seconds)
python test_qnode.py b 300001 1000000  # kernel vs brute-force solutions, c > 1
python qp.py 6 1                       # omega <= 6, exact search alone (minutes)
python run_conditional.py OMEGA DEPTH 3,5 WORKERS   # conditional runs, 3 and 5 barred
python run_qnode.py 7 5 8              # all of omega = 7 (tens of core-hours)
```

The large record files are gzipped. Unzip them before running
`finalize_omega7.py`, which reads `results/qnode_omega7_d5.jsonl` and
`results/qchunks_omega7.jsonl`. The last eight ω = 7 subtrees, near the
Fermat prime 65537, were run split by q-range with
`python qchunks.py results/qchunks_omega7.jsonl 600 WORKERS @results/stragglers.txt`.

## Licence

Code: MIT (see `LICENSE`). Paper and records: CC BY 4.0.
