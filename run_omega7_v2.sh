#!/bin/sh
# omega = 7 re-run on qnode with magnitude-scaled slack (19 Oct 2026)
python -u run_qnode.py 7 5 4 @results/stragglers.txt > results/qnode_omega7.log 2>&1
python -u qchunks.py results/qchunks_omega7.jsonl 600 4 @results/stragglers.txt > results/qchunks_omega7.log 2>&1
python -u finalize_omega7.py > results/omega7_final.log 2>&1
