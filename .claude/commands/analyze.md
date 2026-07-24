---
description: Run the standard analysis suite for a completed run
argument-hint: [run_id]
---
Run the analysis for `$1` (a completed run with plotfiles under runs/$1/diags):

1. `python scripts/make_inputs.py runs/$1 --verify` — confirm the run matched its config.
2. `python scripts/run_checks.py runs/$1` — scales vs Table I + conservation.
3. `python scripts/make_figures.py runs/$1` — paper diagnostics A–D.
Summarize the key numbers (compression, front, reflected fraction) and flag anything
off vs the paper / vs RESULTS.md. Then append a dated entry to RESULTS.md.
