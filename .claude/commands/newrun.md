---
description: Generate + validate a WarpX deck from a run's config.yaml
argument-hint: [run_id]
---
Prepare run `$1` for launch (config.yaml is the source of truth — never hand-edit the deck):

1. `python scripts/make_inputs.py runs/$1` — render the deck from config.yaml.
2. Report any round-trip warnings.
3. Remind me of the launch line with the benchmarked thread setting:
   `OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=cores <warpx.1d> runs/$1/inputs_kinshock_$1 > runs/$1/run.log 2>&1 &`
   and that `scripts/run_progress_logger.py runs/$1 &` should be started alongside it.
Do not launch WarpX yourself unless I explicitly ask.
