# `studies/` — numerical experiments and compute studies

Heavier, multi-run experiments that *launch WarpX* (as opposed to `scripts/`, which
analyze an existing run, or `tests/`, which are fast pytest unit checks). Each study
has its own directory with a `README.md`, a runner, and an analysis script. WarpX
scratch output goes to a gitignored `scratch/` inside the study (regenerable from the
runner), and figures land in `media/testing/`.

| Study | Question | Runner | Analysis |
|---|---|---|---|
| `bfield_convergence/` | Are the B-field oscillations physical or a numerical (grid) artifact? | `run_variants.sh` | `analyze.py`, `scripts/bfield_diagnostic.py` |

Findings are recorded in the top-level `RESULTS.md`.
