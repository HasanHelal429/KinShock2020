# KinShock2020

Reproduction and verification of the piston-driven magnetized collisionless shocks of

> D. B. Schaeffer, W. Fox, J. Matteucci, K. V. Lezhnin, A. Bhattacharjee, K. Germaschewski,
> *"Kinetic simulations of piston-driven collisionless shock formation in magnetized
> laboratory plasmas,"* Phys. Plasmas **27**, 042901 (2020), doi:10.1063/1.5123229

using **WarpX**'s ported PSC laser-ablation surrogate (`ParticleHeater` + `TargetInjector`).
The goal is to verify that these modules can generate Schaeffer-2020-class perpendicular
collisionless shocks and reproduce the paper's quantitative shock-formation signatures.

The paper PDF is **not** included (copyrighted); get it from
[doi:10.1063/1.5123229](https://doi.org/10.1063/1.5123229) or
[OSTI 1643946](https://www.osti.gov/pages/biblio/1643946).

## Documents
- **`OVERVIEW.md`** — physics summary of the paper (setup, Table I, the 7 shock-formation
  criteria, the three timescales, parameter scans).
- **`REPLICATION_PLAN.md`** — the reproduction plan: parameter mapping, run matrix, analysis
  spec, runtime estimates.
- **`RESULTS.md`** — running validation log (structure tests, R0 smoke, heater-calibration
  cross-check, …).

## Layout
```
src/kinshock/     reusable analysis library (config, units, io, metrics, plotting)
scripts/          thin drivers: make_config, run_checks, make_figures, make_movies
runs/RXX/         per-run config.yaml (single source of truth) + WarpX deck (inputs_*)
tests/            structure tests (config/units/round-trip/metrics), no WarpX needed
media/            figures & movies (testing/ = bring-up, RXX/ = per-run results)
```
Simulation output (`runs/*/diags/`) and run logs are git-ignored — regenerate by running the
decks. Per-run parameters live in `runs/RXX/config.yaml`; **no physical constants are
hard-coded in the scripts** — `src/kinshock/units.py` derives all scales from the config
primaries, and `scripts/make_config.py` can regenerate/verify a config from a run's
`warpx_used_inputs`.

## Design principle
Each run's `config.yaml` is the single source of truth. `units.py` derives every physical /
normalized scale (ω_pe, d_e, d_i, ω_ci0, ρ_i0, C_s,ab, v_A, M_A, M_ms, …) from it, and the
driver scripts read the config rather than embedding constants.

## Running
Requires a WarpX build with the `ParticleHeater` + `TargetInjector` modules, and a Python
env with `yt`, `numpy`, `matplotlib`, `pyyaml`.

```bash
# structure tests (no WarpX needed)
python tests/test_structures.py

# a run (example: R0 smoke)
cd runs/R0
MPICH_GPU_SUPPORT_ENABLED=0 mpirun -np <N> <warpx.1d> inputs_kinshock_R0

# verify the config matches what actually ran, then analyze
python scripts/make_config.py runs/R0            # verify vs warpx_used_inputs
python scripts/run_checks.py  runs/R0            # bring-up / progress figures -> media/testing
python scripts/make_figures.py runs/R0           # shock diagnostics -> media/R0
python scripts/make_movies.py  runs/R0           # movies -> media/R0
```

## Runs
| id | tier | purpose |
|----|------|---------|
| `R0` | smoke | structure / pipeline verification (tiny grid) |
| `R1_core` | core | full R1 physics on a moderate domain (→ t*₂), calibration check |
| `R1` | full | Schaeffer 2020 Table I representative run (M_A ≈ 14), → t*₃ |
| `xcheck_flatfoil_1d` | — | 1D-vs-2D heater-saturation cross-check |

See `REPLICATION_PLAN.md` §3 for the full run matrix (controls, Mach scan, multi-species).

## Status
Work in progress — see `RESULTS.md` for the current validation state.
