# KinShock2020 — WarpX replication of Schaeffer et al. 2020

WarpX 1D3V PIC replication of a **piston-driven perpendicular magnetized collisionless
shock** (Schaeffer et al., *Phys. Plasmas* **27**, 042901, 2020). Goal: reproduce the
paper's shock-formation signatures (compression, reflected ions, B pileup) with a
config-driven WarpX deck + a custom kinetic heating operator.

## Read these for context (don't duplicate them here)
- `OVERVIEW.md` — physics, goals, the paper's method (the physics reference).
- `REPLICATION_PLAN.md` — run/analysis plan, Table I parameters, resolution targets.
- `RESULTS.md` — **the running lab notebook.** Dated entries for every run & finding.
  Read it first to learn current state (half-domain validation, benchmarks, B-field study).
- `scripts/README.md`, `studies/README.md` — tool docs.

## The one rule: `config.yaml` is the single source of truth
`runs/<ID>/config.yaml` holds the intuitive primaries (densities as fractions of n0,
θ = kT/mₑc², lengths in d_e/d_i, speeds/c). `scripts/make_inputs.py` renders the WarpX
deck from it. **Never hand-edit a deck** (`inputs_kinshock_*`) — edit config.yaml and
regenerate; `--verify` checks `warpx_used_inputs` against the config after a run.

## Layout
- `src/kinshock/` — the package: `units` (scale derivation), `deck` (config→deck), `io`, `metrics`.
- `scripts/` — analysis/driver CLIs (see table below). All take a `run_dir`, are config-driven.
- `studies/` — heavier experiments that *launch* WarpX (benchmarks, convergence tests).
- `tests/` — fast pytest checks (`test_structures.py`).
- `runs/<ID>/` — config + deck + `warpx_used_inputs` (tracked); `diags/`, `*.log` gitignored.
- `media/<ID>/`, `media/testing/` — figures/movies (gitignored, regenerable).

## Runs
- `R0`, `R0_half` — smoke tests (short, pre-shock).
- `R1_core`, `R1_core_half` — physics runs (dz=0.3 d_e; `_half` = one-sided z≥0 domain).
- Scales: mᵢ/mₑ=100, d_i0=100 d_e, ρ_i0≈1040 d_e, B0 along x (perpendicular shock).

## Typical workflow
```bash
python scripts/make_inputs.py runs/<ID>                 # config.yaml -> deck
OMP_NUM_THREADS=8 <warpx.1d> runs/<ID>/inputs_* > runs/<ID>/run.log 2>&1 &
python scripts/run_progress_logger.py runs/<ID> &       # wall-clock progress/ETA log
python scripts/make_inputs.py runs/<ID> --verify        # deck == config?
python scripts/run_checks.py   runs/<ID>                # scales vs Table I, conservation
python scripts/make_figures.py runs/<ID>                # A–D diagnostics
```
| Script | Purpose |
|---|---|
| `make_inputs.py` | config.yaml → WarpX deck (+ `--verify`, `--check`) |
| `run_checks.py` | bring-up / conservation checks |
| `make_figures.py` | paper diagnostics A–D |
| `make_movies.py` | density + phase-space movies |
| `run_progress_logger.py` | sidecar wall-clock progress/ETA log (`<run_dir>/progress.log`) |
| `bfield_diagnostic.py` | B-field fluctuation: physical vs numerical (spectra/polarization/particle-response) |

## Hard-won conventions & gotchas
- **Half-domain (`layout: one_sided`).** Keep `foil.lo`/`injector.lo` at `−slab` even
  one-sided — the domain clips the heated region to `[0, slab]`; moving them to 0
  **doubles the PSC heating rate** (rate ∝ 1/width). This is why foil geometry is *not*
  rewritten for the half domain. BC token map lives in `kinshock.deck._BC_MAP`.
  Half-domain reproduces full-domain z≥0 to ~3–13% (see RESULTS); bulk energy = 0.5× exactly.
- **Performance.** Launch with `OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=cores`
  — near-linear to 8 cores (~1.8× vs 4), memory-bandwidth-bound beyond. `max_grid_size`,
  tiling, and `sort_intervals` were **benchmarked as neutral-to-negative** here — don't bother.
- **Physics caveat.** dz/λ_D ≈ 7 (Debye under-resolved) → grid heating. Near-shock foot
  turbulence is physical; small-scale far-upstream B fluctuations are (being confirmed as)
  numerical. See `studies/bfield_convergence/` and the RESULTS B-field entry.
- **Env.** conda env at `/opt/anaconda3/envs/physics`; WarpX binary
  `/home/hhelal/warpx-cda/build/bin/warpx.1d` (OMP/CPU, double precision).

## Working preferences
- Work in the **regular repo folders** (not git worktrees). Commit to `main`.
- Keep `RESULTS.md` updated with a dated entry per substantive run/finding — that is how
  context survives between sessions. Anything worth keeping goes in the repo (scratch under
  a job's tmp does not persist).
