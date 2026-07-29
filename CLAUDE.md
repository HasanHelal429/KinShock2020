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
- `R1_warm` — the current full reference run (warm ablative piston ions, settled M_A≈14).
- `R1_coll` — R1_warm's collisional twin: same dimensionless setup, ambient pinned to
  10¹⁸ cm⁻³, pairwise Coulomb collisions (see the collisions gotcha below).
- Scales: mᵢ/mₑ=100, d_i0=100 d_e, ρ_i0≈1040 d_e, B0 along x (perpendicular shock).

## Typical workflow
```bash
python scripts/make_inputs.py runs/<ID>                 # config.yaml -> deck
scripts/launch.sh -b -L runs/<ID>                       # launch WarpX (+ progress logger)
python scripts/make_inputs.py runs/<ID> --verify        # deck == config?
python scripts/run_checks.py   runs/<ID>                # scales vs Table I, conservation
python scripts/tune_shock.py   runs/<ID>                # fit v_sh + front BY EYE -> shock_fit.yaml
python scripts/make_figures.py runs/<ID>                # A–D diagnostics (reads shock_fit.yaml)
```
| Script | Purpose |
|---|---|
| `make_inputs.py` | config.yaml → WarpX deck (+ `--verify`, `--check`) |
| `launch.sh` | **the** way to start a run: cd's into the run dir, benchmarked OMP settings |
| `run_checks.py` | bring-up / conservation checks |
| `tune_shock.py` | fit v_sh + front trajectory BY EYE vs the B/n_e streaks → `runs/<ID>/shock_fit.yaml` |
| `make_figures.py` | paper diagnostics A–D |
| `make_movies.py` | density + phase-space movies |
| `run_progress_logger.py` | sidecar wall-clock progress/ETA log (`<run_dir>/progress.log`) |
| `bfield_diagnostic.py` | B-field fluctuation: physical vs numerical (spectra/polarization/particle-response) |

## Hard-won conventions & gotchas
- **Shock kinematics come from `runs/<ID>/shock_fit.yaml`, fit BY EYE** (`scripts/tune_shock.py`),
  not auto-detection. It is the single source of truth for `v_sh` and `z_front(t)` (linear
  `z0 + v_sh*t`, plus optional per-time overrides); `make_figures` reads it via
  `metrics.load_shock_fit` and falls back to auto `track_front`+`speed_from_trajectory`
  (with a warning) only when it is absent. Automatic `v_sh` had drifted inconsistently
  between scripts — always tune once per run so every diagnostic shares one speed/front.
- **Half-domain (`layout: one_sided`).** Keep `foil.lo`/`injector.lo` at `−slab` even
  one-sided — the domain clips the heated region to `[0, slab]`; moving them to 0
  **doubles the PSC heating rate** (rate ∝ 1/width). This is why foil geometry is *not*
  rewritten for the half domain. BC token map lives in `kinshock.deck._BC_MAP`.
  Half-domain reproduces full-domain z≥0 to ~3–13% (see RESULTS); bulk energy = 0.5× exactly.
- **Launch with `scripts/launch.sh runs/<ID>`, never by hand.** The generated deck sets **no
  `diag*.file_prefix`**, so WarpX writes plotfiles to `diags/` *relative to the launch CWD*.
  Launching two runs from the repo root makes them **share `./diags/` and clobber each other**
  (WarpX leaves `.old.NNNN` rename files as the tell) — cost a rerun on the R2/R3 controls
  (RESULTS 2026-07-26). `launch.sh` exists to make that unrepeatable: it cd's into the run dir,
  applies the benchmarked OMP settings, picks the single deck, logs to `<run_dir>/run.log`, and
  **refuses to start when `diags/` already holds output** (`--force` to override). `-b` detaches,
  `-L` also starts the progress logger, `-n` dry-runs, and anything after `--` is passed to WarpX
  as ParmParse overrides (smoke tests only — they will trip `make_inputs.py --verify`).
- **Collisions (`collisions:` block, `runs/R1_coll`).** Three traps. (1) **lnΛ is a knob, not a
  physical value**: at real c, θ_e = 0.078 ⇒ T_e,ab = 39.9 keV, so ν_ei ∝ n T^(−3/2) keeps the
  plasma collisionless (physical lnΛ → mfp ≈ 4×10⁵ d_e,ab) at *any* attainable density. Set
  `collisions.target` to the physics you want and let `units.coulomb_log_for` invert it.
  (2) **λ_ab ≠ mfp/d_e**: the paper's λ_ab ≡ ω_ce,ab/ν_ei,ab is mfp/ρ_e,ab, and ρ_e,ab ≈ 28 d_e,
  so Table I's λ_ab = 20 means mfp = 559 d_e, *not* 20 d_e. Pick `quantity: lambda_ab` vs
  `mfp_over_de` deliberately. Since n0 is a pure scale factor, a collisional run is the only
  reason to pin an absolute density — the collisionless physics is identical either way.
  (3) **λ_ab and the upstream λ_ii/d_i0 cannot both match the paper.** Table I quotes λ_ab = 20
  *and* λ_mfp/d_i0 = 350; one lnΛ can only hit one of them at real c. R1_coll targets λ_ab = 20
  ⇒ lnΛ = 7713 = 667× physical ⇒ **upstream λ_ii/d_i0 = 0.52, so criterion 2 ("collisionless")
  FAILS in all 51 frames** (RESULTS 2026-07-29). Criterion 2 reads `Scales.mfp_ii_amb`
  (`units.nu_ii`, NRL ion-ion, μ from the real ion mass) — it used to be a hard-coded 350.0,
  which happens to equal the physical-lnΛ answer, so it silently masked this. For a run that is
  collisional *and* still a collisionless shock, target the upstream number: λ_ii/d_i0 ≈
  347×(11.567/lnΛ), so lnΛ ≈ 400 ⇒ λ_ii/d_i0 ≈ 10. R1_coll still reproduces R1_warm's shock on
  all six other criteria (peak B⊥ 12% lower) — the paper's Fig. 13 claim, tested harder.
- **Performance.** Launch with `OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=cores`
  — near-linear to 8 cores (~1.8× vs 4), memory-bandwidth-bound beyond. `max_grid_size`,
  tiling, and `sort_intervals` were **benchmarked as neutral-to-negative** here — don't bother.
- **Physics caveat.** dz/λ_D ≈ 7 (Debye under-resolved) → grid heating. Near-shock foot
  turbulence is physical (converged); far-upstream small-scale (lambda<~2-3 d_e) hash is
  numerical grid noise (filter_npass=8 cuts it 31%). See `studies/bfield_convergence/` + RESULTS.
- **`B_compression` in `criteria.json` is a GLOBAL max — cut the outer 2 d_i0 before quoting it.**
  From t*ω_ci0 ≈ 5.5 the open hi boundary throws an ~80× B⊥/B0 spike (vs ~15–19 for the real
  ramp) in every run, so late-time values are the artifact, not the shock (RESULTS 2026-07-29).
- **Env.** conda env at `/opt/anaconda3/envs/physics`; WarpX binary
  `/home/hhelal/warpx-cda/build/bin/warpx.1d` (OMP/CPU, double precision).

## Working preferences
- Work in the **regular repo folders** (not git worktrees). Commit to `main`.
- Keep `RESULTS.md` updated with a dated entry per substantive run/finding — that is how
  context survives between sessions. Anything worth keeping goes in the repo (scratch under
  a job's tmp does not persist).
