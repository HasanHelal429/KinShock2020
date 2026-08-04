# Speedup study — R1_paper_470eV

Staged and **run** 2026-08-04. Three levers, tested individually, each against the
production deck `runs/R1_phase/R1_paper_470eV` so the numbers apply to the run we actually
care about.

> **Results are in [`runs/opt_phase/SUMMARY.md`](../../runs/opt_phase/SUMMARY.md)** (generated
> by `make_summary.py`, never hand-edited) with figures in `media/opt_phase/` and the
> narrative in `RESULTS.md`. Headline: **GPU 7.89×** with `amr.max_grid_size=30000`,
> CPU threads **1.82×** at 20, θ-implicit **13.5× slower** and therefore not a speedup.
> This file documents the *harness and reasoning*; the predictions below are kept as
> written so they can be checked against what actually happened.

Baseline to beat: **0.1113 s/step** (WarpX's own timer) at 8 OMP threads on an otherwise
idle chablis, from `runs/R1_phase/R1_paper_470eV_pilot`. The full run is 3,224,046 steps →
**5.33 d idle / 7.39 d at load 24**.

## Why these three and not mesh refinement

The pilot's TinyProfiler (`R1_paper_470eV_pilot/run.log:202640`) settles it:

| region | % of 5567 s |
|---|---|
| `GatherAndPush` | 49.9 |
| `CurrentDeposition` | 20.9 |
| `BinaryCollision::LoopOverCollisions` | 9.9 |
| `Redistribute_partition` | 8.5 |
| `ApplyBoundaryConditions` | 4.2 |
| **all field/grid work** | **~0.9** |

98% particles, 6.0e6 macroparticles of which 99.3% are ambient at a flat 100 ppc across
all 30,000 cells. Any lever that does not touch particle-steps cannot matter, which is
what rules out mesh refinement (ceiling ~2×, and dominated by plain coarsening).

## The tests

Run them one at a time. Each script prints a summary table at the end; per-point
artifacts land in `out/<label>/` (`run.log`, `meta.txt`, `result.txt`, `diags/`).

```
studies/speedup/lever1_threads.sh          # ~19 min expected, 55 min hard cap
studies/speedup/lever2_implicit.sh A       # ~1-5 min   go/no-go
studies/speedup/lever2_implicit.sh B       # ~8-30 min  per-step cost at same dt
studies/speedup/lever2_implicit.sh C       # ~30-60 min large-dt Newton, 2 cfl points
studies/speedup/lever3_gpu.sh              # ~6-10 min, IF no rebuild needed
```

`bench.sh` suppresses plotfiles (277 MiB each; disk is at 94%) but keeps the EP/PN
reduced diags, which cost 11 s per 50k steps and are what lever 3's agreement check
compares. It reports the mean of WarpX's per-step `This step = X s` with the first 20%
of steps dropped — the pilot's step 1 took 0.241 s against a 0.111 asymptote, so the
warm-up has to come out.

Every point has a wall-clock cap, so a pathological configuration costs the cap and not
an afternoon. A capped point still reports a valid s/step from the steps it completed.

### Lever 1 — OMP threads (8 → ?)

RESULTS 2026-07-23 measured near-linear scaling to 8 cores but never went above it; the
standing note is a ~20× collapse somewhere above 12. This finds the knee. Tile supply is
not the constraint — the deck runs **235 grids** over 30,000 cells (112–128 cells each,
`run.log:28`), so even 24 threads gets ~10 boxes apiece. Two `max_grid_size=64` variants
test whether granularity rather than count is what bites.

Upside if it scales: **1.4–1.6×**, free, zero physics change.

### Lever 2 — `algo.evolve_scheme = theta_implicit_em`

The only lever that attacks the *blocker* rather than the cost. `parameters.rst:246-249`
claims exactly the two properties the pilot failed on:

> Robust to finite-grid instability (does not require cells that resolve the plasma Debye length)
>
> Numerically stable for large Δt (does not require resolving the plasma period or satisfying the CFL condition for light waves)

and θ = 0.5 is *exact* energy conservation — no grid heating by construction, rather than
by resolution.

Compatibility was checked in source, not assumed:

- **Custom operators are safe.** `TargetInjector` and `ParticleHeater` are applied in the
  outer Evolve loop (`WarpXEvolve.cpp:286,291`), before the `if (m_implicit_solver)`
  branch at `:299` — scheme-agnostic. Contrast `TargetInjector.cpp:225`, which *does*
  hard-abort under mesh refinement; there is no equivalent guard for implicit.
- **Collisions are supported.** `doCollisions()` is called inside the implicit branch
  (`:301`), and the scheme's first cited reference is Angus et al., *"…implicit
  particle-in-cell method coupled with a binary Monte-Carlo algorithm for Coulomb
  collisions"* — this pairing is the published use case.
- **Deposition.** The deck uses Esirkepov (`run.log:39`, the FDTD default), which is
  allowed — but **not** with `use_mass_matrices_jacobian`, the option documented to give
  "large speed ups for simulations with many particles". Stage C therefore switches to
  `villasenor` (charge-conserving *and* mass-matrix compatible).
- **1D works.** `Examples/Tests/implicit/inputs_test_1d_theta_implicit_picard`, also at
  `particle_shape = 2`.

Expect implicit to be **slower per step** — each Picard/Newton iteration redoes the
push+deposit that is 71% of runtime. It wins only via larger dt, which is what stage C
measures. Stage A's real output is the iteration count per step, since that sets the
multiplier.

### Lever 3 — GPU

Best cost/benefit, worst logistics. 98% particle work over 6.0e6 particles in 1D is close
to ideal GPU shape and fits one 4070's 12 GiB easily. **Both GPUs were idle at 15 MiB /
0% when this was staged.**

The binary at `build_cuda1d/bin/warpx.1d` is dated **Jul 28**, and commit `9f981dea2`
landed **Jul 31** with the message *"nvcc rejects an extended `__device__` lambda inside a
private or protected member function … which broke the CUDA build"* — i.e. the binary
predates a CUDA-specific fix to the very operator this run depends on. The change is
behaviour-neutral by its own description, so the stale binary is probably fine, which is
why stage 3C is an **agreement check against the CPU binary** rather than another timing
run.

A rebuild is the fallback, not the first move: `ParticleHeater.H` is included by
`Source/WarpX.H`, which nearly all 391 translation units include, so the edit invalidates
almost the entire CUDA build and ccache cannot help (every TU's preprocessed output
changed). **30–50 min at -j24, up to ~90.**

3A uses the deck's decomposition, which is CPU-tuned and bad for a GPU (~25k particles
per kernel launch); 3B/3B2 test 1 and 8 boxes. GPU 1 is left free for other users
(`CUDA_VISIBLE_DEVICES=0`).

## Follow-on work these tests do NOT cover

- **Whether implicit actually kills the grid heating.** Lever 2 measures cost and
  convergence only. Confirming it removes the pilot's +3.9%/3.41 t_ab needs a run matched
  to the pilot's physical time and analysed with `scripts/grid_heating.py` — and needs
  `numerics.evolve_scheme` plumbed through `scripts/make_inputs.py` first. No passthrough
  block exists in config.yaml today, and the repo rule is never to hand-edit a deck. The
  ParmParse overrides used here are legitimate for benchmarking, not for physics.
- **Noticed in passing, unrelated:** the pilot log reports
  `Unused ParmParse Variables: [TOP]::boundary.reflect_symmetry_axis`. The deck sets a key
  WarpX ignores. `boundary.field_lo = pec` and `particle_lo = reflecting` are doing the
  z=0 wall, so this looks like a harmless leftover — but it is unverified and worth its
  own look.
