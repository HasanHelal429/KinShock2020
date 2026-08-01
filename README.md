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
src/kinshock/     reusable analysis library (config, units, deck, io, metrics, plotting)
scripts/          thin drivers: make_inputs, run_checks, make_figures, make_movies
runs/RXX/         per-run config.yaml (single source of truth) + generated WarpX deck (inputs_*)
tests/            structure tests (config/units/deck-generation/metrics), no WarpX needed
media/            figures & movies (testing/ = bring-up, RXX/ = per-run results)
```
Simulation output (`runs/*/diags/`) and run logs are git-ignored — regenerate by running the
decks. Per-run parameters live in `runs/RXX/config.yaml`; **no physical constants are
hard-coded in the scripts** — `src/kinshock/units.py` derives all scales from the config
primaries, and `scripts/make_inputs.py` generates the WarpX deck from the config (and, with
`--verify`, confirms `warpx_used_inputs` matches the config after a run).

## Design principle
Each run's `config.yaml` is the single source of truth: you author the intuitive primaries
(densities/n0, θ = kT/m_e c², lengths in d_e/d_i, speeds/c) and everything downstream is
derived from them. `make_inputs.py` (via `kinshock.deck`) **generates the WarpX input deck**
from the config, and `units.py` derives every physical / normalized scale (ω_pe, d_e, d_i,
ω_ci0, ρ_i0, C_s,ab, v_A, M_A, M_ms, …). The deck is a build artifact — never hand-edit it;
edit `config.yaml` and regenerate. All WarpX-deck details can therefore be ignored while
reasoning about a run.

## Modeling ablation: heating operator vs. laser deposition
The WarpX fork carries two ways to drive the piston. This project uses the first; the second
is available but unused, and the gap between them is worth understanding before reaching for it.

| | `ParticleHeater` (used here) | `LaserDeposition` (unused) |
|---|---|---|
| model | prescribed momentum diffusion in a fixed slab | geometric-optics ray trace + inverse-bremsstrahlung (IB) absorption |
| energy delivery | `du_i = sqrt(H dt) N(0,1)` | **the same Gaussian kick**, `sqrt((2/3) H dt) N(0,1)` |
| what sets `H` | `fac = 8 θ^1.5 / (√(m_i/m_e)·width/d_e)` — constant, dimensionless | `H = K·I/(n_e m_e)`, `K ~ n_e² T_e^{-3/2}/√(1−n_e/n_cr)` |

**Neither operator's kick constrains the timestep.** Both are drag-free Wiener processes, and
Wiener increments are self-similar under time-splitting — N kicks of variance `H dt` are exactly
distributed as one kick of variance `H·N dt`. That is why `heater.intervals: 20` is *exact*, not
an approximation. Note there is no `-(u−u_target)/τ` relaxation term: `theta` is a heating *rate*
parameter, not a thermostat setpoint; T_e,ab emerges from heating balanced against ablation loss
and `TargetInjector` replenishment.

**The laser's cost is in `H`, not in the laser.** The wavelength is never resolved (the ray march
is instantaneous on a frozen density snapshot; `ray_cfl` is an arc-length knob, a fraction of
`min_dx`). Instead:
1. `H` depends on the plasma state, so freezing it over `dt_dep` is O(dt) error — subcycling
   drops from exact to first-order.
2. T_e feedback with exponent −3/2: deposit → T_e rises → K falls. Needs `(3/2)ΔT_e/T_e ≪ 1`.
   The heater already adds ~2.5% of thermal energy per application, so this is marginal at the
   same cadence.
3. `1/√(1−n_e/n_cr)` puts the whole beam in a thin layer at critical; the heater spreads its
   power over `width = 40 d_e` ≈ 133 cells. Same power, ~100× the local ΔT_e/T_e.
4. The critical surface *moves*, adding a hydrodynamic CFL `dt_dep ≲ dz/v_front`.
5. `n_cr = ε₀m_eω²/e²` and `I₀` [W/m²] are **absolute** scales. The heater has none — which is
   exactly why `n0` is a free scale factor here (see the collisions gotcha in `CLAUDE.md`).

`LaserDeposition.cpp` enforces **no** dt limit; all of the above is on the user.

### Root cause: θ_e = 0.078 means T_e,ab ≈ 39.9 keV
We run the paper's dimensionless problem at the **real** speed of light, so M_A ≈ 14 has to come
from an inflated temperature rather than the paper's reduced c. That one number breaks the laser
three ways: IB absorption ∝ T^{−3/2} (mfp ≈ 4×10⁵ d_e ⇒ transparent), required I₀ ~10¹⁷ W/cm²
(~100× a real HED laser), and it is the same trap as the lnΛ knob, in a third place.

### Fix: emulate reduced c in config, with s ≡ √(39900/470) = 9.21
| primary | change |
|---|---|
| `theta_e_heat`, `theta_e_init`, `theta_i_init`, `theta_0` | ÷ s² = 84.8 |
| `B0_tesla` | ÷ s (3.2075e-3 → 3.483e-4 T) — v_A/c then falls 0.01 → 1.09e-3 on its own |
| `max_step` | × s (250k → 2.3M) |
| lengths in d_e, `dz_over_de`, `cfl`, `mass_ratio` | unchanged |

Every dimensionless target is invariant: M_A (v_p and v_A both ÷s), β (T ÷s², B² ÷s²), ρ_i0/d_e,
d_i0/d_e. The heater stays self-consistent — `fac` ∝ θ^1.5 ÷ s³, so its e-fold time θ/fac grows
by s in step with everything else. Only ω_ci0 ÷ s, hence s× more steps. This buys T_e,ab = 470 eV,
I₀ ≈ 1.3×10¹⁴ W/cm² (real HED range), and lnΛ ≈ 80 (≈7× physical) for 9× the compute.

### If `LaserDeposition` is ever wired in
Free today: pin `n0` to a real ablation density (step-count-neutral) and pick λ₀ so n_cr/n_e,ab
matches; the module never reads `mass_ratio`, so 100 is fine; start in `temperature_mode = fixed`
to kill the stiff feedback while calibrating. The top-hat target `nt*(abs(z)<slab)` must gain a
preplasma ramp — with n_target > n_cr the critical surface is a discontinuity and the near-critical
treatment (`L_eff = 1/drds` on a locally linear profile) degenerates.

Worth changing in C++, ranked: (1) accept `critical_density` directly — `m_n_cr` is welded to
`wavelength` at `LaserDeposition.cpp:448`, with no override; (2) a named `K_scale` instead of
overloading `coulomb_log` (already defaulting to a non-physical **2**) — the R1_coll lesson is that
hiding a knob behind a physical name costs you a silently-wrong criterion; (3) a `laser.target`
block mirroring `collisions.target`, inverted in `units.py` (Python-side only); (4) semi-implicit
T_e, only if `local` mode proves stiff.

## Mapping to the 2019 OMEGA experiment
Schaeffer et al., *"Direct observations of particle dynamics in magnetized collisionless shock
precursors in laser-produced plasmas,"* PRL **122**, 245001 (2019),
[doi:10.1103/PhysRevLett.122.245001](https://doi.org/10.1103/PhysRevLett.122.245001) — the
experiment the 2020 simulation paper is written to support.

**Quoted in the paper.** OMEGA. Ambient: one beam, 351 nm, 100 J, 1 ns on a planar CH target.
Piston (12 ns later, t₀): two drive beams, 351 nm, 350 J, 2 ns on a second CH target. Thomson
probe: 527 nm, 30–50 J, 2 ns, scattering volume 50×50×70 µm³, α ≈ 1.5. Field: By, peak **10 T**
near the piston target, ∝1/x, perpendicular to the expansion. Ambient: **n_e0 = 0.9±0.2×10¹⁸ cm⁻³,
T_e0 = 40±10 eV**. Observations at **x = 3–4 mm**, streaked 2 ns starting 3–4.5 ns after t₀.
Results: v_sh ≈ **750 km/s**, M_s ≈ **15**, n_e/n_e0 ≈ **10**, T_ex/T_e0 ≈ **10**, ramp τ_n ~ 200 ps.
Focal-spot size and on-target intensity are **not quoted**; 350 J / 2 ns over a nominal OMEGA
phase-plate spot gives ~10¹³–10¹⁴ W/cm².

**Derived here** (CH, A/Z = 13/7; B ≈ 3 T at the probe, since 10 T is at the target and falls ∝1/x
— this is the dominant uncertainty, and M_A ∝ 1/B, β ∝ 1/B²):

| quantity | value | note |
|---|---|---|
| d_e0 | 5.60 µm | |
| d_i0 | 327 µm | d_i0/d_e0 = 58.4 (real CH mass ratio) |
| λ_D0 | 0.050 µm | |
| v_A | 50.6 km/s | |
| **M_A** | **14.8** | self-consistent with the quoted M_s ≈ 15 — the reason B ≈ 3 T |
| β_e0 | 1.6 | vs 0.4 in R1_warm (same convention) |
| ω_ci0⁻¹ | 6.5 ns | |
| observed domain | **12 d_i0** | 4 mm |
| observed window | **t·ω_ci0 ≈ 0.5–1.0** | the entire experiment is inside one gyroperiod |

**R1_coll's ambient is already 10¹⁸ cm⁻³ — the experiment's density**, so its collisionality
question is directly testable against the paper's τ_pa/τ_s ≫ 1 claim.

### Cost to reach it
The 2019 experiment is a *small* box (12 d_i0) over a *short* time (≲1 gyroperiod) — far less than
the 2020 representative run (80 d_i0, 6.5 gyroperiods). Matching it is **cheaper**, not harder:

| | R1_warm | experiment-matched |
|---|---|---|
| domain | 7500 d_e (75 d_i0) | 1200 d_e (12 d_i0) |
| cells | 25 000 | ~4 000 |
| steps | 250 000 | ~67 000 (t·ω_ci0 = 1.5) |
| dz/λ_D | 6.7 | 3.3 (better) |
| relative cost | 1× | **~1/23×** |

Experiment-matched keeps `B0_tesla = 3.2075e-3` (v_A/c = 0.01) and `theta_e_heat = 0.078` (M_A ≈ 14–15 ✓) and only
raises `theta_0` 0.002 → 8.05e-3 to match β_e0 ≈ 1.6, which *improves* Debye resolution.

### Adding the laser: n_cr is a dimensionless parameter, not 351 nm
The instinct is that the laser forces a huge density contrast, because
n_cr(351 nm) = 9.05×10²¹ cm⁻³ is 10⁴× the experiment's ambient and solid CH is another 10× above
that. **Neither of those is a requirement of the physics.** IB absorption
(`K ∝ n_e² lnΛ T^{-3/2}/n_cr`) works at any density; what needs n_cr is the *critical surface* —
the turning point, the reflection, the self-regulating front — and that is a statement about
**n_target/n_cr**, not about absolute density.

- **Solid density is unnecessary.** The ray cannot see past n_cr; everything above it is optically
  inaccessible and matters only as a mass reservoir, which is `TargetInjector`'s job. Truncating
  the target at ~2 n_cr is worth ~5.6×.
- **351 nm is unnecessary.** λ₀ enters *only* through n_cr. R1_warm's existing contrast is already
  2.5 : 0.01 = **250×** — ample room to host a critical surface. Put n_cr at ≈1.25 n0 (half the
  target density, 125× the ambient): critical sits inside the target, the ambient is at
  n/n_cr = 8e-3 so refractive index 0.996 (no bending, no absorption), and dz is set by d_e at the
  target density, which the run already resolves.

| | 351 nm + solid CH | truncate at 2 n_cr | **n_cr ≈ 1.25 n0** |
|---|---|---|---|
| density contrast | 10⁵ | 2×10⁴ | **250 (already there)** |
| relative cost | 4×10³ … 5×10⁴× | ~8×10³× | **~1×** |

At R1_coll's pinned n0 = 10²⁰ cm⁻³ that is n_cr = 1.25×10²⁰ ⇒ **λ₀ ≈ 3 µm** (mid-IR, not exotic).
Two further easings in 1D: there is exactly **one ray**, so the ray-march cost the header warns
about (80 ms/application, dominant in 2D) is negligible; and with a *resolved* preplasma ramp of
~30 d_e (~100 cells) the deposition is no more localized than the heater's 133-cell slab, so the
T^{−3/2} stiffness is comparable and `intervals` in the 5–20 range should hold.

What is given up: λ₀ becomes a derived diagnostic, not 351 nm. Every dimensionless ratio governing
deposition geometry is preserved (n_cr/n_target, n_cr/n_ambient, ray optics, the moving critical
surface, the T^{−3/2} self-regulation). Absorption *magnitude* still needs tuning, but scaling n
and n_cr down together by f gives τ ∝ √f, so a 100× scale-down costs only ~10× in lnΛ — on top of
the ~80 from the reduced-c emulation, giving lnΛ ≈ 800, well inside R1_coll's precedent of 7713.

This is what makes **C++ change (1) — accepting `critical_density` directly — load-bearing rather
than a convenience.** Without it the deck must say λ₀ = 3 µm where it means "put n_cr at 1.25 n0",
which is exactly the naming trap R1_coll already paid for once.

If a *faithful* 351 nm / solid-CH calculation is ever wanted, the way to afford it is a small,
high-resolution laser run near the target used purely to *compute* T_e,ab and the ablation rate,
feeding `theta_e_heat` to an R1-class heater run (~23× R1_warm, one-time) — replacing the by-eye
recalibration history (0.092 → 0.062 → 0.078) with a first-principles number.

## Running
Requires a WarpX build with the `ParticleHeater` + `TargetInjector` modules, and a Python
env with `yt`, `numpy`, `matplotlib`, `pyyaml`.

```bash
# structure tests (no WarpX needed)
python tests/test_structures.py

# a run (example: R0 smoke)
python scripts/make_inputs.py runs/R0            # config.yaml -> runs/R0/inputs_kinshock_R0
cd runs/R0
MPICH_GPU_SUPPORT_ENABLED=0 mpirun -np <N> <warpx.1d> inputs_kinshock_R0

# verify what actually ran matches the config, then analyze
python scripts/make_inputs.py runs/R0 --verify   # warpx_used_inputs == config?
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
