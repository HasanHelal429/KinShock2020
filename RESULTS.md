# KinShock2020 — Results / validation log

Status log for the WarpX reproduction of Schaeffer et al. 2020 (Phys. Plasmas 27,
042901). See `OVERVIEW.md` (paper) and `REPLICATION_PLAN.md` (plan).

---

## 2026-07-22 — Scaffold, structure tests, and R0 smoke test

**Environment.** Prebuilt `warpx-cda/build/bin/warpx.1d` (built 2026-07-22, includes the
merged `ParticleHeater` + `TargetInjector` modules). Analysis via the `omegashock` conda
env (yt 4.4.2, numpy, matplotlib, PyYAML).

### Structure tests — 9/9 PASS (`tests/test_structures.py`, no WarpX needed)
- `config.load` + required sections; `units.derive` reproduces Table I from config
  primaries: **M_A = 13.95** (target 14), **M_ms = 12.74** (target 13), C_s,ab/c = 0.0303,
  v_sh/c = 0.1395, dt·ω_pe = 0.225, n_cell = 50000, d_i0/d_e = 100, ρ_i0/d_e ≈ 1040.
- `validate(R1)` clean; R0 shares every physics primary with R1.
- **`make_config.py` round-trip**: parsing both decks and rebuilding the config matches the
  authored `config.yaml` exactly (the "config = what was simulated" guarantee).
- metrics: Eq. 1 expansion speed, Eq. 2 RH ratio (→ ~3.9, strong-shock ~4×), reflected
  fraction G, front tracking, and the 7-criteria evaluator all correct.

### R0 smoke run — PASS (structure) (`runs/R0/`, 2000 cells, 2000 steps, ~exit 0)
Verified the deck + operators + full analysis pipeline:
- **No NaN / no recorded warnings**; all fields finite (checked Bx,By,Bz,Ex,Ey,Ez on the
  last plotfile). The "114 error matches" in a naive grep were all `INFO` lines (`inf`).
- **dt = 6.647e-12 s ⇒ dt·ω_pe = 0.375** — exactly as derived for dz = 0.5 d_e, CFL 0.75.
- **Heater active**: "particle_heater: generated heating rate for species" (foil preset).
- **Injector active**: replenishing the slab each interval (12 → ~80 particles/step).
- **B pile-up**: Bx reaches ~0.037 T (~11× B₀ = 3.2 mT) at the piston edge — the snowplow
  begins, as expected.
- **Analysis pipeline** runs end-to-end on real data and writes all artifacts:
  - `media/testing/`: `R0_config_summary.png`, `R0_loaded_state.png`, `R0_operator_balance.png`
  - `media/R0/`: `shock_streak.png`, `shock_trajectory.png`, `shock_lineouts.png`,
    `shock_phase.png`, `shock_reflected.png`, `criteria.json`, `shock_ni.mp4`, `shock_phase.mp4`
  - Criteria correctly report **no shock** at t·ω_ci0 ∈ [0, 0.08] (`is_shock: None`); the
    trajectory's "M_A = 31" is the meaningless early piston-edge transient (expected: R0 is
    a structure test, not a physics run). `validate` correctly warns R0's domain is too small
    for a real shock.
- Fixed one bug found here: numpy `bool_`/`float64` not JSON-serializable in
  `criteria.json` → added a native-cast `default=` in `make_figures.py`.

### OPEN PHYSICS ITEM — heater over-heats piston electrons (not yet resolved)
Measured piston-electron temperature (drift-subtracted, |z| < 20 d_e) vs the heater target
θ_e = 0.092:

| t·ω_pe | θ_thermal | ratio to θ_e |
|---|---|---|
| 0   | 0.001 | 0.01 |
| 150 | 0.073 | 0.79 |
| 375 | 0.144 | 1.56 |
| 675 | 0.190 | 2.07 (still rising, rate ~halving each ~250 ω_pe⁻¹) |

The overshoot is **genuine thermal** (bulk drift <u_z> ≈ 0.005 ≪ thermal). The ported
operator's fac calibration relaxes the slab toward θ_e *via free-streaming escape of hot
electrons from the slab* — and R0's ±500 d_e domain is almost certainly **too small for that
escape**, so heat accumulates. The validated flatfoil case (which saturated at θ vs PSC,
`warpx-cda/heating_operator/run_flatfoil_compare/`) expanded into a large low-density region.

**Hypothesis (testable):** R1's ±7500 d_e domain will let hot electrons stream away and the
slab will saturate near θ_e. **Next check before trusting R1's absolute piston/shock speed:**
1. Run a **Core-tier R1** (±3600 d_e, ~125k steps) and re-measure the piston-electron θ(t) —
   expect saturation closer to 0.092.
2. If it still overshoots, compare directly to the flatfoil validation deck and test
   `foil.n0 = n0` vs `= nt` (`foil.n0 = nt` gives ~2× larger H → worse; so n0 is the better
   choice, consistent with the PSC normalization-density convention already in the deck).
This matters because C_s,ab = √(θ_e/μ)·c sets the piston speed and hence M_A; a 2× θ error
→ ~1.4× speed error.

### RESOLVED — flatfoil `foil.n0` cross-check + dimensionality (supersedes the R0 hypothesis)
Ran while R1_core was in progress. Two things settled:

1. **`foil.n0 = n0` is the validated convention.** The PSC-cross-validated deck
   `warpx-cda/heating_operator/run_flatfoil_compare/inputs_flatfoil_compare` uses
   `particle_heater.foil.n0 = n0` (ablation reference) — exactly our R1/R0 choice. **The
   `foil.n0` concern is closed: n0 is correct.** (`= nt` would give ~2× larger H → worse.)

2. **The over-heating is dimensionality of the heating geometry, NOT domain size.** Measured
   heated-electron θ_thermal (drift-subtracted) at fixed validated params (θ_heat=0.04,
   slab 10 d_e, bg 0.002 n₀, `foil.n0=n0`):
   - **2D flatfoil (validated vs PSC), Gaussian spot:** saturates at **0.031 = 0.77 θ_heat** (stable).
   - **1D analog (`runs/xcheck_flatfoil_1d/`), same params, uniform slab:** saturates at
     **0.096 = 2.4 θ_heat** (stable by t·ω_pe≈1000).
   - R0 (1D, θ_heat=0.092) heads to the same ~2.4× (0.19 = 2.07× and rising at t·ω_pe=675).

   → In 2D the heater is a *localized spot*, so hot electrons escape **transversely** into cold
   plasma → cools below θ_heat. In 1D / transversely-uniform slabs (`spot_radius=0`) there is
   **no transverse escape** — only z-streaming out of the slab — so it saturates at ~2.4 θ_heat.
   The R0 "domain too small" hypothesis above was **wrong**: R0's ±500 d_e is ample; geometry is
   the driver. Figure: `media/testing/heater_dimensionality_xcheck.png`.

**Quasi-1D context (from the paper — confirms the calibration direction).** Schaeffer's PSC
runs are "gridded in the x–z plane ... quasi-1D, with only a few cells and **uniform driving
conditions in the transverse direction**" (30000×12 cells; uniform in x, NOT a localized spot).
So (a) our 1D uniform slab is the *faithful* geometry — it matches their transverse-uniform
driving, where electrons escape only along z; the 2D Gaussian-spot flatfoil was a different
(plume) test. And (b) Table I's T_e,ab = 0.092 is the *physical* plume temperature that sets
C_s,ab = 0.030 c (√θ/μ·c), so the calibration target is **θ_actual (measured, saturated) =
0.092**, achieved by lowering θ_heat. The acceptance test remains the measured **M_A ≈ 14**
(and v_p ≈ 0.104 c, v_sh ≈ 0.138 c), not the θ bookkeeping.
(Caveat: our port matched PSC to <1% in the 2D-spot validation; a uniform-1D PSC-vs-WarpX
heater comparison would be the check if bit-level 1D agreement is ever required.)

**Calibration implication for the 1D reproduction.** In our 1D (transversely-uniform, like
Schaeffer's quasi-1D uniform driving), θ_heat is a *rate* knob, not the saturation
temperature: T_e,actual ≈ 2.4 θ_heat. To hit Schaeffer's C_s,ab = 0.030 c (⇒ θ_actual ≈ 0.09,
since C_s,ab=√(θ/μ)c) and hence M_A ≈ 14, the piston-electron temperature — not the nominal
θ_heat — must equal ~0.09. With the ~2.4× factor that means **θ_heat ≈ 0.037–0.04**, not 0.092.
**R1_core is running at θ_heat = 0.092 and will therefore over-drive the piston** (θ_actual ~0.22
⇒ C_s ~1.55× ⇒ M_A ~ 20–22 instead of 14). R1_core still validates *shock formation* and the
analysis pipeline; the absolute M_A needs the recalibration below.

**Recommended before the Full R1:** a quick 1D θ-scan at the R1 slab width (|z|<20 d_e) —
2–3 short runs at θ_heat ∈ {0.03, 0.04, 0.05} — measuring the saturated piston-electron θ, to
pick the θ_heat giving θ_actual ≈ 0.09 (C_s,ab = 0.030 c). Then set that in `config.yaml`
(`plasma.piston.theta_e_heat`) for R1_core/R1. This keeps the reproduction faithful in the
speed/Mach that defines the shock.

### Verdict
- **Code structure: verified** (tests 9/9; deck runs clean; config round-trips; pipeline
  produces every planned artifact).
- **`foil.n0`: resolved** (= n0, matches PSC-validated deck).
- **Heater calibration: understood** — 1D uniform heating saturates at ~2.4 θ_heat (vs 0.77×
  for the 2D localized spot); recalibrate θ_heat so θ_actual ≈ 0.09 before the faithful M_A=14 run.
- **Physics calibration: one open item** (piston-electron heating saturation), with a clear,
  cheap next step (Core-tier R1) to close it before committing to the Full run.

### Next
- [ ] Core-tier R1 run → confirm piston-electron θ saturates near θ_e in a large domain
- [ ] then Full R1 (§7 runtime ~730 core-h) → acceptance criteria (plan §6)
- [ ] R2 (B₀=0), R3 (n_e0=0) negative controls
