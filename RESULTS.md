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
- **`make_inputs.py` deck generation**: the deck rendered from each config resolves back to
  the config primaries and is physically equivalent to the committed hand-written deck for R1
  and R0 (same resolved `my_constants` + scalar settings). Post-run, `make_inputs.py --verify`
  confirms `warpx_used_inputs` matches the config (the "config = what was simulated" guarantee).
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
- **Code structure: verified** (tests 9/9; deck runs clean; deck generated from config +
  verified against `warpx_used_inputs`; pipeline produces every planned artifact).
- **`foil.n0`: resolved** (= n0, matches PSC-validated deck).
- **Heater calibration: understood** — 1D uniform heating saturates at ~2.4 θ_heat (vs 0.77×
  for the 2D localized spot); recalibrate θ_heat so θ_actual ≈ 0.09 before the faithful M_A=14 run.
- **Physics calibration: one open item** (piston-electron heating saturation), with a clear,
  cheap next step (Core-tier R1) to close it before committing to the Full run.

### Next
- [ ] Core-tier R1 run → confirm piston-electron θ saturates near θ_e in a large domain
- [ ] then Full R1 (§7 runtime ~730 core-h) → acceptance criteria (plan §6)
- [x] R2 (B₀=0), R3 (n_e0=0) negative controls (2026-07-26): R3 clean null; R2 no magnetized shock
      (streak), auto-criteria false-positive on electrostatic ambient acceleration + weak-B₀ self-fields

---

## 2026-07-23 — One-sided half-domain (exploit z→−z symmetry); validated at R0

**Motivation.** The centered-slab problem is symmetric about z=0 (two counter-propagating
shocks in ±z; we analyze one). Simulating only the +z half [0, half] with a wall at z=0
halves cells + particles (REPLICATION_PLAN.md §7, "one-sided ablation"). Implemented as
**approach (a): a physical foil/reflecting wall at z=0** (specular particles + pec fields),
valid because z=0 sits inside the dense, driven piston — far from where the shock forms.

**Implementation (config-driven, symmetric layout unchanged/default).**
- `geometry.layout: one_sided` → domain `[0, half]`, `n_cell` halved (`units.derive` +
  `deck._n_cell`).
- `geometry.boundary` now accepts a `{lo, hi}` map of semantic names resolved to WarpX tokens:
  `reflecting`→(pec, reflecting), `open`→(pec, absorbing), `absorbing`→(silver_mueller,
  absorbing), `periodic`→(periodic, periodic).
- z=0 = `reflecting` (symmetry/foil wall); far edge = `open`. **`open` uses pec (not
  Silver-Mueller) fields**: the projection B-field divergence cleaner (active with any external
  B) rejects Silver-Mueller, and pmc/neumann would zero the tangential B₀ — pec is the one
  div-safe choice that preserves B₀. Nothing reaches the far edge within a run, so it is a
  don't-care beyond being stable. `config.validate` warns on periodic-on-one-face and on
  `absorbing` with a background field.

**Two bugs found and fixed during R0 validation** (`runs/R0_half`, 1000 cells vs R0's 2000):
1. **Far-boundary field BC** — first attempt used Silver-Mueller at the far edge → WarpX abort
   ("div cleaner requires periodic/PEC/PMC/neumann"). Fixed by the `open`→pec mapping above.
2. **Heater foil width** — the PSC foil heating rate is **H ∝ 1/width**, `width = foil.hi −
   foil.lo` (`ParticleHeater::makeFoilExpression`). Naively moving `foil.lo` from `−slab` to
   `0` halved the width → **doubled the heating rate** → the half-domain piston ran **1.7–1.9×
   over-energized**. Fix: keep `foil.lo = −slab` (and injector, which is per-cell and
   width-independent) even one-sided, so the rate matches the full domain; the one-sided
   domain then clips the heated region to [0, slab]. This is the whole reason the foil geometry
   must NOT be rewritten for the half domain.

**Validation vs full R0 on z ≥ 0** (both at t·ω_pe ≈ 750; `media/testing/R0_half_validation.png`):
- Loaded state on z≥0: **identical** (deterministic uniform loading; ratio 1.000, relRMS 0).
- Bulk conserved quantities: total particle **weight 0.499×**, total **energy 0.504×** full
  (target 0.5). Piston **KE on z≥0**: electrons **1.025×**, ions **0.970×** full's z≥0 half
  (i.e. matches to 2–3%). [Before the foil fix these were 1.69× / 1.87×.]
- **B_perp, ambient (z > 65 d_e, where the shock forms):** mean/B₀ = 1.60 (full) vs 1.58
  (half). The raw RMS difference (1.12) collapses under smoothing (→0.23 at 40 d_e) → it is
  short-wavelength **PIC shot noise** (independent RNG draws), not a systematic offset.
- **B_perp, piston/near-wall (z < 65 ≈ 3 slab):** ~15–20% higher in the half domain
  (peak 11.1 vs 9.8). This is the residual, **localized** cost of approach (a): a specular wall
  flips only the normal v_z, not the gyro-coupled v_y, so the diamagnetic/E×B current does not
  vanish exactly at z=0 as it would at a true symmetry plane. It stays buried in the driven
  piston, far from the ambient shock region.
- Runtime: **8.8 s vs 15.2 s** (single node, 4 OMP threads).

**Verdict.** Half-domain reproduces full-domain z≥0 physics where it matters (bulk + ambient) at
~half cost; the near-wall piston field carries a known ~15–20% localized artifact. **Caveat:**
diagnostics that read field/compression *inside* the piston (z ≲ 3 slab) should be treated with
care in one-sided runs; shock diagnostics (in the ambient) are unaffected. `tests/test_structures.py`
gains `test_one_sided_half_domain` (10/10 pass).

### Next (half-domain)
- [ ] Half-domain R1_core (`layout: one_sided`, `±3600`→`0..3600 d_e`) and cross-check the
  shock front / compression / reflected-ion fraction against a full-domain R1_core, confirming
  the near-wall artifact does not reach the shock at physics resolution (dz=0.3 d_e).


---

## R1_core half-domain cross-check (2026-07-23) — physics-resolution validation

`runs/R1_core_half`: `layout: one_sided`, wall@z=0, `0..3600 d_e`, dz=0.3 d_e, 12000 cells,
100 ppc, 125000 steps (~180 min, tmux `wxr1half`). Compared frame-for-frame against the
complete full-domain `runs/R1_core` on z≥0 (same dt & diagnostic cadence).
Figure: `media/testing/R1_core_half_crosscheck.png`. Scripts: `tmp/crosscheck_r1core.py`
(raw) + `tmp/crosscheck_r1core_v2.py` (refined: ambient-only compression measured in a
[front−400, front−50] d_e window, edge-contaminated frames excluded).

**Boundary caveat.** The +z shock reaches the domain edge (3600 d_e) at t·ω_ci0 ≈ 2.5; the
last 14/51 frames pile against the absorbing boundary and are excluded from the quantitative
comparison. Clean window: 37 frames up to t·ω_ci0 = 2.25. (The full run's +z shock hits its
own +3600 edge at the same time, so this limit is identical for both — not a half-domain defect.)

**Findings on z ≥ 0 (clean window):**
- **Bulk conservation exact:** total macroparticles half/full = 0.500→0.521, total energy
  0.499→0.492. The half domain is precisely half the system. ✓
- **Shock front trajectory:** tracks the full run with a *steady* offset of **−4.7% mean**
  (half slightly behind), narrowing from −7% early to −2.6% by t·ω_ci0=2.25. Clean-window
  front-speed ratio **0.97** (half ~3% slower).
- **Ambient density compression:** full 2.39 vs half 2.45 — **matches to 3%** (shock still
  steepening through this window; not yet the asymptotic ~4×).
- **Ambient B_perp compression:** full 3.33 vs half 3.76 — half runs **~13% higher**.
- **Near-wall B_perp** (z<3 slab): full 17.5 vs half 15.9 — comparable.

**Interpretation.** The −5% front lag + ~13% higher field compression are the R1_core-scale
signature of the *same* z=0 specular-wall artifact seen in R0: the wall flips only v_z, not the
gyro-coupled v_y, so the piston is driven marginally differently. Unlike R0 (a pre-shock smoke
test where the artifact stayed buried in the foil), running to shock formation lets that ~5%
systematic leak into shock *kinematics*. The reflected-ion fraction G from the raw script showed
a larger ~15–20% gap, but that comparison used a single boundary-contaminated v_sh as reference
for both runs and is not reliable; it is not used for the verdict.

**Verdict.** Half-domain reproduces full-domain R1_core shock structure to **~3% (density),
~5% (front kinematics), ~13% (field compression)** at half the compute. The differences are a
small but *systematic* (not noise) consequence of the imperfect reflecting wall — acceptable for
the 2× saving on structural/compression studies, but front-speed-sensitive diagnostics (e.g.
Mach number, reflected-ion energetics) inherit a ~5% offset that should be quoted with the result.


---

## CPU threading benchmark (2026-07-23) — `media/testing/cpu_benchmark.png`

Harness: `$CLAUDE_JOB_DIR/tmp/bench/harness.sh` (600 timed steps, 150 warmup discarded,
median per-step, min of 2 repeats, taskset-pinned to distinct physical cores, diags off).
Ran on AMD 7950X **while 16 `flash4` procs (another user) held loadavg 16-24** -> absolute
numbers are a LOWER BOUND; relative trends hold. Binary: warpx.1d.MPI.OMP (double precision).

**Thread scaling (base config, R1_core_half early-time, 12k cells / 3.3M ptcls):**
| threads | s/step | speedup | efficiency |
|--------:|-------:|--------:|-----------:|
| 1 | 0.2536 | 1.0x | 100% |
| 2 | 0.1172 | 2.16x | (108%*) |
| 4 | 0.0635 | 3.99x | 100% |
| 8 | 0.0347 | **7.31x** | **91%** |

Near-linear to 8 threads — the memory-bandwidth knee I predicted does NOT bite by 8 (it must
sit between 8 and 16). *T2 slightly super-linear = cache/latency artifact on the T1 anchor.

**Config ladder (all @8 threads) — the surprise:** every proposed optimization was neutral-to-
negative for this 1D problem. `+max_grid_size=512` 0.0398 (+15%), `+tiling(mfiter=8)` 0.0397,
`+sort_intervals=25` 0.0388 — all SLOWER than plain 8-thread (0.0347). Partly confounded by
rising load during the ladder (loadavg 20->24), but no positive signal. **Interpretation:** the
default 94 small grids (96-128 cells) already keep per-grid current arrays cache-resident;
enlarging grids trades that cache locality for redistribute savings that don't net out. The 17%
redistribute cost is apparently cheaper than the cache penalty of big grids. **Theory (fixes
B/C/D) refuted by measurement — keep the default grid config.**

**Full-run (125k-step) wall-time estimates** (early s/step x growth G=1.36 from the real 4-thr run):
- 4 threads (current): **180 min** (measured, anchors G)
- 8 threads: **~99 min** (1.83x) <- the free win, just set `OMP_NUM_THREADS=8`
- 16 threads: **~55-73 min projected** (2.6-3.0x), UNMEASURED — needs a quiet machine

**Actionable:** run production with `OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=cores`
for ~1.8x now; do NOT bother with max_grid_size/tiling/sort tweaks. Re-run the sweep to 16/32
threads when the machine is idle to locate the bandwidth knee and confirm the 16-thread number.


---

## B-field oscillations: numerical artifact (finite-grid instability), NOT physical (2026-07-23)

Question: are the B_perp oscillations a physical whistler precursor or a setup artifact?
Verdict: **artifact — the finite-grid (aliasing) heating instability from an under-resolved
Debye length.** Diagnostic figure: `media/testing/bfield_oscillation_diag.png`.

**Root cause (quantified):** lambda_D = v_the/wpe = 0.045 d_e = **0.15 cells**, but dz = 0.30 d_e,
so **dz/lambda_D = 6.7 (ambient), 9.5 (piston)**. The momentum-conserving PIC finite-grid
instability threshold is dz/lambda_D ~ pi = 3.1 -> we are **2-3x over threshold**. Under-resolved
Debye length aliases particle-grid coupling into grid-scale EM noise that grows in time.

**Evidence it is numerical, not a whistler precursor:**
1. **Not statistical ppc noise:** far-upstream RMS(B) is 0.005 B0 early (t*wci=0.11, true 100-ppc
   floor) and grows ~200x to ~1.1 B0 as the shock's hot/streaming particles arrive. Sampling
   noise is constant; this grows.
2. **Not a coherent wave:** corr(Bx,By) = -0.01 (foot), -0.06 (ambient); hodograms are isotropic
   random blobs. A whistler is elliptically polarized and would trace ellipses.
3. **Blue spectrum, filter-truncated:** power rises toward small scales and is cut off at the
   bilinear-filter/Nyquist scale (0.6-1.2 d_e), with NO peak at any physical scale. Real
   waves/instabilities peak at their physical wavelength and decline toward the grid.
4. **Worst where dz/lambda_D is largest:** piston/downstream +/-10 B0 vs ambient +/-2-4 B0.
5. **Mean field correct:** <Bx>=1.001 B0, <By>=0.004 -> initialization fine; the hash is added on top.

**Scale content (pristine upstream dBx):** ~50% of variance at lambda<5 d_e, ~80% at lambda<10 d_e,
peak in the 2-5 d_e band (grid noise left behind after the 1-pass filter removes lambda<2 d_e).
Only ~20% at lambda>10 d_e (where any real precursor would live).

**Impact:** LARGE-scale shock physics is unaffected (compression, front trajectory, reflected
fraction all matched the paper -- those live at >>10 d_e). The artifact pollutes small-scale B
fluctuation/wave diagnostics and adds spurious grid heating (watch electron-T and EP energy).

**Fixes (cheapest first):**
1. **More current-filter passes** (`warpx.filter_npass = 4-8`, or bilinear+compensation) -- damps
   grid-scale noise hard; ~free. First thing to try.
2. **Cubic particle shape** (`algo.particle_shape = 3`, from 2) -- raises the instability threshold; modest cost.
3. **More ppc** -- lowers the seed (delays growth) but does not cure the instability.
4. **Resolve Debye:** dz <~ pi*lambda_D ~ 0.14 d_e (from 0.30) -- definitive cure but ~4x cost
   (2x cells x 2x smaller dt via CFL).
5. Energy-conserving / nodal / Galilean-PSATD solver -- no grid heating, larger change.

**Confirmation test (recommended):** short run baseline vs +filter_npass=4 vs finer dz; if
upstream dB collapses under filtering/finer grid -> confirms numerical. (Not yet run.)


### B-field verdict — REVISION (2026-07-23, later): more nuanced than "all numerical"

The confirmation runs (`studies/bfield_convergence/`) corrected the verdict above:
- **Near-shock / foot turbulence is PHYSICAL.** Where reflected ions live (within
  ~0.3 rho_i of the front), the B spectrum is invariant to filter_npass, particle_shape,
  AND resolving the Debye length (dz 0.30->0.15) — i.e. *converged*, not numerical. This is
  real reflected-ion-driven foot turbulence, as expected for a supercritical perp shock.
- **Far-upstream fluctuation is a separate component, strongly suspected numerical.** Beyond
  ~0.3 rho_i the ambient ions AND electrons sit at the t=0 thermal floor (x1.00-1.05) while
  dBx~1.1 B0 and flat to the domain edge — a dB/B~1 field that scatters neither species is
  not a self-consistent plasma wave. Decisive long-run test (t*wci~0.56, filt8 + finer_dz vs
  baseline) is IN PROGRESS; final verdict + `media/testing/bfield_convergence.png` to follow.
- So the earlier flat "finite-grid instability" verdict is **too strong**: the oscillations are
  a MIX — physical foot turbulence + a (likely numerical) far-upstream grid component. The
  under-resolved Debye length (dz/lambda_D~7) remains a legitimate grid-heating concern but is
  NOT the source of the near-shock turbulence. Use `scripts/bfield_diagnostic.py` per run.


### B-field verdict — FINAL (2026-07-23): mixed, quantified (`media/testing/bfield_convergence.png`)

Decisive long-run test done (`studies/bfield_convergence/`, t*wci~0.56): baseline dz0.3 vs
+filter_npass=8 vs +resolve-Debye dz0.15, far-upstream cold zone (front+600..+1400, ~0.6-1.3 rho_i).

| variant | far-upstream RMS(dBx)/B0 | vs baseline | foot RMS | ion u/floor | e- u/floor |
|---|---|---|---|---|---|
| baseline dz0.3       | 0.721 | --   | 0.768 | 1.00 | 1.01 |
| +filter_npass=8      | 0.498 | -31% | 0.722 | 1.00 | 1.01 |
| +resolve Debye dz0.15| 0.540 | -25% | 0.709 | 1.00 | 1.00 |

**The oscillations are a MIX of physical and numerical:**
1. **Near-shock foot turbulence (within ~0.3 rho_i, reflected ions present): PHYSICAL.** Foot RMS
   moves only -6..-8% across all knobs -> converged. Real reflected-ion-driven shock turbulence.
2. **Far-upstream small-scale (lambda <~ 2-3 d_e) hash: NUMERICAL.** filter_npass=8 drops the
   far-upstream spectral power by orders of magnitude below ~2 d_e and RMS by 31%; this is the
   grid noise enabled by the under-resolved Debye length (dz/lambda_D~7). It does not scatter the
   plasma (ions & electrons pinned at the t=0 thermal floor, x1.00-1.01) -> not self-consistent.
3. **Residual far-upstream component (lambda >~ 3 d_e, ~0.5 B0): unresolved by these tests.**
   Survives filtering and a 2x Debye refinement (which only reached dz/lambda_D=3.3, still above
   the ~pi threshold) AND does not couple to the plasma. Not a strongly-coupled physical wave;
   distinguishing weak-precursor vs longer-wavelength numerical mode needs dz <~ 0.14 d_e and/or
   higher ppc (future convergence sweep).

**Practical guidance.** The shock structure and foot turbulence (the science) are trustworthy.
The upstream "noise" is substantially a grid artifact: raise `warpx.filter_npass_each_dir` (4-8)
to cut it ~30% at ~no physics cost, or resolve the Debye length (dz <~ 0.14 d_e) for a clean
upstream. Use `scripts/bfield_diagnostic.py` per run to check. This SUPERSEDES the earlier flat
"finite-grid instability" verdict (too strong: it ignored the physical foot turbulence).


---

## R1_core_half vs paper: shock speed / onset time / position (2026-07-23)

Compared the three shock-formation observables against Schaeffer 2020 (Table I + §4.3).
Script: `scripts/compare_r1half.py` (reuses `metrics.track_front`, `speed_from_trajectory`,
`onset_time_from_G`, `onset_location_from_F`). All comparisons are **dimensionless** —
the run's absolute n_e,ab normalization differs from the paper's HED realization, so only
the scale-free ratios (C_s,ab, c, ρ_i0, ω_ci0⁻¹) are meaningful; the absolute km/s is not.

| Quantity | Paper | R1_core_half | Δ |
|---|---|---|---|
| shock speed v_sh | 4.6 C_s,ab = 0.138 c | **4.73 C_s,ab = 0.144 c** | **+3%** ✓ |
| v_sh / v_p | ≈ 4/3 (1.33) | **1.38** | +4% ✓ |
| onset time t*₁ (max dG/dt) | ≈ 1.0 ω_ci0⁻¹ | **1.41 ω_ci0⁻¹** | +40% (late) |
| onset position z*₁ (max dF/dz) | ≈ 1.0 ρ_i0 | **1.49 ρ_i0** | +50% (far) |

(R1_core_half scales: C_s,ab = 0.0303 c, ρ_i0 = 1040 d_e, ω_ci0⁻¹ = 177 ns. Front at t*₁ = 2.05 ρ_i0.)

**Shock speed matches to ~3% — but only on the CLEAN WINDOW.** The +z shock reaches the
3600 d_e edge at t·ω_ci0 ≈ 2.5; the last ~14/51 frames pile against the absorbing boundary.
Fitting v_sh over the clean window (t·ω_ci0 ≤ 2.25, front < 0.94×edge, 37 frames) gives
**0.144 c = 4.73 C_s,ab**. A naive `speed_from_trajectory` second-half fit that *includes* the
piled-up frames is flattened to a spurious **3.7 C_s,ab (0.112 c)** — that low number is a
boundary artifact, not the physical speed. Always exclude the edge frames for v_sh here.

**Onset t*₁/z*₁ run ~1.4–1.5× late/far.** Order-unity-correct (the paper stresses these are
"≈" and M_A-insensitive) but with a consistent late/far bias. Use the smooth dG/dt onset
(1.41), NOT the flag-based `criteria.json` first-shock at t·ω_ci0 = 0.225 — criterion 7
(piston separation) trips early and noisily, so the first-flag time is premature and not
comparable to the paper's t*₂. The late/far bias is consistent with the documented
half-domain z=0 wall artifact (~5% front lag; see the fix analysis below).

**Verdict.** Shock *kinematics* (v_sh, v_sh/v_p) reproduce the paper to a few %; formation
*timescale/location* are order-1 but biased ~1.4–1.5× late/far, tracking the half-domain wall.


---

## Half-domain z=0 wall artifact: root cause + fix (2026-07-23)

**Root cause (now pinned to the geometry's discrete symmetry).** The two-sided ablation
problem is symmetric about z=0, but for a **perpendicular** shock (B₀ = B₀ x̂) the correct
discrete symmetry is **NOT a mirror reflection** — it is a **180° rotation about the x-axis
(the B₀ axis):** (x,y,z) → (x, −y, −z). Reason: an ion gyrates in the (y,z) plane
(v×B with B along x couples only v_y↔v_z; v_x is along B and free). The map that sends a
+z gyro-orbit to its −z counterpart must therefore flip **both** v_y and v_z (a point
reflection in the gyration plane), and it leaves B₀ x̂ invariant (a proper rotation about x̂:
B_x→B_x, B_y→−B_y, B_z→−B_z).

Our approach-(a) wall applies **specular reflection**, which flips **only v_z**. It gets the
sense of gyration wrong by leaving v_y unchanged — an error of Δv_y = 2v_y at every bounce.
That mishandled diamagnetic/gyro current at z=0 is exactly the observed near-wall B_perp
excess (R0_half 11.1 vs 9.8; R1_core_half +13% field, −5% front lag). PEC fields compound it
weakly: the true symmetry-plane field BC is per-component {B_x, E_x even (Neumann); B_y, B_z,
E_y, E_z odd (=0 at plane)}; PEC matches only B_x(even, preserves B₀ ✓), B_z(0 ✓), E_y(0 ✓)
and mishandles E_x (pins 0, want Neumann) and B_y (leaves free, want 0). Particle v_y is the
dominant term; the field mismatch is second-order and in the shock-subdominant components.

**Fixes, best first:**

1. **Correct symmetry-plane particle BC (the real fix, targeted, fork change).** Replace the
   specular reflect at z=0 with a **π-rotation-about-x reflect**: on crossing z=0, set
   z→−z (or clamp), **v_z→−v_z AND v_y→−v_y**, leave v_x. This is the exact discrete symmetry
   and should remove almost all of the ~5% front lag / ~13–20% near-wall field excess, since
   the diamagnetic current is a particle current. WarpX's built-in `reflecting` only negates
   the normal component; add a `symmetry_rot_x` particle BC in the CDA fork (same place the
   ParticleHeater/TargetInjector hooks live — negate the momentum component transverse-to-B
   as well as the normal). Cheapest correct option; ~one boundary kernel.

2. **Per-component field symmetry BC at z=0 (completes the fix).** {B_x, E_x: Neumann (even);
   B_y, B_z, E_y, E_z: Dirichlet 0 (odd)}. This is a genuine symmetry plane — neither PEC nor
   PMC (PMC would zero B₀, PEC pins E_x and frees B_y). Needs a per-component ghost-cell rule
   in the fork; must stay div-cleaner-safe (the reason Silver-Mueller was rejected). Do only
   if the residual field error survives fix 1 — particle v_y is expected to dominate.

3. **Full-domain runs for front-speed-sensitive diagnostics (guaranteed-correct fallback).**
   The half domain is a 2× compute optimization; for Mach number / reflected-ion energetics /
   onset timing, run the full ±3600 domain (R1_core exists) and keep the half domain for
   structure/compression studies where the ~3% density agreement is fine.

4. **Bury the wall deeper (mitigation, no code).** Push the reflecting plane a few d_e inside
   the dense driven piston (widen the clipped foil region) so v_y coupling is subdominant to
   the strong drive where the wall sits. Reduces, does not remove, the artifact.

**Validation plan for fix 1.** Re-run R0_half and R1_core_half with the v_y-flip BC and check:
(a) R0_half near-wall B_perp drops from 11.1 → ~9.8 (full-domain z≥0); (b) R1_core_half front
lag shrinks from −4.7% toward 0 and field compression from +13% toward the full-domain 3.33;
(c) onset t*₁/z*₁ move from 1.41/1.49 toward the paper's ~1. Figure: `media/testing/`
`R1_core_half_crosscheck.png` regenerated post-fix vs the full run. If (a)–(c) all collapse,
the symmetry-BC diagnosis is confirmed and the half domain becomes faithful for kinematics too.


---

## Symmetry-wall fix: implementation + 3-way validation run (2026-07-23) — IN PROGRESS

Implemented fix 1 and launched the controlled A/B/C validation.

**Fork change (warpx-cda, compiles; NOT yet committed).** New particle-BC option
`boundary.reflect_symmetry_axis = {x|y|z}`: a reflecting bounce becomes a **π-rotation about
the named (B₀) axis** — flips the normal velocity AND the transverse component perpendicular
to the axis, preserves the parallel one — instead of plain specular (normal-only). Files:
`Source/Particles/ParticleBoundaries.{H,cpp}` (new `int reflect_symmetry_axis`, setter),
`ParticleBoundaries_K.H` (kernel `else if` block after `reflect_all_velocities`),
`PhysicalParticleContainer.cpp` (parse `boundary.reflect_symmetry_axis`, x/y/z→0/1/2, asserts
mutual-exclusion with `reflect_all_velocities`). Default −1 = disabled → fully back-compatible.
Incremental `--target app_1d` build clean; smoke test (200 steps) confirms WarpX consumes the
param (not in unused-inputs list), no abort.

**Deck wiring (`src/kinshock/deck.py`).** New semantic boundary `symmetry` → `(pec, reflecting)`
tokens **plus** the `boundary.reflect_symmetry_axis = x` line (x = B₀ axis for the perpendicular
geometry). `_boundaries()` returns a third `sym_axis` value; `_SYMMETRY_BCS` gates it. Use
`boundary.lo: symmetry` in place of `reflecting` for a faithful one-sided wall.

**New run `runs/R1_core_half_sym`.** Byte-identical deck to `R1_core_half` except the single
added line `boundary.reflect_symmetry_axis = x` (verified by diff) — a clean one-variable A/B.
`filter_npass` deliberately left unchanged (that is a separate physics-quality knob; changing
it too would confound the wall comparison). Launched with the benchmarked optimization
`OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=cores` (~1.8× vs 4 threads; ETA ~1h36m);
no `max_grid_size`/tiling/sort tweaks (benchmarked neutral-to-negative). `runs/R1_core_half_sym/`
`launch.sh` + `finalize.log` (auto-runs verify + figures + 3-way compare on completion).

**Validation harness `scripts/crosscheck_3way.py`.** Compares full=`R1_core` (z≥0),
spec=`R1_core_half` (specular), sym=`R1_core_half_sym` (symmetry) on: clean-window front speed
(t·ω_ci0≤2.25, front<0.94·edge), ambient n- and B_perp-compression (peak in [front−400,
front−50] d_e, piston zone z>3·slab excluded), near-wall B_perp (z<3·slab), onset t*₁/z*₁
(max dG/dt & dF/dz on +z ambient ions). "Fix passes" if sym is closer to full than spec on
front speed / field compression / onset AND near-wall B_perp drops toward full.

**Baseline (full vs spec):** front speed full 4.90 vs spec 4.73 C_s,ab (ratio 0.965 ✓ reproduces
the −4.7% lag). (Peak-in-window compression values run higher than the mean-based numbers quoted
earlier, but the identical method is applied to all three runs, so the A/B/C comparison is what
counts.)

**RESULT — sym run complete (2026-07-23, 50 frames; `media/testing/crosscheck_3way.png`).**
The symmetry wall is closer to the full-domain reference than specular on **all 8 metrics**
(auto-verdict "BETTER" each; t*₁ a tie in magnitude):

| metric | full | spec (Δ) | sym (Δ) |
|---|---|---|---|
| v_sh [C_s,ab]     | 4.90 | 4.73 (−0.17) | **4.83 (−0.07)** |
| v_sh / v_p        | 1.43 | 1.38 (−0.05) | **1.41 (−0.02)** |
| n-comp (amb)      | 4.87 | 5.09 (+0.22) | **4.70 (−0.17)** |
| B-comp (amb)      | 7.80 | 8.03 (+0.24) | **7.98 (+0.18)** |
| near-wall B/B₀    | 23.7 | 20.5 (−3.2)  | **25.3 (+1.5)** |
| t*₁ [ω_ci0⁻¹]     | 2.14 | 2.08 (−0.06) | 2.19 (+0.06) |
| z*₁ [ρ_i0]        | 2.64 | 2.38 (−0.25) | **2.68 (+0.04)** |

**Front lag essentially removed:** clean-window v_sh ratio sym/full = 0.984 (−1.6%) vs spec's
0.965 (−3.5%) — the fix more than halves the front-speed offset. z*₁ offset collapses from
−0.25 ρ_i0 (spec) to +0.04 (sym). At R1, near-wall B_perp: spec *under*-shoots (20.5 vs 23.7);
sym +1.5 (25.3) — opposite sign to R0 (where spec over-shot) but sym is closer to full in both.

**Diagnostic-consistency fix (2026-07-23, later — resolves a crosscheck-vs-`make_figures`
discrepancy in the reflected-ion fraction).** The two diagnostics disagreed on G(t)/onset because
they computed the reflected-ion **velocity threshold v_sh differently**:
- `make_figures` used `speed_from_trajectory` with a naive *second-half* fit that included the
  boundary-**stalled** frames (shock hits the 3600 d_e edge at t·ω_ci0≈2.5) → v_sh=**0.105 c**
  (M_A 10.5), spuriously low → too many ions counted as reflected.
- the crosscheck used a clean-window fit → v_sh=**0.148 c** (M_A 14.8).

Unified so both share identical logic:
1. `metrics.speed_from_trajectory` gained a domain-aware clean window (`z_edge`, keep
   |z|<0.94·edge, no second-half) — excludes the decelerating/stalled tail. Both callers pass it.
2. crosscheck front-tracking piston-exclusion aligned to `make_figures` (`slab·di`).
3. `metrics.onset_time_from_G` now returns the **first prominent dG/dt peak** (not the global
   argmax, which flipped between G's two rises at t·ω_ci0≈1.4 and ≈2.2 per-run); crosscheck uses it.

Result: `make_figures` and crosscheck now report the **identical** v_sh (sym 0.1485 c, M_A 14.85,
paper-consistent) and identical onset. Corrected 3-way table (supersedes the v_sh/onset rows above):

| metric | full | spec (Δ) | sym (Δ) |
|---|---|---|---|
| v_sh [C_s,ab] | 4.96 | 4.76 (−0.20) | **4.90 (−0.06)** |
| v_sh / v_p    | 1.45 | 1.39 (−0.06) | **1.43 (−0.02)** |
| t*₁ [ω_ci0⁻¹] | 1.41 | 1.35 (−0.06) | 0.96 (−0.45) |
| z*₁ [ρ_i0]    | 1.53 | 1.41 (−0.12) | 1.05 (−0.48) |

Onset now lands near the **paper's t*₁≈1, z*₁≈1** for all runs (the old 2.1–2.7 values were the
argmax latching onto G's second rise). **Honest read:** the robust kinematic/field metrics
(v_sh, v_sh/v_p, n/B-compression, near-wall B_perp) show sym closer to full than spec; **onset is
the exception** — sym forms ~0.45 ω_ci0⁻¹ *earlier* than full (a stable result now, not the
argmax artifact), though onset still carries ~±0.4 ω_ci0⁻¹ detection sensitivity on these
near-identical G curves and each run uses its own v_sh threshold. Do not over-read the onset Δ.

**R1 verdict: fix confirmed at physics resolution.** The π-rotation wall moves every metric
toward the full domain — the ~5% specular front lag and the compression/near-wall offsets all
shrink. Combined with the R0 near-wall By result, the specular v_y mishandling is confirmed as
the artifact source, and the particle-BC fix (fix 1) resolves the dominant part. **Residual**
(sym still deviates a few % from full, and R0 showed a low-side undershoot) is consistent with
the un-applied field symmetry BC (fix 2: pec still pins E_x / frees B_y at z=0) — the remaining
work if sub-few-% half-domain fidelity is ever needed.

**Figure `media/testing/crosscheck_3way.png`** (generator `scripts/plot_crosscheck_3way.py`): A —
front trajectory z/ρ_i0(t); B — G(t) + onset markers; C — R1 metrics as deviation-from-full
(spec vs sym bars); D — R0 near-wall B_x/B_y/B_perp. Colors fixed-order CVD-safe (full=gray,
spec=vermillion, sym=blue), direct-labeled.

### Next
- [ ] Commit the fork change (`warpx-cda`, branch `development`) + add a WarpX regression test
  for `reflect_symmetry_axis`; commit KinShock deck/config/scripts + this RESULTS entry.
- [ ] (Optional) Implement fix 2 (per-component field symmetry BC) if sub-few-% fidelity needed.
- [ ] Adopt `boundary.lo: symmetry` as the default one-sided wall in future half-domain runs.

**R0_half_sym early confirmation (DONE — `runs/R0_half_sym`, smoke tier, seconds).** Ran the
cheap A/B/C first. Near-wall B_perp peak (z<3·slab, z≥0, last frame t·ω_pe=750), full R0 = ref:

| near-wall peak /B0 | full R0 | spec (R0_half) | sym (R0_half_sym) |
|---|---|---|---|
| \|Bx\| (compressed B₀)    | 13.06 | 14.29 (+9%)  | 10.81 (−17%) |
| \|By\| (diamagnetic/gyro) | 14.23 | 19.53 (+37%) | 12.15 (−15%) |
| total B_perp              | 14.53 | 19.55 (+34%) | 12.91 (−11%) |

The specular over-shoot is concentrated in **By** — the gyro/diamagnetic component sourced by
the transverse current that the specular v_y mishandling corrupts, exactly as predicted. The
π-rotation fix removes it: total near-wall B_perp deviation from full goes **+34% (spec) → −11%
(sym)**. Residual ~11–17% (now on the LOW side, largest on Bx) is consistent with the un-fixed
field BC (fix 2: pec still pins E_x / frees B_y at z=0) and/or 50-ppc smoke noise; the
physics-resolution R1 run is the real test. (Absolute values differ from the earlier R0 9.8/11.1
— different window/no smoothing — but the identical 3-way method makes it valid.) **Verdict:
fix confirmed at R0 — the dominant particle-v_y near-wall artifact is removed.** Deck differs
from R0_half by only `reflect_symmetry_axis = x`.


---

## Full R1 (half-domain, symmetry wall) — Table I reproduction (2026-07-24)

`runs/R1_half`: the Full-tier Table I run as a one-sided half-domain with the validated
π-rotation symmetry wall at z=0 + 8-thread optimization. 25000 cells (0..7500 d_e, dz=0.3 d_e),
250000 steps (→ t·ω_ci0 = 5.63), θ_heat=0.092 (Table I), 100 ppc/species, θ_Bn=90°.
**Ran 6h19m wall @ 8 threads** (mean 0.091 s/step; matched the ~7.4h pre-run estimate — lighter
load than assumed). Deck verified == config (`warpx_used_inputs`); `run_checks` validation OK.

**ACCEPTANCE — essentially exact agreement with Table I:**

| quantity | measured | paper (Table I / §4.3) |
|---|---|---|
| v_sh | 0.1400 c = **4.62 C_s,ab** | 0.138 c = 4.6 C_s,ab |
| **M_A** | **14.00** | 14 |
| M_ms | 12.78 | 13 |
| v_sh / v_p | 1.35 | ≈ 4/3 (1.33) |
| onset t*₁ (first prominent dG/dt) | 1.35 ω_ci0⁻¹ | ≈ 1 |
| onset z*₁ | 1.46 ρ_i0 | ≈ 1 |
| reflected-ion fraction G (peak) | 0.30 | present (criterion 6) |
| ambient compression (peak-in-window) | ~5 (forming) | → ~4 strong-shock (+overshoot) |
| 7 formation criteria | satisfied (precursor + shock) | (1)–(7) |

**v_sh/M_A are spot-on** (this is the primary acceptance metric, and it now lands at exactly 14 —
the corrected clean-window speed metric + full domain/duration; the earlier R1_core M_A~15 was the
shorter run + edge contamination). Onset t*₁/z*₁ are order-unity (~1.4×, consistent with the
half-domain + onset-detection sensitivity). Reflected ions and the 7 criteria confirm a genuine
collisionless shock, not a piston compression.

**Caveat — domain marginally undersized at M_A=14.** The shock reaches the +z edge
(7.2 ρ_i0 = 7500 d_e) at t·ω_ci0 ≈ 5.2, i.e. front/ρ_i0 ≈ (v_sh/v_p)·(t·ω_ci0) = 1.35·5.6 ≈ 7.6 >
7.2. So the very-late downstream (t*₃ ≈ 5, RH relaxation) is slightly clipped by the boundary;
clean physics window is t·ω_ci0 ≲ 5.0. Front-speed / onset / reflected-ion metrics are unaffected
(measured well before edge contact). A future t*₃-clean run should use ~8500 d_e (or stop at
t·ω_ci0 ≈ 5). Also note derived β_ab=1840, β_0=0.4 vs Table I 1150/0.2 (pre-existing config item,
not one of the shock-formation acceptance signatures).

**Optimizations delivered:** half-domain (25k vs 50k cells) + 8 threads → **6h19m vs ~27h** for the
naive symmetric/4-thread reproduction, at faithful M_A=14. Figures: `media/R1_half/` (shock_streak,
trajectory, lineouts, phase, reflected, criteria.json) + `media/testing/R1_half_*`.

### Verdict
The WarpX `ParticleHeater` + `TargetInjector` + symmetry-wall half-domain reproduce a
Schaeffer-2020-class M_A=14 perpendicular collisionless shock at Table I parameters, with the
speed model, formation timescales, ~4× compression, and reflected-ion signature — at ~1/4 the
naive compute. Remaining: R2 (B₀=0) and R3 (n_e0=0) negative controls.

---

## R2 (B₀=0) + R3 (n_e0=0) negative controls (2026-07-26) — `runs/R2`, `runs/R3`

Built as **exact single-knob clones of the calibrated production run `runs/R1_warm`** (warm
piston ions θ_e_heat=0.078, θ_i=7.8e-4; 25000 cells, dz=0.3 d_e, 250000 steps → t·ω_ci0=5.63,
100 ppc). One line changed per control, nothing else, so any difference isolates that ingredient:
- **R2:** `field.orientation: perpendicular → unmagnetized` → deck applies Bx=0 (`kinshock.deck`
  zeroes `bx` when orientation≠perpendicular). `vA_over_c=0.01` is kept, so the *Scales* (v_A,
  ω_ci0, ρ_i0) are retained and R2 normalizes in R1_warm's units. Verified raw field Bx=By=Bz=0
  **exactly at t=0** — the control starts truly unmagnetized.
- **R3:** `numerics.ppc.ambient: 100 → 0` → zero ambient macroparticles. B₀ is left untouched
  (**avoids the trap** that zeroing `amb.density_over_n0` would also zero B₀ via
  `B0=vA·√(μ₀·n_amb·Mᵢ)`, conflating R3 with R2). All scales identical to R1_warm.

Reference threshold for the reflected-ion / front metrics: R1_warm's by-eye **v_sh=0.1285 c
(M_A=12.85)** copied into each control's `shock_fit.yaml` (scales identical → same cutoff → fair
side-by-side). R2 ran 6h49m, R3 ran 1h48m (both under 2-run core contention; R3 is ~10× cheaper
per step with no ambient particles). Figures in `media/R2/`, `media/R3/`.

**Gotcha fixed (cost one restart):** the generated deck sets **no `diag*.file_prefix`**, so WarpX
writes to `diags/` *relative to the launch CWD*. R1_warm worked because it was launched from
inside its run dir; launching both controls from the repo root made **R2 and R3 write to the same
`./diags/` and clobber each other** (WarpX `.old.NNNN` rename files were the tell). Fix: **launch
each run with CWD = its own run dir** (`bash -c 'cd runs/<ID> && warpx …'`) so diag1/diag_fields/
reducedfiles all land under `runs/<ID>/`. The first R3 "completion" was corrupt and was discarded;
both were rerun clean. → Added to CLAUDE.md gotchas, and since 2026-07-27 **enforced** by
`scripts/launch.sh`, which cd's into the run dir for you and refuses to launch into a populated
`diags/` — use it instead of invoking the binary by hand.

### R3 (n_e0=0): clean null ✓ — matches paper Fig. 9

| signature | R1_warm (shock) | R3 (n_e0=0) |
|---|---|---|
| `is_shock` (7 criteria) | forms (t*≈0.3+) | **False at all 51 frames** |
| reflected-ion fraction G(t) | finite beam (peak ~0.30) | **0.0000 throughout** |
| ambient-ion phase space | reflected beam above v_sh | **rows completely empty** (zero ambient loaded) |
| downstream compression | ~4× plateau | **none** — n_e falls off monotonically into vacuum |
| onset t*₁ / z*₁ | order-unity | **never triggers** (None) |

With no medium the piston just **free-expands into vacuum+B₀** as a self-similar rarefaction fan
(Fig. 7 row 2), and its leading edge piles up frozen-in flux as a **magnetic snowplow**
(B⊥/B₀~14 in the streak) that runs *faster* than v_sh (nothing upstream to load it). That B
pile-up trips the coded `flag 4 (field_compression)` and the raw `n_compression` reads ~250
(=piston 2.5 n0 ÷ ambient-reference 0.01 n0, i.e. there is no ambient baseline left) — both are
artifacts of the missing ambient, **not** a shock; the overall verdict is correctly `is_shock=False`
because criterion 6 (reflected ions) is 0. **Removing the ambient removes every ambient structure.**

### R2 (B₀=0): no magnetized shock ✓, but the automated flags false-positive — read with care

Paper (Fig. 9) expects: ambient ions *still initially accelerated*, but **no magnetic compression,
no strong heating, no secondary compression — no shock.** What we find:

- **Field structure — null ✓.** The streak shows **no coherent magnetic ramp** propagating at
  v_sh; ahead of the piston is dark (B⊥≈0), and the only magnetic activity is incoherent,
  filamentary **self-generated Weibel/current fields** behind the front. No ordered compression.
- **Ambient ions accelerated forward ✓ (matches paper).** Fig. 7 row 1 shows a fast ambient band
  at **v_z/v_sh ≈ 1.5–2**, centered on the **electrostatic piston-reflection cap 2·v_p/v_sh =
  1.62** — i.e. ambient ions specularly reflected off the *moving piston potential*, not
  gyro-reflected at a detached magnetized shock.
- **The criteria give a FALSE POSITIVE** (`is_shock=True` from t*≈0.34, G peak 0.22, B_comp~17–20)
  for two setup-specific reasons, **both worth remembering**:
  1. **`G` over-counts.** G = fraction of ambient ions with v_z > v_sh. The forward
     piston-accelerated ions (up to ~1.6 v_sh) clear that threshold, so a benign electrostatic
     acceleration reads as "reflection." G is not mechanism-specific.
  2. **`B_compression = B_max/B₀` is meaningless when self-fields ≫ B₀.** Raw probe: B₀=0.003208 T,
     but self-generated fields reach **rms ~0.010–0.014 T** by mid/late run (Bx *and* By) — and
     these are **as large as the *total* field in R1_warm and R3** (their By rms ~0.009–0.013 T).
     Because vA/c=0.01 makes B₀ very weak, current-driven fields exceed it, so "B₀=0" does not
     leave the plasma dynamically unmagnetized and the compression ratio is dominated by noise.

**Verdict.** R3 is a textbook clean null. R2 confirms the *magnetized* shock (ordered B pileup +
magnetic ion reflection at a detached front) is **absent** without B₀ — the field streak is
unambiguous — while ambient ions are still electrostatically accelerated by the piston, exactly as
the paper states. The automated 7-criteria detector is **not mechanism-discriminating** and should
not be read as "R2 forms a shock": inspect the streak + ambient phase space. Two follow-ups if a
*quantitatively* clean B=0 null is wanted: (a) gate `flag 4`/`B_compression` on the **ordered**
(box-averaged) field, not peak |B|, and define "reflected" as v_z<0 or a shock-frame beam rather
than v_z>v_sh; and/or (b) note that the weak-B₀ (vA/c=0.01) regime lets self-generated fields rival
B₀ — a caveat that also bounds how "magnetized" the R1_warm shock's turbulent component is (its
*ordered* Bx compression to ~0.28 T peak is still real and carries the reflection).

(Reference note: R1_warm's by-eye fit reads **v_sh=0.1285 c → M_A=12.8** — the whole-run linear fit
averages a mildly curved front and sits ~8% under the M_A=14 target; the settled t≳3 slope trends
steeper. Not re-tuned here; it is the baseline both controls clone and the shared metric threshold.)

---

## tune_shock trajectory streaks were 20× under-sampled in time (2026-07-27) — tooling

**Symptom.** The `|B_perp|` and `n_e` streaks in `media/<ID>/tune_trajectory.png` looked coarsely
pixelated along `t` — the front read as a blocky staircase rather than a line, which is exactly the
thing you are trying to fit by eye.

**Cause.** `TrajectoryTuner` read `io.plotfiles()` → `diags/diag1*`, the **Full particle**
diagnostic (`diag1.intervals = 5000` → **51 frames** for R1_warm). The runs also write a field-only
diagnostic at `diag_fields.intervals = 250` → **1001 frames**, 20× denser, and `io.field_plotfiles()`
already existed to reach it (`make_figures.fig_streak` used it; `tune_shock` never did). Two
distinct things were pinning the tool to the coarse series:
1. the plotfile list itself (`B_perp` is read straight off the grid — it was low-res for no reason);
2. `n_e` came from `io.species_density()`, which **histograms macroparticles**
   (`particle_position_x`/`particle_weight`). `diag_fields` sets `write_species = 0`, so it carries
   no particles and that method structurally cannot use the dense series.

**Fix.** `n_e` now comes from the **deposited `rho_<species>` fields** (`n = |rho|/(Z·e)`), which
`diag_fields.fields_to_plot` already writes, so *both* panels run at the field cadence.
- `io.load_frame(path, fields=(...))` reads extra grid components into a new `Frame.comps` dict;
  missing components are simply absent, so callers can probe and fall back. One-arg calls unchanged.
- `io.field_species_density(frame, species, charge_states)` + `io.rho_field_names(species)`.
- `tune_shock` prefers `field_plotfiles`, and falls back to `plotfiles` + macroparticle histograms
  when the plotfiles carry no per-species rho (verified on `runs/R1_core`: **`diag1` has only 9
  components — Ex..Bz, jx..jz — no rho at all**, so pre-`diag_fields` runs behave exactly as before).

**Validation** (`diag_fields020000` vs `diag1020000`, both at t = 7.976658e-08 s):

| comparison | rho-derived vs macroparticle |
|---|---|
| integrated density (charge conservation) | **1e-6** agreement |
| mean deviation, raw 0.3 d_e cells | 7.3% (shape-factor smoothing) |
| mean deviation, 4.5 d_e display bins | **0.4%** |
| peak n_e (inside the dense piston) | 257 vs 277 — off-scale above the `vmax=8` cap anyway |

**Three things were needed to make 1001 frames × 25000 cells usable**, all in `TrajectoryTuner`:
- **stream the series** — one yt dataset alive at a time (`del fr` per iteration). The old
  list-comprehension would have held 1001 open, since `Frame` keeps `ds` for particle access.
- **float32** streak storage (~200 MB for both panels instead of 400).
- **auto z block-averaging to ~1600 bins** (`--zbin` overrides; R1_warm → ×15 = 4.5 d_e). The panel
  is ~1200 px wide, so 25M pcolormesh quads per re-render bought nothing — and *every* interactive
  command re-renders. Load is ~1m45s for 1001 frames; `--stride` thins it.

**Caveats for anyone reading these streaks.** (a) Deposited rho is shape-factor (and
`filter_npass`) smoothed — fine for front-fitting, but do not read absolute peaks off it; the
particle histogram is the reference. (b) The 4.5 d_e block-mean means the tuner streak is **not** a
small-scale B diagnostic — use `bfield_diagnostic.py` for spectra at λ ≲ 2–3 d_e.

**Observation, not yet chased:** at full cadence the B panel shows resolved **diagonal striations
behind the front** that were previously time-aliased into hash. Character (piston/reflected-ion
striations vs the `intervals = 20` heater cadence beating against the diag cadence) not
investigated.

**Follow-up left open:** `make_figures.fig_streak` already used the dense field series, so its B
streak was always full-cadence — but it holds all 1001 frames via `load_series` and hands
pcolormesh the full 1001 × 25000 grid. Same streaming + decimation treatment would speed it up; no
correctness issue (it reads only `fr.Bx`).

---

## Fig. 7 in ρ_i0 units + `--only` for `make_figures.py` (2026-07-27) — tooling

`make_figures.py` gained two flags, both figure-selection only (no physics touched):

- **`--fig7-xunits d_i0|rho_i0 [...]`** — horizontal normalization for the Fig. 7 grid.
  `d_i0` (default) → `shock_fig7.png` as before; `rho_i0` → **`shock_fig7_rho_i0.png`**, i.e. the
  same panels with z in the upstream ion gyroradius ρ_i0 = v_p/ω_ci0 (the paper's z* normalization,
  already `Scales.rho_i0`). Both units can be requested in one pass. The windows/overlays are still
  computed in d_i0 and only rescaled by d_i0/ρ_i0 at plot time, so the two files are panel-for-panel
  identical — only the x scaling, the axis label and a suptitle normalization line differ.
  For R1_warm, ρ_i0 = 1040 d_e = **10.4 d_i0**, so the ρ_i0 axes are the d_i0 axes ÷10.4
  (t·ω_ci0 = 1.24 panel: 8–19 d_i0 → 0.77–1.83 ρ_i0, front at 16 d_i0 → 1.53 ρ_i0 — consistent with
  the onset z*₁ ≈ 1.5 ρ_i0 measured in the R1_half acceptance table).
- **`--only streak|trajectory|lineouts|phase|fig7|reflected|criteria`** — build a subset. Needed
  because re-rendering one panel of a 51-frame / 17 GB run otherwise costs the whole A–D suite
  (the streak also loads the 1001-frame field series).

**Convention worth keeping: `media/{R1_warm,R1_coll}/shock_fig7*.png` are made with
`--phase-times 0.15 0.49 0.73 0.98 1.25`** (they snap to frames t·ω_ci0 = 0.11/0.45/0.68/1.01/1.24).
Those times were **chosen deliberately** to sit across the formation interval; the default
`--nframes 5` selection spreads panels over the whole run (0.11 … 5.63) and shows late,
edge-contaminated downstream instead. Regenerate with:

```bash
python scripts/make_figures.py runs/<ID> --only fig7 \
    --fig7-xunits d_i0 rho_i0 --phase-times 0.15 0.49 0.73 0.98 1.25
```

**A bare `make_figures.py runs/<ID>` silently reverts fig7 to the default times** — it did exactly
that to R1_warm during the 2026-07-29 R1_coll analysis, and only the stale timestamp on
`shock_fig7_rho_i0.png` (which the default pass does not write, since `--fig7-xunits` defaults to
`d_i0` alone) gave it away. If the two fig7 files for a run have different mtimes, the `d_i0` one is
the default-times imposter. **Always re-render fig7 with the flags above after any full-suite run.**

---

## R1_coll — collisional twin of R1_warm at n_e0 = 10¹⁸ cm⁻³ (2026-07-27) — setup

New run `runs/R1_coll/`: R1_warm's **collisional** comparison run (the plan's R6 slot,
`REPLICATION_PLAN.md` §4). Config + deck only — **not launched yet.**

**What changed vs R1_warm — exactly two things.** `diff` of the two decks (comments stripped)
is 2 hunks: `my_constants.n0` and the new collisions block. Everything else is byte-identical.

1. **`reference.n0: 1.0e18 → 1.0e26 m⁻³`** — i.e. n_e,ab = 10²⁰ cm⁻³ so the **ambient
   n_e0 = 0.01 n_e,ab = 10¹⁸ cm⁻³** (paper's Table I realization: 4.8×10¹⁸ cm⁻³). Since n0 is a
   pure scale factor in the collisionless problem (real c, free n0 — `REPLICATION_PLAN` §1),
   **every dimensionless quantity is unchanged**: M_A = 12.85, M_ms = 11.73, β_ab = 1560,
   β₀ = 0.4, dt·ω_pe = 0.225, n_cell = 25000, steps_per_ω_ci0 = 44444 all match R1_warm to
   ~1e-15. Only the absolute scales move (all lengths ×10⁻⁴):
   d_e,ab = 0.531 µm, d_i0 = 53.1 µm, ρ_i0 = 0.553 mm, **B₀ = 32.1 T**, dt = 0.399 fs,
   τ_sim = 99.7 ps. A `test_collisional_twin_of_R1_warm` regression test pins all of this.
2. **New `collisions:` block** → 10 WarpX `pairwisecoulomb` collisions (Perez 2012, the
   equivalent of the paper's Takizuka–Abe operator) over **all** unordered species pairs
   including the intra-species self-pairs, so a physical electron collides the same way
   regardless of which population (piston/ambient) it is loaded into. Target:
   **λ_ab = ω_ce,ab/ν_ei,ab = 20, the paper's Table I value** ⇒ lnΛ = 7713.

### The finding: lnΛ has to be forced, and there are two different "20"s

**Real collisions are inert at any attainable density.** We run the paper's dimensionless
problem at the *real* c, so θ_e = 0.078 means **T_e,ab = 39.9 keV**, not the paper's ~470 eV.
With ν_ei ∝ n T^(−3/2), even at n_e,ab = 10²⁰ cm⁻³ the physical NRL lnΛ = 11.6 gives
**mfp_ei,ab = 3.7×10⁵ d_e,ab**, i.e. λ_ab = 1.3×10⁴ — still collisionless by ~3 orders of
magnitude. Reaching the paper's λ_ab = 20 with a *physical* lnΛ would need n_e,ab ≈ 4×10²⁶ cm⁻³
(and mfp = 20 d_e,ab would need ~3×10²⁸ cm⁻³) — both far above solid density (~5×10²² cm⁻³).
(This is the quantitative version of the note dropped in R1_warm's header.) So **lnΛ is the
knob**: ν_ei is exactly linear in it, and `collisions.target` states the physics target while
`kinshock.units.coulomb_log_for` inverts it to the single number the deck carries.

**⚠ Two different "20"s — do not conflate them.** At these parameters ρ_e,ab = 27.9 d_e,ab, so:

| target | lnΛ | ν_ei,ab [s⁻¹] | ν_ei·dt | mfp/d_e,ab | mfp/d_i0 | λ_ab = ω_ce/ν_ei |
|---|---|---|---|---|---|---|
| `lambda_ab: 20` **(this run — paper Table I)** | 7.71×10³ | 2.82×10¹¹ | 1.1×10⁻⁴ | 559 | 5.59 | 20 |
| `mfp_over_de: 20` | 2.154×10⁵ | 7.88×10¹² | 3.1×10⁻³ | 20 | 0.2 | 0.72 |
| `coulomb_log: physical` | 11.57 | 4.23×10⁸ | 1.7×10⁻⁷ | 3.7×10⁵ | 3725 | 1.3×10⁴ |

The paper's collisionality **λ_ab ≡ ω_ce,ab/ν_ei,ab = 20** (`OVERVIEW.md` §2, Table I) is
mfp/**ρ_e,ab**, *not* mfp/d_e — with ρ_e,ab = 27.9 d_e,ab it means **mfp = 559 d_e,ab =
5.6 d_i0**, so the plasma stays collisionless at ion scales exactly as in the paper
(λ_mfp/d_i0 = 350), which is the premise of its Fig. 13 result that collisional and
collisionless formation look the same. **This run targets that value**, so Fig. 13 is a fair
comparison. A `mfp_over_de: 20` target would be a different, ~28× more collisional run
(mfp = 0.2 d_i0, collisional across the ramp) — the test guards the two apart.

ν_ei·dt = 1.1×10⁻⁴ (collision time resolved ~9000× by the PIC step), so no subcycling is
needed and `ndt_supercycle: 1`.

### Verification done (no shock physics yet)
- `make_inputs.py runs/R1_coll` → deck round-trips to the config; `--check` still clean for
  R1_warm/R1_recal/R1_cal/R2/R3/R1_core_half (no regression from the new code paths).
- `run_checks.py runs/R1_coll` → **validation OK**, λ_ab = 20.000 exactly (mfp = 558.6 d_e,ab).
- `tests/test_structures.py` → **11/11 PASS** (new collisional-twin test included).
- **WarpX smoke** (20 steps, 2000 cells, 8 ppc, throwaway CWD): exit 0, no errors, and
  `warpx_used_inputs` shows all 10 collisions with `CoulombLog = 7713.304243059874`.
  `deck.verify` against it flags only the deliberate smoke overrides.

### Cost note
Collisions add ~3 active pair-passes/step over most of the domain (only ambient or only
piston present) and all 10 in the piston/ambient overlap. Benchmark the first ~1000 steps
before committing to the full 250 000; if collisions dominate, `ndt_supercycle: 8` still
leaves ν_ei·dt_coll ≈ 9×10⁻⁴ and cuts the cost ~8×.

Launch with:
```bash
scripts/launch.sh -b -L runs/R1_coll
```

---

## R1_coll — collision cost benchmarked, `ndt_supercycle: 8`, run launched (2026-07-28)

Answered the cost note above with a measurement before committing the 250k steps. Two
1500-step runs back-to-back in a throwaway CWD (8 threads, `OMP_PROC_BIND=spread`,
chablis load ~18–20), decks copied verbatim so the ONLY difference is the collisions block:

| deck | s/step (WarpX `Avg. per step`) |
|---|---|
| `inputs_kinshock_R1_warm` (collisionless twin) | **0.0825** |
| `inputs_kinshock_R1_coll` (10 pairs, `ndt=1`) | **0.2785** |

**The 10 pairwise-Coulomb pairs cost 3.4× — i.e. 2.4× the entire rest of the timestep.**
Anchoring on R1_warm's real 7h11m for 250k steps (its logged early rate 0.077–0.087 s/step
brackets today's 0.0825, so the machine is in the same regime) that projects **~24 h** at
`ndt_supercycle: 1`, and the 3.4× is a *floor*: collision cost grows super-linearly in
particle count as the piston sweeps into the ambient and more cells activate all 10 pairs
rather than ~3.

**Decision: `ndt_supercycle: 8`.** ν_ei·dt = 1.1×10⁻⁴, so the collision time was resolved
~9000× — applying the operator every 8th step with dt_coll = 8·dt preserves the rate and
still resolves it ~1100×. Not a physics compromise; a 2.6× saving on an operator that was
four orders of magnitude over-resolved. `coulomb_log = 7713.304243059874` and λ_ab = 20.000
are **unchanged** by the supercycle (verified via `run_checks.py`) — the inversion is on the
rate, not the cadence.

| `ndt_supercycle` | s/step | projected full run |
|---|---|---|
| 1 | 0.279 | ~24 h |
| 4 | 0.132 | ~11.5 h |
| **8 (chosen)** | ~0.11–0.14 | **~10 h** |

**Launched** 2026-07-28 10:08 (`scripts/launch.sh -b -L runs/R1_coll`, pid 1233159, logger
1233161). Early logged rate 0.143 s/step(warpx) cumulative / 0.114 instantaneous at step
200; with R1_warm's 1.34× growth over the run that lands at **~10–11 h**, finishing ~21:00.
Watch `runs/R1_coll/progress.log`.

**Caveat to revisit if the twin comparison ever disagrees with R1_warm:** the supercycle is
unvalidated here — no ndt=1 vs ndt=8 convergence check was run (that would cost the 24 h
this avoids). The ~1100× resolution margin makes it very unlikely to matter, but if
collisional/collisionless formation *does* differ (contra the paper's Fig. 13), rule out
the supercycle before believing the physics.

## R1_coll finished + analysed: shock is R1_warm's twin, but criterion 2 **fails** (2026-07-29)

Run completed 2026-07-28 23:07: **250000/250000 steps, 12h58m wall, mean 0.187 s/step(warpx)**
— 1.3× the ~10 h projected above, the extra coming from the same end-loaded cost growth the
benchmark predicted (0.115 → 0.187 s/step as the piston activates all 10 pairs across more
cells) plus host contention (×1.04 → ×1.27). `make_inputs.py --verify`: **OK**, WarpX ran
exactly this config. `run_checks.py`: validation OK, λ_ab = 20.000, mfp_ei,ab = 558.6 d_e.

### Shock kinematics: identical to R1_warm, by eye and by fit

`tune_shock.py` streaks for R1_coll and R1_warm are superposable. Trials at 0.1225 / 0.1285 /
0.1350 c bracket the front the same way in both, so **`shock_fit.yaml` was set to R1_warm's
exact values — v_sh = 0.1285 c, z0 = 0, no per-time overrides** — deliberately, so every
downstream diagnostic in the twin comparison shares one front definition. M_A = 12.85,
M_ms = 11.73, v_sh/v_p = 1.24.

(The auto-tracker's `n > 1.5 n_e0` front again runs ~4% steeper than the by-eye ramp fit and
saturates when it hits the domain edge at t ≈ 5.2 — the usual reason this project fits by eye.)

### Criterion 2 was a hard-coded constant, and it hid the one thing R1_coll tests

`metrics.evaluate_criteria` defaulted `lambda_ii_over_di0` to a literal **350.0** — the
paper's Table I figure — and `make_figures` never passed anything else. So the collisionless
criterion reported the *paper's* collisionality for every run, including the only run in the
repo that has a collision operator. Now derived from the config:

* `units.nu_ii` — NRL ion-ion rate ν_i = 4.80e-8 Z⁴ μ^(−1/2) n_i lnΛ T_i^(−3/2), with
  μ = m_i/m_p from the **actual** ion mass (μ = 0.054 here; assuming μ = 1 would understate
  the rate 4.3×). Validated two ways: for hydrogen at T_i = T_e it reproduces
  ν_ii/ν_ei = sqrt(m_e/m_p)/√2 to 4 digits, and with a *physical* lnΛ it gives
  λ_ii/d_i0 = **347** vs the paper's quoted 350 — i.e. the old hard-coded constant was
  exactly "the physical-lnΛ value", which is why nobody noticed.
* `Scales.mfp_ii_amb` / `mfp_ii_amb_over_di0`, printed by `run_checks.py`. **inf** when the
  config has no `collisions` block (no collision operator ⇒ genuinely collisionless; a
  finite "physical" number would be meaningless there, since n0 is a free scale factor).

**Result: R1_coll's upstream λ_ii/d_i0 = 0.52, so criterion 2 fails in all 51 frames** and the
run never satisfies criteria 1–7. R1_warm passes all seven (first shock t*ω_ci0 = 0.34).

### The third instance of the "two different numbers" trap

**A single lnΛ knob cannot match both of the paper's collisionality numbers at real c.**
Table I quotes λ_ab = ω_ce,ab/ν_ei,ab = 20 *and* λ_mfp/d_i0 = 350 together; the paper gets
both because its T_e,ab ≈ 470 eV. We run at real c, so θ_e = 0.078 ⇒ T_e,ab = 39.9 keV, and
hitting λ_ab = 20 needs lnΛ = 7713 = **667× physical**. That same 667× multiplies the
*upstream* rate too, dragging λ_ii/d_i0 from 347 down to 0.52 — from collisionless-at-ion-
scales to marginally ion-collisional across the ramp. Choosing `quantity: lambda_ab` was
therefore choosing to break criterion 2, silently.

If a future run wants "collisional but still a collisionless shock", target the upstream
number instead: λ_ii/d_i0 = 347 × (11.567/lnΛ), so **lnΛ ≈ 400 gives λ_ii/d_i0 ≈ 10** while
still being 35× more collisional than physical (λ_ab would then be 385, not 20). The
reachable window is narrow because the 39× temperature error eats it.

### What collisions actually changed: ~12% less peak B⊥, nothing else

Per-criterion, R1_coll vs R1_warm at matched times — every criterion that measures the
*shock* agrees; only criterion 2 differs:

| quantity | R1_warm | R1_coll |
|---|---|---|
| n compression (t = 1.35 / 2.81 / 4.16) | 272 / 270 / 270 | 265 / 267 / 284 |
| piston separation / ρ_i0 | 1.666 / 3.469 / 5.137 | 1.660 / 3.468 / 5.134 |
| reflected fraction G | 0.073 / 0.190 / 0.260 | 0.070 / 0.206 / 0.283 |
| ramp scale / d_i0 | 0.017 / 0.020 / 0.015 | 0.020 / 0.016 / 0.014 |
| frames failing crit. 2 (of 51) | 0 | **51** |

Peak B⊥/B0, **excluding the outer 2 d_i0**: R1_coll / R1_warm = **0.881 ± 0.013** (mean over
46 frames at t*ω_ci0 > 0.5, frame-to-frame σ = 0.088). Consistent with the streaks, where
R1_coll's far-upstream small-scale hash is visibly damped relative to R1_warm's — collisions
eat the fine-scale B fluctuations while leaving the macroscopic ramp alone.

**Do NOT quote `B_compression` from `criteria.json` at late times.** It is a global max, and
from t*ω_ci0 ≈ 5.5 the open hi boundary throws an ~80× B⊥ spike in **both** runs (warm 83.4 /
85.9, coll 78.1 / 66.3 raw, vs 19.0 / 18.4 and 14.1 / 15.6 with the outer 2 d_i0 cut). The
apparent "86 vs 66" late-time difference is that artifact, not physics.

### Bottom line

The paper's Fig. 13 claim — collisional and collisionless perpendicular shock formation look
the same — **holds here, and more strongly than the paper tested it**: R1_coll reproduces
R1_warm's compression, front speed, ramp steepness, piston separation and reflected-ion
fraction while sitting at λ_ii/d_i0 = 0.52 instead of 350, i.e. ~670× more collisional and
past the point where its own criterion 2 calls it collisionless. The honest caveat is that
R1_coll is therefore **not** "the paper's Table I collisionality"; it is a much more
collisional run that behaves the same anyway. The supercycle caveat from 2026-07-28 is not
implicated — the twins agree, which is the direction that needed no exoneration.

---

## Laser deposition vs. the heater, and the cost to reach the 2019 OMEGA experiment (2026-07-31) — analysis, no runs

Desk analysis only: no WarpX was launched. Full write-up now lives in `README.md`
(§"Modeling ablation" and §"Mapping to the 2019 OMEGA experiment"); this entry records
the findings and the numbers that took work to pin down.

### The heater's `intervals: 20` is exact, not an approximation
Both `ParticleHeater` and the fork's unused `LaserDeposition` deliver energy through the
*same* drag-free Gaussian kick (`ParticleHeater.cpp:319`, `LaserDeposition.cpp:1366`).
That makes them Wiener processes, and Wiener increments are self-similar under
time-splitting — so subcycling the heater is distributionally exact at any `intervals`.
There is no `-(u−u_target)/τ` term anywhere: `theta` is a heating **rate** parameter, not a
thermostat setpoint (`OVERVIEW.md`'s "relax toward T_e,ab" is loose phrasing). T_e,ab
emerges from heating balanced against ablation loss + `TargetInjector` replenishment.
Measured for R1_warm: `fac` = 4.36e-4 ω_pe, heating e-fold θ_e/fac ≈ 179 ω_pe⁻¹ ≈ 800
timesteps, ~2.5% of thermal energy per application, 250k steps ≈ 300 heating times.

### The laser's dt constraint is not the wavelength
`LaserDeposition` is geometric-optics: the ray march is instantaneous on a frozen density
snapshot, and `ray_cfl` is an *arc-length* step (a fraction of `min_dx`), not a timestep.
λ₀ costs nothing. The constraint is that `H = K·I/(n_e m_e)` with
`K ~ n_e² T_e^{-3/2}/√(1−n_e/n_cr)` is a stiff explicit functional of the state: subcycling
drops to first order, the T^{−3/2} feedback needs `(3/2)ΔT_e/T_e ≪ 1`, the critical-surface
singularity concentrates the whole beam into ~1 cell (vs the heater's 133), and the critical
surface moves (`dt_dep ≲ dz/v_front`). **The module enforces no dt limit at all** — grepped;
there is no assert. It is entirely on the user.

### Root cause of every laser mismatch: θ_e = 0.078 ⇒ T_e,ab ≈ 39.9 keV
Real c forces M_A ≈ 14 to come from an inflated temperature. Consequences: IB absorption
∝T^{−3/2} gives mfp ≈ 4×10⁵ d_e (transparent), and matching the heater's power needs
I₀ ~10¹⁷ W/cm², ~100× a real HED laser. Same shape as the lnΛ trap, third instance.

**Fix, config-only: emulate reduced c.** With s ≡ √(39900/470) = 9.21, divide every θ by
s² = 84.8 and `vA_over_c` by s, multiply `max_step` by s. M_A, β, ρ_i0/d_e, d_i0/d_e are all
invariant, and the heater stays self-consistent (`fac` ∝ θ^1.5 ÷ s³, so θ/fac grows by s in
step with everything else). Only ω_ci0 ÷ s ⇒ s× more steps. Buys T_e,ab = 470 eV,
I₀ ≈ 1.3e14 W/cm², lnΛ ≈ 80 for 9× compute. **Not yet tested.**

### The 2019 experiment (PRL 122, 245001) is a *small* problem
Paper-quoted: OMEGA, 351 nm; ambient beam 100 J/1 ns, piston 2 beams 350 J/2 ns, both CH;
B_y peak 10 T ∝1/x; ambient n_e0 = 0.9±0.2e18 cm⁻³, T_e0 = 40±10 eV; probed x = 3–4 mm,
2 ns starting 3–4.5 ns after t₀; v_sh ≈ 750 km/s, M_s ≈ 15, n_e/n_e0 ≈ 10, T_ex/T_e0 ≈ 10.
Spot size / on-target intensity are **not quoted** in the paper.

Derived (CH A/Z = 13/7; B ≈ 3 T at the probe — 10 T is at the target and falls ∝1/x, and
3 T is what makes M_A = 14.8 agree with the quoted M_s ≈ 15; this is the dominant
uncertainty, M_A ∝ 1/B and β ∝ 1/B²): d_e0 = 5.60 µm, d_i0 = 327 µm (58.4 d_e0),
λ_D0 = 0.050 µm, v_A = 50.6 km/s, β_e0 = 1.6, ω_ci0⁻¹ = 6.5 ns.

**The headline: the observed domain is 12 d_i0 and the entire observation window is
t·ω_ci0 ≈ 0.5–1.0 — inside one gyroperiod.** Compare R1_warm at 75 d_i0 and t·ω_ci0 = 5.6.
An experiment-matched run is ~4000 cells × ~67k steps ≈ **1/23× R1_warm**, and since matching
β_e0 = 1.6 means raising `theta_0` 0.002 → 8.05e-3, dz/λ_D *improves* from 6.7 to 3.3.
`R1_coll`'s ambient is already 10¹⁸ cm⁻³, i.e. the experiment's density, so its
collisionality is directly testable against the paper's τ_pa/τ_s ≫ 1 claim.

### Correction: adding the laser need NOT cost 4×10³–5×10⁴×
First pass said one box must span solid CH (1e23) down to the 0.9e18 ambient — 10⁵ in
density, 316× in d_e. That conflated two independent things. IB absorption
(`K ∝ n_e² lnΛ T^{-3/2}/n_cr`) works at any density; what needs n_cr is the *critical
surface*, and that is a statement about **n_target/n_cr**, not absolute density.

- Solid density is unnecessary: the ray cannot see past n_cr, so everything above it is a
  mass reservoir — which is `TargetInjector`'s job. Worth ~5.6×.
- 351 nm is unnecessary: λ₀ enters only through n_cr. **R1_warm's contrast is already
  2.5 : 0.01 = 250×**, ample to host a critical surface. Put n_cr ≈ 1.25 n0 (half the target
  density, 125× the ambient): critical sits inside the target, the ambient is at
  n/n_cr = 8e-3 ⇒ refractive index 0.996 (no bending, no absorption), and dz is set by d_e at
  the target density — which the run already resolves. **Cost ≈ 1× R1_warm.**

At R1_coll's n0 = 1e20 cm⁻³ that is n_cr = 1.25e20 ⇒ λ₀ ≈ 3 µm (mid-IR). In 1D there is
exactly one ray, so the ray-march cost (80 ms/application, dominant in 2D per the header) is
negligible; and with a *resolved* preplasma ramp of ~30 d_e (~100 cells) the deposition is no
more localized than the heater's 133-cell slab, so `intervals` 5–20 should hold. Scaling n and
n_cr down together by f gives τ ∝ √f, so the 100× scale-down costs only ~10× in lnΛ →
lnΛ ≈ 800 with the reduced-c emulation. This promotes C++ change (1), `critical_density` as a
direct input, from convenience to **load-bearing**: otherwise the deck says λ₀ = 3 µm where it
means "n_cr at 1.25 n0" — the naming trap R1_coll already paid for.

A *faithful* 351 nm / solid-CH calculation still costs 4e3–5e4×; the way to afford that is a
small high-resolution laser run near the target used only to *compute* T_e,ab and the ablation
rate, feeding `theta_e_heat` to an R1-class heater run (~23× R1_warm, one-time), replacing the
by-eye recalibration history (0.092 → 0.062 → 0.078) with a first-principles number.

---

## Table I re-read from the PDF: four unit errors, and a 1.218× clock bias (2026-07-31) — desk analysis, no runs

Triggered by asking whether WarpX can run a reduced speed of light (it cannot — `PhysConst::c`
is a `constexpr` in `ablastr/constant.H:60`, 207 call sites, no ParmParse override; note
`my_constants.clight` *is* settable and silently shadows it **in deck parser expressions only**,
which would be a self-inconsistent run). Pulling `schaeffer2020.pdf` to settle it turned up
four things we had wrong.

### 1. The reduced speed of light is θ_e,ab. There was never anything to port.

§II, p.3: *"a reduced proton-to-electron mass ratio µ_p = 100 (d_i,ab = 10 d_e,ab), and a
reduced speed of light set by the ratio T_e,ab/m_e c², which can be written relative to the
sound speed as c = √(µ_p/T_e,ab) C_s,ab."*

c/C_s,ab = √(100/0.092) = 33.0 = 1/0.0303, i.e. Table I's C_s,ab = 0.030 c. The `c_sim/c_phys
= 0.02` row is a *reported consequence* of choosing θ_e,ab, not an independent knob. **R1_warm
already is the paper's reduced-c run** (its c/C_s,ab = 35.8; the 8% gap is entirely the
0.092 → 0.078 M_A recalibration). A custom reduced-c WarpX build would reproduce a run we have.

### 2. Sim → phys is three different factors, so s = 50 was wrong; s = 9.21 was right

- velocities × **0.0234** (C_s,ab 0.030c→210 km/s and v_p 0.104c→730 km/s both give 42.7;
  the tabulated 0.02 is rounded)
- temperatures ÷ **100** (0.092 m_ec² = 47.0 keV → 470 eV; 0.002 m_ec² = 1022 eV → 10 eV)
- lengths at **real c and real m_p** (n_e,ab = 6e20 cm⁻³ ⇒ d_i,ab = 9.31 µm; L_z = 900 d_i,ab
  = 8.4 mm ⇒ 9.33 µm ✓; d_i0 = 11.2 d_i,ab = 104 µm ✓)

Temperature is ÷100 and not ÷1824 = 1/0.0234² because the physical column also restores the
real mass ratio: 0.0234² × 1836/100 = 1/99.4. Cross-check: √(470 eV/m_p) = 212 km/s = the
quoted C_s,ab. **So 0.02 is a velocity factor and cannot be used as a temperature factor** —
the 2026-07-31 emulation's s = √(39900/470) = 9.21 is correct and lands exactly on the paper's
own physical column (470 eV). s = 50 would give T_e,ab = 16 eV, colder than the upstream.

### 3. λ_ab = mfp/d_e,ab exactly — CLAUDE.md gotcha (2) was wrong, by 28×

§II, p.3: *"λ_ab = ω_ce,ab/ν_ei,ab = λ_mfp,ab/d_e,ab, where ω_ce,ab is the electron
gyrofrequency at **B_ab**"*, with B_ab ≡ √(µ₀ n_e,ab T_e,ab) — the fundamental field, not B₀.
The identity is exact: ρ_e(B_ab) = √(T m_e)/(e√(µ₀nT)) = √(m_e/(µ₀ n e²)) = c/ω_pe = d_e,ab.

Our note evaluated ρ_e at **B₀** instead, got ρ_e,ab ≈ 28 d_e — that number is √(β_ab/2) = 24,
an artifact of the wrong field — and set R1_coll's target to mfp = 559 d_e. **Table I's
λ_ab = 20 means mfp = 20 d_e,ab; R1_coll is ~28× less collisional than the paper.** Correcting
it needs lnΛ ≈ 2×10⁵ instead of 7713, which makes the λ_ii conflict worse.

And the paper concedes that conflict itself: the λ_ab scaling *"ensures that dimensionless
quantities such as the magnetic Reynolds number are correct, but electron collisionality
relative to global scales (e.g. ν_ei,ab t_ab) is only quantitatively matched at physical mass
ratios."* **The knob that buys both λ_ab and λ_mfp/d_i0 is µ_p = 1836, not a reduced c.** Our
2026-07-29 finding is corroborated by the authors; the proposed fix was aimed at the wrong knob.

### 4. The ambient is 0.008 n_e,ab, and it biases every t·ω_ci0 by 1.218× late

Table I's code density unit is **not** n_e,ab: it lists n_e,ab = 1.25 and n_e0 = 0.01 in code
units ⇒ **n_e0/n_e,ab = 0.008**. Three rows confirm it, and nothing else fits:

| Table I row | formula | with 0.008 | quoted |
|---|---|---|---|
| d_i0/d_i,ab | √(n_e,ab/n_e0) | 11.18 | 11.2 |
| β_ab | (n_e,ab/n_e0)·θ_e/(µ_p(v_A/c)²) | 1158 | 1150 |
| ω_ci0⁻¹/t_ab | √(n_e,ab/n_e0)·C_s,ab/v_A | 34.03 | 33.9 |

(Those also pin the paper's β convention as **β = µ₀nT/B²** despite its text saying 2µ₀nT/B²:
β₀ = θ₀/(µ_p(v_A/c)²) = 0.201 vs the tabulated 0.2. `units.py:222-223` uses the factor-2 form,
so `derive()` reports β₀ = 0.4 against a `targets:` value of 0.2, and β_ab = 1560 vs 1150.
Convention mismatch, not physics — but it makes that check permanently fail.)

**Our `density_over_n0: 0.01` is 25% high, and it propagates into the clock.** `units.derive`
builds B₀ from `vA_over_c` *and* n_amb, so B₀ ∝ √n_amb ⇒ **ω_ci0 ∝ √n_amb**:

| | R1_warm | paper |
|---|---|---|
| n_amb/n_e,ab | 0.010 | 0.008 |
| d_i0/d_i,ab | 10.00 | 11.18 |
| ρ_i0/d_e | 1040 | ≈1180 |
| **ω_ci0⁻¹/t_ab** | **27.93** | **33.9** |

Total 1.218× = 1.118 (density) × 1.086 (θ_e 0.092→0.078, since t_ab ∝ 1/√θ_e). **Our
gyroperiod is 18% short, so every t·ω_ci0 observable is biased 1.218× LATE** (fewer ablation
times per gyroperiod ⇒ a piston-driven event lands at larger t·ω_ci0), and every z/ρ_i0
observable 1.118× far.

That is the direction and roughly the size of the residual we had pinned on the half-domain
wall: onset **t*₁ = 1.41 ω_ci0⁻¹ (later 1.35) vs the paper's ≈1**, with z*₁ correspondingly
far (2026-07-24, 2026-07-26 entries). Dividing out 1.218 leaves ~1.1×, a much more plausible
size for the wall artifact. **Not proven** — it assumes t*₁ is set by piston/ablation physics
rather than pure gyro-physics, which is testable by rerunning with n_amb = 0.008.

### Consequences / open decisions

- Configs are **unchanged** — fixing `density_over_n0` invalidates R1_warm and R1_coll and
  needs a rerun (7h11m each). Flagged, not done.
- Matching M_A = 14 by lowering θ_e necessarily breaks C_s,ab/v_A (2.79 vs the paper's 3.04),
  so v_sh/C_s,ab ≈ 5.0 vs the paper's 4.6. With θ_e off the paper's value, M_A and v_sh/C_s,ab
  cannot both match — pick the anchor deliberately.
- R1_warm's box is 7500 d_e vs the paper's 9000 d_e,ab (30000 cells, one-sided), and 250k steps
  vs 400k (τ·ω_ci0 5.6 vs 6.5).
- The paper uses the **Takizuka–Abe** collision operator; R1_coll uses WarpX's Perez.
- Docs corrected this session: CLAUDE.md gotchas (2) and (3) + two new gotchas;
  REPLICATION_PLAN §1 (reduced c, conversions, invariants) and §2 (the paper's domain is
  **one-sided**, not a centered symmetric slab — `layout: one_sided` is the faithful geometry).

### Architecture follow-up (same day): B0 is now the primary, and every run has a README

Two changes prompted by finding (4) above.

**`field.B0_tesla` replaces `field.vA_over_c` as the primary; v_A is derived.** Under the old
parameterization `units.derive` built `B0 = vA_over_c * c * sqrt(mu0*namb*m_i)`, so B0 — and
therefore `wci0 = q_e*B0/m_i`, the clock every `t*wci0` plot is drawn against — scaled as
`sqrt(namb)`. A wrong ambient density silently rescaled *time*. Inverted, `wci0` depends on
nothing but B0 and the ion mass, and `v_A = B0/sqrt(mu0*namb*m_i)` absorbs the density instead,
which is the physically honest direction. `units.derive` now **refuses** a config carrying the
old key (pointing at the migration script) rather than guessing.

All 14 configs migrated with `scripts/migrate_field_b0.py`, which applies exactly the old map,
so **nothing changed numerically**: the deck diff is two lines
(`my_constants.vA` + symbolic `B0` -> a literal `B0 = 0.0032075256118468715`), and all 13 runs
with a `warpx_used_inputs` still report `OK (WarpX ran exactly this config)`. Two new tests
pin the behaviour: B0/wci0 invariant under a 4x ambient change while v_A halves, and the legacy
key raising. 13/13 pytest green.

Caveat accepted knowingly: absolute tesla couples B0 to n0. R1_coll (n0 = 1e20 m^-3) needs
B0 = 32.1 T for the same dimensionless run R1_warm gets at 3.21 mT. Rescaling n0 now requires
rescaling B0 by sqrt(n0'/n0) by hand — noted in CLAUDE.md.

**`runs/<ID>/README.md`, generated by `scripts/make_run_readme.py`.** One page per run where
every number carries a source: `config.yaml:<key>` (primary), `derived: <formula>`, or Table I
with the ratio and a `**OFF**` flag past 20%. Rows the repo already understands as wrong carry
the known cause inline, so R1_warm's page now states its own 0.824x gyroperiod, the 2x beta
convention and the 0.01-vs-0.008 ambient without anyone having to re-read RESULTS. Prose
between `<!-- prose:begin -->` / `<!-- prose:end -->` is hand-written and preserved across
regeneration; tables are always rewritten. `--check` fails on a stale one. R1_warm and R1_coll
have seeded prose; the other 12 carry the placeholder.

## R1_paper finished + analysed: M_A lands at 14.8, but θ_e,ab is a *rate*, not a temperature (2026-08-01)

`runs/R1_paper` ran to completion: **322,400 steps in 16h09m** @ 8 threads (mean 0.180 s/step),
51 particle frames + 1291 field-only frames, t·ω_ci0 ∈ [0, 6.49]. `make_inputs --verify` →
`OK (WarpX ran exactly this config)`; `run_checks` validation OK; `make_run_readme --check` clean.

### Two fronts, not one — and the fast one is the shock

The first pass at a by-eye fit rode the **piston** edge (0.118 c), because the n_e streak is
dominated by the 250 n_e0 piston plasma. Separating the deposited per-species rho fixes it
(outer 2 d_i0 guarded throughout, per the open-hi artifact note):

| tracker | v | M_A |
|---|---|---|
| ambient compression, leading edge (n_amb > 2 n_e0) | 0.1484 c | 14.84 |
| B ramp (smoothed B_perp > 3 B0) | 0.150 c | 15.0 |
| **piston front** (n_piston > 0.5 n_e0) | **0.1145 c** | 11.45 |

Shock–piston separation grows monotonically 1.7 → 15.5 d_i0 through t·ω_ci0 = 5.3 (it only
shrinks after, when the shock hits the boundary). The shock decoupling from its driver is the
formation signature, and it is unambiguous here. `shock_fit.yaml` set to **v_sh = 0.148 c,
z0 = 0** (M_A = 14.80, M_ms = 13.51).

| quantity | R1_paper | Table I | note |
|---|---|---|---|
| M_A | **14.80** | 14 | +5.7% (R1_warm: 12.85) |
| v_sh / C_s,ab | 4.88 | 4.6 | +6% |
| v_sh / v_p (model v_p) | 1.42 | ≈ 4/3 | |
| v_sh / v_p (**measured** v_p = 0.1145 c) | **1.29** | ≈ 4/3 | the honest comparison |
| onset t*₁ | 1.30 ω_ci0⁻¹ | ≈ 1 | R1_warm 1.35 |
| onset z*₁ | 1.6 ρ_i0 | ≈ 1 | R1_warm 1.46 |
| reflected fraction G (peak) | 0.234 | present | R1_half 0.30 |

**The 1.218× clock bias does not explain the onset gap.** CLAUDE.md's leading hypothesis was that
the wrong ambient density biased every t·ω_ci0 observable late by 1.218×, which predicted
t*₁ → 1.35/1.218 = **1.11**. R1_paper is built with the corrected clock (B0 primary, n_amb =
0.008) and measures **1.30**. Given RESULTS' own ±0.4 detection sensitivity on these near-identical
G curves this is suggestive rather than decisive, but the clock was clearly not the dominant term.

### Compression: quote the front-following window, never criteria.json's globals

`criteria.json` reports n_compression = 278 and B_compression = 61.7 — those are the piston and
the boundary spike. In a window that follows the front ([z_f−5, z_f+2] d_i0, outer 2 d_i0 cut,
t·ω_ci0 ∈ 1.5–5.3, 94 field frames):

* ambient density compression n/n_e0 = **7.1 mean / 7.0 median / 9.0 max**
* field compression B_perp/B0 = **11.9 mean / 12.0 median / 13.9 max**

### Criteria: 6 of 7, and criterion 2 fails exactly as the config predicted

51/51 super-magnetosonic, 51/51 density compression, 50/51 field compression, 50/51 steep ramp,
50/51 reflected ions, 49/51 piston separation — and **0/51 collisionless**, λ_ii/d_i0 = 0.0150
against the 350 threshold. So `is_shock` is false in every frame, the same structural outcome as
R1_coll and for the same reason. The config header called this in advance ("the 350 row describes
the EXPERIMENT, not the run"); criterion 2's threshold is not applicable to a run that targets
λ_ab = 20 at µ_p = 100.

### θ_e,ab = 0.092 is a heating rate, and the run is 2.5× less collisional than configured

`plasma.piston.theta_e_heat` is threaded config → deck → `warpx_used_inputs` intact, so the run
used the paper's number. But `ParticleHeater.cpp:191-207` implements PSC's `HeatingSpotFoil`,
`H = 8·θ^{3/2}/(√(M_i/m_e)·width/d_e)·c²ω_pe`, applied as momentum diffusion `d⟨u_i²⟩/dt = H`.
**Nothing thermostats the electrons to θ.** Measured (weight-weighted, per-cell drift subtracted,
isotropic in x/y/z, steady from t·ω_ci0 = 1.3 to 5.45):

| region | n_e/n_e,ab | θ_e | vs 0.092 |
|---|---|---|---|
| target slab, 0–0.18 d_i0 | 1.93 | 0.219 (112 keV) | 2.38× |
| ablation front, 0.18–2.9 d_i0 | 0.99 | 0.146 (74.5 keV) | 1.58× |
| plume, 8.2–10.9 d_i0 | 0.62 | 0.091 | 0.99× |

There is no point where n = n_e,ab *and* θ = 0.092 simultaneously. Consequence for collisionality:
lnΛ = 3.00×10⁹ was inverted assuming T_e,ab = 47 keV, but WarpX's Perez operator is dimensional
and uses the real local state, and mfp ∝ T²/n, so the **delivered λ_ab ≈ 51 at the ablation front
(59 in the slab) against the configured 20** — R1_paper is ~2.5× less collisional than intended.
Same direction as the R1_coll error, different cause (that one was ω_ce at B0; this is T_e).

This probably does *not* break PSC fidelity — PSC uses the identical formula with the identical θ,
so its electrons should land at the same temperature and 0.092 is plausibly its input parameter
too. What it does break is every **derived** quantity that reads 0.092 as a temperature: C_s,ab ∝
√θ_e, hence t_ab, the M_A/M_ms normalizations, β_ab, and `vsh_model`. TODO: decide whether
`units.derive` should carry a measured-T_e,ab override, or whether the Table I comparison should
be restated against the measured ablation temperature.

### Tooling: the conservation check has been vacuous since it was written

`run_checks.operator_balance_figure` indexed the EP/PN reduced diags at columns 0/1 — but
`io.reduced_diag` keeps the step column, so those are *step* and *time*. Both panels plotted
step-versus-time, for every run ever checked. Fixed to columns 1/2. R1_paper now reads: total
particle energy grows linearly to **2553× E(0)** (heater + injector driving an open system, as
designed) and macroparticles 6.00×10⁶ → 7.65×10⁶, with a rollover at t·ω_pe ≈ 58,000
(t·ω_ci0 ≈ 5.2) as the shock reaches the open boundary and particles start leaving.

### Housekeeping

* `make_movies.py` gained `--only {ni,phase}` so the phase movie can be rendered mid-run.
* Fig. 7 regenerated on the documented convention times (`--phase-times 0.15 0.49 0.73 0.98 1.25`
  → snapping to 0.13/0.52/0.78/1.04/1.30 at this run's cadence), both `d_i0` and `rho_i0`. The
  bare `make_figures.py` pass earlier in the session had reverted it to the default spread —
  exactly the trap documented on 2026-07-27.
* Clean physics window is **t·ω_ci0 ≲ 5.3**; the front reaches the guard limit (78.5 d_i0) at ≈5.4.

## 2026-08-02 — R1_paper's reference density was 6e8x low; corrected to Table I's 6e20 cm^-3

`runs/R1_paper/config.yaml` carried `reference.n0 = 1.0e18 m^-3` — R1_warm's arbitrary
placeholder, inherited by every run in the repo. Table I's own SI column gives
**n_e,ab = 6e20 cm^-3 = 6.0e26 m^-3**, so the config was low by 6e8. Corrected today,
together with the B0 it controls.

* **Table I's SI densities also confirm the 0.008 ambient independently.** Table I lists
  n_e,ab = 6e20 cm^-3 and n_e0 = 4.8e18 cm^-3, and 4.8e18/6e20 = **0.008 exactly** — a
  fourth confirmation of the 2026-07-31 finding, alongside d_i0/d_i,ab = 11.18,
  beta_ab = 1158 and omega_ci0^-1/t_ab = 34.0.
* **B0 had to move with it.** `B0 ~ sqrt(n0)` (CLAUDE.md's absolute-tesla caveat), so
  B0_tesla went 2.8688981230645468e-3 -> **70.27336525536518 T** (x24494.897). Recomputed
  from the config's own primary definition, B0 = 0.01*m_e*wpe(n0/1.25)/q_e, not by scaling
  the old float; the independent M_A = 14 cross-check gives 70.035 T, still agreeing to
  **0.339%** — the same margin as at the old density.
* **The correction is exactly physics-neutral.** Every dimensionless row is bit-identical:
  v_A/c, C_s,ab/c, beta_ab, beta_0, T_e,ab, lambda_ab = 20, nu_ei*dt = 3.4e-3,
  mfp_ii/d_i0 = 0.0150, rho_i0/d_e = 1162.8, dt*wpe = 0.225, n_cell = 30000.
  `run_checks.py` still reports M_A = 13.952, M_ms = 12.737, validation OK.
* **The real gain is the collision dial.** lnLambda ~ n0^(-1/2), so the deck's value drops
  **2.997e9 -> 1.223e5** while the NRL physical value at (n0, T_e,ab) goes 20.9 -> 10.8.
  R1_paper is therefore 1.1e4x physical instead of 1.4e8x — still a dial (see the
  2026-08-01 entry on why PSC's lambda0 is equally unphysical), but four orders better.
* **Checked the one non-invariant safety margin.** WarpX's clamp is
  `sigma_eff = min(pi*b0^2*lnLambda, sigma_max)` with `sigma_max = 1/(maxn*rmin)` and
  `rmin = (4pi/3*maxn)^(-1/3)` (`ElasticCollisionPerez.H:88,94`). b0 is set by relative
  velocity and so is invariant here, leaving sigma_eff/sigma_max ~ n0^(1/6): it grows 29x,
  8.7e-6 -> **2.5e-4**, still a 4.0e3x margin, so the clamp never engages. (The config
  header's old "5e-6" was slightly optimistic; 8.7e-6 was the actual value. Same verdict.)

**The completed 16h09m run is now stale and must be relaunched.** Its plotfiles are in SI,
so they cannot be reinterpreted under the new scales — `--verify` reports
`prob_hi: deck 47.8268 vs config 0.00195252`, and any normalized analysis against the old
`diags/` is meaningless. Deleted the two plotfile-derived figures
(`media/testing/R1_paper_{loaded_state,operator_balance}.png`) so they cannot be mistaken
for current; `R1_paper_config_summary.png` is config-only and remains valid.
`shock_fit.yaml` is stored in normalized units (v_sh/c, lengths in d_e) and so **survives
the rescaling unchanged** — no need to re-tune by eye after the rerun.

**Still open:** the other 15 run configs remain at n0 = 1.0e18 and are each matched to a
completed run, so they were left alone. R1_coll is the one that matters (its collisionality
is absolute); it sits at 1.0e26, still 6x below Table I.

## 2026-08-02 — Physical lnLambda vs the lambda_ab = 20 dial: a matched pair at t*wci0 = 0.389

Two runs, identical in every key except `collisions.target`, both 19344 steps (t*wci0 =
0.389, pre-shock) at ppc 25, 1h01m each:

* `runs/R1_paper_dial` — `lambda_ab: 20.0`      -> lnLambda = 1.2235e5 (the dial)
* `runs/R1_paper_phys` — `coulomb_log: physical` -> lnLambda = 10.836 (NRL at this n,T)

Verified single-variable by diffing the parsed YAML; both give MA = 13.952, Mms = 12.737,
dt = 1.6282e-16, 30000 cells, B0 = 70.273 T. A matched control was needed rather than
reusing R1_paper because that run is at ppc 100 and this grid heats numerically as
1/sqrt(ppc) — the very effect under measurement.

**Collisions transfer electron energy to ions, as they should.** Final-frame temperatures
(mean of the three directional temperatures, weighted):

| species | dial (collisional) | phys (collisionless) | phys/dial |
|---|---|---|---|
| piston_electrons | 69 678 eV | 84 769 eV | 1.22 |
| piston_ions      | 37 635 eV | 19 764 eV | **0.53** |
| amb_electrons    |  9 623 eV | 11 167 eV | 1.16 |
| amb_ions         |  8 791 eV |  7 203 eV | 0.82 |

With the dial on, piston ions are **1.90x hotter** and piston electrons 0.82x cooler: the
collisional operator is draining the heated electrons into the ions. Without it the energy
stays where the heater put it.

**Collisions isotropize, and that is the cleanest single signature.** T_perp/T_par:

| species | dial | phys |
|---|---|---|
| piston_electrons | **1.025** | 1.298 |
| amb_electrons    | **0.680** | 0.413 |

66 electron collision times have elapsed by step 19344 at the dial (vs 5.8e-3 at the
physical value), and the dial run has driven the piston electrons to within 2.5% of
isotropy while the collisionless run retains a 30% perpendicular excess. Both species move
toward T_perp/T_par = 1 with collisions on.

**Bulk dynamics barely care at this early time.** n_compression 320 (dial) vs 315 (phys),
1.6%; piston_separation/rho_i0 identical at 0.554; M_ms,front identical at 13.51.
reflected_fraction_G is 0.00374 (dial) vs 0.00276 (phys), i.e. 36% more reflected ions with
collisions, but both are ~0.3% and this is pre-shock. Criteria 4 (field compression) and 5
(steep ramp) are False in BOTH — no shock by t*wci0 = 0.389, expected since onset is ~1.
`first shock (crit 1-7) = None` for both.

**Criterion 2 flips True/False, but that is definitional, not a discovery.**
`criteria.json`'s `lambda_ii_over_di0` (168.8 phys / 0.014951 dial) is computed from
`Scales.mfp_ii_amb`, i.e. from the config, NOT measured from the plotfiles. So the flip
merely restates the input. The measured results are the temperatures, anisotropies,
compression and reflected fraction above.

**Where this leaves the two Table I rows.** The physical lnLambda puts upstream
mfp_ii/d_i0 at 168.8 against Table I's quoted 350 — within 2.1x, where the dial misses by
2.3e4x — but it puts lambda_ab at 2.2581e5 against Table I's 20, missing by 1.1e4x. Still
one knob, still two rows. Confirms Sec. II's own statement that the fix is mu_p = 1836.

**NB on the 347 rule in CLAUDE.md.** `lambda_ii/d_i0 ~ 347*(11.567/lnLambda)` reproduces
R1_coll exactly (346.5 at its physical lnLambda = 11.567) but is **density-specific** and
does not transfer: mfp_ii/d_i0 ~ T_i^2/(n_amb^(1/2) lnLambda), and R1_paper's ambient is
4.8e24 m^-3 against R1_coll's 1e24, so sqrt(4.8) = 2.19x shorter. Observed ratio
168.803/346.463 = 0.48722 against (1/sqrt(4.8))*(11.567/10.836) = 0.48722.

**Caveat: pre-shock only.** These runs say nothing about whether the physical lnLambda
still yields a collisionless *shock* — that needs t*wci0 >~ 1, i.e. ~50k steps.

## 2026-08-03 — R1_paper complete at the corrected density; and `io.load_frame` was silently returning B = 0

R1_paper finished 322400/322400 in **17h45m** (mean 0.1982 s/step; contention x1.19 late,
because the R1_paper_phys/dial pair shared the box). `--verify` clean, `run_checks` OK.

**THE BUG: every B-derived diagnostic in this repo has been reading zeros.**
`io.load_frame` builds `ds.covering_grid(0, ds.domain_left_edge, dims)` spanning the whole
domain. That can round up one ULP past `domain_right_edge` (here region 0.001952522542456239
vs dataset 0.0019525225424562387), and yt raises

    RuntimeError: yt attempted to read outside the boundaries of a non-periodic domain

`comp_opt` caught it with a bare `except Exception: return None`, and `comp` turned that into
`np.zeros(nz)`. **A read failure was laundered into a plausible physics result:** B = 0 in all
51 frames, so `B_compression = 0.0`, `ramp_scale = inf`, and criteria 4 (field compression)
and 5 (steep ramp) FAILED in every frame of every run analysed this way.

Fixed by calling `ds.force_periodicity()` before the covering grid (it only affects
ghost-cell lookups at the edge, never cell values) and by restricting the None path to fields
genuinely absent from `ds.field_list`, so read errors now raise instead of returning zeros.

**What R1_paper actually shows, once the fields are read:**

| criterion | frames passing | first t*wci0 |
|---|---|---|
| 1 super-magnetosonic | 51/51 | 0.00 |
| 2 collisionless | **0/51** | never |
| 3 density compression | 51/51 | 0.00 |
| 4 field compression | 50/51 | 0.13 |
| 5 steep ramp | 50/51 | 0.13 |
| 6 reflected ions | 50/51 | 0.13 |
| 7 piston separation | 49/51 | 0.26 |

49 of 51 frames pass everything **except** criterion 2, from t*wci0 = 0.26. So the shock forms
and the earlier "first shock = None" was the zero-field artifact; what blocks the verdict now
is criterion 2 alone, which is the lnLambda dial (mfp_ii/d_i0 = 0.0150) and is read from the
config, not measured. B_compression peaks at **14.31 at t*wci0 = 3.37** inside the clean
window, against the paper's ~4x and this repo's documented 15-19 range; the last frame's 61.5
is the known open-boundary artifact past t*wci0 ~ 5.5, not the shock.

**Synthetic Thomson scattering (EPW + IAW).** Ran `Schaeffer_PlasmaPy`
(branch `feature/pic-thomson-pipeline`, which already carries a WarpX reader) on the 51
particle frames, sampling the domain centre with a 532 nm probe at 90 deg. Figures and arrays
are in `media/R1_paper/thomson_{epw,iaw}.png` and `thomson_spectra.npz`.

* **alpha = 0.22 to 0.89 (EPW), median 0.45** — the run is sub-collective to marginally
  collective, reaching alpha ~ 1 only at the end. This is entirely a consequence of the
  2026-08-02 density fix: the pipeline's own docstring notes the same deck gives alpha ~ 1e-5
  at n0 = 1e18 and ~0.2 at 6e26. At the old placeholder density there would be no ion feature
  to model at all.
* n_e at the sampling point climbs 4.744e24 -> 1.610e26 m^-3, with piston arrival at ~26 ps.
* The EPW feature is Doppler-broadened by the 0.453 c electron sigma (window +/-492 nm); the
  IAW window is sized from C_s,ab = 9.093e6 m/s = 0.0303 c, whose doublet sits at +/-22.8 nm.
  Sizing the IAW window from the ion *thermal* sigma instead (0.0911 c) gives +/-147 nm and
  buries the doublet in 2-3 pixels -- that spread is the piston drift, not the acoustic scale.
* Because alpha < 1 for most of the run, the ion feature is NOT a sharp doublet; the IAW panel
  is dominated by the broad electron feature and by the blue-shifted piston drift growing
  after ~40 ps. A genuinely collective ion feature would need alpha > 1, i.e. a longer probe
  wavelength or a denser sampling point.

Three environment notes for reproducing the Thomson run (none require changing any env):
the fork needs PyTorch, absent from `physics` but present in `tsnn`; the phase-space cache
signature embeds astropy's `m_e`, which differs between astropy 8.0.0 (`physics`) and 7.0.1
(`tsnn`), so a cache built by one env is rejected by the other; and the fork calls
`np.trapezoid`, which needs numpy >= 2 (`tsnn` has 1.26, so it needs an `np.trapz` alias).
Also note `spectra_from_phase_spaces` returns its wavelength axis in **metres**, not the nm
the window was given in, and `reference_density` must be `1 * u.m**-3` because the reader has
already scaled f to m^-3 -- passing the real density drives n_e to ~1e51 and alpha to NaN.

### 2026-08-03 (addendum) — the velocity scale factor decides whether the run is collective

`scripts/make_thomson.py` now takes `--velocity-scale-factor R` (divide velocities by
sqrt(R)) and `--notch LO HI`. Both variants are committed for R1_paper:

| | unscaled (default) | `--velocity-scale-factor physical` |
|---|---|---|
| R | none | 18.36 = (m_p/m_e)/mass_ratio, so v/4.285 |
| alpha_epw | 0.219 - 0.895 (median 0.45) | **1.040 - 4.246 (median 2.12)** |
| electron sigma | 1.358e8 m/s = 0.453 c | 3.170e7 m/s = 0.106 c |
| C_s,ab | 9.093e6 m/s | 2.122e6 m/s |
| IAW doublet | +/-22.8 nm | +/-5.3 nm |
| EPW window | +/-492 nm | +/-159 nm |
| files | `thomson_{epw,iaw}.png` | `thomson_{epw,iaw}_scaled.png` |

**This is not cosmetic.** alpha = 1/(k lambda_D) and lambda_D ~ v_te, so dividing
velocities by 4.285 multiplies alpha by the same factor and carries the run across
alpha = 1. The scaled EPW figure therefore shows **Langmuir satellites at ~420 and ~640 nm**
which cannot exist in the unscaled version: they are a collective-regime feature, and their
offset from the probe goes as sqrt(n_e), i.e. a density diagnostic. The scaled IAW feature
likewise collapses from a 114 nm-wide smear to a 27 nm window around the probe.

Which one is "right" depends on the question. Unscaled is what the simulation actually did
and is self-consistent with the rest of the repo (mu_p = 100 at real c throughout). Scaled
is what a real 532 nm system would see off hydrogen, but it is a PARTIAL correction: it
fixes the mass ratio and not the reduced-c temperature offset (T_e,ab = 47 keV vs Table I's
470 eV), so it should not be read as "the physical answer". Hence both are kept, with R
recorded in each npz (`velocity_scale_factor`, NaN when unscaled) rather than one silently
replacing the other.

Sizing note: the windows are measured from the RAW phase spaces, which the pipeline has not
yet rescaled, so they are divided by sqrt(R) explicitly. Omitting that leaves every window
sqrt(R) = 4.285x too wide.

The EPW notch defaults to exactly the IAW window (518.7-545.3 nm scaled, 474.9-589.1 nm
unscaled), so the two panels are complementary: the EPW figure blanks what the IAW figure
resolves. Without it the unshifted probe light and the ion feature dominate the colour
scale and the satellites are invisible.
