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

### R0 smoke run — PASS (structure) (`runs/R0_phase/R0/`, 2000 cells, 2000 steps, ~exit 0)
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
  - `media/R0_phase/R0/`: `shock_streak.png`, `shock_trajectory.png`, `shock_lineouts.png`,
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
   - **1D analog (`runs/xcheck_phase/xcheck_flatfoil_1d/`), same params, uniform slab:** saturates at
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

**Two bugs found and fixed during R0 validation** (`runs/R0_phase/R0_half`, 1000 cells vs R0's 2000):
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

`runs/R1_phase/R1_core_half`: `layout: one_sided`, wall@z=0, `0..3600 d_e`, dz=0.3 d_e, 12000 cells,
100 ppc, 125000 steps (~180 min, tmux `wxr1half`). Compared frame-for-frame against the
complete full-domain `runs/R1_phase/R1_core` on z≥0 (same dt & diagnostic cadence).
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

**New run `runs/R1_phase/R1_core_half_sym`.** Byte-identical deck to `R1_core_half` except the single
added line `boundary.reflect_symmetry_axis = x` (verified by diff) — a clean one-variable A/B.
`filter_npass` deliberately left unchanged (that is a separate physics-quality knob; changing
it too would confound the wall comparison). Launched with the benchmarked optimization
`OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=cores` (~1.8× vs 4 threads; ETA ~1h36m);
no `max_grid_size`/tiling/sort tweaks (benchmarked neutral-to-negative). `runs/R1_phase/R1_core_half_sym/`
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

**R0_half_sym early confirmation (DONE — `runs/R0_phase/R0_half_sym`, smoke tier, seconds).** Ran the
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

`runs/R1_phase/R1_half`: the Full-tier Table I run as a one-sided half-domain with the validated
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
naive symmetric/4-thread reproduction, at faithful M_A=14. Figures: `media/R1_phase/R1_half/` (shock_streak,
trajectory, lineouts, phase, reflected, criteria.json) + `media/testing/R1_half_*`.

### Verdict
The WarpX `ParticleHeater` + `TargetInjector` + symmetry-wall half-domain reproduce a
Schaeffer-2020-class M_A=14 perpendicular collisionless shock at Table I parameters, with the
speed model, formation timescales, ~4× compression, and reflected-ion signature — at ~1/4 the
naive compute. Remaining: R2 (B₀=0) and R3 (n_e0=0) negative controls.

---

## R2 (B₀=0) + R3 (n_e0=0) negative controls (2026-07-26) — `runs/R2_phase/R2`, `runs/R3_phase/R3`

Built as **exact single-knob clones of the calibrated production run `runs/R1_phase/R1_warm`** (warm
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
per step with no ambient particles). Figures in `media/R2_phase/R2/`, `media/R3_phase/R3/`.

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
  when the plotfiles carry no per-species rho (verified on `runs/R1_phase/R1_core`: **`diag1` has only 9
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

New run `runs/R1_phase/R1_coll/`: R1_warm's **collisional** comparison run (the plan's R6 slot,
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
- `make_inputs.py runs/R1_phase/R1_coll` → deck round-trips to the config; `--check` still clean for
  R1_warm/R1_recal/R1_cal/R2/R3/R1_core_half (no regression from the new code paths).
- `run_checks.py runs/R1_phase/R1_coll` → **validation OK**, λ_ab = 20.000 exactly (mfp = 558.6 d_e,ab).
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
scripts/launch.sh -b -L runs/R1_phase/R1_coll
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

**Launched** 2026-07-28 10:08 (`scripts/launch.sh -b -L runs/R1_phase/R1_coll`, pid 1233159, logger
1233161). Early logged rate 0.143 s/step(warpx) cumulative / 0.114 instantaneous at step
200; with R1_warm's 1.34× growth over the run that lands at **~10–11 h**, finishing ~21:00.
Watch `runs/R1_phase/R1_coll/progress.log`.

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

`runs/R1_phase/R1_paper` ran to completion: **322,400 steps in 16h09m** @ 8 threads (mean 0.180 s/step),
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

`runs/R1_phase/R1_paper/config.yaml` carried `reference.n0 = 1.0e18 m^-3` — R1_warm's arbitrary
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

* `runs/R1_phase/R1_paper_dial` — `lambda_ab: 20.0`      -> lnLambda = 1.2235e5 (the dial)
* `runs/R1_phase/R1_paper_phys` — `coulomb_log: physical` -> lnLambda = 10.836 (NRL at this n,T)

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
are in `media/R1_phase/R1_paper/thomson_{epw,iaw}.png` and `thomson_spectra.npz`.

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

### 2026-08-03 (addendum 2) — the fig7 default-times trap caught R1_paper too

The bare `make_figures.py runs/R1_phase/R1_paper` passes run earlier today silently rebuilt
`shock_fig7.png` on the **default** `--nframes 5` spread instead of the convention times,
exactly as documented on 2026-07-27 and as happened to R1_warm on 2026-07-29. The tell was
again the mtime mismatch: `shock_fig7.png` freshly written while `shock_fig7_rho_i0.png`
still carried an 08-01 date, because the default pass does not write the rho_i0 variant.

Rebuilt both on the convention times:

    python scripts/make_figures.py runs/R1_phase/R1_paper --only fig7 \
        --fig7-xunits d_i0 rho_i0 --phase-times 0.15 0.49 0.73 0.98 1.25

which snap to t*wci0 = **0.13 / 0.52 / 0.78 / 1.04 / 1.30** at this run's 51-frame cadence
(R1_warm's snapped to 0.11/0.45/0.68/1.01/1.24 -- the requested times are the same, the
available frames differ). The panels show the reflected-ion branch above v_z/v_sh = 1
building through the formation interval, with B_x/B_0 reaching ~6 and the n_e/n_e0 ramp
co-located, consistent with the criteria table above (49/51 frames pass everything but
criterion 2).

Also regenerated `shock_phase.mp4`, which was dated 08-01 13:19 -- older than even the
superseded run, and rendered while `io.load_frame` was still returning B = 0. Rebuilt from
the new plotfiles with `--vsh-c 0.148` so its reflected-ion threshold uses `shock_fit.yaml`'s
by-eye v_sh rather than the config model value of 0.1395 c.

Still stale in `media/R1_phase/R1_paper/`: `tune_trajectory.png` (08-01 21:54). That one comes from
`tune_shock.py`, which is a BY-EYE fit and deliberately not auto-run.

### 2026-08-03 (addendum 3) — by-eye shock fit for R1_paper: v_sh = 0.143 c

`runs/R1_phase/R1_paper/shock_fit.yaml` re-fitted by eye against the new run's streaks:
**v_sh 0.148 -> 0.143 c** (-3.4%), z0 = 0, no per-time overrides. Derived Mach numbers
move with it: **M_A 14.80 -> 14.30**, **M_ms 13.51 -> 13.05**.

That lands closer to Table I's M_A = 14 and M_ms = 13 than the previous fit (now +2.1% and
+0.4%, against +5.7% and +3.9% before), and closer than the config's own model value
(vsh_over_Csab = 4.6 -> 0.1395 c, M_A = 13.95).

Every v_sh consumer regenerated: `shock_{streak,trajectory,lineouts,phase,reflected}.png`,
`criteria.json`, both `shock_fig7*` and `shock_phase.mp4` (`--vsh-c 0.143`).

The criteria verdict is unchanged by the refit -- 49/51 frames still pass everything except
criterion 2, first at t*wci0 = 0.26, and peak B_compression is still 14.31 at t*wci0 = 3.37
(a field quantity, independent of v_sh). What does move is `reflected_fraction_G`, 0.1132 ->
0.1243 at the last frame, since the reflected-ion threshold is v_z > v_sh and a slower shock
admits more particles.

Two conventions had to be honoured in the same pass, so it was split in two:
`args.phase_times` feeds BOTH `fig_phase` (shock_phase.png) and `fig_fig7`
(`make_figures.py:669,672`), but the documented `--phase-times 0.15 0.49 0.73 0.98 1.25`
convention applies only to fig7. Running the suite with those times would have silently
retimed shock_phase.png; running it without them would have re-triggered the fig7
default-times trap. So: `--only streak trajectory lineouts phase reflected criteria` first,
then `--only fig7` with the convention times.

### 2026-08-03 (addendum 4) — Fig. 7 velocity axes pinned per row

`make_figures.py` gained `--v-ambient`, `--v-piston` and `--v-electron` (each LO HI, in
v_z/v_sh). Previously the two ion rows shared fig_phase's band (-1 ... 3) and the electron
row was auto-sized SYMMETRICALLY from the 99.5th percentile of |v_z|/v_sh, clipped to
[3, 10] — so it varied from run to run and no two runs' Fig. 7 electron rows were directly
comparable. R1_paper's Fig. 7 is now drawn on:

    --v-ambient -0.5 2.0   --v-piston -0.2 1.7   --v-electron -6 12

The ion bands are tightened to the populated region (the shared -1 ... 3 wasted the top
third on both rows), and the electron range is deliberately ASYMMETRIC, which the old
auto-sizing could not express. The twin axis carrying B_x/B_0 and n_e/n_e0 still aligns its
zero with v_z = 0 at any range: the alignment is computed from the actual v-axis endpoints,
so at -6 ... 12 the zero sits a third of the way up and the right-hand scale runs -4 ... 8.

Defaults are unchanged when the flags are omitted.

---

## 2026-08-03 — The three free parameters: Table I reproduced in PSC / physical / WarpX units

Re-derived the whole setup from first principles at the user's direction, and built
`scripts/table1.py` to render Table I in all three unit systems at once with the paper's own
values in a check column. **Every Table I row now reproduces**, several of them exactly, and
three long-standing "cannot match" items in this log turn out to have been comparisons
between different quantities.

### The framework: exactly three free choices
A run of this problem is a set of dimensionless numbers; it corresponds to a *family* of real
plasmas. One member is picked by choosing three things — everything else follows:

1. **ablation density** `n_e,ab` (real 6e20 cm⁻³; code O(1), PSC uses 1.25)
2. **ablation temperature** `T_e,ab` (real 470 eV; code 0.092, O(0.1) to stay
   non-relativistic). This **is** the reduced speed of light: §II sets
   `c = sqrt(mu_p/T_e,ab)·C_s,ab`, so θ_e,ab and c_sim are the same choice.
3. **collisionality** `lambda_ab = mfp_e,ab/d_e,ab = 20`

Order of operations: pick the three real values → pick the code values → derive the rest of
the code column → derive the rest of the physical column, where **β_ab sets B₀**.

### Table I, checked (physical column, real c and real proton mass)
| row | derived | Table I |
|---|---|---|
| d_i,ab | 9.30 µm | 9.31 µm |
| C_s,ab | 212 km/s | 210 km/s |
| v_p / v_sh | 728 / 976 km/s | 730 / 980 km/s |
| **B₀ (from β_ab)** | **7.03 T** | **7 T** |
| T_0 | 10.2 eV | 10 eV |
| d_i0 | 104 µm | 104 µm |
| 1/ω_ci0 | 1.49 ns | 1.5 ns |
| M_A / M_ms | 13.95 / 12.74 | 14 / 13 |
| **β_ab / β_0** | **1150.0000 / 0.20000** | **1150 / 0.2** |
| τ_ei,ab | 477 fs = 0.0109 t_ab | 0.43 ps = 0.009 t_ab |
| c_sim/c_phys | 0.0233 | 0.02 |

### Three items previously logged as unmatchable — all resolved
**(1) The β convention is `mu0·n·T/B²`, no factor of 2 — now pinned, not a taste call.**
Two exact identities: (a) in PSC's normalization `mu0·n·kT/B² = θ·n_code/B_code²`, and Table
I's own code primaries give `0.092·1.25/0.01² = 1150` and `0.002·0.01/0.01² = 0.2`, both
*exactly*; (b) §II's `1/ω_ci0 = (Z_ab/Z_0)·sqrt(β_ab)·t_ab` needs `sqrt(1150) = 33.9`, which
is Table I's own ω_ci0 row (a factor-2 β gives 48.0). `units.py` now drops the 2 and reports
1150.0000 / 0.20000 for R1_paper. β is diagnostic only (reported + `check_factor`), never fed
to a deck, so this is inert for physics. CLAUDE.md's earlier "don't fix that 2×" note is
superseded. This also independently re-confirms `density_over_n0: 0.008` and `n_e,ab = 1.25`.

**(2) `λ_mfp/d_i0 = 350` is the DIRECTED-ion mfp, not the thermal one.** Rutherford
momentum transfer, ∝v⁴: `4πε₀²m_i²v⁴/(nZ⁴e⁴lnΛ)` gives **261 d_i0 at v_p and 845 d_i0 at
v_sh**, bracketing 350 (≈3.6 cm). The thermal upstream ion-ion mfp at T_0 = 10 eV is
**8.9e-4 d_i0** — six orders down, and no lnΛ bridges that gap. So `λ_ab = 20` and 350 were
never competing for one lnΛ; the 350 row is the statement "the experiment is globally
collisionless". ⚠ **Criterion 2 still compares the thermal `mfp_ii_amb` against a 350
threshold, which is a category error — its verdicts should not be trusted until reworked.**
That also explains the 2026-07-29 R1_coll finding that criterion 2 "FAILS in all 51 frames".

**(3) The physical Coulomb logarithm is ≈10, and λ_ab = 20 needs no dial physically.**
Deriving through `lambda_ab` (a pure electron-scale ratio, so mass-ratio-safe):
v_te,ab = 9.09 Mm/s, d_e,ab = 217 nm, mfp = 20·d_e,ab = 4.34 µm, ν_ei = 2.095e12 s⁻¹,
τ_ei = 477 fs = 0.0109 t_ab — against Table I's own **0.009 t_ab / 0.43 ps**. Inverting
`ν_ei = C·n·lnΛ·T^-3/2` gives **lnΛ = 9.0** at C = 3.95e-6, or 12.2 at NRL's 2.91e-6;
NRL's `24 - ln(sqrt(n)/T)` gives 6.2 at the same (n, T). Table I's 0.43 ps is reproduced
exactly by C = 3.95e-6 with lnΛ = 10, which is strong evidence the paper used that
coefficient and a physical lnΛ of 10.

### The 18.4× trap (why an earlier hand-derivation gave lnΛ = 0.5)
`ν_ei,ab·t_ab = (v_te/mfp)(d_i/C_s) = mu/lambda_ab`, which carries the ion mass explicitly:
**5.0** at mu_p = 100 but **91.8** at mu = 1836. §II flags exactly this quantity: the λ_ab
scaling "ensures that dimensionless quantities such as the magnetic Reynolds number are
correct, but electron collisionality relative to global scales (e.g. ν_ei,ab t_ab) is only
quantitatively matched at physical mass ratios." Routing the physical lnΛ through
ν_ei t_ab = 5 (the *code* value) with the *real* t_ab gives ν_ei = 1.14e11 s⁻¹ and
**lnΛ = 0.49**, low by exactly mu_phys/mu_p = 18.36. A lnΛ below 1 is also outside the
validity of the Coulomb logarithm itself. Reproduced both routes in `--show-work` so the
distinction is checkable rather than asserted.

### What did NOT change
The **deck is untouched and remains correct.** WarpX has no reduced-c option, so it must run
PSC's dimensionless problem at real c: θ_e,ab = 0.092 → T_e,ab = 47 keV, B₀ = 70.27 T
(= β_ab 1150 at *that* T, verified), and lnΛ = 1.22e5 as an unavoidable dial (at 47 keV the
plasma is genuinely collisionless; making lnΛ physical would need n0 ≈ 9e34 m⁻³). The
physical column's 7 T and the WarpX column's 70.27 T differ by the reduced-c factor
√100.8 ≈ 10, exactly as `field.B0_tesla`'s comment already said. Confirmed the deck's clock
too: 220 t_ab = 52.6 ps = 6.50/ω_ci0, matching `max_step`'s comment. 11/11 tests pass.

### 2026-08-03 addendum — the physical column must use the run's OWN mass ratio

The entry above built its physical column at µ = 1836 so it would land on Table I's printed SI
values. **That was the wrong framing.** µ = 100 is a *physical* choice, not merely a code
convenience: the plasma the run represents really is a light-ion plasma, and mixing µ = 1836
scales with µ = 100 relations is what generated the "18.4× trap" narrative in the first place.
Using µ = 100 throughout, the whole thing is self-consistent and no trap exists.

**The setup, restated.** Choose, in physical units: n_e,ab = 6e26 m⁻³, T_e,ab = 470 eV,
λ_ab = 20 d_e, µ = 100, β_ab = 1150, T_0 = 10 eV, n_e0 = 0.008 n_e,ab. Everything else follows;
β_ab sets B₀ = **7.0264 T**. `scripts/table1.py` was rewritten around this.

| | PSC (code) | Physical / WarpX | Table I |
|---|---|---|---|
| T_e,ab | 0.092 m_e c² | 470 eV | 470 eV |
| B₀ | 0.01 √(m_e c²) | **7.026 T** | 7 T |
| T_0 | 0.001957 m_e c² | 10 eV | 0.002 / 10 eV |
| C_s,ab | 0.0303 c | 909.2 km/s | 210 km/s (µ=1836) |
| d_i,ab | 10 d_e,ab | 2.169 µm | 9.31 µm (µ=1836) |
| d_i0 | 11.18 d_i,ab | 24.26 µm | 104 µm (µ=1836) |
| 1/ω_ci0 | **33.91 t_ab** | 80.92 ps | 33.9 t_ab / 1.5 ns |
| β_ab / β_0 | 1150 / 0.1957 | 1150 / 0.1957 | 1150 / 0.2 |
| M_A / M_ms | 13.95 / 12.76 | 13.95 / 12.76 | 14 / 13 |
| c_sim/c_phys | **0.100** | 1 (real c) | 0.02 |

Every dimensionless row matches. The ion rows differ from Table I's SI column by √18.36 (speeds)
or 18.36 (gyro-times) purely because that column is real hydrogen — the caption calls it "one
possible set of experimentally-relevant physical values", an illustration, not a unit map.

**β_0 is over-determined.** Given n_e0, T_0 and B₀ it is fixed at 0.1957 (Table I's 0.2). A
requested 0.02 would need n_e0 = 0.00082 n_e,ab or T_0 = 1.02 eV, contradicting the other
inputs; `--beta-0` now reports this rather than silently accepting it.

**Bonus: at µ = 100 the reduced c becomes self-consistent.** c_sim/c_phys = 0.100 fits *both*
the temperature rows (470/0.092) and the velocity rows (v_p = 0.104 c_sim = 3117 km/s). Table I's
printed 0.02 only fits its velocity rows; the 4.3× = √(1836/100) gap that earlier entries
recorded as an internal inconsistency is just an artifact of its µ = 1836 physical column.

**A real bug this exposed.** `1/ω_ci0 = √β_ab·t_ab` held in the real-c column but came out
339 t_ab (10× off) in the reduced-c one, because ω_ci0 was taken as qB₀/m_i (physical) against
a reduced-c t_ab. ω_ci0 must be taken in the *same* normalization as that column's lengths and
times, i.e. from B_code. Now 33.91 in both, asserted in the test along with exact c-scaling:
n, T, all speeds, all frequencies, λ_D, β and dt are c-independent and must be bit-identical
between columns, while d_e, d_i, t_ab, dz, mfp, ρ_i0 and 1/ω_ci0 must scale precisely as c.

**The deck this implies, and its price.** θ_e_heat 0.092 → 9.19767e-4, θ_0 0.002 → 1.95695e-5,
B0_tesla 70.273 → 7.0264468095, max_step 322400 → 3224046. The grid is **unchanged**
(d_e,ab = c/ω_pe depends only on n0 and real c): dz = 65.08 nm, 30000 cells, 1.953 mm.
The collisions block is **also unchanged** — and that is the payoff: `quantity: lambda_ab,
value: 20` now resolves to **lnΛ = 12.23** on its own, a physical value, where at 47 keV it
needed 1.22e5. (Not `value: physical`, which is `24 − ln(√n/T)` = 6.23, a different quantity.)

Two costs, both from real c, neither avoidable in WarpX:
1. **10× the timesteps** for the same 220 t_ab — dt is CFL-locked to dz/c while t_ab ∝ c.
2. **dz/λ_D,ab goes 0.99 → 9.89.** λ_D is c-independent but d_e ∝ c, so PSC's d_e,ab is 3.3 λ_D
   against 33 λ_D at real c. The paper's 0.3 d_e,ab cell is then 9.9 λ_D wide and **will
   grid-heat**; resolving it needs dz ~10× smaller and dt with it, ~100× on top of the 10×.

So this deck is physically faithful where the current one is dimensionally faithful, and the
choice between them is a real trade, not a bug fix. **Not applied to any config** — the cost
and the Debye-resolution consequence are the user's call. 15/15 tests pass.

---

## 2026-08-04 — Pilot of R1_paper_470eV: grid heating is REAL. Full run NOT cleared.

`runs/R1_phase/R1_paper_470eV_pilot`, 50,000 steps = 3.41 t_ab = 0.101/ω_ci0, completed clean in
**1h33m** at 0.1117 s/step (load ~7, 8 threads), zero failure signatures. Physics primaries
bit-identical to the parent. **Verdict: do not launch the 6-day run as configured.**

### The global ParticleEnergy diagnostic is useless for this — measure spatially
`EP.txt`'s whole-domain ambient-electron mean rises **5.9×** (10 → 59 eV) in 3.41 t_ab, which
looks catastrophic and is almost entirely an artifact of averaging. The profile at the final
frame shows why (12 bins across the domain, matched against R1_paper at 4.40 t_ab):

| z/L | pilot T [eV] | pilot T/T_0 | R1_paper T/T_0 |
|---|---|---|---|
| 0.00–0.08 | 604.1 | **60.4** | **21.1** |
| 0.08–0.17 | 11.36 | 1.136 | 1.148 |
| 0.42–0.50 | 10.80 | 1.080 | 0.999 |
| 0.75–0.83 | 10.45 | 1.045 | 0.997 |
| 0.92–1.00 | 10.39 | **1.039** | **0.991** |
| GLOBAL | 59.01 | 5.90 | 2.65 |

The inner bin is **physical** piston heating and is present in both runs (60× vs 21× T_0), so
the global mean is dominated by it. `scripts/grid_heating.py` restricts to the outer 25%.

### The real signal: +3.9% vs −0.9% in the far upstream
Pilot far upstream: **+3.9% in 3.41 t_ab** (1.14 %T_0/t_ab). R1_paper over the same window:
**−0.9%**, i.e. statistically zero. With ~62,500 macroparticles/bin the noise floor is ~0.3%,
so the 4.8% difference is real, not sampling. **R1_paper has no grid heating; this run does.**

Note `grid_heating.py`'s "RATIO" headline is meaningless here (it printed −103×) because the
baseline slope is ~0 and slightly negative — divide-by-near-zero. Use the absolute fractional
rate, not the ratio, whenever the baseline is flat.

### Why R1_paper is immune — same √T relation, one line
The relevant number is the *ambient* Debye resolution, not the ablation one:

| | T_0 | λ_D,amb | dz/λ_D,amb |
|---|---|---|---|
| R1_paper | 1022 eV | 108.5 nm | **0.60** |
| R1_paper_470eV | 10 eV | 10.73 nm | **6.07** |

dz = 65.08 nm is pinned by the density and identical in both. λ_D ∝ √T, so R1_paper's 100×
hotter upstream already resolves it and sits below threshold; the 470 eV run is 6× under.

### Where it ends up is the open question — and the range spans "fine" to "fatal"
Grid heating stops when λ_D,amb grows to dz, i.e. at **T_amb = 368 eV = 36.8× T_0**. The shock
reaches the far boundary at 195.7 t_ab, so the upstream ahead of it heats for that long:

| | T_0 | C_s0 | M_ms | β_0 |
|---|---|---|---|---|
| initial | 10.0 eV | 132.6 km/s | 12.76 | 0.196 |
| linear extrapolation | 32.4 eV (3.2×) | 238.6 km/s | 10.92 | 0.633 |
| saturation ceiling | 367.9 eV (36.8×) | 804.4 km/s | **4.87** | **7.20** |

Linear is a 14% M_ms degradation — survivable, though β_0 leaves Table I's 0.2. Saturation
destroys the regime (M_ms 12.8 → 4.9). **3.41 t_ab is 1.5% of the run and cannot distinguish
these**; the far-upstream points (10, 9.983, 9.964, 10.02, 10.41) even dip before rising, so the
trajectory is genuinely unconstrained. The t^2.27 power law I first fitted came from the
*contaminated* global signal and should be ignored.

### Refined cost, and the levers
0.1117 s/step measured × R1_paper's 1.279 particle-growth drift = 0.1429 s/step mean:
**5.33 d idle, 7.39 d at R1_paper's load ~24.**

~~No current filtering is enabled in the deck at all~~ **WRONG, corrected 2026-08-05.**
`warpx.use_filter` defaults to **ON** for the explicit evolve scheme (`parameters.rst:3842`),
so the deck's silence means 1-pass bilinear filtering is already active: the completed run's
profiler shows `Filter::ApplyStencil` called 8,358,775 times for **0.16%** of runtime. The
lever is therefore not "enable filtering" but "raise `filter_npass_each_dir` from its default
1", which at 0.16% per pass is nearly free. RESULTS 2026-07-2x measured `filter_npass = 8`
cutting far-upstream noise 31%.
Raising ppc cuts heating ~1/ppc but multiplies an already-10× cost. Halving dz is ~100×.

**Recommended next step, not taken:** a ~250k-step (17 t_ab, ~8 h) pilot with filtering on,
which both pins the trajectory past the ambiguous window and tests the cheap fix at once.
Committing 5–7 days before knowing which of the two columns above applies would be premature.

---

## 2026-08-04 — `runs/` and `media/` regrouped by phase; run discovery centralised

`runs/` had grown to 20 flat directories (14 of them R1 offshoots) and `media/` mirrored it.
Both are now grouped by phase, with a run dir always at **depth 2**:

```
runs/R0_phase/{R0, R0_half, R0_half_sym}
runs/R1_phase/{R1, R1_cal, R1_coll, R1_core, R1_core_half, R1_core_half_sym, R1_half,
               R1_paper, R1_paper_470eV, R1_paper_470eV_pilot, R1_paper_dial,
               R1_paper_phys, R1_recal, R1_warm}
runs/R2_phase/R2   runs/R3_phase/R3   runs/xcheck_phase/xcheck_flatfoil_1d
runs/opt_phase/    <- performance sweeps, no config.yaml
```

The `_phase` suffix is not cosmetic: `R1` is itself a run id, so a bare `runs/R1/` container
would be simultaneously a run and a folder of runs, which is exactly what forces tooling to
handle two depths. `media/` mirrors this; `media/testing/` stays put as cross-study figures.

### The trap this created, and the fix
`make_run_readme.py --all` and `migrate_field_b0.py --all` both hard-coded
`glob.glob(runs/*/config.yaml)`. After the regrouping that pattern matches **nothing**, and
because both feed the result straight into a `for` loop the failure mode is a silent
"processed 0 runs" — no error, no output, looks like success. Replaced with one definition,
`kinshock.find_runs(root)` (searches depth 2 *and* depth 1), plus `kinshock.unphased_runs(root)`
so a run not yet filed into a phase is reported rather than quietly bypassing the convention.
`.gitignore` needed no change — its patterns are all `runs/**/…`, already depth-agnostic.

131 lines across 72 tracked files plus 43 lines across the run dirs' own README/config
comments were repathed by exact-id regex (longest id first, so `R1_core_half_sym` is not
mangled by an `R1` rule; idempotent, so `runs/R1_phase/R1` cannot become
`runs/R1_phase/R1_phase/R1`). Placeholders (`runs/<ID>`) and gitignore globs are untouched.

Verified: 15/15 tests pass; `make_inputs.py --check` clean on R1_warm, R1_paper, R0, R2;
`make_run_readme.py --all` discovers all 19 runs.

### Incidental: the β-convention regeneration finally landed
`--all --check` flagged 17 run READMEs stale. The diff is **entirely** β lines — the
2026-08-03 convention fix (`2*mu0*n*T/B^2` → `mu0*n*T/B^2`) had never been regenerated, as
CLAUDE.md had noted. Values halve, targets halve, every ratio is unchanged (R1_warm stays
0.800×, R1_recal 0.678×), so this is inert for physics and only the printed labels move.

### One media wrinkle worth knowing
8 tracked files under `media/R0/` read as pure deletions after the move, because `media/` is
gitignored and so `git add -A` cannot stage the destination. The files were on disk the whole
time; `git add -f` at the new path preserved tracking. Anything force-added under `media/`
will do this on every future move.

---

## 2026-08-04 — optimization sweep: GPU is 7.89x, threads 1.82x, implicit is not a speedup

Three levers tested individually against `runs/R1_phase/R1_paper_470eV`, chosen because the
pilot's TinyProfiler says **~98% of the run is particle work** (GatherAndPush 49.9%,
CurrentDeposition 20.9%, collisions 9.9%, Redistribute 8.5%) over 6.0e6 macroparticles,
99.3% of them ambient at a flat 100 ppc — and **all field/grid work is ~0.9%**. Harness in
`studies/speedup/`, generated tables in `runs/opt_phase/SUMMARY.md`, figures in
`media/opt_phase/`.

Baseline 0.11169 s/step at 8 threads → **5.33 d** for the full 3,224,046 steps (incl. the
1.279 particle-growth drift).

### Lever 3, the winner — GPU, and the decomposition matters more than the hardware

| configuration | s/step (mean) | vs 8 thr | full run |
|---|---|---|---|
| GPU, deck default (235 boxes) | 0.09268 | 1.21x | 4.42 d |
| GPU, 8 boxes | 0.01626 | 6.87x | 0.78 d |
| **GPU, 1 box (`amr.max_grid_size=30000`)** | **0.01415** | **7.89x** | **0.68 d** |

A 6.5x swing from one ParmParse line. AMReX picks a CPU-friendly 235 grids of 112-128 cells,
which starves a GPU at ~25k particles per kernel launch. **A GPU benchmark left at the deck
default reads as "barely worth it" (1.21x)** — the wrong conclusion about the fastest lever
available. The same knob is neutral-to-negative on CPU. See memory `warpx-gpu-max-grid-size`.

### Lever 1 — OMP threads, and the "cliff" is conditional

4 → 0.53x, 8 → 1.00x, 12 → 1.28x, 16 → 1.58x, **20 → 1.82x**, 24 → 1.91x. Monotonic, no
cliff, flattening after 20. This does **not** contradict the 2026-07-27 measurement of a ~20x
collapse above 12 threads; that note's own mechanism (oversubscription against other users)
is the explanation. Today the box was idle (load 0.79) and the deck is 5x larger (30000 cells
/ 6.0e6 particles vs 6400 / 1.2e6), so even 24 threads gets ~10 of the 235 boxes each instead
of spinning on barriers. Rule is now conditional: idle box + production deck → 20 threads;
busy box or smoke deck → stay at 8-12. `max_grid_size=64` is worse at every thread count.

### Lever 2 — theta_implicit_em: correct, compatible, and 13.5x slower

It **runs**, which was the real risk: `TargetInjector` and `ParticleHeater` are applied in the
outer Evolve loop (`WarpXEvolve.cpp:286,291`) before the `if (m_implicit_solver)` branch at
`:299`, so they are scheme-agnostic; `doCollisions()` is called inside the implicit branch
(`:301`) and that pairing is the scheme's first cited reference (Angus et al. on implicit PIC
with binary Monte-Carlo Coulomb collisions); and all four assertions at `WarpX.cpp:1372-1395`
pass for this deck — including field gathering, which must not be momentum-conserving and
already is not.

Picard: **1.33961 s/step = 13.5x explicit**, 21-23 iterations, every step exiting on relative
tolerance (zero hit the ceiling). Break-even needs dt past **cfl ~10**. Not a speedup. Its
value is exact energy conservation at theta = 0.5 — no grid heating *by construction*, which
is the thing actually blocking this run.

Newton was **misconfigured by me**, twice over, and both are recorded because they are the
kind of mistake that looks like a result:
  1. `ImplicitSolver.cpp:811` gates the cheap particle path on
     `use_mass_matrices_jacobian && skip_particle_picard_init`. I set only the first, so every
     Newton iteration ran a full 21-iteration particle Picard update — 23.95 s/step at cfl
     0.75 (241x explicit), with `max particle iterations: 21` printed in the log and ~1680
     GMRES iterations per step. The retune's N1 point deliberately reproduces this.
  2. `newton.verbose=0` and `gmres.verbose_int=0` meant the first attempt logged no iteration
     counts at all, so a 568x cost could not be attributed. Timing verbosity and diagnostic
     verbosity are different requirements.

### Three harness bugs, all of the "looks like it worked" kind
- **`diag*.intervals=0` does not suppress plotfiles.** `dump_last_timestep` defaults to 1 and
  `parameters.rst:4702` says the last step is written "regardless of this parameter". Every
  point dumped a 277 MiB plotfile: 1.7 GiB on a disk at 94% before it was caught.
- **The first GPU-vs-CPU agreement check compared nothing.** The deck writes EP/PN every 5000
  steps and the benchmark is 1500, so it reported "1 common rows, steps 0..0" — which reads
  like a pass. bench.sh now rescales the reduced-diag interval to the run length (~30 rows).
- **Benchmark means were contaminated by my own filesystem work.** thr20/thr24 ran during the
  runs//media/ regrouping (a `find` over 130 GB, `git add -A` across the tree); mean/median
  was 2.69 against 1.13 for the clean re-run. Diagnostic: external I/O inflates the mean and
  leaves the median alone, whereas a real scaling collapse moves the median. Both kept.

Corollary on statistics: for *clean* runs prefer the **mean** — this deck's collision
supercycle (`ndt_supercycle: 10`) makes every 10th step genuinely expensive, and the median
hides work you pay for. The GPU points show it plainly (median 0.01094 vs mean 0.01415).

### What this changes for R1_paper_470eV
The 2026-08-04 pilot verdict ("full run NOT cleared") was a judgement about risk per day. At
**0.68 d instead of 5.33 d** the calculus differs: the grid-heating trajectory that 3.41 t_ab
could not distinguish — linear (T_0 → 32 eV, survivable) vs saturation (T_0 → 368 eV,
M_ms 12.8 → 4.9) — becomes cheap to settle by just running far enough to see it, rather than
by extrapolating. Still gated on the GPU agreement check; the CUDA binary is dated Jul 28 and
commit `9f981dea2` (Jul 31) fixed an nvcc rejection in `ParticleHeater`, so it predates a fix
to an operator this deck needs. It runs without aborting, which is necessary, not sufficient.

---

## 2026-08-04 (addendum) — GPU validated; implicit at large dt; the full lever matrix

Completes the sweep. Generated tables in `runs/opt_phase/SUMMARY.md`, figures in
`media/opt_phase/`, 30 benchmark points. **All costs below are `steps(cfl) x s/step x 1.279`
with `steps(cfl) = 3,224,046 x 0.75/cfl`** — dt is proportional to cfl, so a large-dt point
reaches 220 t_ab in fewer steps. Using a fixed step count (as the first version of
`make_summary.py` did) overstates the cfl 3.0 and 7.5 rows by 4x and 10x.

| configuration | s/step | full run | vs best |
|---|---|---|---|
| **explicit GPU, 1 box** | 0.01415 | **0.68 d** | 1.00x |
| **implicit GPU, cfl 3.0** | 0.10876 | **1.30 d** | 1.92x |
| implicit GPU, cfl 7.5 | 0.29374 | 1.40 d | 2.08x |
| explicit CPU, 20 thr | 0.06143 | 2.93 d | 4.34x |
| implicit GPU, cfl 0.75 | 0.07296 | 3.48 d | 5.16x |
| implicit CPU, cfl 7.5 | 0.84305 | 4.02 d | 5.96x |
| explicit GPU, 235 boxes | 0.09268 | 4.42 d | 6.50x |
| explicit CPU, 8 thr | 0.11169 | 5.33 d | 7.84x |
| implicit CPU, Picard cfl 0.75 | 1.33961 | 64.77 d | 95.3x |

### The GPU binary is validated
30 diagnostic rows over steps 0-1479, GPU vs CPU on an identical deck: no diverging column,
**ambient electron count bit-identical at every step**, weights agreeing to 1e-13, and the
piston species' relative energy differences *shrinking* 5-14x over the run as injection grows
their populations — initial RNG sampling noise averaging out, not divergence. Per-species
energy differences scale as 1/sqrt(N) at step 0 (amb_electrons 3.0e-4 with N=2.99e6 against
1/sqrt(N)=5.8e-4; piston_electrons 2.6e-2 with N=6667 against 1.2e-2), which is what
independent RNG streams look like. So the Jul 28 CUDA build is usable despite predating
commit `9f981dea2`'s nvcc fix to `ParticleHeater`. Caveat: 1479 steps is ~0.1 t_ab.

### Implicit converges at 10x dt but the cost rise beats the step saving
`hit_max_iter = 0` in **every** implicit point on both devices, including cfl 7.5 — the regime
the docs flag as where Picard "often" fails. Newton handles it. But per-step cost rises with
dt, and faster on GPU than CPU:

|  | cfl 0.75 -> 3.0 | cfl 3.0 -> 7.5 |
|---|---|---|
| CPU | 1.15x | 1.69x |
| GPU | 1.49x | **2.70x** |

At 2.70x against a 2.5x step saving, **cfl 7.5 is past the optimum** — implicit GPU is
cheaper at cfl 3.0 (1.30 d) than at 7.5 (1.40 d). Larger dt is not monotonically better, and
`SUMMARY.md` now orders rows by projected cost rather than s/step because with cfl varying the
two orderings disagree (implicit GPU cfl 0.75 has the fastest step of the implicit points and
is the most expensive run).

### The GPU helps implicit less than explicit
5.95x vs 7.89x, and the implicit penalty is *worse* on GPU (5.16x) than CPU (3.88x). GMRES and
mass-matrix assembly are grid work, and grid work is under 1% of the explicit run — the part
the GPU is least busy with. So the implicit scheme's extra cost lands where acceleration helps
least. **My cfl 7.5 GPU projections were 0.51 d, then 0.68 d, then 0.77-0.88 d; the answer is
1.40 d.** Each extrapolated a ratio measured at small dt into large dt, and the nonlinear
solve got harder faster than any of them assumed. Do not project implicit cost across a dt
change without measuring the interval.

### Recommendation
**Explicit GPU, 0.68 d**, for throughput: `amr.max_grid_size=30000` plus the CUDA binary.
**Implicit GPU at cfl 3.0, 1.30 d**, if the finite-grid instability should be made
structurally impossible instead of out-run — 1.92x the cost buys exact energy conservation
(theta = 0.5), which removes the uncertainty the pilot could not resolve (T_0 -> 32 eV
survivable vs T_0 -> 368 eV and M_ms 12.8 -> 4.87). Both are now cheap enough to run and
compare, which was not true at 5.33 d.

Two things are still required before either is a physics run, and neither is done:
`numerics.evolve_scheme` (plus `implicit_evolve.*`) plumbed through `scripts/make_inputs.py`
— no config passthrough exists and hand-editing decks is against the repo rule, so every
implicit number here came from ParmParse overrides that are legitimate for benchmarking only
— and a `grid_heating.py` measurement at matched physical time to confirm implicit actually
removes what the pilot found rather than merely being entitled to.

---

## 2026-08-04 (decision) — `R1_paper_470eV` configured for GPU: 5.33 d -> 0.68 d

`numerics.max_grid_size: 30000` added to the config (= `n_cell`, one box) and
`scripts/launch.sh` grew a `-g/--gpu` flag. The run is staged, **not launched**.

`max_grid_size` is plumbed through `kinshock.deck` properly rather than passed as a launch
override: it is rendered from config, it round-trips through `key_params` (sentinel 0 for
"absent", because absent means "AMReX decides" and that is a *different run* on GPU by 6.5x),
and so `--verify` will catch a deck whose decomposition silently changed. The key is opt-in,
not defaulted, because it is neutral-to-negative on CPU — `mgs=64` measured worse at every
thread count, and a CPU run should simply not carry it.

`launch.sh -g [N]` selects the CUDA build, drops to one thread, pins
`CUDA_VISIBLE_DEVICES=N` so the second card stays free for other users, prints the device,
and **refuses to start if the config lacks `numerics.max_grid_size`**. That last part is a
hard failure rather than a warning on purpose: forgetting it turns a 0.68 d run into a 4.42 d
run with no error message, which is exactly the kind of silent 4-day mistake this session
produced three other examples of.

Launch with:

    scripts/launch.sh -b -g runs/R1_phase/R1_paper_470eV

Preconditions verified: deck `--check` clean, README up to date, `diags/` empty (so no
`--force`), 111 GB free against the ~14 GB of 50 plotfiles, both GPUs idle, and the other
six run configs still `--check` clean after the deck.py change (15/15 tests pass).

**The upstream temperature is a monitored OUTPUT of this run, not a controlled input.** The
pilot's grid heating is real and unresolved: 3.41 t_ab could not distinguish linear
(T_0 -> 32 eV, M_ms 12.76 -> 10.92, survivable) from saturation (T_0 -> 368 eV,
M_ms -> 4.87, regime destroyed). The reason to take the GPU path is precisely that at 0.68 d
the question can be settled by running far enough to see it. Analyse with
`scripts/grid_heating.py` against `runs/R1_phase/R1_paper` at matched t/t_ab, restricted to
the outer 25% — the global `ParticleEnergy` mean is dominated by real piston heating and is
useless for this. The 50 plotfiles give a frame every 4.4 t_ab, against the pilot's entire
3.41 t_ab baseline, so the trajectory should be well resolved.

If it saturates, the fix is `theta_implicit_em` (1.30 d on GPU at cfl 3.0, 1.92x this run),
which conserves energy exactly and cannot grid-heat — not more resolution, which is ~100x.

**max_step truncated 3224046 -> 2784400 (220.0 -> 190.0 t_ab), so the run is 0.58 d / 14.0 h.**
`field_hi = pec` reflects fields (only particles are absorbed), so every step after the shock
reaches the far wall is contaminated rather than clean outflow, and the Table I duration spent
13.6% of the run there. Arrival is 190.9 t_ab at R1_paper's *fitted* v_sh = 4.72 C_s,ab and
195.7 t_ab at the model 4.6, so 190.0 sits just under the earlier estimate. This costs zero
science: the discarded steps were all post-arrival.

Deliberately NOT cut to the 186.5 t_ab B-spike onset (RESULTS 2026-07-29, ~80x B_perp/B0 at
the open boundary from t*wci0 ~ 5.5). That artifact is *localised* and already handled by
cutting the outer 2 d_i0 when quoting `B_compression`, whereas clean post-formation shock
propagation is only ~4.2-4.4/wci0 against the 5 wanted. Trading 4.5 t_ab of clean propagation
for 0.3 h of compute is the wrong direction.

`plotfile_intervals` rescaled 64481 -> 55688 so the frame COUNT stays 50 rather than dropping
to 43 — one frame every 3.80 t_ab, against the pilot's entire 3.41 t_ab baseline for the
heating measurement that is this run's key monitored output. NB R1_paper ran the full
220 t_ab, so compare over the overlapping window.

### Two GPUs: 1.77x, and the agreement check that first said no

Measured 2026-08-04, 1500 steps per point, `runs/R1_phase/R1_paper_470eV`:

| configuration | s/step | vs 1 GPU / 1 box |
|---|---|---|
| 1 GPU, 1 box (`mgs` 30000) | 0.01415 | 1.00x |
| 1 GPU, 2 boxes (`mgs` 15000) | 0.01462 | 0.97x — splitting alone costs 3.3% |
| **2 GPUs, 2 boxes** | **0.00801** | **1.77x** |
| 2 GPUs, 4 boxes (`mgs` 7500) | 0.00849 | 1.67x |

Against the matched 1-GPU/2-box control that is **1.83x, i.e. 91% parallel efficiency** — in
1D the MPI cost is a single guard-cell interface, and 3e6 particles still saturates a 4070.
Two boxes beats four, same direction as the single-GPU finding. Adopted:
`numerics.max_grid_size: 15000` and `launch.sh -g 0,1`, giving **7.9 h** (from 14.0 h on one
card, and 4.60 d on 8 CPU threads at this max_step).

Both custom operators are MPI-safe, checked in source rather than assumed:
`TargetInjector.cpp:277` calls `n_meas.SumBoundary(geom.periodicity())` after its CIC
deposition, so per-cell density accumulates correctly across rank boundaries, and
`ParticleHeater` makes no MPI calls at all (a per-particle kick from a position-evaluated
parser). The injector's only other MPI use is a `ReduceLongSum` for its log line.

**The agreement check initially FAILED and the check was the thing that was wrong.**
1-GPU vs 2-GPU flagged `piston_ions` at 1.6%, `total` at 0.5%. Two runs differing only in
`warpx.random_seed` then gave 1.03% and 0.31% on those columns — so the 2-rank differences
are ~1.6x pure stochastic scatter, while `amb_ions` (the bulk plasma) agrees to **0.1%**
(7.566e-4 vs 7.555e-4). Nothing structural.

Two defects in `compare_diags.py`, now fixed. Its rule was `last > max(10*first, 1e-3)`,
which (a) degenerates when the two runs start identical — `first = 0`, so anything above
1e-3 fires regardless of cause — and (b) had no notion of how much scatter is inherent to the
quantity. It therefore failed a run against *itself with a different seed*, which no correct
check can do. It now accepts `--floor C D` (two seed-only runs) and requires 3x that measured
floor. This also **strengthens** the earlier CPU-vs-GPU pass: those differences (total
2.79e-3, piston_ions 1.41e-3) are *below* this noise floor.

Note for reuse: the seed-only floor is the right tolerance for any comparison of stochastic
PIC runs here — thread count, rank count, device, or binary. Comparing to `first` only works
when the runs start from different RNG draws AND the column is not RNG-dominated.

---

## 2026-08-05 — R1_paper_470eV COMPLETED on 2 GPUs. Grid heating: 11.2x T_0, beta_0 -> 2.20

**2,784,400 / 2,784,400 steps in 10h47m, zero failure signatures**, 51 plotfiles, 20 GB.
First full-length run of the 470 eV physical-units rebuild.

### Cost: I under-projected by 36%
Projected 7.9 h from a 1500-step benchmark; actual **10h47m** (0.0139 s/step vs the
benchmarked 0.00801). Two causes, neither visible in a benchmark at t = 0:
  * **Particle growth 1.29x** — total macroparticles 6.00e6 -> 7.75e6, with
    `piston_electrons` going 6,667 -> **1,014,414** (152x). The 1.279 drift factor I applied
    was right, and is already included in the 7.9 h.
  * **Load imbalance, ~1.35x** — the remaining gap. Rank 0 owns the piston half of the
    domain, so it carries ~2e6 extra particles by late time and sets the pace for both
    ranks. A 2-rank split is balanced at t = 0 and progressively is not.
Lesson for future GPU projections here: benchmark from a *restart* partway through, or
apply a load-imbalance factor on top of the particle-growth drift.

### The grid heating: between the two branches, and closer to the bad one

**First measurement was wrong and nearly got reported.** `grid_heating.py`'s default outer-25%
window reads 3899 eV at t = 190 t_ab — but its own `piston reach` column shows the piston at
91% of the domain, past the window's inner edge from ~135 t_ab. That number is *shocked
plasma*, not grid heating. The docstring warns about exactly this ("If `piston reach` ever
approaches `window lo`, the measurement is contaminated"). Re-measured with
`--upstream-frac 0.95`:

| t/t_ab | t*wci0 | T_upstream | xT_0 | piston reach |
|---|---|---|---|---|
| 0.0 | 0.000 | 9.99 eV | 1.00 | 0.2% |
| 34.2 | 1.008 | 43.8 eV | 4.38 | 26% |
| 68.4 | 2.017 | 70.6 eV | 7.06 | 40% |
| 102.6 | 3.026 | 92.2 eV | 9.22 | 59% |
| **136.8** | **4.034** | **112.2 eV** | **11.2** | 76% |
| 148.2 onward | | 238 -> 1343 eV | | 81-91% — CONTAMINATED, do not quote |

Clean window ends ~137 t_ab. The rate **decelerates** (1.2 eV/t_ab early, 0.55 by t = 137),
consistent with approaching saturation without reaching it: at 112 eV, `dz/lambda_D,amb` has
fallen 6.07 -> 1.81.

### Consequences for the shock parameters

| T_0 | xT_0 | C_s0 | M_ms | beta_0 | |
|---|---|---|---|---|---|
| 10.0 eV | 1.0 | 132.6 km/s | 12.76 | 0.196 | Table I target |
| 32.4 eV | 3.2 | 238.7 | 10.91 | 0.634 | the pilot's linear extrapolation |
| **112.2 eV** | **11.2** | **444.2** | **7.80** | **2.196** | **MEASURED at 137 t_ab** |
| 367.9 eV | 36.8 | 804.4 | 4.87 | 7.201 | saturation ceiling |

**M_A is untouched at 13.95** — v_A depends on B0 and n_amb, not T — so the Alfvenic Mach
number and everything keyed to it survive. But **M_ms falls 12.76 -> 7.80 and beta_0 rises
0.196 -> 2.196**, an order of magnitude off Table I, and the upstream stops being
magnetically dominated. The pilot's linear extrapolation was optimistic by 3.5x; the
saturation ceiling was pessimistic by 3.3x.

**How to use this run.** The shock forms and propagates and M_A is preserved, so it is not
worthless -- but after the first ~30 t_ab it is NOT a Table I beta_0 = 0.2 replication.
Anything depending on upstream beta or the fast-magnetosonic Mach number must be quoted
against the measured T_0(t) above, not the nominal 10 eV. Treat T_0 as a time-dependent
output of the run.

**The fix, already measured:** `theta_implicit_em` at cfl 3.0 conserves energy exactly and
cannot grid-heat -- 1.30 d on one card, so ~0.75 d on both. Blocked only on plumbing
`numerics.evolve_scheme` (and `implicit_evolve.*`) through `scripts/make_inputs.py`.

### Two latent bugs this run exposed — one of which silently scrambled a figure

**1. Plotfile ordering broke past 1,000,000 steps.** `kinshock.io.plotfiles` and
`field_plotfiles` used `sorted()` on the *filename*. WarpX pads the step to SIX digits, so a
run crossing 1e6 steps emits mixed 6- and 7-digit names and string order interleaves them:

    diag_fields100000  <  diag_fields1000000  <  diag_fields997500

Every earlier run was under 1e6 steps (R1_paper 322,400; the 470 eV pilot 50,000), so string
sorting happened to be correct and the bug could not manifest. This run is 2,784,400 steps and
it did: the first `tune_shock` streak came out in blocky vertical bands with the whole domain
saturating after t*wci0 ~ 3.1. Fixed with a numeric `_step_key` on the trailing digits;
verified monotonic across the boundary (995000, 997500, 1000000, 1002500, ...).

Blast radius: anything that **assumed frame order** was silently wrong — streaks, movies,
front trajectories. Anything that **reads its own time per frame** was correct in value and
only mis-ordered on print, which covers `grid_heating.py`, so the heating table above stands
(its rows were sorted before quoting). Re-rendered streak is clean and physical.

**2. `media/` output escaped the mirrored tree.** `plotting.media_dir(run_id=...)` built
`media/<run_id>` from the id alone, so after the 2026-08-04 phase regrouping it wrote to a
flat `media/R1_paper_470eV/` instead of `media/R1_phase/R1_paper_470eV/`. Figures still
appeared, just outside the mirror — noticed only because `tune_shock` prints its output path.
`media_dir` now resolves the phase from the filesystem (`runs/*/<run_id>/config.yaml`), which
keeps all six call sites unchanged and falls back to the flat layout for an unknown id. Stray
flat directories migrated.

Both are the same shape as the harness bugs from the optimisation sweep: **the failure mode
was a plausible-looking output, not an error.** A blocky streak could be read as physics if
you were not expecting it, and a figure in the wrong folder is invisible until someone looks
for it.

### Diagnostics: 6 of 7 criteria pass; "no shock" is the known criterion-2 defect

`make_figures.py` produced the full A-D set plus movies into
`media/R1_phase/R1_paper_470eV/` (streak, trajectory, lineouts, phase, fig7, reflected,
`shock_ni.mp4`, `shock_phase.mp4`). It reported *"first shock (crit 1-7) at t*wci0 = None"*.
That headline is a **tooling artifact, not a physics result**:

| criterion | frames passing (of 50) |
|---|---|
| 1_super_magnetosonic | 50/50 |
| **2_collisionless** | **0/50** |
| 3_density_compression | 50/50 |
| 4_field_compression | 50/50 |
| 5_steep_ramp | 50/50 |
| 6_reflected_ions | 50/50 |
| 7_piston_separation | 49/50 |

Criterion 2 fails in **every** frame at `lambda_ii_over_di0 = 0.014321` against a threshold of
350 — exactly the category error CLAUDE.md already flags: 350 is Table I's **directed**-ion
mfp, while this measures the **thermal** upstream ion-ion mfp, six orders down. It cannot pass
by construction, so `is_shock` can never be True and the composite verdict is meaningless
until criterion 2 is reworked. The six physical criteria pass from t*wci0 = 0.22 onward.

**Two criteria values are CONSTANT across all 50 frames** — `M_ms_front = 15.5589` and
`lambda_ii_over_di0 = 0.014321`. Both are derived from `config.yaml`'s nominal T_0 = 10 eV, not
measured per frame. Given this run heats the upstream to 112 eV, criteria 1 and 2 are being
evaluated against an upstream state that stops being true after ~30 t_ab. Reworking criterion
2 should therefore also make both read the measured T_0(t).

Also worth not trusting yet: `n_compression ~ 270` and `B_compression ~ 15-19` in every frame
from t*wci0 = 0.11, i.e. before any shock could exist. Those look like global maxima over the
domain (as CLAUDE.md already warns for `B_compression`) rather than jumps across the shock
ramp, so criteria 3 and 4 passing is not independent evidence of a shock either.

**The auto v_sh fallback is 21% high and should not be used.** With no `shock_fit.yaml`,
`make_figures` fell back to `track_front` and got v_sh = 0.0170c (M_A = 17.01) — but the
re-rendered streak shows the bright B_perp ridge sitting *below* the model 0.0140c trial line
at late time (ridge ~70 d_i0 vs line ~72 at t*wci0 = 5.4), so the true v_sh is a little
*under* 0.0140c, not 21% over. This is the inter-script drift CLAUDE.md warns about, now
quantified. **Every figure in this set is provisional until `tune_shock.py` is run by eye and
the set regenerated.**

Front reached ~71 of 80.5 d_i0, so the shock never touched the wall and the boundary artifact
appears only in the last ~0.3/wci0 — the truncation to 190 t_ab was well placed.

**Movies regenerated after the ordering fix.** `shock_ni.mp4` and `shock_phase.mp4` were
rebuilt with the numeric `_step_key` in place. Their sizes changed (346,779 -> 359,619 and
364,361 -> 417,368 bytes), and since the encode is deterministic for fixed frames, that size
change is proof the frame content differed — i.e. the originals *were* mis-ordered. Note the
timestamps could NOT establish this: the mp4s were written within ~25 s of the `io.py` edit,
and a Python process caches imports at startup, so a job launched just before the edit still
runs the old code. Regenerating beats reasoning about it.

Mismatch to be aware of in the current figure set: `make_movies.py` annotates with the MODEL
v_sh = 0.0140c, while `make_figures.py` fell back to `track_front` and used 0.0170c. The
movies and the PNGs therefore carry different shock speeds right now — another reason the
whole set needs regenerating once `tune_shock.py` has been run by eye.

### By-eye fit adopted: v_sh = 0.0165c; and make_movies was ignoring shock_fit.yaml

`runs/R1_phase/R1_paper_470eV/shock_fit.yaml` written with **v_sh_over_c = 0.0165**
(M_A = 16.50), fitted by eye against `tune_trajectory.png`. Full figure set + both movies
regenerated against it, so the whole set finally shares one speed.

Two caveats recorded in the fit file's header, since they travel with anything derived from
it. (1) 0.0165c sits ABOVE the model 0.0140c and near `track_front`'s 0.0170c fallback; worth
knowing because R1_paper's fit was made the same way, so cross-run comparisons rest on the
ratio. (2) **M_ms = 15.09 is an overestimate for most of the run** — it uses the config's
nominal T_0 = 10 eV, but the upstream grid-heats to ~112 eV by t = 137 t_ab, which raises the
fast-magnetosonic speed ~3.3x. M_A = 16.50 is unaffected, since v_A depends on B0 and n_amb
and not on temperature. **M_A is the robust number for this run; M_ms is not.**

**`make_movies.py` never read the fit.** Line 99 was
`vsh = args.vsh_c * C if args.vsh_c else sc.vsh_model` — by design per its docstring
("measured one is passed with --vsh-c"). So the instant a `shock_fit.yaml` exists, the movies
and the figures are annotated with different speeds unless the flag is remembered: here
0.0140c against the figures' 0.0165c, an **18% disagreement between two outputs of the same
run**. Fixed by giving it `make_figures`' precedence (`--vsh-c` > `shock_fit.yaml` > model)
and printing which source was used, so a silent divergence cannot recur:

    R1_paper_470eV: 51 plotfiles -> movies (v_sh=0.0165 c from by-eye fit (shock_fit.yaml))

That is the third bug this session with the same shape — plausible output, no error — after
the plotfile ordering past 1e6 steps and `media_dir` escaping the phase mirror. Worth a habit:
when two scripts consume the same derived quantity, have each **print its source**, not just
its value.

`is_shock` still reads None, unchanged: criterion 1 shifted slightly with the new v_sh but was
already 50/50, and criterion 2 remains unpassable by construction. When it is reworked,
criteria 1 and 2 should also read the measured T_0(t) rather than the nominal 10 eV.

**Correction: the first two "regenerated" figure sets did NOT use the fit.** `make_figures`
loaded `shock_fit.yaml`, printed `v_sh=0.0165c`, and then two places ignored it:

  * `fig_streak(frames, cfg, sc)` never received the `shock` object at all and read
    `sc.vsh_model` unconditionally, so the white `v_sh` overlay was drawn at the MODEL
    0.0140c. Confirmed arithmetically from the plot: the line hit 78.1 d_i0 at
    t*wci0 = 5.6, which is 0.0140c; 0.0165c reaches 92 (off-scale).
  * `plotting.stamp()` always read `scales.MA / scales.Mms`, which are model-derived, so
    **all four stamped figures** printed `M_A=14.0 M_ms=12.8` instead of 16.5 / 15.1.

Fixed: `fig_streak` and `fig_lineouts` take `shock` and use `shock.v_sh`; all four `stamp`
calls pass it. `v_p` stays a model quantity deliberately — the fit covers the shock front,
not the piston. The stamp now also prints its SOURCE, `M_A=16.5 M_ms=15.1 (fit)` vs
`(model)` / `(auto)`.

The lesson is the one from the `make_movies` bug one step earlier, which I stated and then
failed to apply: **a derived value must carry its provenance into the artifact.** The log line
said "by-eye fit" while the figure body used the model, so the log was never evidence — only
reading the rendered image was. Verified that way this time: the stamp reads `(fit)` and the
`v_sh` line exits the top at t*wci0 = 4.85 as 0.0165c requires.

Final set, all at v_sh = 0.0165c: 7 PNGs + criteria.json + both movies.

### Figure: what the upstream beta saturates to — `media/R1_phase/R1_paper_470eV/upstream_beta.png`

Generated by `scripts/plot_upstream_beta.py <run_dir>` (reusable for any run; re-measures from
the plotfiles, so the figure cannot drift from the data).

**Two saturation numbers, differing by 1.66x, and the figure exists to keep them apart:**

| | T_0 | beta_0 |
|---|---|---|
| initial / Table I target | 10.0 eV | **0.196** |
| measured, end of clean window (141 t_ab) | 114.7 eV | **2.25** |
| **fitted asymptote** | 221.0 eV | **4.33** |
| ceiling, lambda_D = dz | 367.9 eV | 7.20 |

The **ceiling** is where the finite-grid instability must switch off: lambda_D ~ sqrt(T), so
T_ceil = T_0 (dz/lambda_D,0)^2 = 10 x 6.066^2 = 367.9 eV, and beta_0 scales with T at fixed n
and B_0. The **asymptote** is where the measured drive actually runs out. The rate falls off
linearly in sqrt(T) — exactly what a drive weakening as lambda_D -> dz should do —

    dT/dt = 1.718 - 0.1156 sqrt(T)      r = -0.897 over 32 intervals

and reaches zero at 221 eV, i.e. **60% of the ceiling**, where lambda_D/dz is still only 0.78.
So the instability dies before the Debye length quite reaches the cell size. (This reproduces
the independent hand fit of 222.9 eV / 4.36 to within 1%.)

**Within this run beta_0 never gets near either number**: it reached 2.25 at 141 t_ab and was
still climbing at ~0.5 eV/t_ab when the precursor arrived. Propagating the fit to the run's end
gives only ~2.6. The asymptote is where it would land given several hundred more t_ab.

A methodological point the figure forced. The first version used a piston-reach threshold
(<= 0.80) to define the clean window, which admitted the frame at 144.4 t_ab whose local rate
is **5.49 eV/t_ab against a trailing median of 0.62** — the precursor arriving. That one point
flipped the fitted slope positive and produced an "asymptote" of 52 eV, *below* the current
temperature, which is physically impossible for a saturating model. A reach threshold is the
wrong guard: grid heating is smooth and decelerating while the precursor is a step, so the
script now cuts on that contrast (`--break-factor`, default 3x the trailing median) and prints
where it cut. Reach survives only as a backstop.

### Staged: `runs/heat_phase/` — pilot matrix for the grid-heating levers (NOT run)

Five variants derived from `runs/R1_phase/R1_paper_470eV`, each differing from the baseline in
exactly one respect, all real config-driven runs rather than ParmParse overrides. Runner:
`studies/heating/run_matrix.sh` (serial — each variant uses both GPUs, so two at once would
halve each other's cards and make the timings meaningless).

| variant | key under test | cost | time | plotfiles |
|---|---|---|---|---|
| `h0_baseline` | — (control) | 1.00x | 0.98 h | 3.3 GB |
| `h1_filter8` | `filter_npass: 8` | 1.01x | 0.99 h | 3.3 GB |
| `h2_shape3` | `particle_shape: 3` | ~1.40x* | 1.37 h | 3.3 GB |
| `h3_filter8_shape3` | both | ~1.42x* | 1.39 h | 3.3 GB |
| `h4_ppc400` | `ppc: 400` | 4.00x | 3.92 h | 13.3 GB |
| **total** | | | **8.6 h** | **26.6 GB** |

\* estimate, not a measurement — the cost of cubic shape is one of the things being measured.
Disk: 90 GB free (95% used), so the matrix leaves ~63 GB.

**Metric: dT_0 in the far upstream over a fixed 30 t_ab window**, deliberately not the fitted
asymptote — that needs a long clean window and is the least-constrained quantity in the fit.
The baseline rises ~34 eV over 30 t_ab, so a 2x improvement reads as ~17 eV against a ~0.3%
noise floor: one number, no model in the way. 30 t_ab also clears the initial transient
(T passes 25 eV around t = 18).

**What is actually unknown.** Only `h4`'s mechanism is textbook (heating ~ 1/N_ppc, still
unverified here). For `h1` the measured 31% was a NOISE reduction, not a heating one; `h2`'s
effect and its cost are both estimates. `h3` tests whether the two cheap levers compose or
overlap.

**The bar to clear.** Holding beta_0 below 1.0 through 141 t_ab needs ~2.6x less heating;
below 0.40 (2x Table I) needs ~10x. The cheap levers plausibly reach the first, not the
second. `theta_implicit_em` at cfl 3 removes the heating *by construction* for ~2.1x
(~13 h) — cheaper than the ppc-400 package at 24.8 h. So this matrix is really asking
whether the cheap levers get close enough to make the implicit route unnecessary; if the
answer is no, implicit is the better buy and the matrix has cost 8.6 h to establish that.

Also plumbed: `numerics.filter_npass` now renders `warpx.use_filter` +
`warpx.filter_npass_each_dir` and round-trips through `key_params` (absent == 1, the WarpX
default, not off). `queryArrWithParser` reads exactly `AMREX_SPACEDIM` values
(`WarpX.cpp:843`), so in 1D the deck carries ONE integer, not three.

### The real objection to R1_paper_470eV: the magnetic barrier arrives fully formed

Prompted by a visual comparison of the two runs' ion phase spaces. R1_paper goes through the
expected stages -- sweep-up, then a distinct reflected population -- while R1_paper_470eV
accelerates a large fraction of ambient essentially from t = 0 without sweeping it up first.

I first attributed this to the criteria (`"6_reflected_ions": G > 0.0` on a globally counted
quantity, `metrics.py:247`). **That was wrong.** The criterion IS too loose, but the effect is
physical. Coherent magnetic compression, B_perp smoothed over 0.2 d_i0 so this is signal not
noise:

| t*wci0 | R1_paper | R1_paper_470eV |
|---|---|---|
| 0.00 | 1.00 | 1.00 |
| ~0.14 | 3.50 | **10.95** |
| ~0.29 | 4.91 | **10.82** |
| ~0.44 | 6.06 | **11.06** |

**R1_paper builds its barrier 3.5 -> 4.9 -> 6.1; the 470 eV run has ~11x at the first output
frame and stays flat.** The noise-to-coherent ratio is *lower* in the 470 eV run (1.47-1.67 vs
1.66-2.54), so this is not a field-noise artefact. It accompanies 2.5x faster ambient
acceleration (G = 0.021 vs 0.008 at t*wci0 ~ 0.65).

Two hypotheses I checked and discarded:
  * **Pre-heated upstream ions.** Relative heating is similar in both (T_i x5.9 vs x4.8 over
    the same interval) and both upstreams stay cold against their own v_sh (v_th,i/v_sh = 0.027
    vs 0.035). Not it. (My first pass here used the 470 eV v_sh for BOTH runs and produced a
    spurious 10x difference -- discarded.)
  * **Field noise inflating a domain max.** Ruled out by the smoothing above.

The one controlled difference is **dz/lambda_D,amb = 6.07 vs 0.60**; both resolve d_e
identically at dz = 0.3 d_e,ab. HYPOTHESIS, not demonstrated: with the ambient Debye length
unresolved the grid cannot support the charge separation that mediates piston-to-ambient
coupling, so the field is compressed numerically at the piston edge within the first steps and
the over-strong barrier promptly reflects ambient ions.

**This is a more serious objection than the beta_0 excursion.** I had been treating grid
heating as the problem and beta_0 drift as its consequence; if this holds, the Debye
under-resolution is corrupting shock FORMATION from t = 0, not just the upstream state.

### Staged: `runs/implicit_phase/` — the discriminating test (NOT run)

| run | change | time |
|---|---|---|
| `i0_implicit_cfl075` | `evolve_scheme: theta_implicit_em`, Newton, villasenor, no PC | 3.02 h |
| `i1_explicit_villasenor` | villasenor deposition ONLY — the control | 0.33 h |
| | | **3.35 h** |

Both run 149,164 steps = t*wci0 0.30, only far enough to watch the barrier form, with 20 field
frames across it. **cfl stays 0.75 -- the same dt as production -- so the scheme is the only
thing that changes.** theta_implicit_em is robust to the finite-grid instability at any
dz/lambda_D, so if the instant pileup is a Debye artefact it should vanish. The fast implicit
configuration requires villasenor deposition (`ImplicitSolver.cpp:679`), which production does
not use, hence i1 isolating that change on its own.

    i0 ~= i1 ~= 11x                      -> not the scheme; hypothesis wrong
    i0 builds gradually, i1 stays ~11x   -> the implicit scheme fixes it: Debye resolution
    i1 also builds gradually             -> it was the DEPOSITION all along

Note the heat_phase matrix cannot answer this: ppc, shape order and filtering all reduce
NOISE, and none of them resolves lambda_D. If the hypothesis is right the matrix will show
modest heating improvements and no change to the pileup.

Also plumbed: `numerics.evolve_scheme`, `numerics.implicit.*` and `numerics.current_deposition`
now render and round-trip through `key_params`, with the implicit defaults set to the
configuration measured fastest in runs/opt_phase (Newton + mass matrices + skip_particle_
picard_init, and NO preconditioner -- pc_curl_curl_mlmg was 54x slower than none).

## 2026-08-05 — heating matrix h0/h1: 8 filter passes buy 12.5%, not the 2.6x needed

Two of the five pilot variants ran (`runs/heat_phase/`, 440,000 steps = 30.0 t_ab each,
both GPUs). h2/h3/h4 are staged but **paused pending approval**, so the matrix is a pair,
not a matrix, for now.

### The headline: dT_0 over the fixed 30 t_ab window

The deliberately model-free metric — far upstream (outer 5%), no fit in the way:

| variant | dT_0 | T_0 final | beta_0 final |
|---|---|---|---|
| `h0_baseline` | **29.7 eV** | 39.7 eV | 0.78 |
| `h1_filter8` (`filter_npass: 8`) | **26.0 eV** | 36.0 eV | 0.70 |

**12.5% reduction.** Holding beta_0 under 1.0 for a full-length run needs ~2.6x less
heating; under Table I's 0.2 needs ~10x. Filtering alone is not that lever.

### The asymptote looks 3x better and should not be believed

The saturation fit gives 132.2 -> 87.1 eV (beta_0 2.59 -> 1.71), a 34% drop, and h1's fit is
the tighter of the two (r = -0.943 vs -0.879). Do not quote it:
  * a 30 t_ab window ends at T = 36-40 eV while the asymptote sits at 87-132 eV, so the fit
    extrapolates 2-3x beyond its own data;
  * the SAME fit on the same physics over the production run's 141 t_ab window returned
    beta_0 = 4.33, against 2.59 here. The window length, not the physics, moved it.
This is why dT_0 was chosen as the metric when the matrix was designed.

The lambda_D = dz ceiling is **367.9 eV in both**, as it must be — that is set by the grid,
and no amount of current smoothing moves it.

### Cost: NOT measured, do not read the s/step column
h0 ran at load 18.3, h1 at load 40.7 on this shared box. Their means (0.0103 vs 0.0108
s/step) differ by less than the contention between them. The profiler puts filtering at
0.16% of runtime per pass, so h1 *should* be ~1% slower; these runs cannot confirm it.

Figures: `media/heat_phase/{h0_baseline,h1_filter8}/upstream_beta.png`.

## 2026-08-05 — implicit test moved to ONE GPU; AMReX arena must be capped to share a card

`i0_implicit_cfl075` + `i1_explicit_villasenor` reconfigured for a single card, freeing the
other for other users on this shared box.

**`max_grid_size` 15000 -> 30000** in both configs (one box for one rank). Regenerated both
decks, both round-trip verified. Justified by the measurement already in those configs:
two boxes on one GPU costs 3.3% over one box.

### The arena trap
AMReX allocates **3/4 of TOTAL device memory** at init regardless of need. Measured on
`h1_filter8` (`run.log`, "Device Memory Usage"):

    [The Arena] max space allocated   8904 MB     <- 0.75 x 11873 total
    [The Arena] max space USED         189-204 MB <- the actual footprint
    [The Pinned Arena] max used         74 MB

With another user holding 9.1 GB on both cards (2747 MB free), the default allocation
**fails at init** though the run needs ~300 MB. Fix, passed as a runtime ParmParse arg via
`launch.sh -- ` rather than baked into the deck (it is a machine-sharing knob, not physics,
and it still lands in `warpx_used_inputs`):

    scripts/launch.sh -g 1 runs/implicit_phase/i0_implicit_cfl075 \
        -- amrex.the_arena_init_size=1610612736      # 1.5 GB

Result: 1736 MB resident, coexisting with the other job.

### Implicit loses much less to halving the cards than explicit does
Measured **0.0807 s/step on one contended card** -> ~3.3 h for 149,164 steps, against the
~3.0 h projected for two free cards. The explicit scheme gets 1.77x from the second card;
the implicit scheme is Newton-iteration-bound rather than bandwidth-bound and would not.
Newton converges in **2 iterations to rel 4.6e-10** at cfl 0.75, so the runs/opt_phase
solver retune holds at this timestep.

## 2026-08-05 — implicit test NEGATIVE; and the "11x instant barrier" was a wall spike

`i0_implicit_cfl075` (theta_implicit_em) and `i1_explicit_villasenor` (deposition-only
control) both ran to 149,164 steps (t*wci0 = 0.30) on ONE GPU, zero failure signatures,
Newton at 2 iterations to rel 4.6e-10 throughout. i0 3h19m at 0.0802 s/step; i1 0h43m at
0.0173 s/step. i1's progress logger died before writing its DONE marker -- completion was
verified from `run.log` (`STEP 149164 ends`), not the progress log.

### Result: neither the scheme nor the deposition

Pre-registered three-way read (staged 2026-08-05, both configs' headers):

    i0 ~ i1 ~ 11x                     -> not the scheme; hypothesis wrong   <- THIS ONE
    i0 builds gradually, i1 stays 11x -> implicit fixes it: Debye resolution is the cause
    i1 also builds gradually          -> it was the DEPOSITION all along

max coherent |B_perp|/B0 at t*wci0 = 0.29, smoothed over 0.2 d_i0:

| run | wall spike | ramp (z > 0.5 d_i0) |
|---|---|---|
| `R1_paper` | 4.78 | 4.25 |
| `R1_paper_470eV` | 10.70 | 5.98 |
| `i0_implicit_cfl075` | 12.19 | 5.51 |
| `i1_explicit_villasenor` | 10.97 | 6.88 |

Exact energy conservation at theta = 0.5 does NOT remove the effect, and neither does
Villasenor deposition. **The finite-grid instability is not the cause.** Whatever separates
R1_paper from R1_paper_470eV is upstream of both the field solve and the deposition.

### CORRECTION 1: the 11x was a boundary feature, not a shock ramp

The 2026-08-05 entry above reported R1_paper_470eV at ~11x coherent B_perp against
R1_paper's 3.5-6.1, and called it the major finding. That comparison took the max over the
WHOLE domain, and a non-propagating feature at z ~ 0.10 d_i0 owns it in every run. At
t*wci0 = 0.14 the shock should be at M_A * 0.14 = 2.3 d_i0, so the maximum was never the
ramp. Excluding z < 0.5 d_i0, the propagating structure differs by **~1.4x, not ~2.3x**
(5.98 vs 4.25 at t*wci0 = 0.29).

Corroborating: the wall feature's FWHM is 0.25-0.47 d_i0 in the 470 eV runs against
2.5-4.7 d_i0 in R1_paper. At 0.2 d_i0 smoothing, 0.25 d_i0 is the kernel width itself --
an unresolved spike, not a structure. The noise-ratio check that was supposed to catch
artifacts did not catch this one: it tests noise-vs-coherent, not where the max sits.

### CORRECTION 2: it is not "instant" either

With 22 field frames instead of three sampled times, the rise is fully resolved:

    t*wci0    0.00  0.014  0.029  0.043  0.057  0.086  0.114  0.143
    R1_paper  1.00   2.17   2.83   2.77   2.84   3.11   3.57   3.58
    470eV     1.00   2.27   2.80   3.52   5.12   7.82   9.92  11.70
    i0        1.00   2.23   3.08   4.08   5.35   9.01  11.29  12.30

All three share the same initial climb to ~2.2 by t*wci0 = 0.014 and diverge only after
~0.03. "~11x at the first output frame and stays flat" came from sampling three widely
spaced times in a run that had 1115 frames on disk.

### Method now reproducible
`scripts/plot_bperp_pileup.py` replaces the ad-hoc script that produced the retracted
numbers. It still defaults to a domain-wide max -- **add boundary exclusion before reusing
it**; the wall spike is the default answer otherwise.

### What survives
The user's original phase-space observation -- R1_paper passes through the stages of shock
formation while R1_paper_470eV accelerates ambient immediately -- is untouched by any of
this. The B_perp story built to explain it was measuring the wrong thing. The open question
is now the strong, narrow, NON-PROPAGATING field feature at the piston boundary, plus a real
but modest 1.4x difference in the actual ramp.

## 2026-08-06 — E_z at early times; two `make_thomson.py` bugs; 470 eV fig7 was the default-times imposter

### `make_thomson.py`: plotfile sort, output directory, and a cache that hid the fix
1. **Sort.** `pic_thomson._warpx_plotfiles` orders with a plain `sorted()` on the directory
   name. Past 1e6 steps `diag11002384` (step 1002384) sorts BEFORE `diag1111376`
   (step 111376), so every frame after the millionth was read out of order -- a scrambled
   time axis on both spectrograms, no error, a plausible-looking figure. `R1_paper_470eV`
   (2,784,400 steps) is the first run here to trip it; the identical bug hit
   `kinshock.io.plotfiles` on 2026-08-05. Fixed in OUR script, not the vendored fork, by
   passing the numerically-sorted permutation as `timesteps` (documented as "indices into
   the sorted list", and the reader preserves the order it is given). `numeric_timesteps()`
   returns None when lexicographic order is already correct, so short runs are untouched.
   **The fork is still wrong for any other caller.**
2. **Output directory.** `out_dir` was built by hand as `media/<run_id>`, escaping the
   `media/<phase>/<run_id>` mirror every other figure uses. Now `P.media_dir()`.
3. **The cache would have hidden fix 1.** The reader's cache signature does NOT include
   `timesteps`, so a cache built before the fix is reused verbatim and silently defeats it.
   Reordered reads now get their own `<species>_ordered.npz` rather than deleting anything
   already paid for.

Verified: `thomson_spectra{,_scaled}.npz` now have strictly monotonic `t` over 51 frames.
Regenerated for `R1_paper_470eV`, both unscaled and `--velocity-scale-factor physical`
(R = 18.36, v/4.285): alpha_epw 0.67-4.46 unscaled, 2.88-19.2 scaled -- i.e. the unscaled
spectra are partly SUB-COLLECTIVE (alpha < 1) while the scaled ones are collective
throughout. IAW doublet +/-2.3 nm unscaled, +/-0.5 nm scaled.

### fig7: `R1_paper_470eV` was another default-times imposter
`media/R1_phase/R1_paper_470eV/` had `shock_fig7.png` but no `shock_fig7_rho_i0.png` --
exactly the tell documented at RESULTS 2026-07-27 and hit again on R1_paper 2026-08-03.
Regenerated on the convention times with both x-unit variants; it now picks up the by-eye
fit automatically (v_sh = 0.0165c, M_A = 16.50).

**The 470 eV run's particle cadence cannot resolve early formation.** 51 frames over
t*wci0 = 0..5.60 is 0.11 spacing, so the first fig7 panel snaps 0.15 -> 0.11, and only
THREE particle frames exist below t*wci0 = 0.3 (0.00, 0.11, 0.22). The field diagnostics
have 1115 frames over the same span but carry no particles (`write_species = 0`). Early
phase space needs a rerun with the cadence inverted, not re-plotting.

### `plot_ez.py`: `--tmax`, `--zmax`
Reading STOPS at `--tmax` rather than loading all 1291 frames and slicing, so an early
window is cheap. `--zmax` clips the spatial axis; the two footer rms figures are still
computed over the WHOLE domain, because they are labelled "ALL z" and "far upstream" and a
zoomed number under those labels would misreport what was measured.

**Normalizations are NOT comparable across the two unit systems.** v_sh*B0 is 3.01e9 V/m
for R1_paper against 3.48e7 for R1_paper_470eV -- 87x, from B0 10x and v_sh 8.7x. The
apparent "2.0 vs 23.0" rms E_z is mostly denominator, not physics.

## 2026-08-06 (addendum) — Thomson panel restructured; the convention flags for `R1_paper_470eV`

### The command — a bare re-run reverts the windows, exactly like fig7
```
python scripts/make_thomson.py runs/R1_phase/R1_paper_470eV \
    --velocity-scale-factor physical --notch 525 541 --epw-max 625 --iaw-halfwidths 8
```
(read stage first under `physics`, this is the model stage under `tsnn`.)

Unlike fig7, the artifact now **records its own settings**: `thomson_spectra*.npz` carries
`epw_notch_nm`, `epw_max_nm`, `iaw_halfwidths`, `velocity_scale_factor`, and the
`species_fraction` line-out. A default-window file is now identifiable by its contents
instead of by mtime archaeology.

### Why each flag
* **`--epw-max 625`** (new flag). The derived window is 2 sigma of the electron feature,
  which `--velocity-scale-factor physical` shrinks by sqrt(18.36) = 4.3x -- it was +/-26.4 nm
  and cut the satellites off. Now +/-93 nm (439-625). The satellites run ~508/557 nm early
  to ~487/575 nm at 200 ps, then break into broad structure after 280 ps that the old
  window did not contain at all.
* **`--iaw-halfwidths 8`** (+/-4.3 nm, was the 2.5 default = +/-1.3 nm). The window is sized
  from C_s, but after ~170 ps the bulk drift far exceeds C_s and carried the feature out of
  frame: **30 of 51 frames had their peak PINNED to the window edge**, which read as a flat
  bright band along the bottom and is a clipping artifact, not physics. Now 0/51 pinned,
  median 0.04% of power in the outer 3 bins. What the clipping was hiding: a resolved
  doublet before ~170 ps, a sharp blue shift to ~529 nm at 200-240 ps, then recovery to
  ~530 and a slow red drift back toward line centre.
* **`--notch 525 541`** (default is exactly the IAW window). Verified blanked: mean AND max
  are 0.0 across the band. 541 rather than 535 because residual power sat at 535-540 with
  the brightest bin at 539.2 nm.

### Middle panel: species fraction, not per-timestep normalisation
`--no-species-fraction` restores the old panel. The new one is the piston density fraction
**in the probe volume** (a line-out at z/d_i0 = 40.2, no spatial axis -- the spectra come
from one position, so that position's is the only comparable information).

**It dates the piston arrival, which the spectra alone cannot:** the fraction is ~0.00 until
~230 ps, then rises almost vertically to 1.00 by ~260 ps, crossing 50% at **245 ps**.

That timing changes the reading of the IAW panel. The blue shift BEGINS at ~200 ps while the
probe volume is still pure ambient, and the piston does not arrive for another ~45 ps -- so
the first thing the diagnostic sees is **the shock foot pushing ambient ions ahead of the
contact**, not piston material. The two-stage structure (shift at 200-240, settled band at
~530 nm after) maps onto compressed-ambient-then-piston.

Contact speed from the 0.5 contour: **0.0133 c**, against the by-eye shock fit of 0.0165c.
Different surfaces -- contact discontinuity vs shock front -- and a shock outrunning its
driver by ~25% is the expected ordering, so this is a consistency check that passes.

Cost note: the alpha panel is unchanged, and the run is collective throughout under the
scaling (alpha 2.87-19.2). The UNSCALED spectra dip to alpha = 0.67, i.e. partly
sub-collective, where the IAW panel shows the electron feature and bulk drift rather than a
true ion feature -- read the reported alpha before interpreting either.

## 2026-08-06 — h2 NEGATIVE. The 470 eV upstream noise is the parameter point, not a numerical artifact

`h2_shape3` ran 440,000 steps (30.0 t_ab) in 1h21m on two GPUs, clean, 23 field frames.
It was staged as the decisive test of an aliasing hypothesis. It refuted it.

### Result: cubic splines do essentially nothing

| variant | dT_0 | E_z rms upstream (t*wci0 = 0.08) | grid-scale part |
|---|---|---|---|
| `h0_baseline` (shape 2) | 29.7 eV | 30.51 | 18.65 |
| `h1_filter8` | 26.0 eV | 19.55 | **5.29** |
| `h2_shape3` | **28.0 eV** | **30.00** | **18.27** |
| `R1_paper` (reference) | -- | 2.25 | 0.85 |

E_z rms moves **1.7%**, grid-scale noise 2%. h2 is also WORSE than h1 on dT_0 (5.7% vs
12.5%). Filtering does cut grid-scale noise (72%) yet bought only 12.5% on heating, so the
heating is not dominated by the grid-scale component either -- the aliasing story fails from
both directions.

All E fields here are normalised by **v_A * B0**, configured quantities. Never the fitted
v_sh: that is the thing under investigation (see Retraction 1).

### The actual diagnosis: E_th/(v_A B0) = beta * c/v_ti

The thermal electrostatic scale is E_th = sqrt(n T / eps0), which gives the exact identity
**E_th/(v_A B0) = beta * c/v_ti**.

```
              beta_0   T0[eV]   c/v_ti   beta*c/v_ti   measured   meas/pred   N per lambda_D
R1_paper      0.2000   1022.0    223.6       44.7        2.25       0.050        166.7
h0_baseline   0.1957     10.0   2260.5      442.5       30.51       0.069         16.5
```

Both runs sit at the **same 5-7% of their own thermal noise scale**. Nothing is anomalous
about the 470 eV run's noise physics. beta_0 is identical by construction, but c/v_ti is 10x
larger, so identical relative noise is 10x larger measured against the MHD field scale.

**This is the reduced-c trick being undone.** R1_paper reaches the parameter point with a
reduced c (c/v_ti = 224); the 470 eV rebuild uses real c at real temperatures and gets
c/v_ti = 2261. The upstream electrostatic noise floor relative to v_A B0 rises by exactly
that ratio, and it drives ambient ions before any shock exists.

This is the same trap CLAUDE.md already records for lambda_ab: the paper's §II says the
scaling "ensures that dimensionless quantities are correct, but [collisionality] relative to
global scales is only quantitatively matched **at physical mass ratios**." The knob that
buys both is mu_p = 1836, not a reduced c -- and not a numerical lever.

### No cheap lever closes a 10x
Noise ~ 1/sqrt(N_lambda_D), and there are **16.5 particles per Debye length against
R1_paper's 167**. Recovering 10x needs ~100x ppc, or ~6x finer cells (~36x cost with the
CFL). **This retires h4** (`ppc: 400` buys sqrt(4) = 2x) and h3 (filter+shape, and shape
contributes nothing). The remaining choices are to return to reduced-c for this physics, or
to accept the 470 eV run as upstream-noise-dominated and not comparable to R1_paper on early
shock formation. That is a physics-setup decision, not a numerics one.

### RETRACTION 1 — "there is no false shock" was circular
Earlier today I compared the two runs at equal SHOCK TRAVEL (z_sh = M_A * t*wci0), found the
accelerated-ambient fraction collapsing to a few percent, noted the time offset (1.159)
matched the M_A ratio (1.154) to 0.4%, and concluded the false shock was a normalisation
artifact. **That is circular.** Both runs are configured with the identical dimensionless
driver, v_p/v_A = 10.4. M_A is therefore an OUTPUT that should agree; measuring 16.50 against
14.30 IS the anomaly. I divided by the signal and reported that it vanished.

The reflected-ion fraction is also the wrong diagnostic -- it is a ratio that normalises out
the effect. The right quantities are the ambient phase-space STRUCTURE and the B and E
fields, which is where the 13x upstream E_z shows up immediately.

### RETRACTION 2 — the aliasing hypothesis, and an over-read of i0/i1
Having found the E_z excess, I attributed it to finite-grid aliasing at dz/lambda_D = 6.07
and predicted `particle_shape: 3` would cut it. h2 shows it does not (1.7%).

I also over-read the implicit result. i0/i1 showed theta_implicit_em and Villasenor
deposition change nothing -- correct, and they do rule out energy non-conservation. They do
NOT rule out a Debye-resolution cause, which is what I claimed at the time. Corroborating:
i0 and i1 reach T_up = 12.50 vs 12.46 eV, indistinguishable.

Three hypotheses tested and rejected in order -- energy non-conservation (i0/i1), aliasing
(h2), and "it is only a normalisation" (retraction 1). What survives is the parameter point.

### Cost: still not measured
h2 ran at 0.0110 s/step on an IDLE box against h0's 0.0103 under load 18.3 -- so cubic
splines look ~free here, but the two differ in load as much as in shape order. The cost of
`particle_shape: 3` remains unmeasured on this deck.

---

## 2026-08-11 — Two resolution sweeps staged: `H_phase` (heating) and `S_phase` (early shock)

`R1_paper_470eV` is settled as the baseline. Two defects remain, and both are suspected to
be under-resolution: **(a)** a shock that forms too early, via an upstream E-field pileup,
and **(b)** numerical heating of the upstream. This entry stages a separate sweep for each,
benchmarks every point, and prices the production rerun. **Nothing has been launched.**

### The framing that shapes both sweeps: two knobs, not one

Every previous mitigation attempt varied one knob at the production resolution and failed —
`h1_filter8` (12.5 %), `h2_shape3` (1.7 %), `i0`/`i1` (indistinguishable). What was never
varied is the pair that actually sets the physics of an under-resolved PIC plasma:

- **dz/λ_D,amb = 6.07** — aliasing, which drives the finite-grid instability. `R1_paper`,
  which behaves, sits at 0.60.
- **N_D = ppc·λ_D/dz = 16.5** — particles per Debye length, which sets the *amplitude* of
  the thermal-fluctuation noise (∼1/√N_D). `R1_paper`: 167.

Refining dz improves **both**; raising ppc improves **only** N_D. A one-dimensional dz scan
therefore cannot say which one matters — and the answer sets the price, because at fixed
N_D **ppc is 4× cheaper than dz** (cost ∝ ppc, but ∝ dz⁻² once dt follows the CFL). Both
sweeps are built as grids where lines of constant N_D cross lines of constant dz/λ_D, so
the two effects separate. Designs, grids and readings: `runs/H_phase/README.md`,
`runs/S_phase/README.md`.

### `H_phase` — 8 runs, uniform ambient box, no piston

A periodic 2.06 d_i0 box of ambient plasma and nothing else: no piston, no heater, no
injector, no shock. Everything else byte-identical to `R1_paper_470eV`, including the
collision block. Window 440 000 steps = 30.0 t_ab — **the same window `h0_baseline`
measured (+29.7 eV)**, so `hs_dz1_ppc100` is both the production parameter point and the
control.

With no piston and no injector the domain-wide mean energy per particle *is* the
temperature, so `EP`/`PN` carry the whole measurement — free, ~1000 rows, no plotfile
binning. This is precisely what `grid_heating.py`'s docstring says the reduced diagnostic
cannot do in the full deck. New `scripts/heating_rate.py` reads it. Verified on the
benchmark output: 153 600 macroparticles exactly constant, T = 10.01 eV for both species at
step 0.

### `S_phase` — 6 runs, R1_paper_470eV truncated to 12.0 d_i0 and t·ω_ci0 = 0.30

Identical physics to the production run — same piston, heater, injector, collisions, B₀ —
on 1/6.7 of the domain and 1/19 of the window. By t·ω_ci0 = 0.30 the piston has gone
3.12 d_i0 and the model shock 4.19 d_i0, so there is 2.9× head-room; the anomaly is visible
by 0.39, so it starts inside the window.

⚠ **The truncated domain is not free.** `field_hi = open` is pec, which *reflects* fields.
Light crosses 12 d_i0 in ~6000 steps, so the turn-on precursor makes ~25 round trips here
against ~4 in the full domain. `ss_dz1_ppc100` is therefore a **control, not just the
cheapest point**: it must reproduce the full run over t·ω_ci0 ≤ 0.30 before anything else
in the phase means anything. The full run's `diag_fields` gives ~60 frames in that window,
so the comparison costs no compute.

### Cost, measured (2026-08-11, 1× RTX 4070, 1200-step benchmarks)

`s/step` from WarpX's cumulative `Evolve time` between steps 200 and 1200, so init drops
out. Harness validated against a known number: the production deck measured 0.01533 s/step
against 0.01415 recorded 2026-08-04.

| | serial, 1 GPU | both cards | output |
|---|---|---|---|
| `H_phase`, 8 runs | 5 h 43 m | **2 h 51 m** | 3.2 GB |
| `S_phase`, 6 runs | 3 h 37 m | **1 h 48 m** | 17.7 GB |
| **both sweeps** | 9 h 20 m | **4 h 40 m** | 21 GB |

Cheapest first: `ss_dz1_ppc100` (the control) is **8 min**; the whole constant-N_D aliasing
line in `H_phase` — `dz1_ppc100`, `dz2_ppc50`, `dz4_ppc25` — is **48 min together**.

**Production rerun**, same 9000 d_e,ab domain and 2 784 400-step window, scaled from the
measured anchor with the ×1.279 growth factor. Cost goes as **k²·(ppc/100)**:

| dz/λ_D | ppc | N_D | 2 GPUs |
|---|---|---|---|
| 6.07 | 100 | 16.5 | **8.6 h** ← what already ran |
| 3.03 | 50 | 16.5 | 17.1 h |
| 6.07 | 400 | 65.9 | 34.3 h |
| 3.03 | 100 | 33.0 | 34.3 h |
| 1.52 | 25 | 16.5 | 34.3 h |
| 1.52 | 100 | 65.9 | 5.7 d |
| 0.76 | 50 | 65.9 | 11.4 d |

**The three 34.3 h rows are the decision.** The same budget buys N_D 16.5 → 66 with
aliasing untouched, *or* aliasing 2× better at N_D 33, *or* aliasing 4× better with the
noise untouched. Nothing measured so far distinguishes them. ~9 h of sweep decides how
34 h — or 5.7 d — of production gets spent.

Device memory is never the constraint: **331 MB for the production deck's 6×10⁶ particles**
(55 B/particle), so even dz/8 at ppc 100 is 3.4 GB. Cap `amrex.the_arena_init_size` anyway
when sharing a card — AMReX allocates 3/4 of *total* device memory at init regardless.

### `max_grid_size = n_cell` costs 7.3× on CPU — measured, not inferred

Every sweep config sets `numerics.max_grid_size = n_cell` (one box), which is right for a
single GPU and required by `launch.sh -g`. On CPU it is a trap: AMReX does not tile in 1D,
so one box is one tile and only one OpenMP thread gets work. Measured directly on
`hs_dz1_ppc400` at 8 threads: **0.1073 s/step with one box vs 0.0148 with
`max_grid_size = 64`.** The existing note ("neutral-to-negative on CPU") understates this
badly for a *one-box* decomposition. To run any sweep point on CPU, delete the key from
`config.yaml` and regenerate. Even then the GPU is ~6× faster on these decks.

**And the trap is already live in `R1_paper_470eV`.** That config gained
`max_grid_size: 15000` on 2026-08-04 *for GPU* — two boxes, which on 8 threads measures
**0.3924 s/step against the 0.11169 baseline in `runs/opt_phase/SUMMARY.md`, 3.5×.** The
opt_phase CPU rows predate that key, so **every CPU projection in that table is stale for
the config as it now stands**: the quoted 5.33 d at 8 threads is really ~18 d, and 2.93 d
at 20 threads is not attainable without editing the config back. The GPU rows are
unaffected and remain the ones to use. This was not a hypothesis carried over from the
sweep decks — it was measured on the production deck itself with the same harness.

### Code

`kinshock.deck.render` now treats `operators:` as **optional** — a run with no piston
species emits no heater and no injector. This is not cosmetic: `particle_heater`'s rate is
H ∼ 1/foil_width (`ParticleHeater.cpp:207`), and `slab_halfwidth_di: 0` — which is what
makes the ambient profile uniform — would be a division by zero, not a no-op.
`deck.key_params` guards both operator blocks so `--verify` still catches an operator that
silently vanishes. All 26 pre-existing decks re-render byte-identically (`--check`) and
`tests/` is 15/15.

---

## 2026-08-11 (later) — Controls run. The upstream heating is not grid heating, and both sweeps' premise is falsified

Launched the two cheapest points of the new sweeps as controls before committing to the
rest. The `H_phase` control **failed**, and chasing why produced the result below. Total
cost: 15 minutes of GPU.

### `hs_dz1_ppc100` — the production parameter point, quiescent, heats by nothing

A periodic 2.06 d_i0 ambient box at **identical** dz/λ_D = 6.07, N_D = 16.5, collisions,
B₀, n, T₀, dt, particle_shape and filtering to the production run, over the identical
30 t_ab window:

| | ΔT over 30 t_ab |
|---|---|
| `h0_baseline`, full piston-driven deck, outer 25 %, amb_electrons | **10.0 → 40.7 eV** |
| `hs_dz1_ppc100`, quiescent box, domain-wide | **10.00 → 9.98 eV** |

Macroparticle count exactly constant (153 600), so the mean energy really is a temperature.

**The box is not inert, and this is the key check:** its rms E_z settles at **29.3–30.1
v_A B₀**, which is precisely the production run's far-upstream value (30.2, RESULTS
2026-08-05) and 6.8 % of the exact thermal scale E_th/(v_A B₀) = β·c/v_ti = 442.5. So the
box reproduces the production upstream's field noise **exactly** — and that noise does not
heat the plasma. **The E_z excess and the upstream heating are two different problems.**

### h0's own data says the same thing, independently

Spatially and temporally resolving `h0_baseline` (existing plotfiles, no new compute):

- Heating is **uniform** across the far upstream — 41.9, 41.3, 41.1, 41.0, 41.0, 41.1,
  40.3 eV at z = 28…76 d_i0 at t = 30 t_ab — while the piston has reached only ~16 d_i0.
- At z = 70–80 d_i0 the temperature is **flat at 9.98 eV until t = 2.50 t_ab**, then rises
  linearly at 1.24 eV/t_ab. **Light reaches z = 76.5 d_i0 at t = 2.59 t_ab.**
- **No piston particle ever gets there**: 0 macroparticles beyond 60 d_i0 at t = 30 t_ab,
  and thermal transit at 470 eV would take 85 t_ab.
- The heated distribution is a **perfect undrifting Maxwellian** (v/v_th percentiles match
  to 1 %, mean drift 0.002 v_th).
- Mean |B_perp|/B₀ in that window rises **1.00 → 2.11**, starting at the same t = 2.6.
  ⚠ Checked with a 0.2 d_i0 boxcar, not a raw max — the raw max runs 2.2–2.8× the smoothed
  value, which is exactly the trap behind the retracted "11× instant barrier".

So the far upstream is heated by an **electromagnetic disturbance arriving at c**, not by
the finite-grid instability, and not by any particle that got there.

### And it is the same c/v_ti that explains everything else

`R1_paper` and `R1_paper_470eV` have the same geometry in metres, the same dz, the same
pec boundaries, the same dimensionless piston — and t_ab differing by 10×. So the domain's
light-crossing time differs by 10× **in units of t_ab**, and each run's far upstream is
disturbed exactly at its own:

| run | T₀ | c/v_ti | L/c [t_ab] | measured onset |
|---|---|---|---|---|
| `R1_paper` | 1022 eV | 224 | **25.9** | 22.5–26.3 |
| `R1_paper_470eV` | 10 eV | 2261 | **2.59** | 2.5–3.8 |
| `h0_baseline` | 10 eV | 2261 | **2.59** | 2.50–5.00 |

`R1_paper`'s far upstream sits at |B_perp|/B₀ = 1.0000 and E_z = 2.15 v_A B₀, both dead
flat, until ~25 t_ab. It is not a better-resolved run; it is the same run with the
precursor clocked 10× slower.

**Both defects therefore trace to one quantity — c/v_ti, 10× too large because the
reduced-c trick was undone** (RESULTS 2026-08-05 derived this for E_z; it turns out to set
the precursor crossing time too). **Neither is a resolution problem: no dz and no ppc
changes c/v_ti.**

### Consequences for the sweeps

- **`H_phase` is superseded — do not run the remaining 7 points.** Its premise was that
  grid heating at dz/λ_D = 6.07 warms the upstream. The baseline point measures zero. The
  other points would measure zero more precisely for 5 h 36 m.
- **`S_phase` survives, with its role changed.** Its control passed (below), and it now
  tests a *prediction*: if the early ion acceleration is driven by c/v_ti and not by
  resolution, then refining dz and raising ppc will **not** remove it. That is falsifiable
  and costs 1 h 48 m on two cards.

### `ss_dz1_ppc100` — the domain control, marginal pass

Truncated 12.0 d_i0 domain vs the full 80.5 d_i0 run over t·ω_ci0 ≤ 0.30, on a fixed
8–12 d_i0 upstream window (`scripts/check_domain_control.py` — plot_ez.py's "outer 5 %"
would compare 11.4–12.0 d_i0 against 76.5–80.5 d_i0, i.e. different plasma):

| t·ω_ci0 | rms E_z, trunc / full | coherent \|B_perp\|, trunc / full |
|---|---|---|
| 0.02 | 31.1 / 29.8 = **1.04** | 1.46 / 1.27 = 1.15 |
| 0.10 | 35.8 / 32.4 = **1.11** | 4.16 / 2.99 = 1.39 |
| 0.20 | 44.3 / 34.8 = **1.27** | 6.00 / 5.35 = 1.12 |
| 0.30 | 49.3 / 39.8 = **1.24** | 6.75 / 6.09 = 1.11 |

4 % at t = 0.02, drifting to ~24 % by 0.30 — the direction and growth expected if pec
traps precursor energy in a 6.7× smaller box. **Usable for relative comparison across the
sweep** (every point shares the domain), **not for quoting absolute upstream E_z** against
the full run. Both runs do form a real barrier (|B_perp|/B₀ 1.3 → 6), so the truncation
does not suppress the physics under test.

### Open, and not claimed

Whether the precursor is **physical** (a real ablation plasma does radiate) or an
**artifact** — of the heater operator's abruptness, or of the pec far boundary trapping
radiation that should escape — is **not established here**. Silver-Mueller is unavailable
while a background B is set (`kinshock.deck._BC_MAP`), so the cheap discriminator is domain
size: if it is trapped radiation, the far-upstream heating rate should fall as the domain
grows.

### Code

`resolve_constants` gained `allow_unresolved`, and `verify` uses it for the file under
test. WarpX prunes unused `my_constants` from `warpx_used_inputs` **by name, not by
dependency**, so in the piston-free deck `slab = 0.*di` survives while `di` is pruned —
`--verify` then raised `ValueError: could not resolve my_constants: ['slab']` on every
piston-free run, defeating the "config = what was simulated" guarantee. Decks we generate
are still resolved strictly. Confirmed the check is not weakened: a tampered `max_step` in
`warpx_used_inputs` is still caught. All 40 decks `--check` clean, tests 15/15.

---

## 2026-08-11 (later still) — `boundary.reflect_symmetry_axis` has never done anything

Found while checking what a Perlmutter build would need. **27 runs ask for the π-rotation
symmetry wall and none of them got it.**

`boundary.reflect_symmetry_axis` is a fork-only input that exists **only on
`warpx-cda` branch `feature/reflect-symmetry-axis`, which is not merged into anything we
have ever built.** Neither the OMP nor the CUDA binary contains the string. WarpX parses
the line, never uses it, and says so — every affected `run.log` carries

```
Unused ParmParse Variables:
  [TOP]::boundary.reflect_symmetry_axis(nvals = 1)  :: [x]
```

Confirmed directly: running with `amrex.abort_on_unused_inputs=1` aborts on exactly that
key. So `_BC_MAP`'s `("pec", "reflecting")` is all that took effect and the z=0 wall has
been **plain specular reflection** — which flips only v_z, not the gyro-coupled v_perp.
CLAUDE.md puts that at a ~5 % near-wall artifact (RESULTS 2026-07-23).

**`symmetry` and `reflecting` have therefore been the same simulation.** `R0_half` and
`R0_half_sym` differ in `key_params` by *nothing*; their decks differ by one ignored line.
Any comparison between them measured RNG, not a boundary condition.

Affected: `R1_paper`, `R1_paper_470eV`, `R1_paper_470eV_pilot`, `R1_warm`, `R1_coll`,
`R1_cal`, `R1_recal`, `R1_half`, `R1_core_half_sym`, `R1_paper_phys`, `R1_paper_dial`,
`R0_half_sym`, `R2`, `R3`, all of `heat_phase`, both of `implicit_phase`, and all six
`S_phase` configs.

### Why --verify passed it 27 times, and the fix

`deck.verify` compares config against `warpx_used_inputs` — and **AMReX omits unused
variables from that file**, so an unimplemented key is absent from *both* sides and
matches perfectly. The guarantee "config = what was simulated" had a hole exactly the
shape of a missing feature.

`make_inputs.py --verify` now also reads the run's own `run.log` and reports anything
under "Unused ParmParse Variables" as a MISMATCH (`_unused_parmparse`). `my_constants.*`
are excluded: they are declarations rather than directives, and the piston-free `H_phase`
box deliberately carries unused ones (it keeps `plasma.piston` so `theta_e_heat` can define
t_ab while no species has `role: piston`). Verified precise — clean on `hs_dz1_ppc100`,
MISMATCH on `ss_dz1_ppc100`, `R1_paper_470eV` and `h0_baseline`.

### Not yet decided

The two commits that add it (`d5f2e9917`, `05d74af41`) are self-contained: 7 files,
+220/−0, and they merge cleanly onto `feature/hybrid-laser`. Whether to build **with** them
is a scientific call, not a cleanup: turning the correct wall on makes new runs
**non-comparable to every existing result**, including `R1_paper_470eV` — which is the
reference `ss_dz1_ppc100` exists to reproduce. Deferred to the user.

---

## 2026-08-11 (Perlmutter) — queue characteristics, measured

`sbatch --test-only` on Perlmutter (validates + estimates, submits nothing), for the
single-GPU resource spec `S_phase` needs (`-C gpu -N 1 -n 1 -c 32 -G 1`). Account
`m5032_g`; `$PSCRATCH = /pscratch/sd/h/hhelal`.

| QOS | partition | resources | max wall | max jobs/user | would start |
|---|---|---|---|---|---|
| `shared` | `shared_gpu_ss11` | 1 GPU / 32 cores | 2 d | 5000 | **+11 h** |
| `debug` | `gpu_ss11` | whole node | **30 min** | **5** | **+5 h** |
| `regular` | `gpu_ss11` | whole node | 2 d | 5000 | **+6 days** |

**`shared` works for GPU jobs** — the association carries `gpu_shared`, and `-q shared -C
gpu` resolves to partition `shared_gpu_ss11` with exactly the 1 GPU / 32 cores requested.
That settles the open question in `perlmutter/README.md`; `SWEEP_QOS=shared` is right.

`regular` is **six days** deep and would take a whole 4-GPU node per single-GPU task,
idling three. It is not a fallback worth having for this work.

**The A/B belongs in `debug`**: four ~8-minute runs sits exactly inside debug's 5-job and
30-minute limits and starts ~6 h sooner. The sweep cannot — `ss_dz4_ppc25` (~33 min) and
`ss_dz4_ppc100` (~1 h 42 m) exceed the walltime cap and six tasks exceed the job cap.
`submit.sh` gained `--qos`/`--time` overrides and refuses `--qos debug` above 30 minutes.

⚠ `--test-only` reports "will start no later than", so real starts may be earlier; the
ordering (debug < shared << regular) is the durable part, not the absolute times.

The `reflect_symmetry_axis` cherry-picks are now on `origin/feature/hybrid-laser`
(`acc2d6621..fcb48c9fe`, 8 files, +232/−0), so Perlmutter can build binary B.

**What actually happened (added on merge).** The A/B was submitted to `shared`, not
`debug`, and all four runs started and finished within ~9 minutes each — so the +11 h
`shared` estimate was a loose upper bound in practice, as the ⚠ above anticipated. The
sweep likewise ran in `shared` without a long wait. `debug` remains the right call on
paper; it was simply not needed. The `--qos`/`--time` overrides added to `submit.sh` here
are untested for that reason.
## 2026-08-11 (Perlmutter) — first NERSC execution: the A/B wall test says the wall is invisible

Answers the "Not yet decided" that closes the previous entry. The two cherry-picks were
pushed, both binaries were built on Perlmutter, and the A/B ran. **Result: at this
configuration the pi-rotation wall is not distinguishable from specular reflection above
the GPU noise floor.** Details below, including one ratio that looks like a detection and
is not.

### Bring-up — everything in `perlmutter/` worked on first submission

`perlmutter/README.md` called the first submission a shakeout. It wasn't: `site.conf` →
`build_warpx.sh` → `submit.sh` → `job.sbatch` → `run_warpx` ran end-to-end with no edits
beyond the bug in the next section. Layout is what the README prescribes — both repos on
`$PSCRATCH`, nothing in `$HOME` but the profile.

**The blocker in `perlmutter/README.md` is gone, and its wording is now wrong.** The
cherry-picks are no longer chablis-local: `feature/hybrid-laser` fast-forwarded
`acc2d6621 → fcb48c9fe`, and that tip **is** `WARPX_COMMIT_B` verbatim — the rebuilt
commit reproduced the chablis SHA exactly, so `site.conf.example` needed no change.
`d5f2e9917` also exists on `origin/feature/reflect-symmetry-axis`, but **that branch is
not usable as build B**: it forked before the heater merge and carries zero
`ParticleHeater`/`TargetInjector` files (15 behind `development`). B has to come from
`feature/hybrid-laser`.

Three places asserted the stale "chablis-local" claim — `perlmutter/README.md`, `CLAUDE.md`
(Perlmutter bullet), and `build_warpx.sh`'s closing `echo`. **All corrected later the same
day** (see the addendum at the end of this entry). Remaining "chablis" mentions in
`perlmutter/` are historical provenance (where a measurement was taken), not stale claims.

**Binaries** (A100, `AMREX_CUDA_ARCH=8.0`, 1D CUDA, double, OPENPMD+EB+QED, 535 MB each):

| tag | commit | `strings … reflect_symmetry_axis` | wall |
|---|---|---|---|
| A | `acc2d6621` | 0 hits | specular |
| B | `fcb48c9fe` | 12 hits (4 in `.rodata`) | pi-rotation |

**Environment.** Deps via upstream `install_gpu_dependencies.sh` →
`$PSCRATCH/storage/sw/warpx/perlmutter/gpu` (boost 1.82, c-blosc 1.21.1, adios2 2.10.2,
blaspp/lapackpp 2024.05.31) + venv `venvs/warpx-gpu`, which the profile auto-activates.
Three machine facts that cost time and will cost it again:
- The installer hard-codes `$HOME/src/warpx/requirements.txt` and dies there under
  `set -e`, *after* every expensive C++ build. Fixed with a symlink
  `$HOME/src/warpx → $PSCRATCH/warpx-cda`.
- **Perlmutter's default `python3` is 3.6.15** and cannot parse `from __future__ import
  annotations`. In-job this is fine (the WarpX profile loads `cray-python/3.11.5` before
  `run_warpx` calls anything), but login-node analysis needs the venv explicitly.
- The venv ships yt+numpy+matplotlib+PyYAML but not `astropy`/`plasmapy`; both added.
  Note it carries **astropy 8.0.1**, not chablis `physics`'s 8.0.0 — so a
  `make_thomson.py` cache moved between the two machines is silently rejected by design
  (the CODATA-keying note in CLAUDE.md). Unlike chablis, yt and torch coexist here, so
  the two-env split `make_thomson.py` works around does not apply on Perlmutter.

### `build_warpx.sh`'s feature self-check could never report success (FIXED)

The check that exists precisely to catch a silently-missing fork-only input was inert. It
reported **"does NOT implement reflect_symmetry_axis"** for *both* binaries, including the
one that does.

```bash
set -euo pipefail                                    # top of the script
if strings "$bdir"/bin/warpx.1d* | grep -q reflect_symmetry_axis; then
```

`grep -q` exits on the **first** match and closes the pipe while `strings` is still
streaming a 535 MB binary. `strings` dies of **SIGPIPE (141)**, `pipefail` propagates it,
and the `if` takes the else branch *because the string was found*. Measured both ways on
build B: exit 0 without `pipefail`, exit 141 with it. So the failure modes are
indistinguishable — absent feature → grep reads to EOF, exits 1 → "does NOT" (right answer,
wrong reason); present feature → SIGPIPE → "does NOT" (wrong).

Fixed by draining the stream instead: `grep -c` into a variable, `|| true` to absorb its
exit 1 on zero matches. Now correctly prints `hits=0` for A and `hits=12` for B. **The
binaries were always right; only the report was wrong** — no rebuild was needed.

Worth naming the pattern: this is the second time in two days that a *verification* step,
not the physics, was the thing that was broken (`--verify` missed the unused key for 27
runs; this check could not pass). Both failed silently and in the reassuring direction.

### The runs — 4/4 COMPLETED, and the walls genuinely differed

`submit.sh ab`, jobs `56696895` (A, array 0–1) and `56696896` (B, array 2–3), QOS
`shared` → `gpu_shared`, one A100 each. All four exit `0:0`, elapsed 8:41–8:49 (WarpX
7:58–8:07), 150 field + 30 particle frames each.

**The asymmetry that makes this a real experiment**, from each `run.log`'s
"Unused ParmParse Variables":

| binary | `boundary.reflect_symmetry_axis` | wall actually used |
|---|---|---|
| A (`ss_dz1_ppc100`) | **listed as unused** | specular |
| B (`ss_dz1_ppc100_symwall`) | **absent from the list** → consumed | pi-rotation |

Both decks emit the key — the configs are byte-identical apart from `run_id`/`deck`. The
variable is the *binary*. This is the first time in this project that `symmetry` and
`reflecting` have been different simulations.

### The measurement (`scripts/ab_wall_test.py`)

`|A−B|` over the same-binary replicate floor `max(|A1−A2|, |B1−B2|)`. Ratio ≲ 1 means the
wall is buried in GPU non-determinism (ablastr `RandomSeed.H`: fixed seed does not give
reproducibility on GPU).

| t·ω_ci0 | 0.05 | 0.10 | 0.20 | 0.29 |
|---|---|---|---|---|
| piston front z/d_i0 | 0.19 | 0.38 | 0.56 | 0.49 |
| coherent \|B⊥\|/B0 | 0.50 | 0.08 | 0.30 | 1.99 |
| rms E_z / v_A B0 | 1.09 | 0.17 | **17.59** | 1.84 |

**The 17.59 is a collapsed denominator, not a detection.** Its numerator, |A−B| = 1.263,
is unremarkable — it sits inside the range of |A−B| at every other time (0.224 → 2.029).
What is anomalous is the **floor: 0.0718**, against 0.591 / 1.342 / 1.101 elsewhere. The
floor swings **19×** across four time points while the signal swings 9×; with n = 1
replicate pair per binary each floor is a single draw, and at t·ω_ci0 = 0.20 the two A
replicates happened to land nearly on top of each other. Time-averaged, nothing survives:

| metric | mean \|A−B\| | mean floor | ratio |
|---|---|---|---|
| rms E_z | 1.04 | 0.78 | 1.34 |
| coherent \|B⊥\| | 0.238 | 0.502 | 0.47 |
| piston front | 0.0140 | 0.0290 | 0.48 |

**Verdict: the wall is invisible at this configuration.** Consequently the **~5 % near-wall
artifact CLAUDE.md attributes to the specular approximation (RESULTS 2026-07-23) is not
detectable here** — and this is the first comparison capable of detecting it, every prior
one having been between two identical decks.

Scope, honestly: one resolution, one-sided, t·ω_ci0 ≤ 0.30 (early formation), and a floor
built from a single replicate pair. A gyro-phase effect at the foil could still emerge
later or at finer dz. Four more replicates (~35 GPU-min) would stabilise the denominator
if this needs to be settled rather than indicated.

### Consequence for `SWEEP_BUILD`

Since the two walls are statistically indistinguishable on this problem, **A is the better
choice**: it buys the same physics while keeping the sweep directly comparable to all 27
existing runs and preserving `ss_dz1_ppc100`'s domain control against `R1_paper_470eV`
(itself specular). `perlmutter/site.conf` is gitignored and currently still says `B` —
change it before submitting the sweep, or accept that the control needs re-reading.

### Two traps found in the surrounding tooling

- **Labelled replicates collide in `media/`.** `run_warpx` copies the parent `config.yaml`
  into the labelled work dir unchanged, so `A2`'s config says `run_id: ss_dz1_ppc100` and
  `P.media_dir()` resolves it to the parent's directory. Rendering movies for `A2` silently
  **overwrites A1's**. Same class as the shared-`diags/` clobber `launch.sh` exists to
  prevent, one layer up in the output path. Only the two primary runs were rendered.
- **NERSC's `/usr/bin/ffmpeg` has no H.264 encoder at all** (SUSE omits the
  patent-encumbered ones; only libvpx/vp9/theora/gif/webp are present). `plotting.py:encode`
  hard-requests `-c:v libx264`, so every movie job fails here with a `CalledProcessError`
  that names no cause. `plotting.py:36` already honours an `FFMPEG` env override, so no
  code change is needed:
  `export FFMPEG=$(python -c "import imageio_ffmpeg;print(imageio_ffmpeg.get_ffmpeg_exe())")`
  (`imageio-ffmpeg` installed into the venv; ships a static ffmpeg 7.0.2 with x264). The
  yt pass had already succeeded — only the encode died — so re-encoding the existing frames
  was enough. Phase-space movies for both A and B are in `media/S_phase/*/shock_phase.mp4`
  (gitignored, regenerable), annotated with the **config model** v_sh = 0.0140 c since
  neither run has a `shock_fit.yaml`; re-render if one is ever tuned.

### Next

- ~~Decide `SWEEP_BUILD` and submit the sweep~~ — done, see addendum.
- ~~Correct the three stale "chablis-local" claims~~ — done, see addendum.
- Optional: 2 more replicates per binary to firm up the noise floor.

### Addendum (same day) — sweep submitted, docs corrected

**`SWEEP_BUILD=A`**, on the A/B evidence above, and `SWEEP_TIME` raised `03:00:00 →
04:00:00` (`shared` bills actual elapsed, not requested, so the hour is nearly free).

**Sweep submitted as job `56715249`, array `1-5` — not `0-5`.** Index 0 is
`ss_dz1_ppc100`, which the A/B already produced *on binary A*, i.e. exactly what the sweep
asks for; its `diags/` would in any case have tripped `run_warpx`'s in-place overwrite
guard. Re-running it would have burned ~9 GPU-min to reproduce data already on disk — and
not bit-identically, since GPU runs are not reproducible. The sweep's baseline point is
therefore the A/B's A1 run.

Estimates, anchored to the measured cost unit: 17 min / 33 / 33 / 33–50 / **2 h 09 m**
for indices 1–5; wall-clock ≈ the longest, ~3.9 GPU-h total.

**Per-GPU speed — a claim made here prematurely, and retracted below.** From the cost-1
point alone (479 s/unit vs the ~459 s/point chablis's estimate implied) it looked as though
the A100 gave *no* per-card gain, and `README.md`/`CLAUDE.md` were edited to say so. **That
was wrong**, and the sweep measured it directly a few hours later — see the second addendum.
The error was extrapolating a linear cost model from the smallest, most under-occupied
point. What survives unchanged is the README's thesis that **the win is concurrency**.

Docs corrected in `perlmutter/README.md` (banner, blocker section, the speed/wall-clock
paragraph), `CLAUDE.md` (Perlmutter bullet + the `reflect_symmetry_axis` bullet, which now
carries the A/B verdict), and `build_warpx.sh` (header + closing echo). One new warning
recorded while doing it: **do not build B from `origin/feature/reflect-symmetry-axis`** —
it forked before the heater merge and carries zero `ParticleHeater`/`TargetInjector` files.

---

## 2026-08-11 (Perlmutter, later) — S_phase: the upstream E_z is 1/sqrt(ppc) particle noise. The c/v_ti account is FALSIFIED

The sweep ran and its own premise failed, in the informative direction. `H_phase`
concluded the upstream defects were set by **c/v_ti** — a *physical* thermal-fluctuation
level — and predicted that refining dz and raising ppc **would not** remove them. Raising
ppc removes them, on a clean power law. **It is particle discreteness noise.**

### The runs

Job `56715249`, array `1-5`, binary A, one A100 each, 5/5 COMPLETED exit `0:0`, 183
diagnostic dumps apiece. Index 0 (`ss_dz1_ppc100`) was deliberately not resubmitted — the
A/B had already produced it *on binary A*, i.e. exactly what the sweep asks for.

### The measurement

`scripts/check_domain_control.measure`, upstream window 8–12 d_i0, identical parameters to
the A/B table above, so the two are directly comparable.

**rms E_z / (v_A B0)**

| run | dz/λ_D | ppc | N_D | t*=0.05 | 0.10 | 0.20 | 0.29 |
|---|---|---|---|---|---|---|---|
| `ss_dz4_ppc25` | 1.52 | 25 | 16.5 | 56.79 | 61.62 | 73.02 | 81.90 |
| `ss_dz2_ppc50` | 3.03 | 50 | 16.5 | 43.36 | 48.51 | 57.08 | 66.35 |
| `ss_dz1_ppc100` | 6.07 | 100 | 16.5 | 31.08 | 34.75 | 43.05 | 48.91 |
| `ss_dz2_ppc100` | 3.03 | 100 | 33.0 | 30.06 | 31.51 | 36.87 | 41.87 |
| `ss_dz4_ppc100` | 1.52 | 100 | 66.0 | 26.52 | 27.96 | 30.42 | 32.89 |
| `ss_dz1_ppc400` | 6.07 | 400 | 66.0 | 15.77 | 15.94 | 17.15 | 19.07 |

**It is a 1/sqrt(ppc) law, to within 8% over a 16x range**, anchored on the ppc = 100 mean
(29.22) at t*ω_ci0 = 0.05:

| ppc | 25 | 50 | 100 | 400 |
|---|---|---|---|---|
| observed | 56.79 | 43.36 | 29.22 | 15.77 |
| 1/√ppc predicted | 58.44 | 41.32 | 29.22 | 14.61 |
| obs / pred | 0.97 | 1.05 | 1.00 | 1.08 |

**dz barely matters.** At fixed ppc = 100, refining dz/λ_D from 6.07 → 1.52 (4x) moves
E_z only 31.08 → 26.52, **−15%**. Aliasing is a secondary effect here; the ppc term
dominates completely.

The two controlled comparisons say the same thing from both directions:

- **Fixed N_D = 16.5, dz refined 4x** — E_z *rises* by 1.68–1.83x. That is not a paradox:
  holding N_D = ppc·λ_D/dz fixed while refining dz **forces ppc down** 100 → 25, and
  √(100/25) = 2. It is the ppc law again, wearing the grid's clothes.
- **Fixed dz, ppc raised 4x** — ratio 0.39–0.51 at dz/λ_D = 6.07 and 0.40–0.47 at 1.52,
  against **0.500** for pure 1/√ppc and **~1.0** for a resolution-independent physical
  level. Both pairs land on the noise prediction.

### Why this overturns H_phase

A physical thermal-fluctuation field E_th = √(nT/ε₀) is a property of the plasma state,
so E_th/(v_A B0) = β·c/v_ti is **independent of ppc by construction**. The measured E_z is
not: it falls as 1/√N_particles, which is the signature of finite-particle sampling noise,
not of a physical field. The H_phase argument (RESULTS 2026-08-11 earlier) correctly showed
that both runs sit at 5–7% of *their own* thermal scale, but that consistency check cannot
distinguish "the physical thermal level" from "discreteness noise that happens to scale
with the same n and T". The ppc scan can, and does.

**Consequence:** the early ambient-ion acceleration in `R1_paper_470eV` is a
**resolution artifact after all**, and the fix is ppc, not dz. Criterion-2-style
conclusions that leaned on the c/v_ti account should be revisited.

### The caveat that matters for what to do next

**There is no plateau.** The ppc 100 → 400 ratio is still ~0.4–0.5, i.e. still tracking
1/√ppc with no sign of bottoming out on a physical floor. So E_z at ppc = 400 is **still
discreteness-dominated**, and ppc = 1600 would roughly halve it again. ppc = 400 is
therefore a 2x noise reduction, **not** a converged value, and it must be justified by
whether the early ion acceleration is *gone* — an observable of the ambient phase space —
rather than by E_z having converged. Do not quote ppc = 400 as "resolved" without that check.

### Everything else is flat

`coherent max|B_perp|/B0` and `piston front z/d_i0` show no systematic resolution trend
(front within ±3% across all six points, 0.53–0.56 d_i0 at t* = 0.05 rising to 3.53–4.06 at
0.29). The piston dynamics are converged; it is only the upstream noise that moves.

### Movies

Phase-space animations for all six points plus the symwall twin:
`media/S_phase/*/shock_phase.mp4` (gitignored, regenerable). A light corroboration of the
above, from an unexpected direction: the h264 file sizes fall monotonically with ppc —
187K at ppc = 25, 176K at 50, 166/125/96K at 100, **87K at 400** — because the encoder
compresses a smooth phase space better than a noisy one. Suggestive, not evidence.

### Cost, and the estimate that was wrong by 2.7x

| point | cost | wall | s/unit |
|---|---|---|---|
| `ss_dz1_ppc100` | 1 | 7m59s | **479** |
| `ss_dz2_ppc50` | 2 | 12m39s | 380 |
| `ss_dz2_ppc100` | 4 | 18m25s | 276 |
| `ss_dz1_ppc400` | 4 | 19m54s | 298 |
| `ss_dz4_ppc25` | 4 | 22m46s | 342 |
| `ss_dz4_ppc100` | 16 | 47m44s | **179** |

**GPU cost is not linear in the cost model.** s/unit falls 479 → 179 (2.7x) from the
smallest point to the largest, because `max_grid_size = n_cell` puts one box on one GPU and
the small points leave an A100 under-occupied. A pre-run estimate anchored on the cost-1
point predicted 2 h 09 m for `ss_dz4_ppc100`; it took 48 m. **Never size a job by
extrapolating the cheapest point** — that error is now recorded in CLAUDE.md.

This also **retracts the claim made in the previous entry** that the A100 gives no per-GPU
gain. Measured over the whole sweep: 7 767 GPU-s against chablis's projected 12 960 →
**~1.67x overall**, ~2.3x at the largest point, i.e. inside `perlmutter/README.md`'s
original "1.5–2.5x an RTX 4070" guess. That guess was right and the retraction of it was
wrong; the error was extrapolating from the single most under-occupied point. What survives
is the README's actual thesis, that **the win is concurrency**: wall-clock for the whole
sweep was ~49 min against 2.16 GPU-h of serial work.

### Staged, not yet reported: `R1_paper_470eV_ppc400`

`runs/R1_phase/R1_paper_470eV_ppc400` — `R1_paper_470eV` with **ppc 100 → 400 and nothing
else**. The generated deck differs from its parent in exactly five lines: the four species'
`num_particles_per_cell_each_dim`, plus `target_injector.ppc_reference`, which
`kinshock.deck` propagated on its own. Every derived scale is identical (M_A 13.95,
M_ms 12.76, β_ab 1150, lnΛ 12.23), confirming it is a controlled counterpart.

**It is on FOUR GPUs, not two.** The parent's `max_grid_size = 15000` (two boxes) was tuned
on chablis, a two-card machine, and is an artifact of that hardware; a Perlmutter GPU node
carries four A100s, so two boxes idles half the node. Changed to `7500` = four boxes, one
per GPU. The parent's note *"2 GPU, 4 box — fewer, larger wins"* is **not** an argument
against this: it measured two boxes *per GPU*, i.e. over-decomposition relative to rank
count; four boxes on four GPUs is one box per GPU, the same shape as the 2-box/2-GPU case
that measured 91% efficiency. NB `gpu_shared` caps at `gres/gpu=2`, so >2 GPUs means the
`regular` QOS and a whole node billed (~13% more allocation for ~1.8x less wall-clock).

New `perlmutter/job_multigpu.sbatch` carries this: one job = one run, ranks/GPUs from the
submit line, defaulting to 4. It sources `_common.sh` and calls `run_warpx`, so the
cd-before-launch invariant, the progress logger, the overwrite guard and the post-run
`--verify` are not duplicated. `job.sbatch` stays single-GPU for the one-box S_phase array.

⚠ ~~**QUEUE-DEPTH RISK, surfaced on merge.**~~ **MEASURED AND ACTED ON — see below.** The queue survey two entries above measured
`regular` at **+6 days** against `shared` at +11 h (`sbatch --test-only`, single-GPU spec).
Four GPUs *forces* `regular`, because `gpu_shared` caps at `gres/gpu=2`. So the 2 → 4 GPU
switch buys order 13 h of runtime and may cost days of queue — plausibly a net loss in
wall-clock. Two GPUs on `shared` (`mgs = n_cell/2`) is the slower configuration that may
well finish sooner. This was decided before that measurement was visible in this clone;
re-check with `--test-only` and revert to 2 GPUs if the estimate holds.

Submitted: pilot `56729828` (5 000 steps, 4 GPUs) and the full run `56729831` (48 h,
`gpu_regular`). **Neither has reported** — cost is projected at 3.5–4x the parent's 7.9 h
on two GPUs before any 4-GPU speed-up, and given the 2.7x miss above that projection should
be replaced by the pilot's measured s/step rather than trusted. Output ~59 GiB.

**RESOLVED 2026-08-12 — reverted to two GPUs.** `sbatch --test-only` settled it:

| config | QOS | would start | wait |
|---|---|---|---|
| 4 GPU, 48 h | `regular` | 2026-08-17 20:35 | **~5.8 days** |
| 4 GPU, 24 h | `regular` | 2026-08-17 20:46 | ~5.8 days |
| 2 GPU, 48 h | `shared` | 2026-08-12 06:35 | **~6.5 h** |
| 2 GPU, 24 h | `shared` | 2026-08-12 06:35 | ~6.5 h |

Time-to-result is ~36 h on two GPUs (6.5 h queue + ~30 h run) against ~6.5 days on four
(5.8 d queue + ~17 h run). **The slower configuration finishes about five days sooner**,
and shortening the walltime ask does not help `regular` at all — 24 h and 48 h return the
same date, so that lever is useless.

`max_grid_size` reverted 7500 → 15000, deck regenerated, jobs `56729828`/`56729831`
cancelled (neither had started, no compute consumed) and resubmitted as pilot `56731351`
and full run `56731356`, both `-n 2 -G 2 -q shared`.

### Pilot result 2026-08-12 — ppc x4 costs only 1.32x, not 4x

Pilot `56731351`: 5 000 steps, 2 GPUs, exit 0, 69 s wall.

    Avg. per step = 0.010571 s      (ppc = 400, 2 GPUs)
    parent        = 0.00801  s      (ppc = 100, 2 GPUs, measured 2026-08-04)
    ratio         = 1.32x for 4x the particles

**Projected evolve time: 0.010571 x 2 784 400 = 29 437 s = 8.2 h**, against the 28–32 h
projected from a "particle-dominated, so 3.5–4x" argument. **That estimate was ~3x too
high** — the third time today that predicting instead of measuring produced a wrong number
(after the 2.7x sweep miss and the retracted A100 claim), and the reason is the same one
the sweep already established: **these decks are GPU-under-occupied**, so adding particles
is nearly free until the card is actually fed. The production deck at ppc = 100 was itself
under-occupied — two A100-40GBs running with ~350 MB of pinned arena in use.

⚠ The pilot does NOT include plotfile I/O: `plotfile_intervals = 55688`, so the first
plotfile falls outside a 5 000-step run. Only two `field_intervals` writes are represented.
The parent's quoted 7.9 h against its own 6.2 h of pure evolve implies ~1.7 h of I/O at
~18 GiB; this run writes ~59 GiB, so budget more. Total is still expected well inside the
48 h request.

**This makes a ppc = 1600 point look cheap and worth having.** The sweep found no plateau
in the 1/sqrt(ppc) law, so the physical floor of the upstream E_z has never been located.
If ppc x4 really costs only ~1.3x here, that point is affordable, and it is the measurement
that would say where discreteness noise stops dominating — which is exactly the question
ppc = 400 leaves open.

**The lesson is the general one, not the number.** The 4-GPU switch was argued purely on
node utilisation — four A100s per node, one box per GPU, 91% efficiency precedent — and
every one of those statements is still true. It was wrong anyway, because on a shared
machine *time-to-result is queue + runtime*, and only the second term was optimised. The
queue term was measurable the whole time with `--test-only`, which is free and submits
nothing. **Measure the queue before choosing a resource shape**, not after. The
configuration knowledge is kept in the run's `config.yaml` so the four-box option can be
revived if the queue ever favours it — this is a property of the queue on the day, not of
the deck.

### Next

- Read the pilot's s/step; if 4-GPU scaling disappoints, revert to 2 GPUs on `shared`.
- When `R1_paper_470eV_ppc400` lands: check the **ambient phase space** against `R1_paper`
  for whether the early acceleration is gone. That, not E_z, is the acceptance test.
- Revisit anything that leaned on the c/v_ti account.
- Consider a ppc = 1600 point to find where E_z finally departs from 1/√ppc — that
  departure is the physical floor, and nothing measured so far has located it.

---

## 2026-08-14 — The ambient wedge is set by ELECTRON HEATING, not by resolution

The 470 eV run's "false early shock" is a **broad, multi-valued ambient-ion wedge** where
R1_paper shows a thin laminar arc. Seven numerical knobs and a five-rung resolution ladder
were tested; the wedge survives all of them. A hybrid (fluid-electron) run at the same
parameter point reproduces R1_paper's structure exactly. **The controlling variable is how
freely the electrons can heat.**

### What the scalars hid, and why the phase space was the right diagnostic

`G` (reflected-ion fraction), rms E_z and foot depth each compress the phase space to a
number and lose the structure that distinguishes a real foot from a filled wedge. G moved
*opposite* to visual quality four separate times. The whole investigation only converged
once it was driven by the phase-space movies rather than by scalars.

### The resolution ladder SATURATES (matched t* = 0.26)

| run | dz/lambda_D,amb | sigma(v_e)/v_te,0 in foot | wedge depth | cost |
|---|---|---|---|---|
| dz1 | 6.07 | 19.55 | 3.00 | 1x |
| dz4 | 1.52 | 15.72 | 2.75 | 16x |
| dz8 | 0.76 | 10.34 | 2.50 | 64x |
| **dz8 + shape3 + filter8** | 0.76 | 9.64 | **2.25** | **64x** |
| dz16 | 0.38 | 8.85 | **2.25** | 256x |

**dz16 bought nothing over dz8+s3f8 at 4x the cost.** Electron heating keeps falling (20%,
34%, 14% per rung, decelerating) while the wedge flatlines at 2.25.

**dz16 is FINER than R1_paper itself** (0.38 vs 0.60) and still has a 2.25 d_i0 wedge where
R1_paper has 1.25. At better-than-reference grid quality the structure is still wrong, so
**it is not a resolution artifact.**

### Everything else tested, all negative

| knob | result |
|---|---|
| ppc 100 -> 1600 | rms E_z fell **6x** (1/sqrt(ppc) to 8% over a 16x range) and G did **not move** |
| collisions removed entirely | E_z -30%, \|B_perp\| x2, structure unchanged |
| symmetry wall (pi-rotation) | nothing above the GPU noise floor |
| pec -> Silver-Mueller upstream | E_z **identical** at t* = 0.05; growth only slowed |
| heater/injector intervals 20 -> 1 | **visually identical**; striations are not an impulse train |
| shape3 + filter8 at dz1 (h1/h2/h3) | negative -- but see below |

The ppc result is the sharpest: a 6x reduction in upstream noise left the ion response
untouched. **The 1/sqrt(ppc) law is real and measures something that does not drive the
ions.** Any future "it's particle noise" argument has to explain that.

### The one real win: aliasing suppression, but only near threshold

`shape3 + filter8` was already tried and failed -- at dz1, where dz/lambda_D = 6.07. At dz4
(1.52) the same two knobs remove the bright reflected population, the sharp front spike and
most of the striations. At dz8 they buy a full rung of wedge depth (2.50 -> 2.25) **at zero
extra grid cost**, matching dz16 for a quarter of it. They only work once the grid is within
reach of the Debye length.

### The mechanism, measured

Profiling ambient ions and electrons through the front (matched t*):

- **Peak ion velocity spread is the SAME** in both runs (0.355 vs 0.381) -- the wedge is not
  thicker in velocity, the disturbed layer is **~2x deeper in space**.
- **Far-upstream electrons are NOT grid-heated** in either run (sigma/v_te,0 = 1.58 vs 1.10),
  which excludes global grid heating and is consistent with H_phase.
- **Shocked-layer electrons are ~1.7x hotter** (normalized) in the 470 eV run, and the
  enhancement switches on exactly where the ion precursor deepens.

Hotter electrons behind the front -> higher electron pressure -> deeper ambipolar precursor
-> ambient ions disturbed further upstream -> a filled wedge rather than a thin arc.

⚠ The heating is generated where dz/lambda_D is large (upstream, ramp), **not** where it is
observed: the local dz/lambda_D in the layer measured is 0.02-0.09, resolved by 10x even at
dz1. That is why local reasoning kept failing and why six historical knobs at dz1 did nothing.

### Why R1_paper looks different: its electrons are RELATIVISTIC

| | shocked sigma(v_e) | gamma | sigma(v_e)/v_sh |
|---|---|---|---|
| R1_paper (47 keV) | **0.532 c** | **1.18** | 3.81 |
| 470 eV | 0.091 c | 1.004 | 6.52 |

R1_paper's shocked electrons sit at **half the speed of light**. Its heating is
relativistically capped, which is why its precursor is shallow and its arc thin. That is an
artifact of choosing a 47 keV ambient to hit Table I's dimensionless numbers -- not the
laboratory plasma either run represents.

**WarpX has no "reduced c" capability and none was used.** Both runs use real c; the two
differ by TEMPERATURE (47 keV vs 470 eV). `c_sim/c_phys = 0.100` in the configs describes
PSC's normalization, not a WarpX knob. Corrected here because the loose phrasing was
actively misleading.

### The hybrid run closes it

`runs/hybrid_phase/H3_470eV_dense` (fluid electrons, polytropic gamma = 5/3, 470 eV):

| t* | 0.13 | 0.26 | 0.39 | 0.50 |
|---|---|---|---|---|
| **hybrid (fluid e-)** | 0.60 | **1.25** | 1.80 | 2.25 |
| **R1_paper (kinetic, capped)** | 0.50 | **1.25** | 1.75 | 2.39 |
| dz1 kinetic 470 eV | 0.75 | 2.96 | -- | -- |
| dz8+s3f8 kinetic 470 eV | 1.00 | 2.25 | -- | -- |

**Hybrid and R1_paper agree to ~6% at every time.** Two unrelated brakes on electron heating
-- relativity and a fluid closure -- give the same precursor; unconstrained kinetic electrons
give a different one. Figure: `media/comparison/wedge_ladder_hybrid.png`.

It also quantifies "forms too early": dz8+s3f8 reaches wedge 2.25 at t* = 0.26, hybrid at
t* = 0.50. **The kinetic 470 eV run develops its precursor ~1.9x faster.**

⚠ The hybrid run is **not** a single-variable swap: its piston front moves at 6.7 d_i0 per t*
against R1_paper's 14.8, on an 80 d_i0 domain at dz = 2.5 d_e, from a different project
(`H-PICShock`, different config schema -- scales were read from `warpx_used_inputs`, not
`kinshock.load`). Wedge depth is front-relative so it survives that, but do not present the
row as "the same run with fluid electrons".

### Recommendation

1. **Adopt `dz8 + shape3 + filter8` as production numerics** -- matches dz16 at 1/4 the cost.
2. **Stop refining.** The ladder is saturated and dz32 would cost ~1000x for nothing.
3. **Treat the deeper precursor as the consequence of a genuinely 470 eV plasma with
   unconstrained kinetic electrons**, not as a defect. If R1_paper-like structure is needed,
   hybrid or an intermediate temperature is a modelling choice, not a correction.
4. The temperature trade curve (T_e,ab vs dz/lambda_D vs lnLambda dial) is in this session's
   analysis: 4.7 keV gives dz/lambda_D = 1.92 for a 10x bill instead of 102x, at the cost of
   lnLambda becoming a 143x dial.

### Also fixed: `deck.py` emitted the background B through the parser path

`warpx.B_ext_grid_init_style = parse_B_ext_grid_function` put a **uniform** B0
(divergence-free by construction) on WarpX's general loaded-field path, which enables the
MLMG projection divergence cleaner, which restricts field BCs to periodic/pec/pmc/neumann.
Silver-Mueller aborted at init with exactly that message. The cleaner is enabled only when
`B_ext_grid_type` is neither `default_zero` nor `constant` (`WarpX.cpp:1150`), so emitting
`constant` + `B_external_grid` leaves it off. `_BC_MAP`'s `absorbing` entry had been dead
code since it was written; it now works. 15/15 structure tests pass and existing decks
verify as physically equivalent.

### Next

- Re-render the S_phase figures with the adopted numerics if any are to be published.
- The `argmin` frame matching in the comparison scripts has **no tolerance guard** and will
  silently return the same frame for two different requested times (it did, for the sparse
  hybrid run). Add one before reusing them.

---

## 2026-08-17 — The two runs differ in EXACTLY ONE physics parameter: `eps = v_te,ab/c`

After the numerical branch closed (2026-08-14: the ambient wedge saturates at 2.25 d_i0 and
survives every knob), the remaining question was whether `R1_paper` (47 keV) and
`R1_paper_470eV` are the same *physics*. They were built to be the same dimensionless
problem, so the audit is a check on that construction rather than a guess.

`scripts/dimensionless_audit.py` computes ~50 groups from the config primaries alone.

### Preserved to <1%

`m_i/m_e = 100`; `n_e0/n_e,ab = 0.008`; `n_t/n_e,ab = 2`; `beta_ab = 1150`;
`M_A = 13.952`; `M_ms = 12.74/12.76`; `v_sh/C_s,ab = 4.6`; `v_p/C_s,ab = 3.4288`;
`v_p/v_A = 10.4`; `d_i,ab/d_e,ab = 10`; `d_i0/d_i,ab = 11.18`; `d_i0/d_e0 = 10`;
`rho_i0/d_i0 = 10.40`; `rho_sh/d_i0 = 13.95`; `L/d_i0 = 80.50`; `L_target/d_i0 = 0.1789`;
`(1/w_ci0)/t_ab = sqrt(beta_ab) = 33.912`; `lambda_ab = 20`; **`nu_ei,ab/w_ce = 1.6956`**;
**`rho_e,ab/d_e,ab = 33.912`**.

`beta_0` and `T_0/T_e,ab` differ by 2.1% — Table I's rounding of `T_0` (10 vs 10.217 eV),
already documented in the config header, not a new finding.

### Everything else is an integer power of one number

    eps == v_te,ab/c = sqrt(theta_e,ab)      0.3033 (47 keV)  ->  0.03033 (470 eV)

| group | power | 47 keV | 470 eV |
|---|---|---|---|
| `v_A/c`, `v_sh/c`, `v_p/c`, `C_s,ab/c` | `+1` | 0.0100, 0.1395, 0.1040, 0.0303 | 0.0010, 0.0140, 0.0104, 0.0030 |
| `lambda_D/d_e`, `lambda_D,0/d_i0` | `+1` | 0.3033, 4.47e-3 | 0.0303, 4.42e-4 |
| `nu_ei,ab/w_pe` | `+1` | 0.01517 | 0.001516 |
| `w_pe/w_ce`, `w_pi0/w_ci`, `rho_e/lambda_D` | `-1` | 10.0 upstream, 100 | 100.0, 1000 |
| `theta_e,ab`, `sigma_i`, `sigma_e` | `+2` | 0.092, 1e-4, 1e-2 | 9.2e-4, 1e-6, 1e-4 |
| `lnLambda` (at fixed `lambda_ab = 20`) | `-4` | 1.22e5 | 12.23 |

Every exponent came out to ±0.005 of an integer, which is itself the check that no
independent second parameter is hiding in the set.

**The difference is one-dimensional.** The user's framing was right in substance: both runs
use the real `c`, but `R1_paper`'s temperatures *are* PSC's reduced-`c` values used at real
`c`, so it behaves in every respect like a 10x reduced speed of light — and yes, the
magnetizations differ, `sigma_e = 1e-2` against `1e-4`, `sigma_i = 1e-4` against `1e-6`.

### The two groups that actually matter, and the two that don't

Not on the differing list, and worth stating because it kills two candidates outright:
`rho_e,ab/d_e,ab = 33.9` and `nu_ei/w_ce = 1.696` are **preserved**. Electron magnetization
relative to the *inertial* scales, and collisionality relative to the *field*, are identical.
What differs is electron magnetization relative to the *Debye* scale —
`rho_e/lambda_D = w_pe/w_ce`, 10 in the paper run against 100 at 470 eV.

`v_A/c` is also excluded, and by measurement rather than argument: `H3_470eV_dense` carries
the **470 eV** value (9.999e-4) and still reproduces R1_paper's thin structure. Whatever is
responsible lives in the kinetic-electron sector — which is what the hybrid's fluid closure
removes, and what both surviving candidates are:

- **(a) relativistic capping of electron heating**, `~eps^2`. `gamma(T_e,ab)` = 1.0014,
  1.0044, 1.0138, 1.0436, 1.1380 across a decade in `eps` — flat until the top. R1_paper's
  *shocked* electrons were measured at 0.532 c, `gamma = 1.18`.
- **(b) electron-scale wave regime**, `~eps^-1`. `rho_e/lambda_D = w_pe/w_ce` = 100, 56, 32,
  18, 10 — a smooth decade, and it sets how many Bernstein resonances sit below `w_pe`, i.e.
  how ECDI-like versus Buneman-like the ramp's electron heating is.

Both predict *more* electron heating at 470 eV, which is what was measured (shocked-layer
electrons ~1.7x hotter in normalized units, 2026-08-14). They differ in the **shape** of
`wedge(eps)`: (a) is flat-then-drop, (b) is a straight decade. That is what the ladder reads.

### `runs/E_phase/` — the ladder, staged 2026-08-17, NOT submitted

Five rungs; rung 1 is the existing `S_phase/ss_dz16_ppc100`, rungs 2-5 are new.

| run | T_e,ab | eps | B0 [T] | `w_pe0/w_ce` | `gamma` | cells | steps | wall (1 GPU) |
|---|---|---|---|---|---|---|---|---|
| `ss_dz16_ppc100` | 470 eV | 0.0303 | 7.026 | 100.0 | 1.0014 | 71680 | 2400000 | 7.35 h (done) |
| `es_1p5keV` | 1486 eV | 0.0539 | 12.494 | 56.2 | 1.0044 | 40312 | 759079 | 1.31 h |
| `es_4p7keV` | 4700 eV | 0.0959 | 22.220 | 31.6 | 1.0138 | 22667 | 239998 | 0.23 h |
| `es_15keV` | 14860 eV | 0.1705 | 39.509 | 17.8 | 1.0436 | 12748 | 75909 | 0.04 h |
| `es_47keV` | 47012 eV | 0.3033 | 70.273 | 10.0 | 1.1380 | 7167 | 23994 | 0.01 h |

Verified identical at every rung: `M_A = 13.952`, `M_ms = 12.759`, `beta_ab = 1150`,
`beta_0 = 0.1957`, `rho_i0/d_i0 = 10.40`, `L = 12.02 d_i0`, run `= 0.302/w_ci0`,
**`dz/lambda_D,0 = 0.379`**, `ppc = 100`, **`N_D,PIC = 264`**. The top rung reproduces
R1_paper's `eps`, `gamma` and `w_pe/w_ce` exactly — at better Debye resolution (0.379 vs
0.60) and higher `N_D` (264 vs 167) than R1_paper itself.

**Why `dz/lambda_D` is held and not `dz/d_e`.** `lambda_D ~ eps` while `d_e`, `d_i0` and the
domain are `eps`-INDEPENDENT (they depend on `n` and the real `c` only). A ladder at the
paper's fixed `dz/d_e = 0.3` would drag `dz/lambda_D` along by the same factor of 10 —
reinstating the exact confound the resolution ladder spent a week excluding. Fixing
`dz/lambda_D` and `ppc` pins the finite-grid margin *and* the discreteness noise, and lets
`dz/d_e` float 0.019 -> 0.188 (still >=5 cells per `d_e,ab` at the coarsest rung).

**Cost `~eps^-3`** at fixed `dz/lambda_D` (cells `~1/eps`, steps `~1/eps^2`), so every rung
above 470 eV is *cheaper* than the anchor: **~1.6 GPU-h for all four**. The two top rungs are
small enough that GPU occupancy will be poor and their real cost is likely several times the
scaled estimate — irrelevant at ~1 minute each, and the mistake to avoid is the reverse one
(extrapolating a small point up).

**One honest caveat: `lnLambda`.** Holding `lambda_ab = 20` with `nu_ei ~ n lnL T^-3/2` forces
`lnLambda ~ T^2`, so it rides 12.2 -> 1.22e5 across the ladder — the same unphysical dial
R1_paper always carried. Kept rather than removed, because the collisional state is what both
baselines have and because `ss_dz1_ppc100_nocoll` already showed that deleting collisions
entirely does not touch the wedge. `nu_ei/w_ce = 1.696` is preserved at every rung.

Submit with `perlmutter/submit.sh eps` (4-task array, cheapest-first, `shared` QOS).
Read the wedge from the **visual structure of the ambient-ion phase space** at matched
`t*wci0`, as agreed 2026-08-13 — not from a scalar. All rungs have 30 frames over the same
`0.302/w_ci0`, so frames match by index and the `argmin` tolerance problem does not arise.

15/15 structure tests pass; `make_inputs.py` verified all four decks resolve back to their
config primaries.

### 2026-08-17 (same day) — three E_phase rungs died on `blocking_factor`; grid snapped, guard added

Array `57175431` launched and 3 of 4 tasks aborted ~45 s in:

    amrex::Error::0::domain size not divisible by blocking_factor !!!

**Cause, and it was mine.** The generator sized each rung by `n_cell = round(L/dz)`, which
produces an arbitrary integer; AMReX's default `amr.blocking_factor` is 8.

| rung | n_cell as generated | /8 | |
|---|---|---|---|
| `es_1p5keV` | 40312 | 5039.000 | ran fine — divisible by luck |
| `es_4p7keV` | 22667 | 2833.375 | ABORT |
| `es_15keV` | 12748 | 1593.500 | ABORT |
| `es_47keV` | 7167 | 895.875 | ABORT |

Every pre-existing run in the repo is divisible by 8, which is why this had never surfaced.

**Fix.** `n_cell` snapped to a multiple of 8 (22664 / 12744 / 7168) and `dz`, `max_step`
and the diagnostic intervals re-derived from it. `dz/lambda_D,0` moves by **<0.03%**
(0.3791 -> 0.3792), so the held invariant is untouched; re-verified that `M_A`, `M_ms`,
`beta_ab`, `beta_0`, `rho_i0/d_i0`, `L/d_i0`, `t*wci0` and `N_D` are still identical at
every rung.

**Guard.** `deck._n_cell` now RAISES on a non-multiple, naming the two nearest legal cell
counts and the `dz_over_de` that lands on one. The failure mode this replaces is the
expensive kind — the config generator was perfectly happy, `make_inputs.py --verify`
passed, and the error only appeared after the job had queued, launched and initialised,
with a backtrace naming neither the config nor the offending number.

`perlmutter/submit.sh` now accepts run dirs positionally (`submit.sh eps <dir> ...`) so a
partial re-run is a subset of the same target rather than a hand-rolled sbatch line. Refused
for `ab`, where the array index selects the binary and a subset would silently remap it.

Resubmitted the three as array `57183562`. `es_1p5keV` (task 3 of the original array) was
never affected and ran through — ~1.5 h against the 1.31 h estimate.
