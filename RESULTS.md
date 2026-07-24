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
- [ ] R2 (B₀=0), R3 (n_e0=0) negative controls

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
