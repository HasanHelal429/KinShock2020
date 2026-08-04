# Replication Plan — Schaeffer 2020 piston-driven magnetized shocks in WarpX

**Goal.** Verify that WarpX's ported laser-ablation surrogate (`ParticleHeater` +
`TargetInjector`, in `warpx-cda/`, already operator-validated against PSC) can generate the
*piston-driven magnetized perpendicular collisionless shocks* of Schaeffer et al. 2020
(`OVERVIEW.md`) and reproduce their quantitative formation signatures.

**Scope of this document.** Planning only — parameter mappings, run matrix, analysis
specification, and runtime estimates. No input decks or scripts are written here; §8 lists the
concrete files to author (all under `KinShock2020/`) when execution is approved.

**What already exists (starting point).** `warpx-cda/heating_operator/run_shock_1d/` contains
a working 1D3V piston-shock deck (`inputs_shock_1d`), a plan (`SHOCK_PLAN.md`), and an
analysis script (`scripts/make_shock_figures.py`) that already produces density/field
line-outs, ion phase space, and a shock-front trajectory → Mach number. **That run used
demonstration parameters, not Schaeffer-2020 values** (θ_e = 0.04, n_t/n_amb contrast 25,
M_A ≈ 5, no collisions). This plan's job is to (a) retarget the deck to the paper's Table I
regime and (b) reproduce the paper's specific diagnostics (speed model, 7 criteria, 3
timescales, reflected-ion fraction, streak plots, negative controls) with **new, self-contained
analysis scripts in `KinShock2020/` that are designed after the `warpx-cda` scripts but do not
modify or import them** (§6, §8).

---

## 1. Reproduction strategy: match the dimensionless run, not the µm-scale experiment

**PSC's "reduced speed of light" is not an independent input — it *is* θ_e,ab.** §II, p.3:
"a reduced speed of light set by the ratio T_e,ab/m_e c², which can be written relative to the
sound speed as c = √(µ_p/T_e,ab) C_s,ab." Pick θ_e,ab and µ_p and c/C_s,ab is determined:
√(100/0.092) = 33.0, i.e. Table I's C_s,ab = 0.030 c. **So setting θ_e,ab in real-c WarpX
reproduces the reduced-c run exactly; there is nothing to port.** (R1_warm's θ_e = 0.078 gives
c/C_s,ab = 35.8 vs the paper's 33.0 — an 8% gap that is entirely the M_A recalibration.) The
`c_sim/c_phys = 0.02` row is a *reported consequence* of θ_e,ab, not a knob, and is anyway
rounded: the velocity rows imply 0.0234.

We therefore reproduce the **sim-value / dimensionless columns of Table I directly**, reading
every "× c" as a fraction of the *real* c. This yields hot (≈ 47 keV) piston electrons and a
fast (≈ 0.1 c) piston — i.e. we run PSC's electron-reduced simulation, which is exactly the
object we are validating against.

**Table I sim → phys conversions** (three *different* factors — no single one works):
velocities × 0.0234; temperatures **÷ 100** (= 0.0234² × 1836/100, because the phys column also
restores the real mass ratio — *not* ÷1824); lengths at real c and real m_p, unscaled.
Matching the "physical value" column therefore costs s = √(39900/470) = 9.21 in the
config-only emulation (RESULTS 2026-07-31), **not** s = 50 — 0.02 is a velocity factor.

**Invariants to preserve** (everything else follows): m_i/m_e = 100 · θ_e,ab = 0.092 (heater)
· T₀/T_e,ab = 0.0217 · density contrast n_t : n_e,ab : n_e0 = 2.5 : 1 : **0.008** (Table I's
0.01 is in code units, where n_e,ab = 1.25) · M_A = 14 · M_ms = 13 · β_ab = 1150 · β₀ = 0.2
(both under β = µ₀nT/B²) · collisionality λ_ab = 20 = mfp/d_e,ab.

**Reference density.** Define the WarpX reference n₀ ≡ n_e,ab (the ablation density). Its
absolute value is free (choose for a convenient dt / domain in SI, e.g. n₀ = 1×10²⁴–10²⁶ m⁻³);
d_e ≡ c/ω_pe(n₀), d_i = 10 d_e.

---

## 2. Parameter mapping (Table I → WarpX 1D deck)

Geometry: **1D3V**, z = ablation/propagation axis. **The paper's own domain is one-sided**
(§II, p.3): 30000 cells over 9000 × 5 d_e,ab, target at 0 < z ≤ 2 d_i,ab, ambient plasma and
field at z > 2 d_i,ab — it does *not* ablate in both directions from a centered slab. An
earlier version of this line said otherwise; `layout: one_sided` (R1_warm, R1_coll) is the
faithful geometry, and the symmetric/periodic full-domain runs are the approximation, not the
other way round. Perpendicular field **out of the propagation axis: B_x = B₀** (matches the
paper's B₀ x̂; any transverse component is valid in 1D).

| Paper (Table I) | WarpX deck setting | Notes |
|---|---|---|
| Domain L_z = 9000 d_e,ab (one-sided, 30000 cells) | **tiered** (see below) | R1_warm uses 7500 d_e = 83% of the paper's box |
| dz ≈ 0.3 d_e | `n_cell = L_z / 0.3 d_e` | resolve d_e; `algo.particle_shape = 2–3` |
| Target region 0<z<2 d_i,ab | slab `|z| < 2 d_i` (= 20 d_e) | heater + injector region |
| Target density 2.5 n_e,ab | `n_t = 2.5 n₀` | injector target density |
| Ablation ref n_e,ab | `n₀` | fixes d_e, ω_pe |
| Ambient **n_e0 = 0.008 n_e,ab** (Table I: 0.01 *code units*, n_e,ab = 1.25) | `n_amb = 0.008 n₀` — **runs still use 0.01, 25% high** | contrast 312, not 250. Sets B₀ ∝ √n_amb ⇒ ω_ci0; see CLAUDE.md gotcha |
| mass ratio 100 | `M_i = 100 m_e` | piston & ambient ions |
| T_e,ab = 0.092 m_e c² | `particle_heater.<piston_e>.theta = 0.092` | drives C_s,ab = √(θ/100) c = 0.030 c |
| Heating for whole run | `particle_heater.intervals = 20` | PSC flatfoil cadence |
| Injector maintains 2.5 n_e,ab | `target_injector.density = n_t`, `tau = 40/ω_pe` | replenishment |
| T₀ = 0.002 m_e c² (ambient) | amb electrons `u_std = √0.002`; ions `u_std = √(0.002/100)` | upstream temp |
| B₀ ⇒ M_A = 14, M_ms = 13, β₀ = 0.2 | `B0 = v_A·√(µ₀ n_amb M_i)`, target **v_A ≈ 0.01 c** | tune B₀; verify M_A, β₀ from loaded state |
| collisionality λ_ab = 20 | binary Coulomb collisions on (H+, e⁻) | see §4 — run collisionless first |
| 400 000 steps, τ = 220 t_ab (≈ 6.5 ω_ci0⁻¹) | `max_step` set by target time (see §7) | **≈ 50 000 steps per ω_ci0⁻¹** at dz=0.3 d_e, CFL 0.75 (dt·ω_pe ≈ 0.225) |

Piston/ambient are loaded as **separate species pairs** (`piston_electrons/ions`,
`amb_electrons/ions`) so the injector acts only on the piston, the heater heats only piston
electrons, and reflected *ambient* ions are cleanly identifiable in analysis (criterion 6).

**Self-consistency checks before trusting a run** (measure from the loaded/early state, don't
assume): C_s,ab = √(θ_e/100)·c ≈ 0.030 c; v_A = B₀/√(µ₀ n_amb M_i); M_A = v_sh/v_A ≈ 14;
M_ms = v_sh/√(v_A²+C_s²) ≈ 13; β₀ = 2µ₀ n_e0 k_BT₀/B₀² ≈ 0.2. Adjust B₀ (and, if needed, θ_e)
to land on M_A ≈ 14, β₀ ≈ 0.2.

**Domain tiering (drives cost — see §7).** Symmetric periodic slab → two shocks in ±z, so the
*half*-domain must exceed the shock position at the final analysis time. With
ρ_i0 = v_p/ω_ci0 ≈ 1175 d_e and the shock front at ≈ z\* ρ_i0, the target diagnostic time sets
the box:

| Tier | Target time | Front reaches | Half-domain | `n_cell` (dz=0.3 d_e) | Steps (≈50k/ω_ci0⁻¹) |
|---|---|---|---|---|---|
| **Core** | t\*₂ ≈ 2.5 ω_ci0⁻¹ (onset + separation) | ≈ 2.5 ρ_i0 ≈ 2900 d_e | ±3600 d_e | ~24 000 | ~125 000 |
| **Full** | t\*₃ ≈ 5 ω_ci0⁻¹ (downstream forms) | ≈ 5–6 ρ_i0 ≈ 7000 d_e | ±7500 d_e | ~50 000 | ~250 000 |
| Late | ~10–12 ω_ci0⁻¹ (reformation, RH relax) | paper's 9000 d_e half | ±9000 d_e | ~60 000 | ~500 000–600 000 |

The **Core** tier is enough to confirm shock formation (criteria 1–7, t\*₁, t\*₂); the **Full**
tier is needed for t\*₃, compression-ratio asymptotes, and the negative-control comparison at
matched time; **Late** is optional (cyclic reformation, RH downstream). Alternatively, ablate
one-sided into a 0→9000 d_e box with a reflecting/thermal wall at z=0 to match the paper's
geometry exactly at half the cell count.

---

## 3. Run matrix

Ordered by cost/priority. Start in 1D (the paper is quasi-1D, so 1D is faithful).

| # | Run | Purpose | Key change from R1 |
|---|---|---|---|
| **R0** | Smoke test | Deck sanity: no NaN/warnings, heater+injector active, dt·ω_pe as expected | R1 with `max_step ≈ 2000` |
| **R1** | **Representative (M_A = 14)** | Primary reproduction of Table I | — (baseline) |
| R2 | B₀ = 0 | Negative control: no field → no magnetic compression / reflection / shock (Fig. 9) | `B0 = 0` |
| R3 | n_e0 = 0 | Negative control: no ambient → no ambient structures (Fig. 9) | `n_amb = 0` |
| R4 | M_A ≈ 7 | Mach scan; expect same t*/z* in ω_ci0⁻¹/ρ_i0, more piston–shock separation (Fig. 6) | raise B₀ (larger v_A) |
| R5 | M_A ≈ 28 | Mach scan, high end | lower B₀ |
| R6 | Collisional | Match λ_ab = 20; confirm similar formation (Fig. 13) | enable Coulomb collisions |
| R7 | m_i/m_e = 400 | Mass-ratio convergence (Fig. 14) | `M_i = 400 m_e` |
| R8 | CH multi-species | Multi-species formation, C-dominated profiles (Fig. 12) | add C ions (Z=… ), H+C piston/ambient |

R1–R3 are the core verification (shock + two negative controls). R4–R8 are robustness/scan
extensions to run if the core passes.

Each run writes: full field+particle plotfiles at fine cadence (for phase space and streak
plots — ~100+ frames over the run), and reduced diagnostics (`ParticleEnergy`, `ParticleNumber`
per species) to confirm heater↔injector energy balance and piston replenishment.

---

## 4. Collisions

The paper includes Takizuka–Abe collisions (λ_ab = 20) but shows (Fig. 13) that shock
*formation* is essentially the same collisionless, because ion–ion interactions stay
collisionless at the high flow speeds (λ_mfp/d_i0 ≈ 350). **Plan:** run R1 collisionless first
(cheaper, and matches how the heater/injector were operator-validated); add R6 with WarpX
binary Coulomb collisions on the (H⁺, e⁻) pairs, tuned to λ_ab ≈ 20, to confirm formation is
unchanged and to get realistic electron isotropy/heating in the downstream.

**Status: `runs/R1_phase/R1_coll/` fills the R6 slot** (config + deck done 2026-07-27, not yet launched)
— R1_warm's collisional twin, every dimensionless primary identical, pinned to an absolute
ambient n_e0 = 10¹⁸ cm⁻³ (`reference.n0 = 1e26` m⁻³). Config-driven via a `collisions:` block
(`model`, `pairs`, `target`, `ndt_supercycle`) rendered by `kinshock.deck`; `--verify` covers it.

Two things to know before using it (details + table in `RESULTS.md` 2026-07-27):

- **lnΛ is a knob, not a physical value.** At real c, θ_e = 0.078 ⇒ T_e,ab = 39.9 keV (not the
  paper's ~470 eV), so ν_ei ∝ n T^(−3/2) leaves the plasma collisionless by ~4 orders of
  magnitude even at n_e,ab = 10²⁰ cm⁻³ (physical lnΛ = 11.6 → mfp = 3.7×10⁵ d_e,ab). A
  physically-collisional run would need n_e,ab ~ 3×10²⁸ cm⁻³. `collisions.target` therefore
  states the *physics* target and `units.coulomb_log_for` inverts it to the deck's lnΛ.
- **λ_ab ≠ mfp/d_e.** λ_ab ≡ ω_ce,ab/ν_ei,ab = mfp/**ρ_e,ab**, and ρ_e,ab = 27.9 d_e,ab here.
  Table I's λ_ab = 20 is mfp = 559 d_e,ab = **5.6 d_i0** — still collisionless at ion scales,
  as the paper is. R1_coll uses that value (`quantity: lambda_ab, value: 20` → lnΛ = 7713,
  ν_ei·dt = 1.1×10⁻⁴), so it is a fair comparison against Fig. 13. Beware `mfp_over_de: 20`,
  which sounds the same but is 28× *more* collisional (0.2 d_i0, collisional across the ramp).

---

## 5. Known risks / mitigations

- **Ambient shot noise.** Contrast n_t : n_e0 = 250; if ambient ppc scales with density the
  upstream is starved. **Mitigation:** load the ambient species with its own adequate ppc
  (e.g. 100–400/cell at n_e0), independent of the 1000-at-n_e,ab piston loading — they are
  separate species. Watch upstream ⟨u²⟩ for spurious growth.
- **Numerical grid heating** of the cold (T₀ = 0.002) ambient over a long run. **Mitigation:**
  resolve d_e (dz ≈ 0.3 d_e), `particle_shape = 2→3`, raise ppc; monitor ambient temperature
  vs. time and confirm it is flat far upstream.
- **Boundary wrap.** Periodic BC with two shocks. **Mitigation:** size the domain (±4500 d_e)
  and `max_step` so both shock fronts stay inside through t*₃ + downstream development.
- **Reduced-c relativistic electrons.** θ_e = 0.092 ⇒ v_te ≈ 0.3 c piston electrons. WarpX is
  relativistic so this is fine, but keep the Boris/Vay pusher and CFL appropriate.
- **v_p vs v_sh disambiguation.** The whole point of the paper — do not read piston
  compression / piston-accelerated ambient ions as a shock. The analysis (§6) must implement
  criterion 7 (separation) and the reflected-ion definition, not just "a density bump moved."

---

## 6. Analysis specification (reproduce the paper's figures/criteria)

**Self-contained by design.** Do **not** modify or import the `warpx-cda` scripts. Use them as
a *design reference* — `heating_operator/scripts/make_shock_figures.py` (plotfile loading,
per-species ion-density histogramming, B read, front-trajectory fit) and `make_movies.py`
(yt frame rendering + ffmpeg encode) — and write **new** code under `KinShock2020/`, split into
a reusable library and thin drivers, depending only on the WarpX output and a standard
yt/matplotlib/numpy/PyYAML environment.

**Shared helpers live in `KinShock2020/src/` (a `kinshock` package); driver scripts in
`KinShock2020/scripts/` stay thin.** This avoids the design smell in the reference script,
which hard-codes deck constants at the top with a "keep in sync with the deck" comment — a
drift hazard. Instead:
- `src/kinshock/config.py` — load/validate a per-run config file (§6.0a).
- `src/kinshock/units.py` — compute derived scales (ω_pe, d_e, d_i, ω_ci0, ρ_i0, C_s,ab, v_A,
  M_A, M_ms, …) **from the config's primary quantities**, so there is one source of truth.
- `src/kinshock/io.py` — plotfile discovery/loading, per-species density & (z,u) phase space.
- `src/kinshock/metrics.py` — speed model (Eqs. 1–2), the 7 criteria, reflected-ion F/G,
  t\*/z\* extraction, compression ratios.
- `src/kinshock/plotting.py` — shared styling + figure/movie (ffmpeg) helpers, media paths.
- `src/kinshock/deck.py` — generate a WarpX input deck from a config (`render`), and parse/
  resolve a deck back to numbers for post-run verification (`verify`).
- `scripts/make_inputs.py`, `scripts/run_checks.py`, `scripts/make_figures.py`,
  `scripts/make_movies.py` — thin drivers that take a run id / config path, call `src`
  helpers, and generate the deck or write to `media/`. No physical constants are written in
  the scripts.

**6.0a Per-run config files (single source of truth).** Each run carries a config file
(`KinShock2020/runs/RXX/config.yaml`, matching your OSIRIS `run.yaml` convention) holding the
run's **primary** parameters and metadata — nothing derived is hand-entered:
- metadata: `run_id`, `tier` (Core/Full/Late), `description`, WarpX build, deck path;
- primaries (SI + normalized): `n0`, `mass_ratio`, `theta_e`, `nt`, `namb`, `T0`, `B0` (or the
  target `vA`), `domain`, `n_cell`, `cfl`, heater/injector `intervals`/`tau`, species names,
  `dt` and `max_step` as run.
`units.py` derives everything else from these and can echo a `derived:` block back for
reference; `config.py` validates the config against the Table I targets (M_A≈14, β₀≈0.2, …) and
warns on mismatch. **The config is authored by hand and `scripts/make_inputs.py` generates the
WarpX deck from it** (`kinshock.deck.render`), so the deck is a build artifact that never needs
hand-editing — all WarpX-deck details can be ignored while reasoning about a run. After a run,
`make_inputs.py --verify` parses the WarpX `warpx_used_inputs` file and confirms it matches the
config, so what was simulated provably equals the config (no manual transcription). Every figure
caption records the `run_id` and key parameters pulled from the config, so plots in `media/` are
self-describing.

**All figures and movies — including bring-up/testing figures — are written to
`KinShock2020/media/`** (create `media/testing/` for progress artifacts and `media/RXX/` per
run for final figures) so the project is fully self-contained and progress is visible at a
glance. WarpX plotfiles/reduced diagnostics stay under each run's `KinShock2020/runs/RXX/diags*`.

All quantitative plots are in **paper-normalized units** (z/d_i0, z/ρ_i0, t·ω_ci0,
v/C_s,ab or v/v_sh, n/n_e0, B/B₀):

**0. Progress / testing figures (`media/testing/`).** Emit these continuously during bring-up
so the run matrix's state is legible without re-running anything:
- Loaded-state sanity (R0/short R1): initial n(z) per species, B_x(z), and the derived
  C_s,ab, v_A, M_A, M_ms, β₀ printed on the figure vs. their Table I targets.
- dt·ω_pe and energy-conservation history (should stay within a few %), piston-inventory
  (`ParticleNumber`) and piston-energy (`ParticleEnergy`) histories showing heater↔injector
  balance — the operator sanity check inherited from the flatfoil validation.
- A low-frame-count "quicklook" density+field movie per run as it completes.

**A. Speed model (Fig. 2, 3; Eqs. 1–2).**
- **B_x(z,t) streak plot** with overlaid piston-speed and shock-speed lines.
- Track **v_p** = speed of the peak field tied to the magnetic cavity; **v_sh** = speed of the
  faster magnetic compression in the shocked ambient. Report v_sh/v_p (expect → 4/3 at low B,
  rising with B₀) and compare to the perpendicular RH relation Eq. 2.
- Optional across R1/R4/R5: v_p vs n_e0 and v_sh/v_p vs v_A/C_s,ab against Eqs. 1–2.

**B. Phase-space evolution (Figs. 5, 7).** Panels at several t·ω_ci0 of ambient-ion, piston-ion,
and electron (v_z, z) phase space, with B_x/B₀ and n_e/n_e0 profiles overlaid — the primary
qualitative signature (snowplow → reflection → separation).

**C. Seven formation criteria (Sec. III B).** Per output time, evaluate and log:
1. M_ms > 1 (local); 2. L/λ_ii > 1; 3. n_e/n_e0 > 2; 4. B/B₀ > 2; 5. steep ramp (dB/dz,
dn_e/dz on ~d_i0); 6. reflected ambient ions present; 7. front separated ≥ ¼ ρ_i0 from the
piston peak field. Classify each structure as **precursor** (1–6) vs **shock** (1–7).

**D. Reflected-ion fraction and timescales (Figs. 4, 6).** Define, for ambient ions,
G(t) = N_a,refl/N_a,tot and F(z,t) = f_a,refl/N_a,tot with "reflected" = v_z > v_sh(t).
Extract **t*₁** (max dG/dt), **z*₁** (max dF/dz), **t*₂** (separation ¼ ρ_i0 from piston),
**t*₃** (well-defined downstream). Expect **t* ≈ 1, 2.5, 5 ω_ci0⁻¹** and **z* ≈ 1, 2.5,
5–6 ρ_i0**, roughly independent of M_A (verify with R4/R5). ρ_i0 = v_p/ω_ci0.

**E. Compression ratios & late-time (Figs. 10, 11).** n_e/n_e0 and B/B₀ at the front vs time →
approach ~4 for M_A = 14; look for **cyclic reformation ~1.5 ω_ci0⁻¹** and **RH-consistent
downstream by ~10 ω_ci0⁻¹**.

**F. Negative controls (R2, R3).** Confirm that with B₀ = 0 the magnetic compression, strong
ion heating, secondary compression, and reflected-ion population **disappear**, and with
n_e0 = 0 no ambient structures form — i.e. the signatures are shock-specific, not piston
artifacts.

**Acceptance criteria (module is verified if R1 reproduces):** the v_sh ≈ (4/3–…)·v_p speed
relation; criteria 1–7 satisfied in order; t*₁,₂,₃ ≈ 1/2.5/5 ω_ci0⁻¹ and z*₁,₂,₃ ≈
1/2.5/5–6 ρ_i0 (within the paper's scatter); ~4× compression and a clear reflected-ambient-ion
beam; and R2/R3 controls null. Report tolerances relative to the paper's Fig. 6 scatter, not
exact equality (different code, RNG, real vs reduced c).

---

## 7. Runtime estimates

**Method.** Scaled from the one measured data point — the `run_shock_1d` movie run
(3200 cells, 4 species × 100 ppc ≈ 0.7 M macroparticles, 25 000 steps, 8 MPI ranks, ~35 min
wall on a Perlmutter **CPU** login node) → **≈ 1×10⁶ particle-updates / s / core**. Cost is
dominated by the particle push, so runtime ≈ (macroparticles × steps) / (1e6 × cores).

**Assumptions** (all linear levers — scale the numbers if you change them):
- Ambient loaded at **100 ppc/species** (its own loading, independent of the piston); ambient
  particles dominate the count because they fill the whole box (piston occupies ~40 cells).
- Steps ≈ **50 000 per ω_ci0⁻¹** (dz = 0.3 d_e, CFL 0.75). Grid/steps per tier from §2.
- Throughput anchor is a *login-node CPU* run → treat as a conservative upper bound; batch
  CPU or GPU will be faster. **Overall uncertainty ≈ ±2×.**

**Per-run estimate** (core-hours, and wall time on **256 cores** = 2 Perlmutter CPU nodes):

| # | Run | Tier | ~Macroparticles | Steps | **Core-hours** | Wall @256 cores |
|---|---|---|---|---|---|---|
| R0 | smoke | Core grid, short | 5 M | 2 000 | ~3 | ~1 min |
| **R1** | **M_A=14 baseline** | **Full (→ t\*₃)** | 10 M | 250 000 | **~730** | **~2.9 h** |
| R2 | B₀=0 control | Core (→ t\*₂) | 5 M | 125 000 | ~180 | ~0.7 h |
| R3 | n_e0=0 control | Core (no ambient) | ~1.5 M | 125 000 | ~60 | ~0.25 h |
| R4 | M_A≈7 | Core | 5 M | 125 000 | ~180 | ~0.7 h |
| R5 | M_A≈28 | Core | 5 M | 125 000 | ~180 | ~0.7 h |
| R6 | collisional | Core + collisions (~+30%) | 5 M | 125 000 | ~235 | ~0.9 h |
| R7 | m_i/m_e=400 | to t·ω_ci0≈1.1 (4× gyroperiod, larger box) | 7.4 M | 220 000 | ~450 | ~1.8 h |
| R8 | CH multi-species | Core, +ion species | 7 M | 125 000 | ~350 | ~1.4 h |

**Totals** (serial; independent runs can overlap across more nodes to cut wall time):
- **Minimum verification** (R0 + R1-Full + R2 + R3): **≈ 970 core-hours** ≈ **3.8 wall-hours**
  on 256 cores (R1 dominates; add a ~180 ch Core-tier R1 first for bring-up if desired).
- **Full matrix** (R0–R8): **≈ 2400 core-hours** ≈ **9.4 wall-hours** on 256 cores, or
  **≈ 4.7 h** on 512 cores (4 nodes) — the eight runs are independent and embarrassingly
  parallel across nodes.

**Levers if the budget is tight:** halve ambient ppc (×0.5 cost, watch upstream noise §5);
run R1 at Core tier only (~180 ch) to confirm formation criteria and skip the downstream
asymptote; use the one-sided-wall geometry (§2) to halve `n_cell`; drop R6–R8 (robustness
extras). **Levers for fidelity:** raise ambient ppc to 200–400 and/or extend R1 to the Late
tier (~2000 ch) for cyclic reformation and RH-downstream relaxation.

---

## 8. Files to author at execution time (not part of this planning doc)

Keep everything under **`KinShock2020/`** so the project is self-contained (only external
dependency: the WarpX executable built in `warpx-cda/build`):

```
KinShock2020/
├── OVERVIEW.md, REPLICATION_PLAN.md, schaeffer2020.pdf, RESULTS.md
├── src/kinshock/            # reusable library (§6): config, units, io, metrics, plotting
├── scripts/                 # thin drivers (§6): make_inputs, run_checks, make_figures, make_movies
├── runs/
│   └── RXX/
│       ├── config.yaml           # per-run source of truth (§6.0a); authored by hand
│       ├── inputs_kinshock_RXX   # deck (§2–§3), GENERATED from config.yaml by make_inputs.py
│       └── diags*/               # WarpX output, written alongside the deck
└── media/
    ├── testing/             # §6.0 bring-up / progress artifacts
    └── RXX/                 # per-run final figures & movies
```

- `config.yaml` per run: authored by hand — the single source of truth every script reads;
  **no physical constants are hard-coded in `src/` or `scripts/`.**
- Decks (`inputs_kinshock_RXX`): §2–§3, **generated from `config.yaml` by
  `scripts/make_inputs.py`** (R1 baseline; R2–R8 are config diffs). Verified against
  `warpx_used_inputs` after the run with `make_inputs.py --verify` (§6.0a).
- `src/kinshock/` + `scripts/`: the new, standalone code of §6 (designed after, but not
  importing, the `warpx-cda` scripts).
- `RESULTS.md`: status/validation log (mirroring `run_shock_1d/SHOCK_PLAN.md §8`).

Build/run per `warpx-cda/CLAUDE.md` (1D app or Python build; Perlmutter for the Full/Late
tiers, a short Core-tier or R0 run locally for bring-up).

---

## 9. Sequencing checklist

- [ ] Scaffold `KinShock2020/{src/kinshock,scripts,runs,media/testing}/`
- [ ] `src/kinshock/{config,units}.py` + a first `runs/R1_phase/R1/config.yaml`; unit test that
      `units.py` reproduces Table I targets (M_A≈14, M_ms≈13, β₀≈0.2) from the config primaries
- [ ] R0 smoke test — `scripts/make_inputs.py` generates the deck from `config.yaml`; deck
      loads, heater+injector active, no NaN, dt·ω_pe sane; `make_inputs.py --verify` confirms
      `warpx_used_inputs` matches the config
- [ ] `scripts/run_checks.py` + R1 Core-tier short run — loaded-state sanity figures to
      `media/testing/`; confirm C_s,ab, v_A, M_A, β₀ (read from config) match Table I
- [ ] `src/kinshock/{io,metrics,plotting}.py` + analyses A–F wired into `scripts/make_figures.py`
      and `make_movies.py`, validated on R1 (Core) output, figures to `media/R1_phase/R1/`
- [ ] **R1 Full run (→ t\*₃, ~730 ch / ~2.9 h @256 cores)** — acceptance criteria (§6) evaluated
- [x] R2, R3 negative controls (Core tier), 2026-07-26 → `media/R2_phase/R2`, `media/R3_phase/R3` (RESULTS 2026-07-26):
      R3 (n_e0=0) clean null (G=0, empty ambient, no shock); R2 (B₀=0) no *magnetized* shock (no
      coherent B ramp), though ambient ions are electrostatically piston-accelerated and the coded
      7-criteria false-positive (mechanism-blind G>v_sh + B_comp/B₀ with self-fields ≫ weak B₀)
- [ ] (If core passes) R4–R5 Mach scan → t*/z* invariance (Fig. 6)
- [ ] (Optional) R6 collisional, R7 mass-ratio 400, R8 CH multi-species
- [ ] `KinShock2020/RESULTS.md` written; verdict on module verification
