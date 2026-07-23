# Overview — Schaeffer et al. 2020, *Kinetic Simulations of Piston-Driven Collisionless Shock Formation in Magnetized Laboratory Plasmas*

**Citation.** D. B. Schaeffer, W. Fox, J. Matteucci, K. V. Lezhnin, A. Bhattacharjee, and
K. Germaschewski, *Phys. Plasmas* **27**, 042901 (2020). doi:10.1063/1.5123229.
Full text (DOE accepted manuscript) saved locally as `schaeffer2020.pdf`.

**One-line summary.** Quasi-1D PIC simulations (code **PSC**) of a *laser-ablated piston
plasma* expanding into a pre-magnetized ambient plasma, using a kinetic laser-ablation
model, to establish robust, experimentally-observable **signatures of collisionless
shock formation** in the regime where the shock is still forming and cannot be cleanly
separated from the driving piston.

This document is the physics reference for our WarpX port of the heating operator; the
companion `REPLICATION_PLAN.md` uses it to define the runs and analysis that verify the
WarpX module reproduces these shocks.

---

## 1. Motivation and context

- **Collisionless shocks** thermalize supersonic flows through electromagnetic fields on
  scales far shorter than the collisional mean free path. Above a critical Mach number
  (M_A,crit ≈ 3 for perpendicular geometry), *supercritical* magnetized shocks dissipate
  energy in part by **reflecting incoming ions** at the ramp. Examples span planetary bow
  shocks, the heliospheric termination shock, supernova remnants, and GRBs.
- **Laser laboratory experiments** recreate scaled versions. Two drive types: *piston-driven*
  (the laser-ablated plasma acts as a supersonic piston that sweeps up and accelerates a
  pre-existing magnetized ambient plasma — the focus here) and *obstacle-driven*.
- **The problem this paper addresses.** Experiments are limited to spatio-temporal scales
  *of order the shock-formation scale itself*, so the piston and the shock overlap and are
  hard to distinguish. Classical **Rankine–Hugoniot (RH) jump conditions** describe the
  shock structure only *long after* formation, when a well-defined downstream exists — they
  do **not** apply during formation. New criteria are needed.
- **Goal.** Provide (a) a model for piston and shock speeds when there is no single well-defined
  inflow speed, (b) a set of **criteria to judge shock formation** without RH, and (c) the
  characteristic timescales, profiles, and phase-space signatures observable in experiments.

---

## 2. The PSC laser-ablation (piston) model — what our WarpX module ports

PSC is a fully electromagnetic, relativistic, massively-parallel PIC code. Instead of
simulating the laser and solid target from first principles, it uses a **kinetic ablation
model** applied to a dense "target" slab — this is precisely the surrogate our WarpX
`ParticleHeater` + `TargetInjector` modules reproduce:

1. **Localized heating operator** — mimics laser heating. Slab electrons receive Monte-Carlo
   momentum-space kicks that relax the slab toward the **ablation electron temperature
   T_e,ab**, launching a supersonic ablation/piston plume. (In our port:
   `H = fac·shape`, `fac = 8 θ_e^{3/2} / (√(m_i/m_e)·width/d_e)`, θ_e = k_B T_e,ab/m_e c².)
2. **Target replenishment** — new particles are continuously added inside the target region
   to **maintain the target density** (here 2.5 n_e,ab) as material ablates away, mimicking
   continuous ablation. (In our port: the density-relaxation `TargetInjector`.)

The ablated plume expands, sweeps up magnetic flux (forming a **magnetic cavity**) and the
ambient plasma, and drives the shock. The model is deliberately run with **reduced electron
parameters** — reduced mass ratio (m_i/m_e = 100) and reduced speed of light — while keeping
**ion-scale parameters matched to experiment**, which is what makes the large,
multi-ion-gyroradius HED system sizes affordable.

**Geometry.** Gridded in the x–z plane, ablation along **ẑ**, uniform background field
**B₀ = B₀ x̂** → strictly **perpendicular** geometry (θ_Bn = 90°, most efficient for shock
generation). Quasi-1D: only a few cells transverse, uniform driving, so the ablation is
effectively planar. Target ablated in both ±z directions; only z ≥ 0 analyzed. Coulomb
collisions via the Takizuka–Abe operator, matched to HED collisionality (λ_ab ≡
ω_ce,ab/ν_ei,ab = 20).

---

## 3. Simulation setup and the representative run (Table I)

Fundamental normalization is built on the **ablation** density/skin-depth (n_e,ab, d_e,ab,
t_ab), with the ambient/upstream defined relative to it. Ion inertial length
d_i,ab = 10 d_e,ab (mass ratio 100).

**Domain & numerics (representative run):**
- Box: **30 000 cells in z × 12 cells in x** = 9000 × 5 d_e,ab (= 900 × 0.5 d_i,ab).
- Grid: dz ≈ 0.3 d_e,ab. **1000 macroparticles/cell** at n_e,ab = 1.
- **400 000 timesteps**; heating applied for the whole run. Energy conserved to within 2%.
- Ambient plasma + field fill z > 2 d_i,ab; the **target** (density 2.5 n_e,ab) fills
  0 < z ≤ 2 d_i,ab; x and z boundaries periodic (run stopped before particles reach +z).

**Table I — representative quasi-1D run** (sim / dimensionless values, and one HED-relevant
physical realization):

| Quantity | Symbol | Sim value | Physical value |
|---|---|---|---|
| x size | L_x | 0.5 d_i,ab | 4.7 µm |
| z size | L_z | 900 d_i,ab | 8.4 mm |
| sim time | τ_sim | 220 t_ab | 10.9 ns |
| **Ablation** | | | |
| charge state | Z_ab | 1 | 1 |
| electron density | n_e,ab | 1.25 | 6×10²⁰ cm⁻³ |
| electron temperature | T_e,ab | 0.092 m_e c² | 470 eV |
| collision frequency | ν_ei,ab | 0.009 t_ab | 0.43 ps |
| sound speed | C_s,ab | 0.030 c | 210 km/s |
| piston speed | v_p | 0.104 c | 730 km/s |
| shock speed | v_sh | 4.6 C_s,ab | 980 km/s |
| **Upstream** | | | |
| magnetic field | B₀ | 0.01 √(m_e c²) | 7 T |
| charge state | Z₀ | 1 | 1 |
| electron density | n_e0 | 0.01 n_e,ab | 4.8×10¹⁸ cm⁻³ |
| temperature | T₀ | 0.002 m_e c² | 10 eV |
| ion inertial length | d_i0 | 11.2 d_i,ab | 104 µm |
| ion gyroperiod | ω_ci0⁻¹ | 33.9 t_ab | 1.5 ns |
| **Dimensionless** | | | |
| mass ratio | m_i/m_e | 100 | |
| speed-of-light ratio | c_sim/c_phys | 0.02 | |
| ablation beta | β_ab | 1150 | |
| upstream beta | β₀ | 0.2 | |
| collisionality | λ_ab | 20 | |
| ion mean free path | λ_mfp/d_i0 | 350 | |
| **Alfvén Mach number** | **M_A** | **14** | |
| **magnetosonic Mach number** | **M_ms** | **13** | |

β ≡ 2 µ₀ n_e T_e / B². Upstream is high-Mach (M_A > 1), strongly magnetized ablation
(β_ab ≫ 1), and low-β upstream (β₀ ≲ 1). The parameter scan (Fig. 1) covers
**T_e,ab = 0.03–0.134 m_e c²**, **β_ab = 30–11840**, **β₀ = 0.008–3.2**, **M_A = 3–57**,
**M_ms = 3–27**, plus control runs with **B₀ = 0** and **n_e0 = 0**.

---

## 4. Key results

### 4.1 Model for piston and shock speeds (Sec. III A)
Because the sim is ablation-initialized, there is no single inflow speed. The ablation
density follows a scale-free profile n_e = (n_e,ab − n_e0)·exp[−(z−z₀)/z₀] + n_e0, giving an
expansion speed (**Eq. 1**):

  v = (C_s,ab/2)·[1 − ln((n_e − n_e0)/(n_e,ab − n_e0))].

- **Piston speed v_p** is tracked via the peak field tied to the **magnetic cavity**; it
  corresponds to where the piston density ≈ 1.35 n_e0. Without an ambient, v_p → 6.5 C_s,ab.
- **Shock speed v_sh** is tracked via the faster-moving magnetic compression in the shocked
  ambient. It asymptotes to **≈ (4/3) v_p** at low field and rises with B₀, consistent with
  the perpendicular RH relation (**Eq. 2**):
  v₁/v₂ = −[M_A² + β̃₀] + [(M_A² + β̃₀)² + 8 M_A²]^{1/2}, β̃₀ = (5/2)(1 + β₀).

### 4.2 Seven criteria for piston-driven shock formation (Sec. III B)
Because RH does not apply during formation, the paper defines observationally-motivated
criteria. A structure is a **shock precursor** if it satisfies (1)–(6); it is a **shock**
once it also satisfies (7):

1. Super-magnetosonic: M_ms = v/√(v_A² + C_s²) > 1
2. Collisionless: L/λ_ii > 1
3. Large density compression: n_e/n_e0 > 2
4. Large magnetic-field compression: B/B₀ > 2
5. Steep ramp: dB/dz and dn_e/dz on the scale ~ d_i0
6. Presence of **magnetically-reflected ambient ions** (ambient ions with v_z > v_sh)
7. **Separation from the piston** (shock precursor separated ≥ ¼ ρ_i0 from the piston's peak field)

Criteria (3)–(5) distinguish a real shock from mere interpenetrating flows / piston
compressions; (6) is the key kinetic dissipation signature; (7) is what makes the shock
dynamically independent of the piston. Compression ratios approach the strong-shock limit of
**~4** at large M_A (and can fall below 2 at very low M_A).

### 4.3 Three characteristic timescales (nearly independent of M_A and β_ab; Fig. 6)
Defined via the reflected-ambient-ion functions F(z,t) = f_a,refl/N_a,tot and
G(t) = N_a,refl/N_a,tot (fraction of ambient ions faster than the shock):

- **t*₁ ≈ 1 ω_ci0⁻¹** — **onset of shock formation** (dG/dt is maximum); location z*₁ ≈ 1 ρ_i0.
  First 6 criteria satisfied. Shock still buried in the piston envelope.
- **t*₂ ≈ 2.5 ω_ci0⁻¹** — **clear separation from the piston** (precursor moves ¼ ρ_i0 past the
  piston field); z*₂ ≈ 2.5 ρ_i0. All 7 criteria satisfied → a *shock*.
- **t*₃ ≈ 5 ω_ci0⁻¹** — **downstream develops** on scales ≫ shock width; z*₃ ≈ 5–6 ρ_i0. Only now
  do RH jump conditions begin to apply.

Here ρ_i0 = v_p/ω_ci0 is the upstream-directed gyroradius of an ion moving at the piston speed.
The onset time is insensitive to the heating-operator duration (≳ 1 ω_ci0⁻¹ of heating suffices).

### 4.4 Structure and evolution (Figs. 2, 5, 7, 8)
- **Early "snowplow"** (region I): the plume sweeps up ambient ions and magnetic flux;
  ambient ions are *piston-accelerated* (not yet reflected) via the ambipolar/in-plane E_z
  field, with a comparable "Larmor-coupling" E_y. **Care must be taken not to mistake
  piston-accelerated ions, or piston compressions, for shock signatures.**
- **Onset of reflection** where the swept-up field is compressed → magnetically-reflected
  ambient ions appear upstream; secondary density/field compressions build a proto-foot.
- **Separation** (t*₂): distinct forming-shock (localized steep n_e and B compressions,
  reflected ions) vs. piston regions become visible in phase space.
- **Late time** (Figs. 10, 11): **cyclic shock reformation** roughly every ~1.5 ω_ci0⁻¹; the
  downstream relaxes to RH-predicted compressions (with overshoots) by ~10 ω_ci0⁻¹.

### 4.5 Parameter dependences (Figs. 6, 9)
- Formation is qualitatively similar across M_A; **piston–shock separation is most pronounced
  at lower Mach number** (shock forms closer to the target for lower M_A / slower shocks).
- **Control runs:** with **B₀ = 0**, ambient ions are still initially accelerated but there is
  **no magnetic compression, no strong ion heating, no secondary compression** — no shock.
  With **n_e0 = 0** (no ambient), there are no ambient-ion structures and no compressions.
  These are the negative controls that prove the signatures are shock-specific.

### 4.6 Robustness (Secs. III D–F, Mass ratio)
- **Multi-species (CH):** each species shocks quasi-independently near its own gyroperiod
  (H at ≈ 0.5 ω_ci0⁻¹, C at ≈ 1 ω_ci0⁻¹); global field/density profiles are dominated by the
  species with the largest m/Z (C).
- **Collisions:** collisional vs. collisionless runs show similar shock formation — ion–ion
  interactions stay collisionless due to high flow speeds (λ_mfp/d_i0 ≫ 1); collisions mainly
  smear sub-d_i0 electron phase-space structure and keep electrons isotropic/Maxwellian.
- **Mass ratio:** m_i/m_e = 100 vs. 400 give the same formation times/locations and
  compressions — shock formation is converged at reduced mass ratio.

---

## 5. Implications for the WarpX verification

The WarpX `ParticleHeater` + `TargetInjector` are a direct port of exactly the PSC
heating/replenishment model above, and are already cross-validated operator-for-operator
against PSC (see `warpx-cda/heating_operator/`). This paper supplies the **physics
acceptance test**: if the WarpX modules can drive a piston that reproduces

- the piston/shock **speed model** (v_sh ≈ (4/3–…)·v_p; Eqs. 1–2),
- the **seven formation criteria** in the right order,
- the **three timescales** (t* ≈ 1, 2.5, 5 ω_ci0⁻¹) and locations (z* ≈ 1, 2.5, 5–6 ρ_i0),
- the **~4× density/field compression** and **reflected-ion** phase-space signature, and
- the **B₀ = 0 / n_e0 = 0 negative controls**,

then the module is demonstrated to generate Schaeffer-2020-class magnetized perpendicular
collisionless shocks. `REPLICATION_PLAN.md` turns this into concrete runs and analysis.

---

*Sources:* [Phys. Plasmas 27, 042901 (2020)](https://pubs.aip.org/aip/pop/article-abstract/27/4/042901/319049) ·
[OSTI 1643946 (full text)](https://www.osti.gov/pages/biblio/1643946) ·
companion model paper [Fox et al., Phys. Plasmas 25, 102106 (2018), arXiv:1712.00152](https://arxiv.org/abs/1712.00152).
