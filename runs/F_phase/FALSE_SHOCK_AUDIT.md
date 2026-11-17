# False shock (the "notch"): audit and mechanism candidates — 2026-11-16

Shared working document for the two sessions on this question. All verdicts are the user's, by eye.

**The feature.** In the ion phase space at the front of the C's lower branch there is a plateau at
~0.9 v_sh, then a vertical cliff, which never detaches (RESULTS 2026-09-03). In the fine run it
grows smoothly until the *accelerated upstream* nearly meets the *upstream ions kicked up by the
piston* (user, 2026-10-15).

## 1. What has been varied

| variable | runs | result |
|---|---|---|
| ppc 25 → 1600 at fixed dz | ss_dz1 ×7, ss_dz4 ×4 | no effect |
| collisions off, filter/shape 3, smoothing, int1, symmetric wall | ss_dz1 variants, s3f8 | no effect |
| upstream dz/λ_D: 6.1 → 0.38 at 470 eV | ss_*, R1_470eV_s3f8_t2, fs_470eV_dz0379_t1 | **sets onset only**: by 0.3 if ≥ 1.5, later if ≤ 0.76, still forms at 0.38 |
| ablation dz/λ_D,ab 1.24 → 0.62 | R1_470eV_s3f8_t2 → fs_470eV_dz0379_t1 | still forms |
| coarse grid at 47 keV, dz/λ_D,ab up to 4.9 | fs_47keV_dz1p5, fs_47keV_dz3p0 | clean |
| mass convention: heavy e⁻ (IM), real µ = 1836 (RMR) | im_470eV_t03, um1836_470eV | both form; their 47 keV twins are clean |
| upstream T 10 → 100 eV (β₀ 2, collisionless) | fs_470eV_T100_t1 | **still forms** (2026-11-16) |
| T_e,ab 470 → 26 eV, M_A 3.3 | vA_lowTab | forms |
| B₀ ×0.1 at 47 keV (M_A 139) | mb_r100_t2 | clean |
| hybrid (fluid e⁻, quasineutral) | H3_470eV_dense, hy_470eV_t2 | no false shock (but the piston is ~1.5× slower, with no heat flux) |
| dz/λ_D,ab 0.33 (workbook) | fs_470eV_dz0202_t1 | running (chain 61436281–84) |

## 2. What it is not

- **Particle noise:** ppc is irrelevant.
- **Collisions, filtering, particle shape, the wall:** each tested, no effect.
- **A pure grid artifact:** coarse 47 keV grids stay clean, and the finest 470 eV grids still
  form it.
- **Debye resolution, either the upstream or the ablation λ_D:** the 47 keV runs have
  identical λ_D ratios and are clean; at 470 eV the finest grids still form it.
- **The timestep:** ω_pe dt at 470 eV is 0.014–0.04, 10–80× finer than the clean runs.
- **A cold or collisional upstream:** the 100 eV run still forms it.
- **The mass convention:** heavy-electron and real-µ runs both form it.

## 3. The one variable that separates the runs: ε_ab = v_te,ab/c

| run | ε_ab | ε₀ | c/v_A | verdict |
|---|---|---|---|---|
| vA_lowTab | 0.0071 | 0.0044 | 1000 | forms |
| im_470eV_t03 | 0.0071 | 0.0010 | 4283 | forms |
| 470 eV family (all grids), um1836_470eV, T100 | 0.030 | 0.001–0.014 | 1000–4283 | forms |
| im_47keV_t03 | 0.071 | 0.010 | 428 | clean |
| 47 keV family, um1836_47keV, mb_r100_t2 | 0.303 | 0.044 | 100–1000 | clean |
| PoC twins (CH, Z = 3.5, B₀(z), M_A 6.4) | 0.19 | 0.032 | 59 | forms: **the exception**, different setup |

- **ε₀ does not separate.** T100 (ε₀ 0.014) forms it while im_47keV (ε₀ 0.010) is clean.
- **c/v_A does not separate.** mb_r100_t2 (c/v_A = 1000) is clean.
- **ε_ab does**, within the Schaeffer family: it forms at ≤ 0.030 and is clean at ≥ 0.071.
- The grid changes *when* it forms, and only at low ε.
- **Context:** PSC's reduced c puts the paper at ε_ab = 0.30, which is the clean side. A
  real-c 470 eV plasma (the experiment) sits at 0.03.

## 4. Mechanism candidates

**M1 — physical, electrostatic ion–ion streaming in the low-ε regime (favoured).**
- The ablation front's ambipolar E accelerates ions ahead of the piston (the "accelerated
  upstream"). The magnetic piston reflects upstream ions (the "kicked-up" population).
- Where the two overlap they form counter-streaming cold ion beams, which are unstable to the
  electrostatic ion–ion mode.
- The fastest-growing wavelength is λ ≈ 2π v_rel/ω_pi = 2π (v_rel/c) d_i, which scales ∝ ε:
  - at 470 eV, v_sh/c = 0.014, so λ ≈ 0.09 d_i, and many wavelengths fit inside the overlap
    layer;
  - at 47 keV, v_sh/c = 0.14, so λ ≈ 0.9 d_i, comparable to the layer, and the mode is
    swallowed by the ion-scale magnetic structure.
- The saturated mode is an electrostatic (ion-acoustic) shock that never detaches, which is
  the notch.
- Fits:
  - the onset at the meeting of the two populations;
  - the ε threshold;
  - the 100 eV run (beam energy ≫ 100 eV, so the beams stay cold relative to their drift);
  - the hybrid's absence of the feature (quasi-neutral, no electron Landau physics), though
    that comparison is weak;
  - the grid setting the onset: the ion-acoustic branch reaches kλ_D ~ 1, a few cells at
    dz/λ_D ≥ 1.5, so a coarse grid adds dispersion and aliasing and seeds it earlier.
- If M1 is right, the false shock is **real physics of the real-c plasma** that the paper's
  reduced c suppresses. It is not an artifact to be removed.

**M2 — numerical, accumulated per plasma period.** Explicit PIC heats (and makes noise) at a
rate set per ω_pe⁻¹ and by dz/λ_D. Per ion time 1/ω_ci0 there are ~ω_pe/ω_ci0 plasma periods,
roughly 10× more at 470 eV than at 47 keV, with 10× more steps. That fits the ε dependence and
the dz-sets-onset result. Against it: ppc should matter for heating and does not, and the 0.38
grid is well inside the low-heating regime (H_phase).

**M3 — missing or incorrect electron physics in the ablation source at low ε.** This covers
heater and injector coupling at ω_pe dt ≪ 1 (the heater acts every 20 steps), and cold injected
electrons at λ_D = 0.68 nm. Unlikely: the same ratios hold in the clean runs.

## 5. Discriminating tests, cheapest first

| # | test | cost | separates |
|---|---|---|---|
| A | Measure the notch's wavelength and growth in `fs_470eV_dz0379_t1` against the linear ion–ion dispersion of the two measured beams | analysis only | M1 yes/no |
| B | Heating ledger: upstream and slab T_e, T_i vs t·ω_ci0 from EP/PN, 470 eV vs 47 keV, before onset | analysis only | M2 |
| C | ε threshold: rerun `es_1p5keV` (ε_ab 0.054), `es_4p7keV` (0.096), `es_15keV` (0.17) to t·ω_ci0 ≈ 1 with the phase-space histograms (they stopped at 0.3, where even 470 eV is still clean) | cheap; higher T means coarser cells | M1 predicts a threshold near 0.04–0.07, independent of dz |
| D | 470 eV with B₀ → ~0 (no magnetic reflection) | ~production cost | whether the piston reflection is needed |
| E | `fs_470eV_dz0202_t1` (running) | ~30 node-h | closes ablation-Debye resolution at the workbook's standard |
