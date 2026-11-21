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
| hybrid WITH electron heat flux (Braginskii Spitzer, f 0.06), 470 eV and 47 keV | hc_470eV_t2, hc_47keV_t2 | **no false shock** (user, 2026-11-16). Piston at 86 % of kinetic speed; heat reaches the front (T_e there ×3.9) |
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
- **Electron heat transport or a slow piston:** the conducting hybrids (piston at 86 % of
  kinetic, front T_e ×3.9) are clean at 470 eV. What the kinetic run adds is electron-kinetic:
  finite λ_D and charge separation, non-Maxwellian electrons, electron-scale instabilities, and
  the explicit-PIC numerics that come with them (M1 and M2 both survive; M3 is weakened).

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
- **Caveat on the threshold (sister session, 2026-11-16).** The "clean at 0.071" point is
  `im_47keV_t03`, which stopped at t·ω_ci0 0.30, where the fine-grid 470 eV runs are also
  still clean. It therefore cannot place a threshold, and nor does the T100 vs im_47keV
  argument against ε₀. **Every run judged clean past 0.3 has ε_ab = 0.303 (θ_e,ab = 0.092)
  exactly.** Without the short run, the PoC twins (ε_ab 0.19) are not an exception: the
  ordering is monotone, with formation at ε_ab ≤ 0.19 and clean only at 0.30. The threshold
  could be anywhere in 0.19–0.30. Test C now runs at 0.24 and 0.17 (`ft_30keV_t2`,
  `ft_15keV_t2`; `runs/F_phase/README.md`).
- **`ft_30keV_t2` (ε_ab 0.24, to t·ω_ci0 2.0) is clean** (user, 2026-11-16). It is the first long
  clean run at any θ other than 0.092, and it narrows the threshold to ε_ab 0.19–0.24, if the
  PoC run belongs on the same axis. `ft_15keV_t2` (0.17) is the test of that.
- **`ft_15keV_t2` (ε_ab 0.17, to 2.0) is ALSO clean** (user, 2026-11-16). So the PoC twins (0.19,
  form it) are **not** on the ε_ab axis: their false shock needs its own cause among Z = 3.5, the
  applied B₀(z) gradient, M_A 6.4 and n_target/n_e0 400. **Within the Schaeffer family the
  threshold lies in ε_ab 0.030 (forms) to 0.17 (clean)**, which is consistent with the original
  "clean at ≥ 0.071" placement, though `im_47keV_t03` (0.071, to 0.3 only) still does not pin
  it. Next rungs on the same grid: 4.7 keV (ε_ab 0.096, ~3 GPU-h) and 1.5 keV (0.054, ~16 GPU-h).
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
| C | ε threshold: rerun `es_1p5keV` (ε_ab 0.054), `es_4p7keV` (0.096), `es_15keV` (0.17) to t·ω_ci0 ≈ 1 with the phase-space histograms (they stopped at 0.3, where even 470 eV is still clean) | cheap; higher T means coarser cells | M1 predicts a threshold near 0.04–0.07, independent of dz. **Partly run as `ft_30keV_t2` / `ft_15keV_t2`** (ε_ab 0.24 / 0.17, to 2.0, `ek_47keV_t2`'s grid), 2026-11-16 |
| D | 470 eV with B₀ → ~0 (no magnetic reflection) | ~production cost | whether the piston reflection is needed |
| E | `fs_470eV_dz0202_t1` (running) | ~30 node-h | closes ablation-Debye resolution at the workbook's standard |

## 6. Assignments (user, 2026-11-18: "tackle A B and C across the 2 sessions")

| test | owner | status |
|---|---|---|
| A — linear ion–ion (Penrose) test at the front, `scripts/penrose_ionion.py` | psc-schaeffer2020 session | **done 2026-11-18: does NOT discriminate** (see below) |
| B — heating ledger (particle dumps: upstream, slab and beam-region T_e, T_i) | psc-schaeffer2020 session | **done 2026-11-18: M2 falsified**; the fast-ion front speed orders with ε (see below) |
| C — ε ladder: 30 and 15 keV clean (sister); remaining rungs 4.7 keV (ε_ab 0.096, ~3 GPU-h) and 1.5 keV (0.054, ~16 GPU-h) | **theta-ladder (sister) session** | in progress |

Write results into §1/§3 and RESULTS.md; keep this table current.

### Test A result (2026-11-18): the linear ES ion–ion instability does not discriminate

`scripts/penrose_ionion.py` applies the Penrose criterion with electron shielding to each z column
of the binned phase space. F(v) is the normalized ion v_z distribution and T_e,loc is the
electron u_z variance. S = P·T_e/m_i > 1 means unstable; the unstable band ends at k_c. The table
takes the max over the dense columns, every 25th dump. The max always sits within ~±1 d_i0 of
v_sh t.

| run | ε_ab | verdict | S_max | 2π/k_c [d_i0] |
|---|---|---|---|---|
| fs_470eV_dz0379_t1 | 0.030 | forms | 1.6–6.4 | 0.009–0.025 |
| ft_15keV_t2 | 0.17 | clean | 1.7–4.8 | 0.035–0.08 |
| fs_47keV_dz1p5 | 0.30 | clean | 3.3–8.9 | 0.04–0.08 |

- **Clean runs are at least as unstable.** T_e,loc counts counter-streaming electron drift as
  temperature, most at 470 eV (1–3 keV against T_ab = 470 eV). With T_ab instead, 470 eV drops to
  S ≈ 0.5–1.3 and the clean runs stay above 1.5.
- **So M1 in its linear form is falsified.** The front ion distribution is electrostatically
  unstable in every regime. Whatever selects the false shock at low ε acts after linear onset:
  saturation, trapping, or electron kinetics, consistent with the conducting hybrids being
  clean.
- **The mode's scale tracks ε only weakly.** 2π/k_c is ~3× shorter at 470 eV, not 10×. It is
  ~120 cells at the fine run's dz but ~8 cells at the coarse ss_dz1 grids, which may be why a
  coarse grid only moves the onset.

### Test B result (2026-11-18): no heating signature; the fast-ion front speed orders with ε

The script is `scripts/heat_ledger.py` (plot: `scripts/plot_ledger.py`), run on the particle
dumps of 7 runs. T is the 3-component variance of p/mc. The upstream window is
z ∈ [L−4, L−2.5] d_i0, cut off once ions with v_z > 0.3 v_sh exceed 1% of it. The slab is
[0, 2 d_i,ab]. Figure: `media/F_phase/heating_ledger.png`.

| run | verdict | T_e,up/T₀ @ t·ω_ci0 0.56 | T_i,up/T₀ | slab T_e/T_e,ab | slab T_i/T_e,ab |
|---|---|---|---|---|---|
| fs_470eV_dz0379_t1 | forms | 1.54 | 1.34 | 2.23 | 0.25 |
| fs_470eV_T100_t1 | forms | **1.01** | **0.99** | 2.24 | 0.25 |
| ft_15keV_t2 | clean | **1.99** | 1.38 | 2.35 | 0.23 |
| ft_30keV_t2 | clean | 1.56 | 1.17 | 2.38 | 0.22 |
| ek_47keV_t2 / fs_47keV_dz1p5 | clean | 1.09 / 1.02 | 1.01 / 0.99 | 2.43 / 2.44 | 0.22 |

- **M2 is falsified.** The upstream heating does not order with the verdict: T100 forms the
  notch with no upstream heating, and 15 keV heats the most yet is clean. Per plasma period the
  470 eV upstream heats far *less*: ~1.1 T₀ by ω_pe,ab t = 2×10⁴, where the 47 keV runs are at
  ~6 T₀.
- **The slab is the same within 10%.**
- **Beam region:** slab edge to the fast-ion front, in 0.1 d_i0 bins about the local flow
  (`scripts/beam_te.py`). T_e/T_i,core runs from ~30–45 at t·ω_ci0 0.08 down to 15–18 at 0.56 in
  every run, so it does not discriminate either.
- **What does order with ε is the fast-ion front speed** (99.5th-percentile z of ambient ions
  with v_z > 0.3 v_sh, slope over t·ω_ci0 0.24–0.56):

  | ε_ab | front speed |
  |---|---|
  | 0.030 (470 eV fine / prod) | **1.58 / 1.61 v_sh** |
  | 0.17 (15 keV) | 1.35 v_sh |
  | 0.24 (30 keV) | 1.30 v_sh |
  | 0.30 (47 keV) | 1.24 v_sh |

  The "accelerated upstream" is driven harder at low ε, monotonically. That points at the
  acceleration mechanism at the piston or ablation front (ambipolar or Debye-scale fields)
  rather than at the instability, consistent with test A's "acts after linear onset". T100's
  hot tail crosses the 0.3 v_sh cut, so it is excluded from this metric.

### Correction (2026-11-19): hybrids DO show the feature above a heater threshold; it is a different cutoff

User verdict, relayed by the "hybrid 470eV run diff" session:
- The false-shock-like feature IS present in every higher-heater hybrid: hy_470eV_t2_lo / mid /
  hi, hc_470eV_t2, hc_47keV_t2 and hc_470eV_dense_th63.
- It is absent in the low-heater ones: H3_470eV_dense and hc_470eV_dense.
- This **supersedes** the 2026-11-16 "conducting hybrids: no false shock" line in §1, and the
  inference that the feature needs electron-kinetic physics.

Per that session, heat conduction and box size do not set it; the heater does. The hybrid cutoff
sits at peak slab T_e/(m_i v_A²) between 11.8 and 22. It is not in θ_heat: the 47 keV twin, at
the same ratio, shows the same verdict.

The kinetic runs do not follow that cutoff (`notch/kin_ratio.py`):
- **Along the θ ladder, T_e,ab/(m_i v_A²) = β_ab·n_amb/n_ab = 9.2 is fixed.** The measured slab
  T_e/(m_i v_A²) is 20.5–22.5 both for the forming 470 eV runs and for the clean 15/30/47 keV
  runs: same point on the hybrid axis, opposite verdicts.
- **Off the ladder it runs the wrong way:** vA_lowTab forms at 0.50, mb_r100_t2 is clean at 920.

So the kinetic false shock (beam-front ion-acoustic, set by v_beam/v_te; tests A and B) and the
hybrid feature (set by heater and B) are, on present evidence, **different mechanisms with a
similar phase-space look.**

### Hybrid discriminant (2026-11-20, user verdicts via the hybrid session; commit 62a6354)

At a fixed heater, the upstream B decides the hybrid feature:
- hc_470eV_dense_th63_B15 (B₀ ×1.5, slab T_e unchanged at 1140 eV): ratio 9.9, **clean**.
- hc_470eV_dense_th40 (θ 4e-4, 746 eV): ratio 14.6, **clean**.
- The parent hc_470eV_dense_th63: ratio 22.1, forms.

So the hybrid cutoff is T_e/(m_i v_A²) = 14.6–22.1 (c_s/v_A 4.9–6.1). The kinetic θ ladder sits
fixed at 9.2 (slab ≈ 21) while its verdict flips with ε, so "no common cutoff" stands.
Details: runs/hybrid_phase/README.md.

### Hybrid feature redefined, and mass-ratio independent (2026-11-20/21, user via the hybrid session)

1. **What the hybrid "feature" is** (user, by eye). It is a *transient* common to all four
   B-bracket runs (th63, B11, B12, B15). The ion C at the piston edge gets much harsher corners
   for its reflecting ions, enclosing an empty hole in phase space. All four are identical at
   t·ω_ci0 ≈ 1.2 and settle to identical states by ≈ 3. Only the transient's **size and
   lifetime** differ, clearest at t* ≈ 2.14, shrinking with R = T_e/(m_i v_A²): th63 > B11 >
   B12 > B15 (barely). So the hybrid feature is **not** the kinetic plateau-plus-cliff notch. It
   is a graded transient of the piston-braking / ramp-separation phase, consistent with the
   pressure-balance reading: minimum edge M_A 4.3 / 3.2 / 2.4 / 1.5.
2. **Mass ratio does not matter.** hc_470eV_dense_th63_mp (m_i/m_e 1836, same n, B₀, T and R;
   conduction ~4.3× faster relative to the ions) is identical to th63. So neither m_i/m_e nor
   conduction strength sets it.
3. **Why LaserProdShock P8 shows it unsettled.** Its laser runs accumulate only 0.8–1.7 local
   ion gyro-times by 4.5 ns, because the coil field falls 7.5 → 2 T along the path. They stay in
   this transient phase.
