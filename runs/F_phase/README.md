# F_phase — is the false shock (the notch) a resolution artifact?

The user's by-eye verdicts on 34 runs (RESULTS 2026-10-12) found that at 470 eV the notch
is present by t·ω_ci0 ≈ 0.3 whenever dz/λ_D ≥ 1.5, whatever the ppc, while dz/λ_D ≤ 0.76
only delays it past 0.3. Every 47 keV-class run is clean, but every one is also on a fine
grid. Each run here is a single-variable twin of a judged parent: only dz changes (plus the
box and clock needed to afford it).

| run | parent (verdict) | dz/λ_D | N_D | box | to t·ω_ci0 | cells × steps | cost |
|---|---|---|---|---|---|---|---|
| `fs_470eV_dz0379_t1` | `R1_470eV_s3f8_t2` (0.758; notch after 0.3, onset ~0.36) | **0.379** | 264 | 17.1 d_i0 | 1.0 | 101,984 × 7.95 M | ~33 GPU-h |
| `fs_470eV_T100_t1` | `R1_470eV_s3f8_t2`: only **T₀ 10 → 100 eV** (2026-10-15) | 0.240 (via T₀) | 417 | 17.1 d_i0 | 1.0 | 50,992 × 3.98 M | ~10–13 GPU-h |
| `fs_470eV_dz0202_t1` | `fs_470eV_dz0379_t1` (notch): only dz → the workbook's 0.33 λ_D,ab (2026-11-16) | **0.202** | 495 | 17.1 d_i0 | 1.0 | 191,224 × 14.91 M | ~100–130 GPU-h |
| `fs_47keV_dz1p5` | `ek_47keV_t2` (0.758; clean to 2.0) | **1.516** | 66 | 55 d_i0 | 2.0 | 8,200 × 39,757 | minutes |
| `fs_47keV_dz3p0` | `ek_47keV_t2` | **3.029** | 33 | 55 d_i0 | 2.0 | 4,104 × 19,898 | minutes |

**Caveats.** `fs_470eV_dz0379_t1`'s front reaches the box end near t·ω_ci0 ≈ 1.2, so a
restart can stretch it that far, but not to 2.0. At 47 keV a grid coarse in λ_D is also
coarse in d_e, because λ_D/d_e = ε: dz is 0.75 and 1.5 d_e,ab here, against 0.075–0.30 in
the 470 eV notch runs. dz ×8 (dz/λ_D 6) is omitted because ω_pe,ab·dt = 2.25 > 2. ppc is
held, so N_D moves with dz; the verdicts showed ppc does not matter at fixed dz.

**Reading.**

| exp 1 (fine, 470 eV) | exp 2 (coarse, 47 keV) | conclusion |
|---|---|---|
| clean to 1.0 | notch | resolution artifact: coarse grids make it, fine grids prevent it |
| notch, later | clean | low-ε physics; the grid only sets the onset time |
| notch, later | notch | both |
| clean to 1.0 | clean | a fine grid suppresses a low-ε effect |

Launch with `perlmutter/submit_chain.sh` from the main checkout (all three carry checkpoints).
Judge them by eye from the phase-space movies (`scripts/make_phase_hist_movie.py`).

## Verdicts

| run | status | verdict (by eye) |
|---|---|---|
| `fs_47keV_dz1p5` | done 2026-10-12, 4 min | **no false shock** |
| `fs_47keV_dz3p0` | done 2026-10-12, 2.5 min | **no false shock** |
| `fs_470eV_dz0379_t1` | done 2026-10-14 (36 h); movies `media/F_phase/fs_470eV_dz0379_t1/` | **false shock still forms**. Evolution is smooth until the accelerated upstream almost meets the kicked-up upstream ions (user, 2026-10-15) |
| `fs_470eV_T100_t1` | done 2026-10-17 (10 h 37 m, job 59831119); movies `media/F_phase/fs_470eV_T100_t1/` | **false shock still forms** (user, 2026-11-16): not a cold/collisional upstream |
| `fs_470eV_dz0202_t1` | running, chain 61436281–84 | — |
| `ft_30keV_t2` | done 2026-11-16 (13 m 35 s, debug job 61437037); movie `media/F_phase/ft_30keV_t2/` | **clean** (user, 2026-11-16): ε_ab 0.24 is on the clean side |
| `ft_15keV_t2` | done 2026-11-16 (30 m 46 s, shared job 61437039); movie `media/F_phase/ft_15keV_t2/` | **clean** (user, 2026-11-16): ε_ab 0.17 is also clean, so PoC (0.19, forms) is not on this axis |

## θ_e,ab ladder (`ft_30keV_t2`, `ft_15keV_t2`, 2026-11-16) — audit test C

Twins of `ek_47keV_t2` (θ_e,ab 0.092, clean to t·ω_ci0 2.0) with **only T_e,ab moved**, every
Table I group held the way the E_phase ladder holds them: T₀/T_e,ab, β_ab (B₀ ∝ √T), v_p/C_s,ab,
dz/λ_D 0.758 (dz ∝ √T), the 55 d_i0 box and the 2.0/ω_ci0 clock. `units.derive` gives M_A 13.95,
M_ms 12.76, β₀ 0.196 for both, identical to the parent. The parent's numerics are kept: shape 3,
filter 8, ppc 100, λ_ab = 20. There is one box on one GPU and one checkpoint at the end.

| run | T_e,ab | θ_e,ab | ε_ab | n_cell × steps | × parent |
|---|---|---|---|---|---|
| `ek_47keV_t2` (parent) | 47.0 keV | 0.092 | 0.303 | 16,400 × 79,515 | 1 |
| `ft_30keV_t2` | 30.0 keV | 0.0587 | 0.242 | 20,528 × 124,593 | 2.0 |
| `ft_15keV_t2` | 14.86 keV (= `es_15keV`) | 0.0291 | 0.171 | 29,168 × 251,540 | 5.7 |

**Why.** Sorted by θ_e,ab, every judged run is in order: the false shock forms at θ ≤ 0.037
(PoC included) and is clean only at 0.092. **Every run judged clean past t·ω_ci0 0.3 is at
θ = 0.092 exactly.** The one clean run below that, `im_47keV_t03` (θ 0.005, ε_ab 0.071), stopped
at 0.30, where the fine-grid 470 eV runs are also still clean, so it cannot place a threshold.
These two rungs sit between PoC (0.037) and the parent (0.092) and just below PoC.
**Reading:** if 15 keV forms it and 30 keV is clean, the threshold is a property of the
electron temperature in rest-mass units (ε). That is the condition M1 needs (FALSE_SHOCK_AUDIT
§4), and it gives a ~1 GPU-h model of the false shock. If both are clean, the PoC run's notch
needs its own explanation. If both form it, the paper's θ = 0.092 is close to the edge.

Experiment 2 rules out "pure resolution artifact": coarseness alone does not make it at
47 keV (RESULTS 2026-10-12, later).

## Hot-upstream test (`fs_470eV_T100_t1`, 2026-10-15)

The user's hypothesis is that the notch comes from the cold 10 eV upstream: either it is
under-resolved or an instability grows in it. The test raises only `plasma.ambient.theta_0`
(T_e0 = T_i0) from 10 to 100 eV. dz, dt, ppc, n, B₀, the piston, the heater, collisions,
shape 3 and filter 8 are all the parent's. Box and clock are as `fs_470eV_dz0379_t1`. What moves
with T₀ (`units.derive`):
- β₀: 0.20 → 1.96
- M_ms: 12.8 → 8.1 (M_A is unchanged at 13.95)
- M_s: 31.5 → 10.0
- dz/λ_D,amb: 0.758 → 0.240
- 2 v_sh/v_te0: 6.3 → 2.0
- upstream λ_ii/d_i0: 0.014 → 1.43

The last row is the confound: a clean result cannot separate "warmer" from "less
collisional". The collisionless 10 eV twin was offered and declined. Note also that the clean
47 keV runs share the 470 eV run's T₀/T_e,ab (0.0213), β₀, Mach numbers and upstream mfp.
They differ only in v_te/c (×10). So this test probes temperature sensitivity; it cannot
reproduce the 47 keV condition, which would need T₀ ≈ 1 keV (β₀ ≈ 20). The hot upstream drains
through the open hi boundary at C_s0 ≈ 0.0014 c, about 1.4 d_i0 inward by t·ω_ci0 = 1.0, while
the front is near 14 d_i0. Watch the box end in the last frames.

**Outcome of experiments 1 + 2.** The fine grid still forms the notch, and the coarse 47 keV
grids stay clean. That is the "notch, later / clean" row: **low-ε physics, with the grid
setting at most the onset time.** A pure resolution artifact is ruled out from both sides.
Whether the fine run's onset is actually later than the parent's ~0.36 has not been read
yet. The user's description points at an **interaction between two upstream ion
populations**: nothing happens until the ions accelerated ahead of the piston nearly meet the
upstream ions the shock kicks up. That is the setting of an ion–ion streaming instability. A
cold upstream (v_ti0/v_sh = 0.032 at 10 eV) is the condition it needs, and the 100 eV twin
(v_ti0/v_sh = 0.10, C_s0 ×√10) is the direct test.
