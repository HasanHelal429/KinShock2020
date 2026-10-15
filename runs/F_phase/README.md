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
| `fs_470eV_dz0379_t1` | done 2026-10-14 (36 h); movies `media/F_phase/fs_470eV_dz0379_t1/` | awaiting the user's verdict |
| `fs_470eV_T100_t1` | staged 2026-10-15 | — |

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
