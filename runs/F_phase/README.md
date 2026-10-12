# F_phase — is the false shock (the notch) a resolution artifact?

The user's by-eye verdicts on 34 runs (RESULTS 2026-10-12) found that at 470 eV the notch
is present by t·ω_ci0 ≈ 0.3 whenever dz/λ_D ≥ 1.5, whatever the ppc, while dz/λ_D ≤ 0.76
only delays it past 0.3. Every 47 keV-class run is clean, but every one is also on a fine
grid. Each run here is a single-variable twin of a judged parent: only dz changes (plus the
box and clock needed to afford it).

| run | parent (verdict) | dz/λ_D | N_D | box | to t·ω_ci0 | cells × steps | cost |
|---|---|---|---|---|---|---|---|
| `fs_470eV_dz0379_t1` | `R1_470eV_s3f8_t2` (0.758; notch after 0.3, onset ~0.36) | **0.379** | 264 | 17.1 d_i0 | 1.0 | 101,984 × 7.95 M | ~33 GPU-h |
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
