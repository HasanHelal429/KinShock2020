# E_phase — the `eps = v_te,ab/c` ladder

**The question.** `R1_paper` (47 keV) and `R1_paper_470eV` were built to be the *same*
dimensionless problem — Table I of Schaeffer 2020 — yet they produce visibly different
shock structure: an ambient wedge 2.25 d_i0 deep at 470 eV against 1.25 d_i0 at 47 keV.
S_phase spent a full resolution ladder plus six other numerical knobs failing to remove
that difference. This phase asks the other question: **are the two runs actually the same
physics?**

## The audit: exactly one parameter differs

Computed from the config primaries only (`scripts/dimensionless_audit.py`, folded
into RESULTS 2026-08-17). Preserved to <1%:

| group | value |
|---|---|
| `m_i/m_e` | 100 |
| `n_e0/n_e,ab`, `n_t/n_e,ab` | 0.008, 2 |
| `beta_ab` | 1150 |
| `M_A = v_sh/v_A` | 13.952 |
| `M_ms` | 12.74 / 12.76 |
| `v_sh/C_s,ab`, `v_p/C_s,ab`, `v_p/v_A` | 4.6, 3.4288, 10.4 |
| `d_i,ab/d_e,ab`, `d_i0/d_i,ab`, `d_i0/d_e0` | 10, 11.18, 10 |
| `rho_i0/d_i0`, `rho_sh/d_i0` | 10.40, 13.95 |
| `L/d_i0`, `L_target/d_i0` | 80.50, 0.1789 |
| `1/w_ci0 / t_ab = sqrt(beta_ab)` | 33.912 |
| `lambda_ab = mfp_e/d_e,ab` | 20 |
| **`nu_ei,ab / w_ce`** | **1.6956** |
| **`rho_e,ab / d_e,ab`** | **33.912** |

`beta_0` and `T_0/T_e,ab` differ by 2.1%, which is Table I's rounding of `T_0` (10 vs
10.217 eV) and was already documented in the config header. Everything else that differs
is an **integer power of one number**:

    eps == v_te,ab/c = sqrt(theta_e,ab)      0.3033 (47 keV)  ->  0.03033 (470 eV)

| group | scaling | 47 keV | 470 eV |
|---|---|---|---|
| `v_A/c`, `v_sh/c`, `v_p/c`, `C_s,ab/c` | `eps^+1` | 0.0100, 0.1395, 0.1040, 0.0303 | 0.0010, 0.0140, 0.0104, 0.0030 |
| `lambda_D/d_e`, `lambda_D,0/d_i0` | `eps^+1` | 0.3033, 4.47e-3 | 0.0303, 4.42e-4 |
| `w_pe/w_ce`, `w_pi0/w_ci`, `rho_e/lambda_D` | `eps^-1` | 10.0 (upstream), 100 | 100.0, 1000 |
| `theta_e,ab`, `sigma_i`, `sigma_e` | `eps^+2` | 0.092, 1e-4, 1e-2 | 9.2e-4, 1e-6, 1e-4 |
| `lnLambda` (at fixed `lambda_ab = 20`) | `eps^-4` | 1.22e5 | 12.23 |

So the difference is **one-dimensional**. Either the wedge is a function of `eps` — in
which case it is real physics of the normalization, and the 470 eV answer is the one that
describes an actual 470 eV plasma — or it is not, and the cause is somewhere else entirely.

Note what is *not* on the differing list: `rho_e,ab/d_e,ab = 33.9` and `nu_ei/w_ce = 1.696`
are both preserved. Electron magnetization relative to the *inertial* scales, and
collisionality relative to the *field*, are identical in the two runs. What differs is
electron magnetization relative to the *Debye* scale — `rho_e/lambda_D = w_pe/w_ce`, which
is 10 in the paper run and 100 at 470 eV.

## The ladder

Five rungs, `eps` moving by a decade, **every** Table I group and **every** numerics group
held fixed. Rung 1 is the existing `S_phase/ss_dz16_ppc100`; rungs 2–5 are new.

| run | T_e,ab | eps | B0 [T] | `w_pe0/w_ce` | `gamma(T_e,ab)` | cells | steps | wall (1 GPU) |
|---|---|---|---|---|---|---|---|---|
| `S_phase/ss_dz16_ppc100` | 470 eV | 0.0303 | 7.026 | 100.0 | 1.0014 | 71680 | 2400000 | 7.35 h (**done**) |
| `es_1p5keV` | 1486 eV | 0.0539 | 12.494 | 56.2 | 1.0044 | 40312 | 759079 | 1.31 h pred / **1 h 30 m** |
| `es_4p7keV` | 4700 eV | 0.0959 | 22.220 | 31.6 | 1.0138 | 22664 | 239966 | 0.23 h pred / **20 m** |
| `es_15keV` | 14860 eV | 0.1705 | 39.509 | 17.8 | 1.0436 | 12744 | 75885 | 0.04 h pred / **5 m** |
| `es_47keV` | 47012 eV | 0.3033 | 70.273 | 10.0 | 1.1380 | 7168 | 23997 | 0.01 h pred / **1 m** |

Held identical at every rung, verified numerically:
`M_A = 13.952`, `M_ms = 12.759`, `beta_ab = 1150`, `beta_0 = 0.1957`, `rho_i0/d_i0 = 10.40`,
`L = 12.02 d_i0`, run `= 0.302/w_ci0`, `dz/lambda_D,0 = 0.379`, `ppc = 100`,
`N_D,PIC = 264`.

### Why `dz/lambda_D` is the thing held fixed, not `dz/d_e`

`lambda_D ~ eps`, while `d_e`, `d_i0` and the domain are **eps-independent** (they depend
on `n` and the real `c` only). A ladder at the paper's fixed `dz/d_e = 0.3` would therefore
drag `dz/lambda_D` along by the same factor of 10 — reintroducing exactly the confound the
S_phase ladder spent a week excluding. Fixing `dz/lambda_D` and `ppc` instead pins the
finite-grid margin *and* the discreteness noise, and lets `dz/d_e` float from 0.019 to
0.188 (still ≥5 cells per `d_e,ab` at the coarsest rung).

Every rung sits at `dz/lambda_D,0 = 0.379`, **finer than R1_paper's own 0.60**, so no rung
can be dismissed as under-resolved, and `N_D,PIC = 264` against R1_paper's 167.

### Cost

`cost ~ eps^-3` at fixed `dz/lambda_D` (cells `~1/eps`, steps `~1/eps^2`), so every rung
above 470 eV is *cheaper* than the anchor. All four together are **~1.6 GPU-h**.
The estimate is anchored on `ss_dz16_ppc100`'s measured 7.35 h scaled by `cells x steps`;
the two top rungs are small enough that GPU occupancy will be poor and their true cost is
likely a few times the number quoted — irrelevant at ~1 minute each.
**Measured:** 1 m / 5 m / 20 m / 1 h 30 m, i.e. 1.4–2.0x the prediction, exactly the
small-point under-occupancy penalty and in the safe direction. Cell counts were also
snapped to multiples of AMReX's `blocking_factor = 8` after three rungs aborted at init
(`62e7658`, which added the generation-time guard in `deck.py`); `dz/lambda_D` moves by
<0.03% as a result and `es_1p5keV`, which had already run, was not touched.

## RESULT (2026-08-17): the wedge IS a function of eps — see RESULTS.md

| run | eps | `w_pe/w_ce` | `gamma_sh` (measured) | wedge [d_i0] | `T_e,shock` |
|---|---|---|---|---|---|
| 470 eV anchor | 0.0303 | 100.0 | 1.0012 | **2.062** | 98.7 |
| `es_1p5keV` | 0.0539 | 56.2 | 1.0033 | **2.064** | 81.3 |
| `es_4p7keV` | 0.0959 | 31.6 | 1.0089 | **1.808** | 61.6 |
| `es_15keV` | 0.1705 | 17.8 | 1.0237 | **1.293** | 38.0 |
| `es_47keV` | 0.3033 | 10.0 | 1.0586 | **1.164** | 21.6 |
| `R1_paper` | 0.3033 | 10.0 | 1.0715 | **1.119** | 24.9 |

`es_47keV` reproduces `R1_paper` to **4.0%** at the same eps on a different domain,
`dz/d_e`, `dz/lambda_D` and `N_D`. The wedge is set by eps, not by resolution.

**RETRACTED — the shape does NOT discriminate the mechanism.** This README previously
claimed that a flat-then-drop curve would mean relativistic capping and a straight decade
would mean the electron-scale wave regime. That was wrong: eps is ONE parameter, so both
candidates are monotone functions of it and both are evenly spaced in `log eps`
(`gamma-1 ~ eps^2`, `w_pe/w_ce ~ eps^-1`). Fitted, they are comparable (rms 0.141 vs 0.163).

What does argue against relativistic capping is the **lever arm**: the *measured*
`gamma_shocked` only runs 1.0012 -> 1.0586, a 5.7% change in electron inertia, against a
43.6% change in wedge depth. That leaves the electron-scale wave regime
(`rho_e/lambda_D = w_pe/w_ce`, 100 -> 10) as the surviving explanation — an argument from
lever arm, not a direct measurement of the instability.

## Reading the result

Measure the wedge the same way it was measured for the S_phase ladder and the hybrid
comparison — from the **visual structure of the ambient-ion phase space**, at matched
`t*wci0`, not from a scalar diagnostic. Every rung has 32 plotfiles over the same
`0.302/w_ci0`, so frames match by index and no `argmin` tolerance question arises.

```bash
python scripts/plot_eps_ladder.py     # media/E_phase/wedge_vs_eps.png + the wedge table
```

## Run it

```bash
perlmutter/submit.sh eps --dry     # print the sbatch line, submit nothing
perlmutter/submit.sh eps           # four rungs, one GPU each, cheapest-first
```
