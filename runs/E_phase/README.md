# E_phase — the `eps = v_te,ab/c` ladder

**The question.** `R1_paper` (47 keV) and `R1_paper_470eV` were built to be the *same*
dimensionless problem — Table I of Schaeffer 2020 — yet they produce visibly different
shock structure: an ambient wedge 2.25 d_i0 deep at 470 eV against 1.25 d_i0 at 47 keV.
S_phase spent a full resolution ladder plus six other numerical knobs failing to remove
that difference. This phase asks the other question: **are the two runs actually the same
physics?**

## The audit: exactly one parameter differs

Computed from the config primaries only (`$CLAUDE_JOB_DIR/tmp/dimensionless.py`, folded
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
| `es_1p5keV` | 1486 eV | 0.0539 | 12.494 | 56.2 | 1.0044 | 40312 | 759079 | 1.31 h |
| `es_4p7keV` | 4700 eV | 0.0959 | 22.220 | 31.6 | 1.0138 | 22667 | 239998 | 0.23 h |
| `es_15keV` | 14860 eV | 0.1705 | 39.509 | 17.8 | 1.0436 | 12748 | 75909 | 0.04 h |
| `es_47keV` | 47012 eV | 0.3033 | 70.273 | 10.0 | 1.1380 | 7167 | 23994 | 0.01 h |

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

## What the shape of the answer means

The two candidate mechanisms predict **different shapes**, which is why the ladder needs
interior points and not just the two endpoints already in hand:

- **(a) Relativistic capping of electron heating** (`~eps^2`). `gamma(T_e,ab)` runs
  1.0014, 1.0044, 1.0138, 1.0436, 1.1380 — flat until the top rung. R1_paper's *shocked*
  electrons were measured at 0.532 c, `gamma = 1.18`. Predicts a wedge that stays ~2.25
  d_i0 for the lower three rungs and collapses only at 47 keV.
- **(b) Electron-scale wave regime** (`~eps^-1`). `rho_e/lambda_D = w_pe/w_ce` runs
  100, 56, 32, 18, 10 — a smooth decade, and it sets how many Bernstein resonances sit
  below `w_pe`, i.e. how ECDI-like versus Buneman-like the ramp's electron heating is.
  Predicts a wedge that thins monotonically and roughly log-linearly at every rung.
- **(c) No trend** — the wedge is not an `eps` effect and the difference is elsewhere.

The hybrid run constrains this already: `H3_470eV_dense` carries the **470 eV** value of
`v_A/c` (9.999e-4) and still reproduces R1_paper's thin structure. So no group built purely
from `v_A/c` can be the cause — whatever is responsible must live in the kinetic-electron
sector, which is what the hybrid's fluid closure removes and what both (a) and (b) are.

## Reading the result

Measure the wedge the same way it was measured for the S_phase ladder and the hybrid
comparison — from the **visual structure of the ambient-ion phase space**, at matched
`t*wci0`, not from a scalar diagnostic. Every rung has 30 plotfiles over the same
`0.302/w_ci0`, so frames match by index and no `argmin` tolerance question arises.

## Run it

```bash
perlmutter/submit.sh eps --dry     # print the sbatch line, submit nothing
perlmutter/submit.sh eps           # four rungs, one GPU each, cheapest-first
```
