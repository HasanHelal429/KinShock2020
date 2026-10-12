# U_phase — the physical mass ratio, as a single-variable A/B

**The question.** Every run in this project uses `mu = m_i/m_e = 100`, the paper's value.
Umeda, Kidani, Matsukiyo & Yamazaki (2014, arXiv:1401.5903), on *perpendicular* shocks,
show the foot-instability selector is fixed by the mass ratio once `M_A` and `beta` are
pinned — the identity is `v_te/v_A = sqrt(beta_e * mu)`, so

    v_u/v_te = M_A / sqrt(beta_e * mu)

| | M_A | beta_0 | mu | `v_u/v_te` | regime |
|---|---|---|---|---|---|
| Table I **as we simulate it** | 13.95 | 0.196 | **100** | **3.15** | Buneman (`v_u > v_te`) |
| the **same plasma** at `mu_phys` | 13.95 | 0.196 | **1836** | **0.74** | MTSI (`v_u < v_te`) |

**`mu = 100` puts every run in this project on the wrong side of `v_u = v_te`** — including
`ek_47keV_t2` and `mb_r100_t2`, which look healthy (both are 3.15, because `M_A` and
`sqrt(beta_0)` scale together between them). The Buneman instability needs the drift to
exceed the electron thermal speed; at the real mass ratio it does not, and the MTSI takes
over instead. Different instability, different foot physics.

**The paper's own mass-ratio check could not have caught this.** Fig. 14 compares
`mu = 100` against `400` and finds the ion physics converged — but `mu = 400` gives
`v_u/v_te = 1.58`, still above threshold. Both of its points sit in the Buneman regime.
The crossing happens between 400 and 1836.

## The two runs

Each is a twin of an existing run with **`reference.mass_ratio` as the only choice made**.
`B0 = sqrt(mu0 n T_0/beta_0)` contains no `mu`, so holding `M_A`, `beta_0`, `beta_ab`,
`T_0`, `n_0` and `T_e,ab` leaves B0 and every `targets:` value untouched. Three things
follow from `mu` and nothing else does:

    v_sh, v_p, C_s,ab ~ 1/sqrt(mu)   ->  vp_over_c / 4.285
    d_i0 = sqrt(mu) d_e0             ->  halfwidth_de x 4.285  (box FIXED in d_i0)
    w_ci0 ~ 1/mu                     ->  max_step x 18.36      (t*w_ci0 FIXED)

`dz`, `lambda_D,0`, `d_e,ab` and `ppc` are all `mu`-independent, so `dz/lambda_D` and `N_D`
carry over untouched. `rho_i0/d_i0` and the 24.73 d_i0-per-`t*` precursor rate are
`mu`-invariant too, so the box needs no resizing in `d_i0`. Cost goes as `mu^1.5` = 78.7x.

| run | partner | `dz/lambda_D` | cells | steps | `t*` | wall |
|---|---|---|---|---|---|---|
| `um1836_470eV` | `S_phase/ss_dz1_ppc100_int1` | 6.07 | 19,200 | 5,508,450 | 0.604 | **~4.5 GPU-h** |
| `um1836_47keV` | `E_phase/ek_47keV_t2` | 0.758 | 70,272 | 1,460,015 | 2.000 | **~4.4 GPU-h** |

Verified against the partners: B0, T_0, M_A, M_ms, beta_0, beta_ab, `dz/lambda_D`,
`rho_i0/d_i0`, box in `d_i0` and `t*w_ci0` are all held; only `v_sh/c` and `v_u/v_te` move,
both by exactly `1/sqrt(18.36)` = 0.233.

## Reading it

`ss_dz1_ppc100_int1` shows the notch at **plateau 0.970, drop −52** (RESULTS 2026-09-03).

* **Notch survives in `um1836_470eV`** ⇒ the mass ratio is exonerated and the cause is in
  the `eps` family (`gamma`, `c/v_ti`, `lambda_D,0/d_i0`).
* **Notch gone** ⇒ `mu = 100` was holding the whole campaign in the Buneman regime when the
  real plasma is in the MTSI regime, and every 470 eV result needs re-reading.

`um1836_47keV` is the control: its partner has no notch, so it should not grow one. If it
*does*, the mass ratio is doing something at high `eps` too and the interpretation changes
again.

**The pair is what makes a `dz/lambda_D` = 6.07 run usable.** That is deep in the
grid-heating regime, but grid heating, particle noise, box-filling and the collision
operator are identical in twin and partner, so they are common mode and cannot be
confounders.

⚠ **If the notch vanishes, check the alternative before concluding.** It could be grid
heating scaling differently at higher `mu` rather than the regime change. Refining the pair
costs 16x per 4x in `dz` — `dz/lambda_D` = 1.52 would be ~72 GPU-h for the 470 eV twin — so
the coarse A/B is the affordable version and that check has to be made by argument.

⚠ The 12.02 d_i0 box fills at `t*` ~ 0.49, same as its partner. Fine for the notch, which
is mature by 0.40; the far upstream is not clean past that in **either** run.

## Run it

Both are predicted at ~4.5 GPU-h against `SWEEP_TIME`'s 4 h default, so **pass `--time` or
they die at the wall.** `shared` bills elapsed, not requested, so the extra hours are free.

```bash
perlmutter/submit.sh mu --dry --time 10:00:00
perlmutter/submit.sh mu --time 10:00:00

# 470 eV alone -- it is the one that carries the result
perlmutter/submit.sh mu runs/U_phase/um1836_470eV --time 10:00:00
```
