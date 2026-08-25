# M_phase — the B₀ scan at fixed `eps`

**The question.** E_phase established that the ambient wedge is a function of
`eps = v_te,ab/c` alone — 2.062 d_i0 at `eps` = 0.0303 down to 1.164 at 0.3033, with every
differing row a pure power of `eps` to a worst residual of 3.0e-4. It could not say *which*
consequence of `eps` does the work, and **it never can**, because two exact identities
collapse the candidates into one number:

    w_pe0/w_ce = 1/(sqrt(mu) * v_A/c)          lambda_D,0/d_i0 = sqrt(beta_0) * v_A/c

On that ladder `v_A/c`, `w_pe0/w_ce`, `lambda_D,0/d_i0`, `sigma_i` and `v_sh/c` are the same
parameter renamed. "It's the B field", "it's rho_e/lambda_D", "it's the wave regime" and
"it's relativistic" are one hypothesis with four names. That is why three mechanism tests
(RESULTS 2026-08-19) and the (k,omega) measurement on `ek_47keV_t2` (2026-08-23) all failed
to attribute: there was nothing to attribute *between*.

## The cut

`beta_0` is **held fixed** in E_phase, and that is precisely what locks `lambda_D,0/d_i0` to
`v_A/c`. Let `beta_0` go and they separate. Scanning **B₀ alone**, at fixed `eps`, moves the
magnetization family and freezes everything else — verified numerically with
`scripts/table1_ladder.py` over these five runs:

| | frozen | moving |
|---|---|---|
| **physics** | `eps` 0.3033 · `gamma` 1.138 · T_e,ab 47012 eV · T_0 1000 eV · T_0/T_e,ab · n_0 · d_e0 · d_i0 · `lambda_D,0` · **`lambda_D,0/d_i0` 4.42e-3** · `v_sh/c` · `v_p/c` · `C_s,ab/c` · `v_sh/C_s,ab` · `lambda_ab` 20 · `mfp_dir(v_sh)/d_i0` 2666 | `w_pe0/w_ce` · `v_A/c` · `M_A` · `beta_0` · `beta_ab` · `sigma_i` · `rho_e0/d_i0` · `rho_i0/d_i0` |
| **numerics** | `dz` · `dt` · `dz/lambda_D,0` 0.3791 · `n_cell` 7168 · `ppc` 100 · `N_D,PIC` 264 · `t*w_ci0` 0.3018 | `max_step` (∝ 1/B₀, to hold the clock) |

**`lambda_D,0/d_i0` stays at 4.424e-3 while `w_pe0/w_ce` moves a full decade.** That is the
separation E_phase cannot produce, and it is the whole point of the phase.

## The rungs

Rung 1 **is `runs/E_phase/es_47keV`** — already run, wedge measured at 1.164 d_i0. It is not
copied here: a second `config.yaml` for an identical configuration is exactly what CLAUDE.md's
A/B note forbids. The four new runs share its frame cadence (30 plotfiles / 150 field / 300
reduced over the same `t*w_ci0` = 0.3018), so frames match **by index** and the wedge is read
the same way at `t*w_ci0` = 0.26 with no `argmin` tolerance question.

| run | B₀ [T] | `w_pe0/w_ce` | `v_A/c` | `M_A` | `M_ms` | `beta_0` | `rho_i0/d_i0` | steps | wall (1 GPU) |
|---|---|---|---|---|---|---|---|---|---|
| `E_phase/es_47keV` (rung 1) | 70.273 | 10.0 | 0.01000 | 13.95 | 12.76 | 0.196 | 10.40 | 23 997 | **1 m (done)** |
| `mb_r20` | 35.137 | 20.0 | 0.00500 | 27.9 | 20.9 | 0.783 | 20.80 | 47 994 | ~2 m |
| `mb_r32` | 22.222 | 31.6 | 0.00316 | 44.1 | 25.7 | 1.957 | 32.89 | 75 885 | ~3 m |
| `mb_r56` | 12.497 | 56.2 | 0.00178 | 78.5 | 29.3 | 6.190 | 58.48 | 134 945 | ~6 m |
| `mb_r100` | 7.027 | 100.0 | 0.00100 | 139.5 | 30.8 | 19.57 | 104.0 | 239 970 | ~10 m |

`w_pe0/w_ce` = 100 at rung 5 is the **470 eV end of E_phase**, reached here with `eps` frozen
at the 47 keV value. Predicted walls are `steps` scaled off `es_47keV`'s measured 1 minute;
expect 1.4–2.0× more from the documented small-point GPU under-occupancy penalty, so **~35–45
GPU-min for the whole phase.**

## Reading it

Two-way, and either answer is a result:

* **Wedge climbs back toward 2.06 as `w_pe0/w_ce` → 100** ⇒ the agent is in the magnetization
  family, and `gamma` and `lambda_D,0/d_i0` are exonerated in one shot (both are frozen here).
* **Wedge stays near 1.164** ⇒ the magnetization family is out, and the agent is on the `eps`
  side — `gamma`, `lambda_D,0/d_i0`, or the shocked-electron heating channel. The follow-up is
  then the fixed-upstream control (pin T_0, n_0, B₀; move T_e,ab), ~7.4 GPU-h.

Measure the wedge exactly as E_phase did — from the visual structure of the ambient-ion phase
space at matched `t*w_ci0`, not a scalar diagnostic:

```bash
python scripts/plot_eps_ladder.py     # adapt: the x axis here is w_pe0/w_ce, not eps
```

## Two confounds, stated up front

1. **`M_A` and `w_pe0/w_ce` cannot be separated by this phase.** At fixed `eps`, `v_sh/c` is
   fixed, so `M_A = sqrt(mu) * (v_sh/c) * (w_pe0/w_ce)` **exactly** — they are proportional by
   construction. This phase separates the magnetization family from the `eps` family, *not its
   members*. If the wedge does track it, the next question is whether `M_A` alone accounts for
   it, which is already in tension with the paper's own Fig. 6 Mach scan (t\*/z\* invariant in
   `w_ci0^-1`/`rho_i0` units) and is what REPLICATION_PLAN's R4/R5 were staged for.
2. **`beta_0` swings 0.196 → 19.6, so the far rungs are a different upstream.** By `mb_r100`
   the ambient is thermal-pressure dominated, `C_s0` has overtaken `v_A` (which is why `M_ms`
   saturates near 31 while `M_A` reaches 139), and there may be no magnetized perpendicular
   shock at all. **Read rungs 2–3 first**: together they cost ~5 minutes, reach
   `w_pe0/w_ce` = 31.6 at `beta_0` ≤ 2 and `M_A` ≤ 44, and if the wedge moves there the
   question is answered without the degenerate end.

`targets:` in these configs carries **per-rung check values, not Table I targets** — this phase
deliberately leaves `M_A` = 14 and `beta_0` = 0.2 behind, and `run_checks.py` should be read
with that in mind.

## Run it

```bash
perlmutter/submit.sh mag --dry     # print the sbatch line, submit nothing
perlmutter/submit.sh mag           # four rungs, one GPU each, cheapest-first
perlmutter/submit.sh mag runs/M_phase/mb_r20 runs/M_phase/mb_r32   # just the safe pair
```
