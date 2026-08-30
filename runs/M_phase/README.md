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

## RESULT (2026-08-24): rungs 2–3 ran clean and the measurement is INVALID — the box

`57563903_0/_1`, both COMPLETED, exit 0, full step counts, no NaN or abort, and the
measurement method reproduces `es_47keV`'s recorded wedge of 1.164 d_i0 exactly. The runs
are fine. **The phase design was wrong.**

| run | `w_pe0/w_ce` | z_piston [d_i0] | wedge [d_i0] | `T_e,upstream` | verdict |
|---|---|---|---|---|---|
| `es_47keV` | 10.0 | 3.358 | 1.164 | **2.03** | valid |
| `mb_r20` | 20.0 | 7.301 | 3.122 | **38.76** | INVALID |
| `mb_r32` | 31.6 | 10.547 | `nan` | **33.59** | INVALID |

Two hard failures. The far-upstream control — which RESULTS 2026-08-17 requires to stay
cold or the comparison is void — is 19× hot. And the piston front reaches 7.30 and 10.55
d_i0 inside a 12.02 d_i0 box, so at `mb_r32` no undisturbed ambient layer remains and the
wedge scan finds no crossing at all. **Do not quote 3.122**: with a contaminated upstream
the `SPREAD` threshold is crossed far downstream of the real wedge tip, so it is an
artifact, not an upper bound and not a signal.

**Cause, and it is a design error in this phase, not in the runs.** `w_ci0 ~ B0`, so holding
`t*w_ci0` fixed stretches ABSOLUTE time as `1/B0`, while `v_p` and `d_i0` are frozen. The
piston therefore sweeps `1/B0` further through a box left at a fixed 12.02 d_i0. Predicted
piston positions 6.72 and 10.62 against measured 7.30 and 10.55. The table above lists
`rho_i0/d_i0` as *moving* — the consequence that the DOMAIN must move with it was missed.
`mb_r56` and `mb_r100` would put the piston at 18.9 and 33.6 d_i0 in the same 12.02 box;
they were not submitted, and must not be until the box is fixed.

**FIX APPLIED 2026-08-24 — boxes scale with rho_i0, and rho_i0 is now the unit.**

The domain is held at rung 1's **1.156 rho_i0** at every rung instead of a fixed 12.02
d_i0, i.e. `domain_halfwidth_de ~ 1/B0`. `dz` is UNCHANGED at 0.1875 d_e,ab, so
`dz/lambda_D,0` = 0.3791 and `N_D` = 264 still hold to four figures and only the cell
COUNT grows. Verified across the five runs:

| run | `domain_halfwidth_de` | cells | box [d_i0] | **box [rho_i0]** | `dz/lambda_D,0` | `N_D` | `t*w_ci0` | steps | cost | wall |
|---|---|---|---|---|---|---|---|---|---|---|
| `es_47keV` | 1344.0 | 7 168 | 12.02 | **1.1559** | 0.3791 | 263.8 | 0.3018 | 23 997 | 1× | done |
| `mb_r20` | 2688.0 | 14 336 | 24.04 | **1.1559** | 0.3791 | 263.8 | 0.3018 | 47 994 | 4× | ~6 m |
| `mb_r32` | 4249.5 | 22 664 | 38.01 | **1.1557** | 0.3791 | 263.8 | 0.3018 | 75 885 | 10× | ~15 m |
| `mb_r56` | 7558.5 | 40 312 | 67.61 | **1.1560** | 0.3791 | 263.8 | 0.3018 | 134 945 | 32× | ~47 m |
| `mb_r100` | 13440.0 | 71 680 | 120.21 | **1.1559** | 0.3791 | 263.8 | 0.3018 | 239 970 | 100× | ~2.5 h |

Phase total **~3.5 GPU-h**, against the 35 min originally quoted. Cell counts are snapped
to `blocking_factor` = 8, which moves the box by at most 0.014%.

**The clock/piston dilemma dissolves.** With the box scaling, holding `t*w_ci0` ALSO holds
the piston at a fixed fraction of it — 0.3229 rho_i0 of 1.156, 28%, the anchor's value at
every rung. So `--match clock` and `--match piston` now coincide, and the salvage's second
defect (higher rungs having had proportionally fewer gyroperiods to form a shock) is gone.
Both matched at once is exactly what the fixed box buys.

**rho_i0 is the settled normalization.** What does *not* dissolve is the unit choice, since
`rho_i0/d_i0 = v_p/v_A` still moves a decade — and on the small-box salvage the two
disagreed **in sign** (+23% in d_i0, −38% in rho_i0). It is settled in favour of `rho_i0`:
it is the piston's own sweep scale, it keeps the geometry comparable between rungs, and it
is what the box is now held to. `scripts/plot_mag_ladder.py` reports rho_i0 as the primary
column and keeps d_i0 as a cross-check, printing the disagreement when it occurs rather
than letting a reader pick the flattering one.

**The number to beat: rung 1's wedge is 0.1123 rho_i0** (= 1.168 d_i0, reproducing the
recorded 1.164 to 0.3% on the finer rho_i0 scan step), with the piston at 0.3229 rho_i0 and
`T_e,upstream` = 2.0.

The superseded small-box output is preserved as `diags_smallbox/` in `mb_r20` and `mb_r32`,
with their logs as `*.smallbox`, so the invalid measurement and the valid `t*` = 0.131
salvage point both remain inspectable. `diags/` is clear, so `launch.sh` will start.

## RESULT (2026-08-24): the wedge follows the magnetization family — and one group survives both phases

`57724352_0/_1`, both COMPLETED, full step counts, no NaN, 31+151 frames each, and **all
six panels VALID** — `T_e,upstream` runs 1.2–3.0 against the 5.0 gate, where the small-box
attempt hit 39. The piston lands at 0.323–0.359 rho_i0, essentially the predicted 0.3229,
confirming the box scaling.

| run | `w_pe0/w_ce` | `M_A` | `beta_0` | z_piston [rho_i0] | **wedge [rho_i0]** | wedge [d_i0] | `T_e,up` |
|---|---|---|---|---|---|---|---|
| `es_47keV` | 10.0 | 14 | 0.20 | 0.3229 | **0.1123** | 1.168 | 2.0 |
| `mb_r20` | 20.0 | 28 | 0.78 | 0.3473 | **0.1465** | 3.046 | 2.6 |
| `mb_r32` | 31.6 | 44 | 1.96 | 0.3589 | **0.1544** | 5.077 | 3.0 |

**The wedge grows with `w_pe0/w_ce`: +30% then +37%, at FIXED eps, gamma and
`lambda_D,0/d_i0`.** The same trend appears at `t*` = 0.131 (0.0509 → 0.0586 → 0.0690,
+36%), so it is not a single-frame artifact. This **overturns the small-box salvage**,
which reported −38% in rho_i0; that comparison matched piston position at the cost of
halving the gyroperiods, and the fixed box removes the need for it.

### M_phase reproduces E_phase's curve while moving different things

| `w_pe0/w_ce` | E_phase [rho_i0] | M_phase [rho_i0] | M/E |
|---|---|---|---|
| 10.0 | 0.1119 | 0.1123 | 1.003 |
| 20.0 | 0.1344 (interp) | 0.1465 | 1.090 |
| 31.6 | 0.1738 | 0.1544 | 0.888 |

E_phase varied `eps` with `M_A` and `beta_0` **pinned**; M_phase varies `B0` with `eps` and
`gamma` **pinned** and lets `M_A` run 14 → 44. They agree to 2–11%.

### The two-way elimination

Intersecting the groups each phase moves (`scripts/table1_ladder.py --json`, numerics
excluded) leaves **one dimensionless group**:

| group | E_phase | M_phase | verdict |
|---|---|---|---|
| `eps`, `gamma`, `lambda_D,0/d_i0`, `T_e,ab`, `v_p/c`, `C_s,ab/c`, `t_ab` | moves | FIXED | **excluded by M_phase** |
| `M_A`, `M_ms`, `beta_0`, `beta_ab`, `rho_i0/d_i0`, `rho_e,ab/d_e,ab`, `nu_ei/w_ce` | FIXED | moves | **excluded by E_phase** |
| **`v_A/c`** | moves | moves | **SURVIVES BOTH** |

The six surviving dimensionless rows — `v_A/c`, `w_pe0/w_ce`, `w_pi0/w_ci`,
`w_pe,ab/w_ce`, `sigma_i`, `sigma_e` — are **one group under six names**, verified as exact
powers of `v_A/c` with prefactors constant to 1e-4 in both phases:

    w_pe0/w_ce = 1/(sqrt(mu) v_A/c)    w_pi0/w_ci = 1/(v_A/c)    sigma_i = (v_A/c)^2

**So the wedge is set by `v_A/c`** — equivalently the electron magnetization relative to the
Debye scale, `rho_e/lambda_D = w_pe0/w_ce`. The B-field hypothesis is confirmed, and each
competing family is excluded by the phase that holds it fixed. Note this is also the group
the reduced-c scheme distorts by a decade (RESULTS 2026-08-24): the real experiment sits at
`w_pe0/w_ce` = 100 and PSC's Table I run at 10.

### What is not settled

* **Rungs 4–5 are unrun.** `mb_r56` and `mb_r100` reach `beta_0` = 6.2 and 19.6, where the
  upstream is thermal-pressure dominated and may not support a perpendicular shock at all.
  The trend above spans only 10 → 31.6 of the intended 10 → 100.
* **The 11% shortfall at `w_pe0/w_ce` = 31.6** (0.1544 against E_phase's 0.1738) is larger
  than the 0.3% agreement at the anchor. It could be scatter, or a residual second
  dependence the two phases do not share. Three points cannot tell.
* **Within the magnetization family nothing is separated.** `v_A/c`, `w_pe0/w_ce`,
  `sigma_i` are the same number; this phase cannot say which *name* is the mechanism, only
  that the group is the agent.

## Run it

```bash
perlmutter/submit.sh mag --dry     # print the sbatch line, submit nothing
perlmutter/submit.sh mag           # four rungs, one GPU each, cheapest-first
perlmutter/submit.sh mag runs/M_phase/mb_r20 runs/M_phase/mb_r32   # just the safe pair
```
