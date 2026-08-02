# Run `R1_paper`

> Faithful reproduction of Schaeffer 2020 Table I: theta_e,ab = 0.092, ambient 0.008 n_e,ab, target 2.0 n_e,ab, B0 from Table I's own code units, 30000 cells over 9000 d_e,ab one-sided, 220 t_ab, and PSC's collisionality lambda_ab = 20 translated into WarpX's lnLambda. Every Table I dimensionless row reproduces to <=2% (beta convention aside); only lambda_mfp/d_i0 cannot match, and it cannot in PSC either at mu_p = 100.

**tier** Full · **deck** `inputs_kinshock_R1_paper` · **paper** Schaeffer et al., Phys. Plasmas 27, 042901 (2020), Table I + Sec. II

## Notes

<!-- prose:begin -->
**The faithful Table I run, collisions included.** Built 2026-08-01 after re-reading
`schaeffer2020.pdf` and PSC's own source. Every Table I dimensionless row reproduces
to <=2%: `beta_ab` 2300 vs 2300, `beta_0` 0.400 vs 0.400, `d_i0/d_i,ab` 11.180 vs
11.18, gyroperiod 33.91 vs 33.9 t_ab, `lambda_ab` 20.0 vs 20.

**Collisions are PSC's dial, translated.** `psc_2d_shock.cxx:427` sets
`collision_nu = 3.76*target_Te_heat^2/Zi/lambda0` with `lambda0` a plain free input --
PSC's collision path (`psc_collision_impl.hxx:230` -> Takizuka-Abe in
`binary_collision.hxx:229`) contains no Coulomb logarithm, no SI density and no SI
temperature. WarpX's Perez 2012 operator is the same small-angle physics but fully
dimensional, so the only way in is lnLambda. Inverting gives **lnLambda = 3.0e9**,
which is the correct translation, not a hack: PSC's `lambda0 = 20` is equally
unphysical at the sim's own (n, T). Verified safe -- WarpX's cross-section clamp
(`sigma_eff = min(pi b0^2 lnL, sigma_max)`) sits at sigma/sigma_max = 5.5e-6, and
`nu_ei*dt*ndt = 0.034 << 1`.

**The one row that cannot match, in any code:** Table I's `lambda_mfp/d_i0 = 350`
against this run's upstream `mfp_ii/d_i0 = 0.015`. The ratio of that to `lambda_ab` is
independent of lnLambda (both mfps ~ 1/lnLambda) AND of c (both normalizing lengths
~ c), so the 2.3e4x gap exists inside PSC too at mu_p = 100. Sec. II says why: global-
scale collisionality is "only quantitatively matched at physical mass ratios". **That
row describes the experiment, not the run -- criterion 2's 350 threshold does not
apply here.**

**What changed vs R1_warm** (six primaries, all corrections rather than tuning):
`theta_e_heat` 0.078 -> 0.092, target 2.5 -> 2.0 n_e,ab and ambient 0.01 -> 0.008
n_e,ab (Table I's numbers are CODE units, n_e,ab = 1.25), `theta_e_init` 0.001 ->
0.002 and `theta_i_init` 7.8e-4 -> 2.0e-5 (PSC's shock case injects at
target_Te = target_Ti = 0.002), domain 7500 -> 9000 d_e,ab, plus collisions on.

**The clock is right by construction.** B0 is the primary, so `wci0 = q_e*B0/m_i` no
longer inherits the ambient density's error; the 1.218x short gyroperiod that biased
every `t*wci0` observable late in R1_warm is gone. Headline test: **does onset t*_1
move from 1.35-1.41 to the paper's ~1?**

**Two things to watch, both genuine open questions:**
- `theta_i_init = 2.0e-5` reverts the warm piston-ion change. That change was made
  because cold ions gave a delta-function phase-space ridge instead of the paper's
  broad Fig. 7 row-2 ridge, but it is NOT a PSC input. If the ridge is still narrow,
  the broadening is something our heater/injector does not reproduce -- worth knowing
  separately from the unit errors.
- `psc_2d_shock.cxx` is the QUASI-perpendicular variant (BB_perp/BB_par split,
  lambda0 = 30, target_Te_heat = 0.06) that the paper defers to Ref. 36. The exact
  perpendicular case file is not in the local checkout, so its injection temperature
  and collision_interval are the closest structural guide, not the run itself.

**Residual deviation:** the beta rows read 2x Table I because `units.py` uses
beta = 2*mu0*n*T/B^2 while the paper tabulates mu0*n*T/B^2 (the comparison table
doubles the paper value so like meets like).

Cost: 322400 steps x 30000 cells, plus collisions on 10 species pairs every 10 steps
~ **12-14 h**.
<!-- prose:end -->

## Primaries — `config.yaml` is the only source of truth

Edit these; never the deck. Regenerate with `python scripts/make_inputs.py runs/R1_paper`.

| Quantity | Value | config key |
|---|---|---|
| reference density n0 [m^-3] | 6.0e26 | `reference.n0` |
| mass ratio m_i/m_e | 100 | `reference.mass_ratio` |
| charge state Z | 1 | `reference.charge_state` |
| target density / n0 | 2 | `plasma.piston.density_over_n0` |
| heater theta_e,ab | 0.092 | `plasma.piston.theta_e_heat` |
| piston init theta_e | 0.002 | `plasma.piston.theta_e_init` |
| piston init theta_i | 2e-05 | `plasma.piston.theta_i_init` |
| ambient density / n0 | 0.008 | `plasma.ambient.density_over_n0` |
| ambient theta_0 | 0.002 | `plasma.ambient.theta_0` |
| field orientation | perpendicular | `field.orientation` |
| **B0 [T]** | 70.2734 | `field.B0_tesla` |
| dims | 1 | `geometry.dims` |
| layout | one_sided | `geometry.layout` |
| slab halfwidth [d_i] | 2 | `geometry.slab_halfwidth_di` |
| domain halfwidth [d_e] | 9000 | `geometry.domain_halfwidth_de` |
| dz [d_e] | 0.3 | `geometry.dz_over_de` |
| boundary lo / hi | lo symmetry, hi open | `geometry.boundary` |
| CFL | 0.75 | `numerics.cfl` |
| particle shape | 2 | `numerics.particle_shape` |
| max_step | 322400 | `numerics.max_step` |
| ppc | piston 100, ambient 100 | `numerics.ppc` |
| heater intervals | 20 | `operators.heater.intervals` |
| injector intervals | 20 | `operators.injector.intervals` |
| injector tau [1/wpe] | 40 | `operators.injector.tau_over_wpe_inv` |
| collisions target | quantity lambda_ab, value 20 | `collisions.target` |

## Derived — computed by `units.derive`, not stored

| Quantity | Value | From |
|---|---|---|
| omega_pe [rad/s] | 1.38187e+15 | sqrt(n0 q^2/(eps0 m_e)) |
| d_e [m] | 2.16947e-07 | c/omega_pe |
| d_i,ab [m] | 2.16947e-06 | d_e*sqrt(m_i/m_e) |
| C_s,ab / c | 0.0303315 | sqrt(theta_e,ab/mass_ratio) |
| t_ab [s] | 2.38583e-13 | d_i,ab / C_s,ab |
| n_amb [m^-3] | 4.8e+24 | density_over_n0 * n0 |
| d_i0 [m] | 2.42554e-05 | c/omega_pi(n_amb) |
| d_i0 / d_e | 111.803 | derived |
| **omega_ci0 [rad/s]** | 1.23598e+11 | **q_e*B0/m_i — B0 only, no n_amb** |
| 1/omega_ci0 [s] | 8.09073e-12 | 1/omega_ci0 |
| **v_A / c** | 0.01 | **B0/sqrt(mu0*n_amb*m_i) — DERIVED** |
| rho_i0 / d_e | 1162.76 | v_p/omega_ci0 / d_e |
| v_p / c (model) | 0.104 | config model.vp_over_c |
| v_sh / c (model) | 0.139525 | model.vsh_over_Csab * C_s,ab |
| M_A | 13.9525 | v_sh/v_A |
| M_ms | 12.7368 | v_sh/sqrt(v_A^2+C_s0^2) |
| beta_ab | 2300 | 2*mu0*n0*T_e,ab/B0^2 |
| beta_0 | 0.4 | 2*mu0*n_amb*T_0/B0^2 |
| dz [m] | 6.50841e-08 | dz_over_de * d_e |
| dt [s] | 1.62823e-16 | CFL-limited |
| dt*omega_pe | 0.225 | dt * omega_pe |
| n_cell | 30000 | domain / dz (halved when one_sided) |
| steps per 1/omega_ci0 | 49690.4 | 1/(omega_ci0*dt) |
| run length [1/omega_ci0] | 6.48817 | max_step * dt * omega_ci0 |
| T_e,ab [eV] | 47011.9 | theta_e,ab * m_e c^2 |
| lnLambda (used) | 122348 | units.coulomb_log_for(collisions.target) |
| lnLambda (physical) | 10.8364 | NRL at (n0, T_e,ab) |
| nu_ei,ab [1/s] | 2.09571e+13 | NRL electron-ion |
| nu_ei*dt | 0.00341229 | must be << 1 |
| mfp_ei,ab / d_e | 20 | v_te/nu_ei / d_e |
| lambda_ab | 20 | omega_ce,ab/nu_ei,ab = mfp/d_e,ab |
| mfp_ii,amb / d_i0 | 0.0149509 | upstream ion-ion — what criterion 2 tests |

## vs Schaeffer 2020 Table I

| Quantity | This run | Table I | ratio | known cause if off |
|---|---|---|---|---|
| mass ratio m_i/m_e | 100 | 100 | 1.000x (+0.0%) ok |  |
| C_s,ab / c | 0.0303315 | 0.03 | 1.011x (+1.1%) ok |  |
| piston speed v_p / c | 0.104 | 0.104 | 1.000x (+0.0%) ok |  |
| Alfven Mach v_sh/v_A | 13.9525 | 14 | 0.997x (-0.3%) ok |  |
| magnetosonic Mach | 12.7368 | 13 | 0.980x (-2.0%) ok |  |
| ablation beta (Table I 1150, x2 convention) | 2300 | 2300 | 1.000x (-0.0%) ok |  |
| upstream beta (Table I 0.2, x2 convention) | 0.4 | 0.4 | 1.000x (-0.0%) ok |  |
| d_i0 / d_i,ab | 11.1803 | 11.18 | 1.000x (+0.0%) ok |  |
| gyroperiod in ablation times | 33.9116 | 33.9 | 1.000x (+0.0%) ok |  |
| collisionality mfp/d_e,ab | 20 | 20 | 1.000x (-0.0%) ok |  |

## Files

| Path | What |
|---|---|
| `config.yaml` | the primaries (tracked) |
| `inputs_kinshock_R1_paper` | generated deck (tracked, **never hand-edit**) |
| `warpx_used_inputs` | what WarpX actually ran (tracked) |
| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |
| `diags/`, `*.log` | output (gitignored, regenerable) |

---

_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, not the tables. Prose between the `prose:` markers is preserved._
