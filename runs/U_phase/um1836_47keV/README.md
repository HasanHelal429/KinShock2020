# Run `um1836_47keV`

> ek_47keV_t2 at the PHYSICAL mass ratio mu = 1836.15, single variable. Crosses the foot-instability threshold v_u/v_te = M_A/sqrt(beta_e mu) from 3.15 (Buneman) to 0.74 (MTSI), per Umeda+2014. B0, T_0, T_e,ab, n_0, dz/lambda_D, ppc and every targets: value are unchanged; only the box in d_e, max_step and v_p follow mu. See runs/U_phase/README.md.

**tier** Sweep · **deck** `inputs_kinshock_um1836_47keV` · **paper** Schaeffer et al., Phys. Plasmas 27, 042901 (2020), Table I + Sec. II

## Notes

<!-- prose:begin -->
_No hand-written notes yet. Anything written between the
`prose:begin` / `prose:end` markers survives regeneration — put the run's story here
(why it exists, what it showed, what to distrust)._
<!-- prose:end -->

## Primaries — `config.yaml` is the only source of truth

Edit these; never the deck. Regenerate with `python scripts/make_inputs.py runs/um1836_47keV`.

| Quantity | Value | config key |
|---|---|---|
| reference density n0 [m^-3] | 6.0e26 | `reference.n0` |
| mass ratio m_i/m_e | 1836.15 | `reference.mass_ratio` |
| charge state Z | 1 | `reference.charge_state` |
| target density / n0 | 2 | `plasma.piston.density_over_n0` |
| heater theta_e,ab | 0.0920002 | `plasma.piston.theta_e_heat` |
| piston init theta_e | 0.00195745 | `plasma.piston.theta_e_init` |
| piston init theta_i | 1.95745e-05 | `plasma.piston.theta_i_init` |
| ambient density / n0 | 0.008 | `plasma.ambient.density_over_n0` |
| ambient theta_0 | 0.00195745 | `plasma.ambient.theta_0` |
| field orientation | perpendicular | `field.orientation` |
| **B0 [T]** | 70.2734 | `field.B0_tesla` |
| dims | 1 | `geometry.dims` |
| layout | one_sided | `geometry.layout` |
| slab halfwidth [d_i] | 2 | `geometry.slab_halfwidth_di` |
| domain halfwidth [d_e] | 26348.5 | `geometry.domain_halfwidth_de` |
| dz [d_e] | 0.37495 | `geometry.dz_over_de` |
| boundary lo / hi | lo symmetry, hi open | `geometry.boundary` |
| CFL | 0.75 | `numerics.cfl` |
| particle shape | 3 | `numerics.particle_shape` |
| max_step | 1460015 | `numerics.max_step` |
| ppc | piston 100, ambient 100 | `numerics.ppc` |
| heater intervals | 20 | `operators.heater.intervals` |
| injector intervals | 20 | `operators.injector.intervals` |
| injector tau [1/wpe] | 400 | `operators.injector.tau_over_wpe_inv` |
| collisions target | quantity lambda_ab, value 20 | `collisions.target` |

## Derived — computed by `units.derive`, not stored

| Quantity | Value | From |
|---|---|---|
| omega_pe [rad/s] | 1.38187e+15 | sqrt(n0 q^2/(eps0 m_e)) |
| d_e [m] | 2.16947e-07 | c/omega_pe |
| d_i,ab [m] | 9.29625e-06 | d_e*sqrt(m_i/m_e) |
| C_s,ab / c | 0.00707848 | sqrt(theta_e,ab/mass_ratio) |
| t_ab [s] | 4.38073e-12 | d_i,ab / C_s,ab |
| n_amb [m^-3] | 4.8e+24 | density_over_n0 * n0 |
| d_i0 [m] | 0.000103935 | c/omega_pi(n_amb) |
| d_i0 / d_e | 479.081 | derived |
| **omega_ci0 [rad/s]** | 6.73139e+09 | **q_e*B0/m_i — B0 only, no n_amb** |
| 1/omega_ci0 [s] | 1.48558e-10 | 1/omega_ci0 |
| **v_A / c** | 0.00233371 | **B0/sqrt(mu0*n_amb*m_i) — DERIVED** |
| rho_i0 / d_e | 4982.44 | v_p/omega_ci0 / d_e |
| v_p / c (model) | 0.0242706 | config model.vp_over_c |
| v_sh / c (model) | 0.032561 | model.vsh_over_Csab * C_s,ab |
| M_A | 13.9525 | v_sh/v_A |
| M_ms | 12.7595 | v_sh/sqrt(v_A^2+C_s0^2) |
| beta_ab | 1150 | mu0*n0*T_e,ab/B0^2 |
| beta_0 | 0.195745 | mu0*n_amb*T_0/B0^2 |
| dz [m] | 8.13444e-08 | dz_over_de * d_e |
| dt [s] | 2.03502e-16 | CFL-limited |
| dt*omega_pe | 0.281213 | dt * omega_pe |
| n_cell | 70272 | domain / dz (halved when one_sided) |
| steps per 1/omega_ci0 | 730008 | 1/(omega_ci0*dt) |
| run length [1/omega_ci0] | 2 | max_step * dt * omega_ci0 |
| T_e,ab [eV] | 47012 | theta_e,ab * m_e c^2 |
| lnLambda (used) | 122349 | units.coulomb_log_for(collisions.target) |
| lnLambda (physical) | 10.8364 | NRL at (n0, T_e,ab) |
| nu_ei,ab [1/s] | 2.09571e+13 | NRL electron-ion |
| nu_ei*dt | 0.00426481 | must be << 1 |
| mfp_ei,ab / d_e | 20 | v_te/nu_ei / d_e |
| lambda_ab | 20 | omega_ce,ab/nu_ei,ab = mfp/d_e,ab |
| mfp_ii,amb / d_i0 | 0.0033422 | upstream ion-ion — what criterion 2 tests |

## vs Schaeffer 2020 Table I

| Quantity | This run | Table I | ratio | known cause if off |
|---|---|---|---|---|
| mass ratio m_i/m_e | 1836.15 | 100 | 18.361x (+1736.1%) **OFF** |  |
| C_s,ab / c | 0.00707848 | 0.03 | 0.236x (-76.4%) **OFF** | theta_e_heat recalibrated off the paper's 0.092 |
| piston speed v_p / c | 0.0242706 | 0.104 | 0.233x (-76.7%) **OFF** |  |
| Alfven Mach v_sh/v_A | 13.9525 | 14 | 0.997x (-0.3%) ok |  |
| magnetosonic Mach | 12.7595 | 13 | 0.981x (-1.9%) ok |  |
| ablation beta, mu0*n*T/B^2 (Table I 1150) | 1150 | 1150 | 1.000x (-0.0%) ok |  |
| upstream beta, mu0*n*T/B^2 (Table I 0.2) | 0.195745 | 0.2 | 0.979x (-2.1%) ok |  |
| d_i0 / d_i,ab | 11.1803 | 11.18 | 1.000x (+0.0%) ok |  |
| gyroperiod in ablation times | 33.9116 | 33.9 | 1.000x (+0.0%) ok |  |
| collisionality mfp/d_e,ab | 20 | 20 | 1.000x (-0.0%) ok |  |

## Files

| Path | What |
|---|---|
| `config.yaml` | the primaries (tracked) |
| `inputs_kinshock_um1836_47keV` | generated deck (tracked, **never hand-edit**) |
| `warpx_used_inputs` | what WarpX actually ran (tracked) |
| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |
| `diags/`, `*.log` | output (gitignored, regenerable) |

---

_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, not the tables. Prose between the `prose:` markers is preserved._
