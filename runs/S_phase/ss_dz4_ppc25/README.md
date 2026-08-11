# Run `ss_dz4_ppc25`

> R1_paper_470eV on a truncated domain and window, dz = 0.075 d_e,ab (dz/lambda_D = 1.516), ppc 25 (N_D = 16.5). Early-formation resolution sweep point; see runs/S_phase/README.md.

**tier** Sweep · **deck** `inputs_kinshock_ss_dz4_ppc25` · **paper** Schaeffer et al., Phys. Plasmas 27, 042901 (2020), Table I + Sec. II

## Notes

<!-- prose:begin -->
_No hand-written notes yet. Anything written between the
`prose:begin` / `prose:end` markers survives regeneration — put the run's story here
(why it exists, what it showed, what to distrust)._
<!-- prose:end -->

## Primaries — `config.yaml` is the only source of truth

Edit these; never the deck. Regenerate with `python scripts/make_inputs.py runs/ss_dz4_ppc25`.

| Quantity | Value | config key |
|---|---|---|
| reference density n0 [m^-3] | 6.0e26 | `reference.n0` |
| mass ratio m_i/m_e | 100 | `reference.mass_ratio` |
| charge state Z | 1 | `reference.charge_state` |
| target density / n0 | 2 | `plasma.piston.density_over_n0` |
| heater theta_e,ab | 0.000919767 | `plasma.piston.theta_e_heat` |
| piston init theta_e | 1.95695e-05 | `plasma.piston.theta_e_init` |
| piston init theta_i | 1.95695e-07 | `plasma.piston.theta_i_init` |
| ambient density / n0 | 0.008 | `plasma.ambient.density_over_n0` |
| ambient theta_0 | 1.95695e-05 | `plasma.ambient.theta_0` |
| field orientation | perpendicular | `field.orientation` |
| **B0 [T]** | 7.02645 | `field.B0_tesla` |
| dims | 1 | `geometry.dims` |
| layout | one_sided | `geometry.layout` |
| slab halfwidth [d_i] | 2 | `geometry.slab_halfwidth_di` |
| domain halfwidth [d_e] | 1344 | `geometry.domain_halfwidth_de` |
| dz [d_e] | 0.075 | `geometry.dz_over_de` |
| boundary lo / hi | lo symmetry, hi open | `geometry.boundary` |
| CFL | 0.75 | `numerics.cfl` |
| particle shape | 2 | `numerics.particle_shape` |
| max_step | 600000 | `numerics.max_step` |
| ppc | piston 25, ambient 25 | `numerics.ppc` |
| heater intervals | 20 | `operators.heater.intervals` |
| injector intervals | 20 | `operators.injector.intervals` |
| injector tau [1/wpe] | 400 | `operators.injector.tau_over_wpe_inv` |
| collisions target | quantity lambda_ab, value 20 | `collisions.target` |

## Derived — computed by `units.derive`, not stored

| Quantity | Value | From |
|---|---|---|
| omega_pe [rad/s] | 1.38187e+15 | sqrt(n0 q^2/(eps0 m_e)) |
| d_e [m] | 2.16947e-07 | c/omega_pe |
| d_i,ab [m] | 2.16947e-06 | d_e*sqrt(m_i/m_e) |
| C_s,ab / c | 0.00303277 | sqrt(theta_e,ab/mass_ratio) |
| t_ab [s] | 2.38613e-12 | d_i,ab / C_s,ab |
| n_amb [m^-3] | 4.8e+24 | density_over_n0 * n0 |
| d_i0 [m] | 2.42554e-05 | c/omega_pi(n_amb) |
| d_i0 / d_e | 111.803 | derived |
| **omega_ci0 [rad/s]** | 1.23583e+10 | **q_e*B0/m_i — B0 only, no n_amb** |
| 1/omega_ci0 [s] | 8.09176e-11 | 1/omega_ci0 |
| **v_A / c** | 0.000999873 | **B0/sqrt(mu0*n_amb*m_i) — DERIVED** |
| rho_i0 / d_e | 1162.76 | v_p/omega_ci0 / d_e |
| v_p / c (model) | 0.0103987 | config model.vp_over_c |
| v_sh / c (model) | 0.0139507 | model.vsh_over_Csab * C_s,ab |
| M_A | 13.9525 | v_sh/v_A |
| M_ms | 12.7595 | v_sh/sqrt(v_A^2+C_s0^2) |
| beta_ab | 1150 | mu0*n0*T_e,ab/B0^2 |
| beta_0 | 0.195745 | mu0*n_amb*T_0/B0^2 |
| dz [m] | 1.6271e-08 | dz_over_de * d_e |
| dt [s] | 4.07057e-17 | CFL-limited |
| dt*omega_pe | 0.05625 | dt * omega_pe |
| n_cell | 17920 | domain / dz (halved when one_sided) |
| steps per 1/omega_ci0 | 1.98787e+06 | 1/(omega_ci0*dt) |
| run length [1/omega_ci0] | 0.301831 | max_step * dt * omega_ci0 |
| T_e,ab [eV] | 470 | theta_e,ab * m_e c^2 |
| lnLambda (used) | 12.2287 | units.coulomb_log_for(collisions.target) |
| lnLambda (physical) | 6.231 | NRL at (n0, T_e,ab) |
| nu_ei,ab [1/s] | 2.09544e+12 | NRL electron-ion |
| nu_ei*dt | 8.52965e-05 | must be << 1 |
| mfp_ei,ab / d_e | 20 | v_te/nu_ei / d_e |
| lambda_ab | 20 | omega_ce,ab/nu_ei,ab = mfp/d_e,ab |
| mfp_ii,amb / d_i0 | 0.0143214 | upstream ion-ion — what criterion 2 tests |

## vs Schaeffer 2020 Table I

| Quantity | This run | Table I | ratio | known cause if off |
|---|---|---|---|---|
| mass ratio m_i/m_e | 100 | 100 | 1.000x (+0.0%) ok |  |
| C_s,ab / c | 0.00303277 | 0.03 | 0.101x (-89.9%) **OFF** | theta_e_heat recalibrated off the paper's 0.092 |
| piston speed v_p / c | 0.0103987 | 0.104 | 0.100x (-90.0%) **OFF** |  |
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
| `inputs_kinshock_ss_dz4_ppc25` | generated deck (tracked, **never hand-edit**) |
| `warpx_used_inputs` | what WarpX actually ran (tracked) |
| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |
| `diags/`, `*.log` | output (gitignored, regenerable) |

---

_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, not the tables. Prose between the `prose:` markers is preserved._
