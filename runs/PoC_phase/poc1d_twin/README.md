# Run `poc1d_twin`

> 1D twin of PoC_2d_v2 (FLASH OmegaShock_07-10 surrogate): CH ion Z=3.5, m_i=100 m_e, n_t:n_amb = 400, B decaying 8 T -> 1 T (floored), drive off at 2.1 ns, to 7.05 ns; same dz/dt/shape/filter as the 2D run, 1024 ppc.

**tier** Pilot · **deck** `inputs_kinshock_poc1d_twin` · **paper** Collaborator's PoC_2d_v2 config; FLASH cellar/data/OmegaShock_07-10

## Notes

<!-- prose:begin -->
_No hand-written notes yet. Anything written between the
`prose:begin` / `prose:end` markers survives regeneration — put the run's story here
(why it exists, what it showed, what to distrust)._
<!-- prose:end -->

## Primaries — `config.yaml` is the only source of truth

Edit these; never the deck. Regenerate with `python scripts/make_inputs.py runs/poc1d_twin`.

| Quantity | Value | config key |
|---|---|---|
| reference density n0 [m^-3] | 1.0e18 | `reference.n0` |
| mass ratio m_i/m_e | 100 | `reference.mass_ratio` |
| charge state Z | 3.5 | `reference.charge_state` |
| target density / n0 | 5 | `plasma.piston.density_over_n0` |
| heater theta_e,ab | 0.037 | `plasma.piston.theta_e_heat` |
| piston init theta_e | 0.001 | `plasma.piston.theta_e_init` |
| piston init theta_i | 1e-07 | `plasma.piston.theta_i_init` |
| ambient density / n0 | 0.0125 | `plasma.ambient.density_over_n0` |
| ambient theta_0 | 0.00100559 | `plasma.ambient.theta_0` |
| field orientation | perpendicular | `field.orientation` |
| **B0 [T]** | 0.00325533 | `field.B0_tesla` |
| dims | 1 | `geometry.dims` |
| layout | one_sided | `geometry.layout` |
| slab halfwidth [d_i] | 2 | `geometry.slab_halfwidth_di` |
| domain halfwidth [d_e] | 1235.86 | `geometry.domain_halfwidth_de` |
| dz [d_e] | 0.367816 | `geometry.dz_over_de` |
| boundary lo / hi | lo reflecting, hi open | `geometry.boundary` |
| CFL | 0.353553 | `numerics.cfl` |
| particle shape | 2 | `numerics.particle_shape` |
| max_step | 63000 | `numerics.max_step` |
| ppc | piston 1024, ambient 1024 | `numerics.ppc` |
| heater intervals | 20 | `operators.heater.intervals` |
| injector intervals | 20 | `operators.injector.intervals` |
| injector tau [1/wpe] | 40 | `operators.injector.tau_over_wpe_inv` |
| drive stop [t*w_ci0] | 0.8728 | `operators.stop_t_wci0` |
| B profile | type inverse_linear, decay_length_di0 1.37, floor_over_B0 0.125, apply particle | `field.profile` |

## Derived — computed by `units.derive`, not stored

| Quantity | Value | From |
|---|---|---|
| omega_pe [rad/s] | 5.64146e+10 | sqrt(n0 q^2/(eps0 m_e)) |
| d_e [m] | 0.00531409 | c/omega_pe |
| d_i,ab [m] | 0.028405 | d_e*sqrt(m_i/(Z m_e)) |
| C_s,ab / c | 0.0359861 | sqrt(Z theta_e,ab/mass_ratio) |
| t_ab [s] | 2.63293e-09 | d_i,ab / C_s,ab |
| n_amb [m^-3] | 1.25e+16 | density_over_n0 * n0 |
| d_i0 [m] | 0.254062 | c/omega_pi(n_amb) |
| d_i0 / d_e | 47.8091 | derived |
| **omega_ci0 [rad/s]** | 2.00394e+07 | **Z q_e*B0/m_i — B0 only, no n_amb** |
| 1/omega_ci0 [s] | 4.99017e-08 | 1/omega_ci0 |
| **v_A / c** | 0.0169826 | **B0/sqrt(mu0*(n_amb/Z)*m_i) — DERIVED** |
| rho_i0 / d_e | 306.168 | v_p/omega_ci0 / d_e |
| v_p / c (model) | 0.108756 | config model.vp_over_c |
| v_sh / c (model) | 0.10875 | model.vsh_over_Csab * C_s,ab |
| M_A | 6.40364 | v_sh/v_A |
| M_ms | 6.04538 | v_sh/sqrt(v_A^2+C_s0^2) |
| beta_ab | 359.212 | mu0*n0*T_e,ab/B0^2 |
| beta_0 | 0.122034 | mu0*n_amb*T_0/B0^2 |
| dz [m] | 0.00195461 | dz_over_de * d_e |
| dt [s] | 2.30512e-12 | CFL-limited |
| dt*omega_pe | 0.130043 | dt * omega_pe |
| n_cell | 3360 | domain / dz (halved when one_sided) |
| steps per 1/omega_ci0 | 21648.2 | 1/(omega_ci0*dt) |
| run length [1/omega_ci0] | 2.91018 | max_step * dt * omega_ci0 |
| T_e,ab [eV] | 18907 | theta_e,ab * m_e c^2 |

## vs Schaeffer 2020 Table I

| Quantity | This run | Table I | ratio | known cause if off |
|---|---|---|---|---|
| mass ratio m_i/m_e | 100 | 100 | 1.000x (+0.0%) ok |  |
| C_s,ab / c | 0.0359861 | 0.03 | 1.200x (+20.0%) ~ | theta_e_heat recalibrated off the paper's 0.092 |
| piston speed v_p / c | 0.108756 | 0.104 | 1.046x (+4.6%) ok |  |
| Alfven Mach v_sh/v_A | 6.40364 | 14 | 0.457x (-54.3%) **OFF** | model M_A from model.vsh_over_Csab; the by-eye settled value is in shock_fit.yaml |
| magnetosonic Mach | 6.04538 | 13 | 0.465x (-53.5%) **OFF** |  |
| ablation beta, mu0*n*T/B^2 (Table I 1150) | 359.212 | 1150 | 0.312x (-68.8%) **OFF** | n_amb 0.01 vs Table I's 0.008, and/or theta_e off 0.092 |
| upstream beta, mu0*n*T/B^2 (Table I 0.2) | 0.122034 | 0.2 | 0.610x (-39.0%) **OFF** |  |
| d_i0 / d_i,ab | 8.94427 | 11.18 | 0.800x (-20.0%) ~ | n_amb is 0.01 n0; Table I is 0.008 n0 |
| gyroperiod in ablation times | 18.9529 | 33.9 | 0.559x (-44.1%) **OFF** | n_amb 25% high (1.118x) and/or theta_e recal (1.086x) |

## Files

| Path | What |
|---|---|
| `config.yaml` | the primaries (tracked) |
| `inputs_kinshock_poc1d_twin` | generated deck (tracked, **never hand-edit**) |
| `warpx_used_inputs` | what WarpX actually ran (tracked) |
| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |
| `diags/`, `*.log` | output (gitignored, regenerable) |

---

_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, not the tables. Prose between the `prose:` markers is preserved._
