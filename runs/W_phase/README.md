# W_phase — a second driver: the reflecting wall

**The question.** Every diagnosis of the notch so far has been internal to **one driver** —
the custom kinetic heating operator plus the target injector, ported from PSC and validated
only on a 2D flatfoil case (`warpx-cda/heating_operator/run_flatfoil_compare`). That is the
most WarpX-specific, least standard component in the whole setup, and the notch is a
**formation-stage** feature sitting on the C, which is precisely the ambipolar structure the
heater and injector create.

The paper's own §I makes the driver a physics choice rather than a detail: most
magnetized-shock simulations use counter-streaming flows or an expanding slab, and are
"concerned with timescales long after the shock has formed, **when the dynamics of the shock
is independent of its creation**" — whereas experiments live in the formation stage, where
"observable features may be dominated by the formation mechanisms rather than the shock."

So: drive the same upstream a second, standard way. If the notch is there too, it is the
plasma. If it is not, it is the driver.

## The smoke test comes first

`ws_smoke_470eV` is **not** the comparison — it is the go/no-go. The reflecting boundary has
never been exercised as a *driver* in this deck (only as a domain edge), and
`species.*.drift_uz_over_c` is a new key in `kinshock.deck`. Three yes/no questions:

1. does the deck build and run with **no piston species, no heater, no injector**?
2. does the drifting ambient actually reflect off the lo wall?
3. does a density/field ramp appear and **move away** from the wall?

| | value |
|---|---|
| driver | ambient drifts at −0.01395 c into a **specular** reflecting wall at z = 0 |
| upstream | `M_A` = 13.95, `M_ms` = 12.76, `beta_0` = 0.196, `mu` = 100 — same as the partner |
| grid | 4,480 cells, `dz/lambda_D` = 6.07, ppc 100 — same as `S_phase/ss_dz1_ppc100_int1` |
| clock | 74,545 steps, `t*` = 0.15 |
| cost | **~6 GPU-min** |

`t*` = 0.15 is a quarter of the partner's window: long enough for a ramp to separate from the
wall, far too short for the notch. That is deliberate — this run is not asked to show the
notch, only to show that the driver works.

**`reflecting`, not `symmetry`.** `symmetry` applies the π-rotation about B0 that models a
mirror-image half-domain; here the wall *is* the driver and must reflect specularly. The
generated deck carries no `reflect_symmetry_axis` line, which is the check.

## What cannot be matched, and is not pretended to be

There is no piston population, so `v_p`, `n_t`, the target slab and every
`t_ab`-referenced piston diagnostic have no counterpart — and **`z_piston`, the origin every
notch measurement uses, does not exist**. The wall is the reference instead. What carries
over is the dimensionless upstream, the grid, the clock, and the observable that matters:
the ambient-ion phase space, the C, and the notch on its lower branch.

In the wall frame the upstream flows in at `v_sh` and the shock propagates back out into it,
so the drive speed is set from `v_sh` (Table I's 4.6 `C_s,ab`), not from `v_p`.

## After the smoke test — RAN 2026-09-06, and the plan above was wrong twice

The smoke test passed (drift to 0.01%, reflection 0 → 0.17, both ramps off the wall at
+11.5 d_i0/t*). The pair `ws_470eV` / `ws_47keV` then ran and both COMPLETED. **The driver
is validated; the physics question is not answered.** See RESULTS 2026-09-06.

The plan in the line above was wrong in two ways, both now fixed in the configs:

1. **Boxes.** A wall run is a finite slab with no injector behind it. The front advances at
   11.5 d_i0/t* while the upstream drifts in at 13.952, closing at **25.6 d_i0/t\***
   (predicted 25.5). A 12 d_i0 box exhausts at t* = 0.47, before the notch forms. The runs
   were built at 22 and 80 d_i0 instead.
2. **Grids.** Twinning `ss_dz1_ppc100_int1` would have put the test at N_D = 16.5, where
   onset is noise-seeded (RESULTS 2026-09-05) and a notch would prove nothing. Both runs
   take N_D = 131.9, matching `R1_470eV_s3f8_t2` and `ek_47keV_t2`.

**And the clock was still wrong.** `t*` = 0.70 was chosen to cover the PISTON run's notch
window — but the wall driver is ~2× slower, and the 47 keV control shows compression only
steepening at t* ≈ 1.4 (it is at 2.51 at t* = 0.70, where `ws_470eV` ended at 2.44). The two
drivers must be matched on **shock formation**, not on elapsed t*. Settling the 470 eV case
needs t* ≈ 1.4–1.5, which needs a 40–44 d_i0 box: **~24–28 GPU-h**, not staged.

Costs actually incurred: 3:08, 7:15 and 6:32 wall-clock (~0.05, 0.12 and 6.7 GPU-h).

## Run it

```bash
perlmutter/submit.sh wall --dry
perlmutter/submit.sh wall
```
