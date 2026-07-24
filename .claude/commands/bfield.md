---
description: Diagnose a run's B-field fluctuations — physical vs numerical
argument-hint: [run_id]
---
Assess whether run `$1`'s magnetic-field oscillations are physical shock turbulence
or a numerical (grid) artifact:

1. `python scripts/bfield_diagnostic.py runs/$1 --twci 1.4` (adjust time to a well-formed
   shock frame). Report: dz/lambda_D, per-zone spectra + polarization, and the
   particle-response test (does dB/B~1 scatter ions/electrons, or are they at the floor?).
2. Interpret against the rule: near-shock foot turbulence (reflected ions present) is
   physical; far-upstream fluctuation that doesn't move particles and rises toward the
   grid scale is numerical.
3. If inconclusive, point me at `studies/bfield_convergence/` to run the knob-variant test.
Record the verdict in RESULTS.md.
