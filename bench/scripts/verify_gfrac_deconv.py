#!/usr/bin/env python3
"""SUPERSEDED — use `bench/scripts/gfrac_attribution.py`.

This script printed one number (position-1 G fraction, real vs sim) for one
hard-coded dataset, from a hard-coded corpus path and a `/tmp` sim report. It
was the first probe of the metagenomics g_frac blocker and it is kept only as a
pointer: it answered "is the head deconvolved" and nothing else, and the
answer to the blocker is not one number.

What replaced it, and what each prints:

  gfrac_attribution.py     per (dataset, variant): the real/sim per-position G
                           curve, the diff, its argmax and zone means, the
                           per-mate read lengths, the per-mate terminal G, the
                           LENGTH-MATCHED counterfactual gate, the fitted
                           parameters that drive the G curve, and the N-overlay
                           control. Writes a TSV; refuses to overwrite.
  gfrac_mechanism_summary.py
                           joins that TSV with the fitted run rates and the
                           run-structure agreement (TV) so each candidate
                           mechanism's contribution is one column.
  probe_plane_arms.py      mutates one block of the fitted profile, re-simulates
                           at the same seed and measures the gate -- the
                           counterfactual for each proposed fix.
  probe_frame_shift.py     tests whether the simulated curve is the real one
                           convolved with the head-run distribution.

Findings: docs/GFRAC-MECHANISM.md.
"""
import sys

sys.exit(__doc__)
