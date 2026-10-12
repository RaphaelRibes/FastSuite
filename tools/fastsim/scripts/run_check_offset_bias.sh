#!/bin/bash
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
/usr/bin/python3 scripts/check_offset_bias.py profiles/mandril_E.json
