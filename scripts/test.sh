#!/usr/bin/env bash
# Run every tool's test suite. Prints one line per tool and a grand total.
# Usage: bash scripts/test.sh
set -uo pipefail

cd "$(dirname "$0")/.."

total_pass=0
total_fail=0
status=0

for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper; do
    out=$(cd "tools/$t" && cargo test --offline 2>&1)
    rc=$?
    line=$(echo "$out" | grep -E "^test result" \
        | awk -F'[ ;]' '{p+=$4; f+=$6} END {printf "%d %d", p, f}')
    pass=${line% *}
    fail=${line#* }
    # A non-zero cargo exit, or zero tests parsed, means the suite did not
    # run (build/link failure) — report it as a failure, never a silent 0.
    if [ "$rc" -ne 0 ] || [ "$pass" -eq 0 ]; then
        echo "$out" | grep -iE "^error|error\[|error:" | head -3
        fail=1
    fi
    total_pass=$((total_pass + pass))
    total_fail=$((total_fail + fail))
    printf "%-10s %3d passed, %d failed\n" "$t" "$pass" "$fail"
    [ "$fail" -gt 0 ] && status=1
done

printf "%-10s %3d passed, %d failed\n" "TOTAL" "$total_pass" "$total_fail"

# The pathology corpus: adversarial fixtures for the five defects Phase D found
# only by real-data measurement (none caught by the unit suite). Runs the release
# binaries against committed fixtures in bench/fixtures and asserts the post-fix
# behaviour. Kept out of the per-tool loop above because it is a cross-binary
# integration test, not a cargo test.
echo
if bash bench/tests/pathology_test.sh; then
    total_pass=$((total_pass + 1))
    printf "%-10s %3d passed, %d failed\n" "pathology" 1 0
else
    total_fail=$((total_fail + 1))
    printf "%-10s %3d passed, %d failed\n" "pathology" 0 1
    status=1
fi

printf "%-10s %3d passed, %d failed\n" "GRAND" "$total_pass" "$total_fail"
exit $status
