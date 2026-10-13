#!/usr/bin/env bash
# Run every tool's test suite. Prints one line per tool and a grand total.
# Usage: bash scripts/test.sh
set -uo pipefail

cd "$(dirname "$0")/.."

total_pass=0
total_fail=0
status=0

for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper; do
    line=$(cd "tools/$t" && cargo test --offline 2>&1 | grep -E "^test result" \
        | awk -F'[ ;]' '{p+=$4; f+=$6} END {printf "%d %d", p, f}')
    pass=${line% *}
    fail=${line#* }
    total_pass=$((total_pass + pass))
    total_fail=$((total_fail + fail))
    printf "%-10s %3d passed, %d failed\n" "$t" "$pass" "$fail"
    [ "$fail" -gt 0 ] && status=1
done

printf "%-10s %3d passed, %d failed\n" "TOTAL" "$total_pass" "$total_fail"
exit $status
