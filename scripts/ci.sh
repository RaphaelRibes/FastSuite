#!/usr/bin/env bash
# FastSuite CI gate — the check the suite was missing.
#
# scripts/test.sh prints one line per tool and a GRAND total; on its own it
# cannot tell "471 tests passed" from "463 tests ran because 8 test binaries
# were stale and never recompiled". This script closes that hole: it runs the
# suite, parses the per-tool lines and asserts every count against the
# checked-in manifest scripts/expected-test-counts.txt. Short, long, missing or
# failing rows all exit non-zero with a per-tool expected-vs-actual diff.
#
# It also refuses to run on a rustc/cargo outside the range pixi.toml pins
# (`rust = ">=1.99,<2"`), because the same tree compiled with the ambient 1.95
# toolchain produces a different (1.44x slower) binary and nothing in the log
# records which compiler ran.
#
# Usage: pixi run ci            # sanctioned: pixi puts the env's rust 1.99 first
#        bash scripts/ci.sh     # same check, but warns/fails on a foreign rustc
# No network access, creates no files, safe to run repeatedly.
set -uo pipefail

EXIT_COUNTS=1     # count mismatch, missing row, or a tool reporting failures
EXIT_TOOLCHAIN=2  # rustc/cargo outside the pixi.toml pin (suite not run)
EXIT_SETUP=3      # usage error, unreadable/invalid manifest, unparsable output

COUNTS="scripts/expected-test-counts.txt"
TESTS_OUTPUT=""
RUSTC_CMD="${RUSTC:-rustc}"
CARGO_CMD="${CARGO:-cargo}"
ALLOW_RUSTC=0
PIN_FALLBACK_MIN="1.99"
PIN_FALLBACK_MAX="2"

usage() {
    cat <<'EOF'
usage: bash scripts/ci.sh [options]

Runs scripts/test.sh from the repo root and asserts every tool's test count
against the expected-counts manifest (a staleness detector: a short count means
cargo skipped a rebuild and tests silently did not run).

  --counts <file>         Manifest to compare against
                          (default: scripts/expected-test-counts.txt).
                          A relative path resolves against the repo root.
  --tests-output <file>   Parse an already-saved scripts/test.sh output instead
                          of running the suite (test hook for the parser; the
                          file must hold raw `bash scripts/test.sh` output).
  --rustc <cmd>           rustc to check (default: $RUSTC, else `rustc`).
  --cargo <cmd>           cargo to check (default: $CARGO, else `cargo`).
  --allow-rustc-mismatch  Run anyway when rustc/cargo are outside the pin
                          (deliberate bypass; prints a warning, not for CI).
  -h, --help              This help.

Exit codes: 0 everything matched and nothing failed
            1 count mismatch / missing tool row / failing test
            2 toolchain outside the pin (suite not run)
            3 usage, manifest or parse error
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --counts)               COUNTS="${2:-}"; shift 2 ;;
        --counts=*)             COUNTS="${1#*=}"; shift ;;
        --tests-output)         TESTS_OUTPUT="${2:-}"; shift 2 ;;
        --tests-output=*)       TESTS_OUTPUT="${1#*=}"; shift ;;
        --rustc)                RUSTC_CMD="${2:-}"; shift 2 ;;
        --rustc=*)              RUSTC_CMD="${1#*=}"; shift ;;
        --cargo)                CARGO_CMD="${2:-}"; shift 2 ;;
        --cargo=*)              CARGO_CMD="${1#*=}"; shift ;;
        --allow-rustc-mismatch) ALLOW_RUSTC=1; shift ;;
        -h|--help)              usage; exit 0 ;;
        *) echo "ci.sh: unknown argument: $1" >&2; usage >&2; exit $EXIT_SETUP ;;
    esac
done
if [ -z "$COUNTS" ] || [ -z "$RUSTC_CMD" ] || [ -z "$CARGO_CMD" ]; then
    echo "ci.sh: --counts/--rustc/--cargo need a value" >&2
    exit $EXIT_SETUP
fi

cd "$(dirname "$0")/.." || exit $EXIT_SETUP
REPO=$(pwd -P)

# ---------------------------------------------------------------------------
# Toolchain guard.
#
# The pin is read out of pixi.toml rather than hard-coded here, so the two can
# never drift; the fallback is only used if pixi.toml becomes unparseable, and
# says so out loud.
# ---------------------------------------------------------------------------
PIN_SPEC=""
[ -f pixi.toml ] && PIN_SPEC=$(grep -m1 -E '^[[:space:]]*rust[[:space:]]*=' pixi.toml \
    | sed -E 's/^[^"]*"([^"]*)".*$/\1/')
PIN_MIN="$PIN_FALLBACK_MIN"
PIN_MAX="$PIN_FALLBACK_MAX"
PIN_SOURCE="built-in fallback (pixi.toml rust pin not parsed)"
if [ -n "$PIN_SPEC" ]; then
    pmin=$(printf '%s' "$PIN_SPEC" | sed -n 's/.*>=[[:space:]]*\([0-9][0-9.]*\).*/\1/p')
    pmax=$(printf '%s' "$PIN_SPEC" | sed -n 's/.*<[[:space:]]*\([0-9][0-9.]*\).*/\1/p')
    if [ -n "$pmin" ] && [ -n "$pmax" ]; then
        PIN_MIN="$pmin"
        PIN_MAX="$pmax"
        PIN_SOURCE="pixi.toml rust = \"$PIN_SPEC\""
    fi
fi
PIN_HUMAN=">=$PIN_MIN,<$PIN_MAX"

# version_in_range VERSION MIN MAX — exit 0 iff MIN <= VERSION < MAX.
version_in_range() {
    awk -v v="$1" -v lo="$2" -v hi="$3" '
        function scale(s, p) { split(s, p, "."); return (p[1]+0)*1000000 + (p[2]+0)*1000 + (p[3]+0) }
        BEGIN { exit !(scale(v) >= scale(lo) && scale(v) < scale(hi)) }'
}
# tool_version CMD -> "1.99.0" (empty if the command is missing or silent).
tool_version() {
    "$1" --version 2>/dev/null | awk 'NR==1{print $2}' \
        | sed -E 's/^([0-9]+(\.[0-9]+)*).*$/\1/'
}
majmin() { printf '%s' "$1" | sed -E 's/^([0-9]+\.[0-9]+).*$/\1/'; }

rustc_ver=$(tool_version "$RUSTC_CMD")
rustc_path=$(command -v "$RUSTC_CMD" 2>/dev/null || printf '%s' "$RUSTC_CMD")
cargo_ver=$(tool_version "$CARGO_CMD")
cargo_path=$(command -v "$CARGO_CMD" 2>/dev/null || printf '%s' "$CARGO_CMD")

tc_reasons=""
if [ -z "$rustc_ver" ]; then
    tc_reasons="${tc_reasons}  rustc: not found or produced no version (command: $RUSTC_CMD)
"
elif ! version_in_range "$rustc_ver" "$PIN_MIN" "$PIN_MAX"; then
    tc_reasons="${tc_reasons}  rustc in use : $rustc_ver at $rustc_path  [OUT OF RANGE, needs $PIN_HUMAN]
"
fi
if [ -z "$cargo_ver" ]; then
    tc_reasons="${tc_reasons}  cargo: not found or produced no version (command: $CARGO_CMD)
"
elif ! version_in_range "$cargo_ver" "$PIN_MIN" "$PIN_MAX"; then
    tc_reasons="${tc_reasons}  cargo in use : $cargo_ver at $cargo_path  [OUT OF RANGE, needs $PIN_HUMAN]
"
fi
# A mixed PATH (new cargo driving an old rustc, or the reverse) is its own
# failure mode: cargo and rustc must come from the same toolchain.
if [ -n "$rustc_ver" ] && [ -n "$cargo_ver" ] && [ "$(majmin "$rustc_ver")" != "$(majmin "$cargo_ver")" ]; then
    tc_reasons="${tc_reasons}  MIXED toolchain: rustc $rustc_ver ($rustc_path) does not match cargo $cargo_ver ($cargo_path)
"
fi

echo "==> repo      $REPO"
echo "==> toolchain pin $PIN_HUMAN  ($PIN_SOURCE)"
if [ -z "$tc_reasons" ]; then
    echo "==> toolchain ok  rustc $rustc_ver ($rustc_path), cargo $cargo_ver ($cargo_path)"
else
    echo "==> toolchain BAD rustc ${rustc_ver:-none} ($rustc_path), cargo ${cargo_ver:-none} ($cargo_path)"
    if [ "$ALLOW_RUSTC" -eq 0 ]; then
        printf 'ci.sh: toolchain mismatch — refusing to run the suite (exit %d)\n\n' "$EXIT_TOOLCHAIN" >&2
        printf '%s' "$tc_reasons" >&2
        cat >&2 <<EOF

  Why this is fatal rather than a warning: this repo measured rust 1.99 at 1.44x
  the speed of 1.95 on identical source, and its two machines carry different
  ambient toolchains (1.99 desktop, 1.95 laptop via ~/.cargo/bin). A bare
  \`bash scripts/test.sh\` therefore compiles a different program on each machine
  with nothing in either log recording it. The 463-vs-471 incident was exactly
  this shape — a divergence that stayed invisible until someone compared
  numbers — so a warning here would be ignored just as the count was.

  Sanctioned invocation (pixi puts the env's rust 1.99 first on PATH):
      pixi run ci
      pixi run bash scripts/ci.sh

  Deliberate, non-CI bypass (warns and runs anyway):
      bash scripts/ci.sh --allow-rustc-mismatch
EOF
        exit $EXIT_TOOLCHAIN
    fi
    echo "==> WARNING: --allow-rustc-mismatch given; running the suite on an unpinned toolchain"
fi

# ---------------------------------------------------------------------------
# Read the manifest: `tool count` per line, `#` comments and blanks ignored.
# Validated before the suite runs, so a malformed manifest fails fast.
# ---------------------------------------------------------------------------
if [ ! -f "$COUNTS" ]; then
    echo "ci.sh: manifest not found: $COUNTS  (relative paths resolve against $REPO)" >&2
    exit $EXIT_SETUP
fi

declare -A exp=()
exp_order=()
manifest_bad=0
lineno=0
while IFS= read -r raw || [ -n "$raw" ]; do
    lineno=$((lineno + 1))
    line=$(printf '%s' "${raw%%#*}" | awk '{$1=$1; print}')
    [ -n "$line" ] || continue
    # shellcheck disable=SC2086
    set -- $line
    if [ "$#" -ne 2 ]; then
        echo "ci.sh: $COUNTS:$lineno: expected 'tool count', got: $line" >&2
        manifest_bad=1
        continue
    fi
    name=$1
    count=$2
    case "$name" in ''|*[!A-Za-z0-9_.-]*)
        echo "ci.sh: $COUNTS:$lineno: bad tool name: $name" >&2; manifest_bad=1; continue ;;
    esac
    case "$count" in ''|*[!0-9]*)
        echo "ci.sh: $COUNTS:$lineno: bad count for $name: $count" >&2; manifest_bad=1; continue ;;
    esac
    if [ -n "${exp[$name]+set}" ]; then
        echo "ci.sh: $COUNTS:$lineno: duplicate entry for $name" >&2; manifest_bad=1; continue
    fi
    exp[$name]=$count
    exp_order+=("$name")
done < "$COUNTS"
if [ "$manifest_bad" -ne 0 ]; then
    exit $EXIT_SETUP
fi
if [ "${#exp_order[@]}" -eq 0 ]; then
    echo "ci.sh: manifest $COUNTS has no entries" >&2
    exit $EXIT_SETUP
fi
if [ -z "${exp[GRAND]+set}" ]; then
    echo "ci.sh: manifest $COUNTS has no GRAND line (the gate needs it)" >&2
    exit $EXIT_SETUP
fi
sum_tools=0
for name in "${exp_order[@]}"; do
    [ "$name" = "GRAND" ] && continue
    sum_tools=$((sum_tools + exp[$name]))
done
if [ "${exp[GRAND]}" -ne "$sum_tools" ]; then
    echo "ci.sh: manifest $COUNTS is internally inconsistent: the tool lines sum to" >&2
    echo "       $sum_tools but GRAND says ${exp[GRAND]}. Fix the manifest before trusting" >&2
    echo "       either number (exit $EXIT_SETUP)." >&2
    exit $EXIT_SETUP
fi
echo "==> manifest  $COUNTS (${#exp_order[@]} entries, GRAND ${exp[GRAND]})"

# ---------------------------------------------------------------------------
# Run the suite (or replay a saved run) and parse its per-tool lines.
# ---------------------------------------------------------------------------
if [ -n "$TESTS_OUTPUT" ]; then
    if [ ! -f "$TESTS_OUTPUT" ]; then
        echo "ci.sh: --tests-output file not found: $TESTS_OUTPUT" >&2
        exit $EXIT_SETUP
    fi
    echo "==> reusing saved output: $TESTS_OUTPUT (suite NOT re-run)"
    raw=$(cat -- "$TESTS_OUTPUT")
    run_rc=0
else
    echo "==> running: bash scripts/test.sh"
    echo
    raw=$(bash scripts/test.sh 2>&1)
    run_rc=$?
fi
printf '%s\n' "$raw"

# `fastcheck  108 passed, 0 failed` -> "fastcheck 108 0"
actual_table=$(printf '%s\n' "$raw" \
    | grep -E '^[[:space:]]*[^[:space:]]+[[:space:]]+[0-9]+ passed,[[:space:]]*[0-9]+ failed[[:space:]]*$' \
    | awk '{print $1, $2, $4}')

declare -A act_pass=()
declare -A act_fail=()
while read -r name pass fail; do
    [ -n "${name:-}" ] || continue
    if [ -n "${act_pass[$name]+set}" ]; then
        echo "==> warning: duplicate '$name' row in test output (keeping the first)" >&2
        continue
    fi
    act_pass[$name]=$pass
    act_fail[$name]=$fail
done <<EOF
$actual_table
EOF

# ---------------------------------------------------------------------------
# Compare.
# ---------------------------------------------------------------------------
echo
echo "==> count check (expected vs actual)"
printf '  %-12s %8s %8s  %s\n' "ENTRY" "EXPECTED" "ACTUAL" "VERDICT"
printf '  %-12s %8s %8s  %s\n' "------------" "--------" "--------" "-------"

n_bad=0
problems=""
short_list=""
long_list=""
for name in "${exp_order[@]}"; do
    e=${exp[$name]}
    verdict=""
    if [ -z "${act_pass[$name]+set}" ]; then
        a="-"
        verdict="MISSING (test.sh printed no row)"
    else
        a=${act_pass[$name]}
        f=${act_fail[$name]}
        if [ "$f" -gt 0 ]; then
            verdict="FAILED ($f failed)"
        elif [ "$a" -eq "$e" ]; then
            verdict="ok"
        elif [ "$a" -lt "$e" ]; then
            verdict="SHORT $((a - e))"
        else
            verdict="LONG +$((a - e))"
        fi
    fi
    printf '  %-12s %8s %8s  %s\n' "$name" "$e" "$a" "$verdict"

    case "$verdict" in
        ok) continue ;;
    esac
    n_bad=$((n_bad + 1))
    # A row-level hint can name the tool's own directory; the aggregate GRAND row
    # cannot (there is no tools/GRAND), so it gets the generic wording.
    if [ "$name" = "GRAND" ]; then
        where="the tools whose rows above are short"
    else
        where="tools/$name"
    fi
    case "$verdict" in
        "SHORT -"*)
            short_list="${short_list:+$short_list, }${name} SHORT $((a - e))"
            problems="${problems}  ! $name expected $e, actual $a -> ${verdict}
      cargo considered a STALE test binary fresh, so part of this suite never
      compiled/ran. This is a real failure, not bookkeeping: rebuild cleanly
      through \`pixi run\` (or \`cargo clean\` in $where) and re-run.
      Do NOT lower the expected count.
" ;;
        "LONG +"*)
            long_list="${long_list:+$long_list, }${name} LONG +$((a - e))"
            problems="${problems}  ! $name expected $e, actual $a -> ${verdict}
      Tests were added (or removed) without updating the manifest. That is a
      deliberate, reviewable change: update scripts/expected-test-counts.txt
      (this line and the GRAND line) in the same commit.
" ;;
        "FAILED"*)
            problems="${problems}  ! $name expected $e, actual $a -> ${verdict}
      The tool's own test run reported failures; read its output above.
" ;;
        "MISSING"*)
            problems="${problems}  ! $name expected $e, actual $a -> ${verdict}
      scripts/test.sh never printed this row: the tool loop broke, was renamed,
      or its output changed shape. The gate cannot verify it.
" ;;
    esac
done

if [ "$run_rc" -ne 0 ]; then
    n_bad=$((n_bad + 1))
    problems="${problems}  ! scripts/test.sh exited $run_rc (its rows may match, but the suite failed)
"
fi

cases=""
[ -n "$short_list" ] && cases="${cases}short=[${short_list}] "
[ -n "$long_list" ] && cases="${cases}long=[${long_list}] "

# --- report-schema gate ---------------------------------------------------
# `fasthelper check --report` is the only consumer of REPORT-SCHEMA, and until
# now nothing ran it on the suite's own output: that is how fastdedup shipped
# `metrics.abundance = false` (a boolean in a numeric namespace, schema rule 1)
# while its validator rejected the file. Every tool's report must pass here.
schema_gate() {
    local tmp t rep rc=0
    tmp=$(mktemp -d "${TMPDIR:-/tmp}/fs-schema.XXXXXX") || return 1
    printf '@r1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' > "$tmp/a.fq"
    printf '@r2\nTTGGTTGGTTGGTTGGTTGGTTGGTTGGTTGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' > "$tmp/b.fq"
    local th=tools/fasthelper/target/release/fasthelper
    for t in fastdedup fastprep fastcheck fastsim fastmaps; do
        rep="$tmp/$t.report.json"
        case "$t" in
            fastdedup) tools/fastdedup/target/release/fastdedup -1 "$tmp/a.fq" -o "$tmp/d.fq" --json "$rep" -q --force >/dev/null 2>&1 ;;
            fastprep)  tools/fastprep/target/release/fastprep -1 "$tmp/a.fq" -o "$tmp/p.fq" --json "$rep" --no-adapter -q --force >/dev/null 2>&1 ;;
            fastcheck) tools/fastcheck/target/release/fastcheck "$tmp/a.fq" --json "$rep" -q >/dev/null 2>&1 ;;
            fastsim)   tools/fastsim/target/release/fastsim --profile tools/fastsim/profiles/test_small.json --pairs 50 -o "$tmp/sim" --json "$rep" --force >/dev/null 2>&1 ;;
            fastmaps)  tools/fastmaps/target/release/fastmaps -o "$tmp/maps" "$tmp/a.fq" --json "$rep" -q --force >/dev/null 2>&1 ;;
        esac
        if [ ! -s "$rep" ]; then
            echo "  ! $t wrote no --json report: the schema gate cannot verify it"
            rc=1
            continue
        fi
        if ! "$th" check --report "$rep" >/dev/null 2>"$tmp/$t.err"; then
            echo "  ! $t report failed fasthelper check --report:"
            sed 's/^/      /' "$tmp/$t.err"
            rc=1
        else
            echo "  $t report: schema ok"
        fi
    done
    rm -rf "$tmp"
    return "$rc"
}

echo
echo "==> report-schema gate (fasthelper check --report on each tool's own output)"
if ! schema_gate; then
    n_bad=$((n_bad + 1))
    problems="${problems}  ! a tool's --json report does not satisfy REPORT-SCHEMA (see the gate output above)
"
fi

if [ "$n_bad" -eq 0 ]; then
    printf 'ci: OK — %d/%d entries match (GRAND %s), 0 failing tests, rustc %s, manifest %s\n' \
        "${#exp_order[@]}" "${#exp_order[@]}" "${exp[GRAND]}" "$rustc_ver" "$COUNTS"
    exit 0
fi

echo
printf '%s' "$problems"
printf 'ci: FAIL — %d of %d manifest entries did not match %s(ci exit %d; scripts/test.sh exit %d)\n' \
    "$n_bad" "${#exp_order[@]}" "$cases" "$EXIT_COUNTS" "$run_rc"
exit $EXIT_COUNTS
