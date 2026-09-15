#!/usr/bin/env bash
# FastSuite end-to-end demo: build, then run the whole chain on generated data.
#
#   ./examples/demo.sh            # 20k pairs through fastsim -> fastdedup -> fastprep -> fastmaps
#   ./examples/demo.sh 200000     # bigger
#
# Everything lands in examples/runs/demo/ (deleted first).

set -euo pipefail

cd "$(dirname "$0")/.."
PAIRS="${1:-20000}"
TOOLS="$PWD/tools"

echo "== building (offline) =="
for t in fastdedup fastprep fastsim fastmaps fasthelper; do
    (cd "tools/$t" && cargo build --release --offline 2>&1 | tail -1)
done

echo
echo "== install check =="
"$TOOLS/fasthelper/target/release/fasthelper" doctor --tools-dir "$TOOLS" || true

echo
echo "== config =="
sed -n '1,60p' examples/fastsuite.toml | grep -v '^#' | grep -v '^$'

rm -rf examples/runs

echo
echo "== run ($PAIRS pairs) =="
"$TOOLS/fasthelper/target/release/fasthelper" run -c examples/fastsuite.toml --tools-dir "$TOOLS"

echo
echo "== aggregated metrics =="
cat examples/runs/demo/report.tsv

echo
echo "== artifacts =="
ls -la examples/runs/demo/ examples/runs/demo/*/
echo
echo "open examples/runs/demo/report.html"
