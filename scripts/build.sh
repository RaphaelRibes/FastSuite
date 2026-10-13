#!/usr/bin/env bash
# Build every FastSuite tool in dependency order (fasthelper last).
# Usage: bash scripts/build.sh [--offline]
set -euo pipefail

cd "$(dirname "$0")/.."

OFFLINE=""
if [ "${1:-}" = "--offline" ]; then
    OFFLINE="--offline"
fi

for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper; do
    echo "==> building $t"
    (cd "tools/$t" && cargo build --release $OFFLINE)
done

echo "echo "build OK — run: pixi run doctor"
