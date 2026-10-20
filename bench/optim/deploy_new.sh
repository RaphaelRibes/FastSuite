#!/bin/bash
# Build the tree under test on the desktop and freeze its six binaries into the
# bench's new-bin directory, so the measured arm cannot change under the run.
set -euo pipefail

REPO=/home/raphael/FastSuite
NEW=/mnt/Fast/fastsuite-optim/new-bin

cd "$REPO"
export PATH=/home/raphael/FastSuite/.pixi/envs/default/bin:$PATH
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/cc

echo "== building (HEAD $(git rev-parse --short HEAD))"
cargo build --release --offline --manifest-path tools/fastcheck/Cargo.toml
for t in fastdedup fastprep fastsim fastmaps fasthelper; do
  cargo build --release --offline --manifest-path "tools/$t/Cargo.toml"
done

mkdir -p "$NEW"
for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper; do
  cp -f "$REPO/tools/$t/target/release/$t" "$NEW/$t"
done

echo "== frozen arm"
ls -l "$NEW"
for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper; do
  printf '%-10s ' "$t"; "$NEW/$t" --version
done
echo "== head: $(git rev-parse HEAD)"
