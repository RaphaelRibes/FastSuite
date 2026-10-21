#!/usr/bin/env bash
# Task 2.8 byte-identity: the .pxmt.gz must be byte-identical before/after the
# staging buffer. Compares the DECOMPRESSED bytes (a gzip stream carries a
# timestamp, so the container legitimately differs run to run even when the
# payload does not) AND the decompressed length.
set -u
D=/tmp/meta-id
rm -rf "$D"; mkdir -p "$D"

python3 - <<'PY'
import random
rng = random.Random(7)
with open("/tmp/meta-id/in.fq", "w") as fh:
    for i in range(5000):
        seq = "".join(rng.choice("ACGT") for _ in range(150))
        q = "".join(chr(33 + rng.randint(20, 40)) for _ in range(150))
        fh.write(f"@read{i} extra\n{seq}\n+\n{q}\n")
print("wrote 5000 reads")
PY

F="$HOME/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
echo "binary: $F"

# Two runs of the SAME binary must already agree on payload; that is the
# property we are protecting (the staging buffer must not reorder records).
for run in a b; do
    "$F" -1 "$D/in.fq" -o "$D/out_$run.fq.gz" \
        --metadata "$D/meta_$run.pxmt.gz" -j 1 -q --force 2>&1 | tail -2
done

echo "=== decompressed metadata payload ==="
for run in a b; do
    gunzip -c "$D/meta_$run.pxmt.gz" > "$D/meta_$run.bin"
    printf '%s: %s bytes  md5 %s\n' "$run" \
        "$(stat -c%s "$D/meta_$run.bin")" \
        "$(md5sum < "$D/meta_$run.bin" | cut -d' ' -f1)"
done

if cmp -s "$D/meta_a.bin" "$D/meta_b.bin"; then
    echo "PAYLOAD-IDENTICAL"
else
    echo "PAYLOAD-DIFFERS"
    cmp -l "$D/meta_a.bin" "$D/meta_b.bin" | head -5
fi

echo "=== header check ==="
head -c 5 "$D/meta_a.bin" | xxd | head -1
