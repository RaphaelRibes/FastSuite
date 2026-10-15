#!/usr/bin/env python3
"""Measure fastprep's cut against KNOWN truth on synthetic reads."""
import os, random, struct, subprocess, collections

FP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
D = "/tmp/overtrim_probe"
os.makedirs(D, exist_ok=True)

NO_G = b"ACT"
random.seed(1234)
MODE = os.environ.get("BIO_MODE", "nog")


def make_bio():
    if MODE == "nog":
        return bytes(random.choice(NO_G) for _ in range(80))
    if MODE == "gstart":
        return b"G" + bytes(random.choice(NO_G) for _ in range(79))
    if MODE == "gend":
        return bytes(random.choice(NO_G) for _ in range(79)) + b"G"
    if MODE == "full":
        return bytes(random.choice(b"ACGT") for _ in range(80))
    if MODE == "oneG":
        # biology whose FIRST base is a single G, then no G at all:
        # the run looks exactly one base longer than it is.
        return b"G" + bytes(random.choice(NO_G) for _ in range(60))
    if MODE == "twoG":
        # two G's at the start of biology: the run looks two bases longer
        return b"GG" + bytes(random.choice(NO_G) for _ in range(60))
    raise SystemExit(f"unknown BIO_MODE {MODE}")


def write_fastq(path, reads):
    with open(path, "wb") as fh:
        for name, seq, qual in reads:
            fh.write(name.encode() + b"\n")
            fh.write(seq + b"\n+\n")
            fh.write(qual + b"\n")


def parse_pxmt(path):
    out = {}
    if not os.path.exists(path):
        return out
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        return out
    off = 5
    while off + 2 <= len(buf):
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len].decode("utf-8", "replace"); off += id_len
        if off + 22 > len(buf):
            break
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        (_rl,) = struct.unpack_from("<I", buf, off); off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[rid] = (status, t5a + t5p, t3a + t3p)
    return out


truth, reads = {}, []
for h in (0, 1, 2, 3, 5, 7, 10, 14, 20):
    for t in (0, 1, 2, 3, 5, 7, 10):
        bio = make_bio()
        seq = b"G" * h + bio + b"G" * t
        name = f"@read_{h}_{t}_{random.randrange(1 << 30):010d} 1:N:0:ATCGATCG"
        rid = name[1:].split()[0]
        truth[rid] = (h, t)
        reads.append((name, seq, bytes([40] * len(seq))))

write_fastq(f"{D}/r1.fastq", reads)
for f in (f"{D}/out.pxmt", f"{D}/out_r1.fastq.gz"):
    if os.path.exists(f):
        os.remove(f)
cmd = [FP, "-1", f"{D}/r1.fastq", "-o", f"{D}/out_r1.fastq.gz",
       "--metadata", f"{D}/out.pxmt", "-q", "--no-dedup", "--no-mean-q",
       "--no-auto-phred"]
rc = subprocess.run(cmd, capture_output=True, text=True)
if rc.returncode != 0:
    print("rc", rc.returncode, rc.stderr[-300:])
    raise SystemExit(1)

px = parse_pxmt(f"{D}/out.pxmt")
# PXMT keys carry the whole header line; index by the bare read id
by_id = {}
for key, val in px.items():
    by_id[key.split()[0]] = val

head_d, tail_d = collections.Counter(), collections.Counter()
full_exact = 0
untouched = 0
for rid, (h, t) in truth.items():
    if rid not in by_id:
        untouched += 1
        continue
    status, t5, t3 = by_id[rid]
    if status == 0 and t5 == 0 and t3 == 0:
        untouched += 1
        continue
    head_d[t5 - h] += 1
    tail_d[t3 - t] += 1
    if t5 == h and t3 == t:
        full_exact += 1

n = sum(head_d.values())
print(f"MODE={MODE}  reads={len(truth)} touched={n} untouched={untouched}")
if n:
    print(f"  5' excess: {dict(sorted(head_d.items()))}")
    print(f"  3' excess: {dict(sorted(tail_d.items()))}")
    oh = sum(v for k, v in head_d.items() if k > 0)
    ot = sum(v for k, v in tail_d.items() if k > 0)
    print(f"  5' overtrim {oh}/{n} ({oh/n*100:.1f}%)   3' overtrim {ot}/{n} ({ot/n*100:.1f}%)")
    print(f"  both-exact {full_exact}/{n} ({full_exact/n*100:.1f}%)")
