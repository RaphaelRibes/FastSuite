#!/usr/bin/env python3
"""Can quality separate an INJECTED polyG run from a BIOLOGICAL G?

The only lever left after establishing that base identity cannot fix the
overtrim. Construct reads where the artefact run has real low quality and the
biological G has normal quality (and the reverse), and see what fastprep does.

Cases, all with h = 8 artefact G's then biology:
  A) artefact G's at Q20, biological G at Q40   -> can quality split them?
  B) artefact G's at Q40, biological G at Q40   -> no signal (the control)
  C) artefact G's at Q40, biological G at Q20   -> inverted signal
Also: what does fastsim itself do? Its model gives the first head base normal
quality, so the boundary is deliberately weak.
"""
import os, random, struct, subprocess, collections

FP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
D = "/tmp/overtrim_q"
os.makedirs(D, exist_ok=True)
NO_G = b"ACT"
random.seed(7)

H = 8          # artefact run length (>= min_run 6 so the scan fires)
QBIO = 40      # biology quality
QART = 20      # artefact quality


def build(q_art, q_bio_extra_g):
    """biology = one G at q_bio_extra_g, then 60 non-G at QBIO."""
    reads, truth = [], {}
    for i in range(30):
        seq = b"G" * H + b"G" + bytes(random.choice(NO_G) for _ in range(150))
        quals = [q_art] * (H + 1) + [QBIO] * 150
        name = f"@q{i}_{random.randrange(1<<30):010d} 1:N:0:ATCGATCG"
        reads.append((name, seq, bytes(quals)))
        truth[name[1:].split()[0]] = H       # true head run is exactly 8
    return reads, truth


def write_fastq(path, reads):
    with open(path, "wb") as fh:
        for name, seq, qual in reads:
            fh.write(name.encode() + b"\n" + seq + b"\n+\n" + qual + b"\n")


def parse_pxmt(path):
    out = {}
    if not os.path.exists(path):
        return out
    buf = open(path, "rb").read()
    off = 5
    while off + 2 <= len(buf):
        (idl,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + idl].decode("utf-8", "replace"); off += idl
        if off + 22 > len(buf):
            break
        off += 2
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[rid.split()[0]] = t5a + t5p
    return out


def run(label, q_art, q_bio_g):
    reads, truth = build(q_art, q_bio_g)
    write_fastq(f"{D}/r.fastq", reads)
    for f in (f"{D}/o.pxmt", f"{D}/o.fastq.gz"):
        if os.path.exists(f):
            os.remove(f)
    rc = subprocess.run(
        [FP, "-1", f"{D}/r.fastq", "-o", f"{D}/o.fastq.gz", "--metadata",
         f"{D}/o.pxmt", "-q", "--no-dedup", "--no-mean-q", "--no-auto-phred"],
        capture_output=True, text=True)
    if rc.returncode != 0:
        print(f"{label}: rc={rc.returncode} {rc.stderr[-160:]}")
        return
    px = parse_pxmt(f"{D}/o.pxmt")
    d = collections.Counter()
    for rid, h in truth.items():
        if rid in px:
            d[px[rid] - h] += 1
    n = sum(d.values())
    over = sum(v for k, v in d.items() if k > 0)
    print(f"{label}")
    print(f"    artefact Q={q_art}, biological-G Q={q_bio_g}")
    print(f"    t5 - h dist: {dict(sorted(d.items()))}")
    print(f"    5' overtrim: {over}/{n} ({over/n*100:.1f}%)  [expect t5={H+1} if it trims the biological G, {H} if quality splits them]")


print("=== can quality split artefact G from biology G?")
run("A) artefact LOW Q, biology G HIGH Q", 15, 40)
print()
run("B) both same Q (control, no signal)", 30, 30)
print()
run("C) artefact HIGH Q, biology G LOW Q (inverted)", 40, 15)
