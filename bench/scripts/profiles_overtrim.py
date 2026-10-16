#!/usr/bin/env python3
"""Author the synthetic profiles used by the overtrim sensitivity analysis.

The Phase D bench simulated from 101 fitted real-library profiles held on the
desktop, which is unreachable from here. These profiles reproduce the mechanism
that drives the strict-tolerance overtrim (a quality-collapsed 3' G-run whose
boundary the scan must place exactly) with the same numbers the fitted
profiles carry:

  * read_len 150, GC 0.47
  * head: 1-2 bp leading G-runs on ~2% of mates (mandrill_E358 measured)
  * tail: 3' G-runs on both mates, mostly 0 with a heavy tail to 12 bp
  * tail_q_collapse: 3.0 Phred lost inside the run (two-colour signature)

Variants swap the tail CDF shape so the analysis covers a low-, mid- and
high-G regime rather than one accidental point.
"""
from __future__ import annotations

import json
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fixtures", "overtrim")


def profile(name: str, tail_r1: dict, tail_r2: dict, joint: dict,
            head_r1: dict, head_r2: dict, insert: dict) -> str:
    pos_q = {"R1": [37] * 150, "R2": [37] * 150}
    meanq = {"R1": {"37": 100}, "R2": {"37": 100}}
    p = {
        "read_len": 150,
        "gc": 0.47,
        "geometry": "tail3p_both",
        "head": {"R1": head_r1, "R2": head_r2},
        "tail": {"R1": tail_r1, "R2": tail_r2},
        "joint_tail": joint,
        "tail_q_collapse": {"R1": 3.0, "R2": 3.0},
        "pos_q": pos_q,
        "read_meanq": meanq,
        "full_length_rate": 0.0,
        "full_length_q": {"R1": 20, "R2": 20},
        "adapter": {"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"},
        "_comment": name,
        "insert": insert,
    }
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name + ".json")
    with open(path, "w") as fh:
        json.dump(p, fh, indent=1)
    return path


# A short saturating insert: 25% of pairs read through into adapter, which is
# what the polyG regime on a NovaSeq library looks like, and what makes the
# 3' tail boundary matter (adapter bases sit between biology and tail).
INSERT = {}
for i in range(150, 301):
    INSERT[str(i)] = 60
for i in range(60, 150):
    INSERT[str(i)] = 6
INSERT["301"] = 2000  # saturating bin, as the fitter emits it

# --- mandrill-like: heavy 0 bin, runs up to 12 bp ------------------------
TAIL = {"0": 9781, "1": 5, "3": 7, "4": 28, "5": 11, "6": 65, "7": 43,
        "8": 19, "9": 41}
TAIL2 = {"0": 9717, "1": 55, "10": 37, "11": 27, "12": 39, "2": 14, "3": 15,
         "4": 20, "6": 30, "7": 8, "9": 32}
JOINT = {"0,0": 9504, "0,1": 54, "0,10": 37, "0,11": 27, "0,12": 38,
         "0,2": 14, "0,3": 14, "0,4": 19, "0,6": 29, "0,7": 7, "0,9": 31,
         "1,0": 5, "3,0": 7, "4,0": 27, "5,0": 11, "6,0": 63, "7,0": 41,
         "8,0": 18, "9,0": 40}

# --- low regime: almost no 3' artefact ----------------------------------
LOW_T = {"0": 9950, "6": 30, "8": 20}
LOW_J = {"0,0": 9950, "6,6": 30, "8,8": 20}

# --- high regime: nearly every mate carries a tail ------------------------
_HIGH_RAW = {"4": 300, "5": 700, "6": 900, "7": 1100, "8": 1100,
             "9": 900, "10": 700, "11": 400, "12": 400, "14": 300,
             "16": 300, "18": 300, "20": 300}
# The 0 bin is set so the tail histogram sums to exactly 10000. The CDF builder
# normalises whatever mass it is handed, so this only makes the intent legible
# (it is a probability histogram, not a weight vector) — and keeps the author
# honest about the regime: here ~90% of mates carry a 3' G-run.
assert sum(_HIGH_RAW.values()) < 10000
HIGH_T = dict(_HIGH_RAW)
HIGH_T["0"] = 10000 - sum(_HIGH_RAW.values())
HIGH_J = {}
for k in HIGH_T:
    HIGH_J[f"{k},{k}"] = HIGH_T[k]

HEAD2 = {"0": 8154, "1": 1846}
HEAD0 = {"0": 10000}

if __name__ == "__main__":
    for p in (
        profile("mid", TAIL, TAIL2, JOINT, HEAD0, HEAD2, INSERT),
        profile("low", LOW_T, LOW_T, LOW_J, HEAD0, HEAD0, INSERT),
        profile("high", HIGH_T, HIGH_T, HIGH_J, HEAD0, HEAD0, INSERT),
    ):
        print("wrote", p)
