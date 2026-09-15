#!/usr/bin/env python3
"""Inspect one profiler JSON: full-spike rate, collapse, geometry."""
import json

v = json.load(open("/scratch/users/ribesr/mandril-polyg-2026/stats/E101_seq2.json"))
n = v["r1_stats"]["n_reads"]
t = v["r1_stats"]["tail_len_hist"]
print("R1 tail@149,150:", t.get("149"), t.get("150"))
print("R1 full-rate:", (t.get("149", 0) + t.get("150", 0)) / n)
t2 = v["r2_stats"]["tail_len_hist"]
n2 = v["r2_stats"]["n_reads"]
print("R2 tail@149,150:", t2.get("149"), t2.get("150"), "full-rate:", (t2.get("149", 0) + t2.get("150", 0)) / n2)
h = v["r1_stats"]["head_len_hist"]
print("R1 head@150:", h.get("150"))
print("R1 head mean len:", sum(int(k) * c for k, c in h.items()) / n)
rq = v["r1_stats"]["run_q"]
overall = sum(v["r1_stats"]["qsum"]) / (n * 150)
print("R1 overall mean q:", round(overall, 2), "tail_q:", round(rq["tail_sum"] / rq["tail_n"], 2), "collapse:", round(overall - rq["tail_sum"] / rq["tail_n"], 2))
rq2 = v["r2_stats"]["run_q"]
overall2 = sum(v["r2_stats"]["qsum"]) / (n2 * 150)
print("R2 overall mean q:", round(overall2, 2), "tail_q:", round(rq2["tail_sum"] / rq2["tail_n"], 2), "collapse:", round(overall2 - rq2["tail_sum"] / rq2["tail_n"], 2))
