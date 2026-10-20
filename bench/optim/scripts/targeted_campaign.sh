#!/bin/bash
# Targeted desktop campaign, run to completion while this shell holds the ssh
# session. Serial (--cores 1) so each cell's taskset pin is exclusive.
set -uo pipefail

SSH="ssh -o BatchMode=yes desktop"
CFG=/home/raphael/FastSuite/bench/optim/config-papers.yaml
RES=/mnt/Fast/fastsuite-optim/results-papers

$SSH "bash -c 'export PATH=/home/raphael/FastSuite/.pixi/envs/default/bin:\$PATH; cd /home/raphael/FastSuite/bench/optim; snakemake --unlock --configfile config-papers.yaml >/dev/null 2>&1; rm -rf $RES; mkdir -p $RES; nohup snakemake all --cores 1 --configfile config-papers.yaml --printshellcmds > /mnt/Fast/fastsuite-optim/papers.log 2>&1 & echo launched \$!'"

echo "== launched; waiting for the aggregate to appear"
for i in $(seq 1 110); do
  if $SSH "test -s $RES/results.tsv" 2>/dev/null; then echo "== results.tsv present after $((i*20))s"; break; fi
  if ! $SSH "pgrep -f 'snakemak[e] all' >/dev/null" 2>/dev/null; then
     echo "== snakemake exited without an aggregate"; break
  fi
  sleep 20
done

echo "== tail of the campaign log"
$SSH "tail -6 /mnt/Fast/fastsuite-optim/papers.log"
echo "== summary"
$SSH "python3 -" <<'PY'
import csv, os
p = "/mnt/Fast/fastsuite-optim/results-papers/results.tsv"
if not os.path.exists(p):
    print("NO results.tsv")
    raise SystemExit(0)
rows = list(csv.DictReader(open(p), delimiter="\t"))
print("rows:", len(rows))
for r in sorted(rows, key=lambda r: (r["tool"], r["kind"], r["gz_level"], r["threads"])):
    print(f"{r['tool']:10} {r['kind']:9} reads={r['reads']:>7} j{r['threads']:<2} gz={r['gz_level'] or '-':>3} "
          f"wall={float(r['wall_s']):>8.2f} cpu/wall={r['cpu_over_wall']:>5} rss={int(r['peak_rss_kb'])//1024:>5} MB")
PY
