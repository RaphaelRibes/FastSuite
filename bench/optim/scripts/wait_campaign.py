import subprocess, sys, time, json

TSV = "/mnt/Fast/fastsuite-optim/results/results.tsv"
deadline = time.time() + 150 * 60  # give the campaign 2.5 h
while time.time() < deadline:
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", "desktop",
                        "test -s " + TSV], capture_output=True)
    if r.returncode == 0:
        break
    # bail out if snakemake is gone without producing it
    alive = subprocess.run(["ssh", "-o", "BatchMode=yes", "desktop",
                            "pgrep -c -f 'snakemak[e] all' || true"],
                           capture_output=True, text=True).stdout.strip()
    if alive in ("", "0"):
        print("SNAKEMAKE GONE and no results.tsv -- campaign died")
        sys.exit(2)
    time.sleep(45)
else:
    print("TIMEOUT waiting for results.tsv")
    sys.exit(3)

# summary
code = r'''
import csv
rows = list(csv.DictReader(open("/mnt/Fast/fastsuite-optim/results/results.tsv"), delimiter="\t"))
print("rows:", len(rows))
bad = [r for r in rows if r["kind"] in ("baseline","new")
       and abs(float(r["cpu_over_wall"]) - 1.0) > 0.15
       and abs(float(r["cpu_over_wall"]) - float(r["threads"])) > 0.15]
print("cells whose cpu/wall is neither ~1.00 nor ~threads:", len(bad))
for t in ("fastcheck","fastprep","fastsim","fastmaps","fasthelper","fastdedup"):
    sub = [r for r in rows if r["tool"] == t and r["threads"] == "1" and r["gz_level"] in ("", "6")]
    if not sub:
        print(f"{t:11} no j1 cell"); continue
    big = max(int(r["reads"]) for r in sub)
    for r in sorted([x for x in sub if int(x["reads"]) == big], key=lambda x: x["kind"]):
        print(f"{t:11} {r['kind']:8} reads={r['reads']:>8} wall={r['wall_s']:>8} "
              f"cpu/wall={r['cpu_over_wall']:>5} rss={int(r['peak_rss_kb'])//1024:>5} MB")
'''
out = subprocess.run(["ssh", "-o", "BatchMode=yes", "desktop", "python3 -"],
                     input=code, capture_output=True, text=True)
print(out.stdout or out.stderr[-800:])
