#!/bin/bash
# Cancel stuck sbatch, rebuild profile via srun (worked before).
set -e
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'scancel 9013030 2>/dev/null; cd /scratch/users/$USER/FastSuite/tools/fastsim && timeout 400 srun --partition=cpu-dedicated --account=dedicated-cpu@cirad-normal --time=00:15:00 --mem=4G --cpus-per-task=1 python3 scripts/build_profile.py /scratch/users/$USER/mandril-polyg-2026/stats profiles/mandril_E.json 2>&1 | tail -5; ls -la profiles/mandril_E.json; python3 -c "
import json
p = json.load(open(\"profiles/mandril_E.json\"))
print(\"has pos_q_uncond:\", \"pos_q_uncond\" in p)
if \"pos_q_uncond\" in p:
    print(\"pos_q R1[149]:\", p[\"pos_q\"][\"R1\"][149], \" uncond:\", p[\"pos_q_uncond\"][\"R1\"][149])
    print(\"pos_q R1[0]:\", p[\"pos_q\"][\"R1\"][0], \" uncond:\", p[\"pos_q_uncond\"][\"R1\"][0])
    print(\"run_q_mean R1:\", p[\"run_q_mean\"][\"R1\"])
"' 2>&1 | grep -v 'post-quantum\|store now\|may need'
echo SRUN_DONE
