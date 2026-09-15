#!/bin/bash
# Check profile rebuild status; run srun fallback if still pending/stale.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'cd /scratch/users/ribesr/novaseq-sim && squeue -u ribesr -j 9013030 -h -o "%T" 2>/dev/null | head -1; ls -la profiles/mandril_E.json; python3 -c "
import json
p = json.load(open(\"profiles/mandril_E.json\"))
print(\"has pos_q_uncond:\", \"pos_q_uncond\" in p)
if \"pos_q_uncond\" in p:
    print(\"pos_q R1[149]:\", p[\"pos_q\"][\"R1\"][149], \" uncond:\", p[\"pos_q_uncond\"][\"R1\"][149])
"' 2>&1 | grep -v 'post-quantum\|store now\|may need'
echo CHECK_DONE
