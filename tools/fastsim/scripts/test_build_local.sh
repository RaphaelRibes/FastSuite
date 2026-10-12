#!/bin/bash
# List a few stats files, fetch 3, test build_profile.py locally.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
rm -rf /tmp/stats_test && mkdir -p /tmp/stats_test
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'ls /scratch/users/ribesr/mandril-polyg-2026/stats/ | head -6'
echo "--- fetching ---"
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'cd /scratch/users/ribesr/mandril-polyg-2026/stats && for f in $(ls | head -3); do echo "--- $f ---"; cat "$f"; done' > /tmp/stats_test/combined.txt
python3 - <<'PYEOF'
import json, re
raw = open('/tmp/stats_test/combined.txt').read()
# split on the "--- filename ---" markers
parts = re.split(r'^--- (.+) ---$', raw, flags=re.M)
# parts: [pre, name1, json1, name2, json2, ...]
for i in range(1, len(parts), 2):
    name, content = parts[i], parts[i+1]
    json.loads(content)  # validate
    open(f'/tmp/stats_test/{name}', 'w').write(content)
print("split ok")
PYEOF
python3 scripts/build_profile.py /tmp/stats_test /tmp/profile_test.json
python3 - <<'PYEOF'
import json
p = json.load(open('/tmp/profile_test.json'))
print("has pos_q_uncond:", "pos_q_uncond" in p)
print("pos_q R1[149]:", p["pos_q"]["R1"][149], " uncond:", p["pos_q_uncond"]["R1"][149])
print("pos_q R1[0]:", p["pos_q"]["R1"][0], " uncond:", p["pos_q_uncond"]["R1"][0])
print("run_q_mean R1:", p["run_q_mean"]["R1"])
print("collapse:", p["tail_q_collapse"])
PYEOF
echo LOCAL_TEST_DONE
