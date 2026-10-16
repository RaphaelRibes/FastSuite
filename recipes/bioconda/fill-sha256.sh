#!/usr/bin/env bash
# Fill the sha256 of every recipes/bioconda/*/meta.yaml from the tagged
# release archives. Run this AFTER pushing the tags; before that the URLs 404.
set -euo pipefail
cd "$(dirname "$0")"

fail=0
for d in */; do
  t="${d%/}"
  m="$t/meta.yaml"
  [ -f "$m" ] || continue
  url=$(awk '/^  url: /{print $2; exit}' "$m")
  echo "== $t"
  echo "   $url"
  if ! sha=$(curl -fsSL "$url" | sha256sum | cut -d' ' -f1); then
    echo "   !! could not fetch (tag pushed yet?)" >&2
    fail=1
    continue
  fi
  python3 - "$m" "$sha" <<'PY'
import re, sys
p, sha = sys.argv[1], sys.argv[2]
s = open(p).read()
s2 = re.sub(r"(^\s+sha256: ).*$", r"\g<1>" + sha, s, count=1, flags=re.M)
if s2 == s and sha not in s:
    sys.exit(f"no sha256 line in {p}")
open(p, "w").write(s2)
PY
  echo "   sha256 = $sha"
done
exit $fail
