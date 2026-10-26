#!/usr/bin/env bash
# Move every submodule pin to the head of its branch, and record it in ONE commit.
#
# The umbrella pins each sub-repo at an exact commit (reproducible: a clone in six
# months gets this code). This script is the one-command way to advance them:
# `git submodule update --remote` fetches each sub-repo's configured branch
# (`.gitmodules: branch = ...`) and checks it out, then the new SHAs are staged
# together so the umbrella's history says exactly when it moved and by how much.
#
# Usage: scripts/sync-submodules.sh [--dry-run]
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

dry=0
[[ "${1:-}" == "--dry-run" ]] && dry=1

echo "current pins:"
git submodule status | sed 's/^ /  /'

git submodule update --remote --init --recursive --quiet 2>&1 | tail -5

# What actually moved?
moved=$(git status --porcelain -- tools papers | grep -c '^ M' || true)
if [[ "$moved" -eq 0 ]]; then
    echo "no pin moved — every submodule already sits on its branch head"
    exit 0
fi

echo
git status --short -- tools papers | sed 's/^/  /'
if [[ "$dry" -eq 1 ]]; then
    echo "--dry-run: not committing. Run without --dry-run to record these pins."
    exit 0
fi

git add tools papers
git commit -q -F - <<'MSG'
submodules: bump every pin to its branch head

One commit per advance, so the umbrella records exactly when each sub-repo moved
and by how much: `git diff <this>^ <this> -- tools/ papers/` lists the old and new
SHAs, and each sub-repo's own history explains the change.

Produced by `scripts/sync-submodules.sh` (`pixi run sync-submodules`).
MSG
echo "committed: $(git log --oneline -1)"
echo "push it with: git push origin main"
