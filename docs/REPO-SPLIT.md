# Splitting FastSuite into one repo per tool and per paper

Status: **executed 2026-10-26.** Every tool and both papers live in their own
public repository; this repo pins them as submodules at the same paths, and
`pixi run ci` is green (GRAND 1009) against the pinned trees.

Decisions taken: dev branches (`main` on the new repos is the tidy initial
import; work happens on `dev`), capitalised names matching `FastDedup`, public,
pinned submodules, and one paper repo per paper mounted at `papers/<tool>/`.

FastDedup already had a `dev` branch (the 1.2.2 line, last commit 2026-05-20),
so this suite's version did **not** overwrite it: it lives on `fastsuite-1.4`
and the submodule is pinned to that branch. `main` and `dev` there are untouched.

## Target shape

```
RaphaelRibes/FastDedup        (exists, 1.2.2-era, root-level layout)
RaphaelRibes/FastPrep         FastCheck   FastMaps   FastSim   FastHelper
RaphaelRibes/FastDedup-Paper  (exists, April version)
RaphaelRibes/FastPrep-Paper

RaphaelRibes/FastSuite        umbrella
  docs/ bench/ scripts/ recipes/ examples/ deploy/ pixi.toml .github/
  tools/<tool>      -> submodule, pinned
  papers/<tool>     -> submodule, pinned
```

`tools/<tool>` keeps its current path so nothing in the umbrella has to move:
`pixi run build|ci`, `scripts/*.sh`, `bench/scripts/*` and fasthelper's
sibling-binary lookup (`<manifest>/../<tool>/target/release/<tool>`) all keep
working unchanged.

## What was measured before proposing it

`git subtree split -P <prefix>` per crate, so each new repo carries real history
rather than a squashed import:

| repo | commits | files | standalone `cargo build --release` + `cargo test` |
|---|---|---|---|
| fastdedup | 40 | 36 | builds, tests pass |
| fastprep | 69 | 49 | builds, tests pass |
| fastcheck | 54 | 30 | builds, tests pass |
| fastmaps | 25 | 40 | builds, tests pass |
| fastsim | 40 | 71 | builds, tests pass |
| fasthelper | 35 | 21 | builds; **33 tests need the sibling binaries** (below) |
| fastprep paper | 17 | 5 | `typst compile paper.typ` → exit 0 |
| fastdedup paper | 7 | 1 | `typst compile paper.typ` → exit 0 |

The papers are self-contained: no `#include`, `#image` or `read()` reaches
outside `paper/`, and both compile from a bare clone.

## What the split breaks, and what was done about it

1. **fasthelper's tests resolve siblings at `<manifest>/../<tool>/target/...`.**
   33 call sites (`fastcheck` 10, `fastprep` 8, `fastdedup` 6, `fastmaps` 3,
   `fastsim` 1). A standalone checkout has no siblings, so those tests panic with
   "X is not built under …". Two ways out, and this needs a decision:
   * keep them in the umbrella (fasthelper's repo holds unit + stub-driven CLI
     tests; the integration suite stays where the six tools meet), or
   * have fasthelper's CI check out and build the other five repos into its
     GitHub workspace root, so `<manifest>/..` contains them.
2. **`examples/fastsuite.toml` was read from `../../examples/`.** Fixed here: the
   crate now ships its own `examples/fastsuite.toml` (with a vendored
   `examples/profiles/test_small.json`, because every path in that file is
   relative to the file), and the test validates **every** copy that exists, each
   in place. They cannot be byte-identical — that is why the test validates
   instead of comparing.
3. **`missing_step_binary_fails_the_run_like_check` depended on ambient `PATH`.**
   It passed in the suite tree (nothing on `PATH`) and failed in a standalone
   checkout where `fastcheck` is installed. Fixed with a `PATH`-scrubbed spawn,
   so "missing" means missing on every machine.
4. **Every crate carries a stale `pixi.toml`/`pixi.lock`** pinning
   `rust = ">=1.85,<1.96"` — the pin this project already measured as costing
   ~31 % of suite speed, and the one that reproduces the link failures
   (`__libc_csu_init`) unless `CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER` is
   set. In a split repo that file becomes the authoritative environment. Each new
   repo must get a fresh `pixi.toml` (rust `>=1.99,<2`, c-compiler, cmake,
   pkg-config, the linker env) or none at all.
5. **The sync script** (`fastsuite-local-artifacts/sync-fastsuite.sh`) rsyncs the
   whole tree; it needs `git submodule update --init --recursive` instead of
   copying `tools/`.

## Execution, once the decisions are made

1. For each crate/paper: push `split/<name>` to `<repo>/main`.
2. In FastSuite: `git rm -r --cached tools papers` (history stays), add
   `.gitmodules`, `git submodule add <url> tools/<tool>`, commit the pin.
3. Replace each per-crate `pixi.toml` with the current pin (finding 4) **inside
   the sub-repo**, so the tool is buildable alone.
4. Per-repo CI: build + that repo's tests. Umbrella CI keeps the six-crate
   manifest, the report-schema gate and the closed loop.
5. Re-verify: umbrella `pixi run ci` at GRAND 1009 with submodules checked out;
   each repo green alone; `git submodule status` clean on both machines.
6. Update `docs/`, the papers' build instructions, and
   `references/` in the `fastsuite` skill (the per-crate `pixi.toml` pitfall
   becomes "the sub-repos pin rust themselves").

## Open decisions

1. `FastDedup` / `FastDedup-Paper` already exist (1.2.2-era, older paper). Reuse
   them (tag the old head `v1.2.2`, then push the new history), freeze them and
   name the new repos differently, or force-push over them?
2. Names for the five new tool repos: `FastPrep`-style (matching `FastDedup`) or
   lowercase (matching the binaries)?
3. Visibility: public (as `FastDedup`), private, or mixed?
4. Containment: pinned submodules (reproducible), or tracking `main`?
5. Papers: `papers/<tool>/` in the umbrella, or left inside each tool repo?

Question 1 is the only one that can destroy something; the other four are cheap
to change later.
