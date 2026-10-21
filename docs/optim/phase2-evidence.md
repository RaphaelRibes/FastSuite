# Phase 2 — the speed levers, measured

Running evidence log for `docs/plans/2026-10-21_implementation-round2-speed-undertrim.md`
Phase 2. Every entry names the artifact it came from. Claim tiers follow the
usual convention: **[M]** measured with a named artifact, **[E]** estimated,
**[V]** verified against current source or binary.

The reference machine for all speed numbers is the **desktop**
(`cachyos-x8664`, 16 threads, 46 GB RAM, NVMe), per `bench/profile/README.md`.
The laptop is used only for functional verification (tests, byte-identity).

---

## Wave A — free wins

### Task 2.14 — fastmaps `lto = "fat"`, `codegen-units = 1`

**[V] The premise held.** fastmaps' `[profile.release]` carried only
`opt-level = 3`; the other five tools all carried both flags. This was the only
tool left without them.

**Byte-identity [M]** — 2 M-pair corpus (`bench/profile/fastmaps_ab.py`), all
five grid files compared by md5:

```
GRID-IDENTICAL
GRID-MD5-IDENTICAL        (5 files: heatmap.data, heatmap.preview,
                           heatmap.json, heatmap.png, heatmap_3d.html)
```

The `fastmaps.json` report differs **only** in `command` (the binary path) and
`elapsed_s` (a timing field); normalised, the two reports are byte-identical:

```
normalised identical: True
data fields: 971 bytes vs 971 bytes
```

**Binary size [M]:** 1 692 288 B → 1 308 368 B, a **23 % reduction**.

**2 M-pair timing [M]**, 5 reps, wall/cpu `min` of the run set:

| build | wall min | wall median | cpu min | peak RSS |
|---|---:|---:|---:|---:|
| baseline (no lto) | 1.859 s | 1.868 s | 3.610 s | 21.4 MB |
| `lto=fat, cgu=1` | 1.847 s | 1.852 s | 3.591 s | **16.4 MB** |

Wall is **not** meaningfully different (0.6 %, inside run-to-run noise);
**peak RSS falls 23 %** (21.4 → 16.4 MB). This is a memory win and a
code-size win, not a wall-time win, and it is reported as such.

---

### Task 2.13 — fastmaps: cheap hasher, `Vec` for per-file counts

**What changed [V]:**

* the four pixel maps (`counts`, `counts_r1`, `counts_r2`, and the worker's
  local) move from `std::collections::HashMap` (SipHash-1-3) to
  `rustc_hash::FxHashMap` — the suite's existing crate, the same version
  fastcheck pins, so no second hashing dependency was introduced;
* `per_file_counts` moves from `HashMap<u16, u64>` to `Vec<u64>`. The key is a
  dense file id, which is what an array is for.

**Why it is safe [V]:** the keys are `(u16, u16)` pairs derived from read
content, never untrusted network input, so SipHash's DoS resistance buys
nothing here.

**Byte-identity [M]:** `per_file_counts` is serialised by **index order**
(`format.rs` loops `0..files.len()`), not hash order, so the `Vec` is a
drop-in — and an out-of-range id still serialises as `0`, exactly as the old
`.get().unwrap_or(0)` did. Same md5-identical grid as above.

**The merge stays exact [V]:** `spiller.rs` accumulates by index over the
partial's `Vec`, which is exact in any `Done` arrival order (FM-12) and grows
to fit an id beyond the current length rather than dropping it.

**Tests [M]:** fastmaps 96 pass, 0 fail; full `pixi run ci` green at 765.

---

### Task 2.20 — fasthelper: "parse each report once"

**[V] The premise did not hold on current source — no change was made.**

The claim was that `check → aggregate → concatenate → html` re-reads every
report. It does not:

* `cmd_concatenate` calls `collect(&a.inputs, …)` **once** and hands the
  resulting `Vec<ToolReport>` to `write_reports`;
* `write_reports` builds the matrix at most once (`owned_matrix`) and reuses it
  for the TSV, and renders `json`, `md` and `html` from the same slice;
* `runner::run` likewise calls `collect_named` once, then `write_reports`.

Measured rather than assumed, with a parse-dominated corpus
(`bench/scripts/fh_parse_once_probe.py`: 40 reports, 9.3 MB, 40 000 counts
each) — if each format re-parsed, four formats would cost about four times one:

```
1 format  (json)            : 0.048 s
4 formats (json,tsv,md,html): 0.048 s
ratio 4/1                   : 1.00x

VERDICT: parsed ONCE — adding three output formats costs render only.
```

A re-parse would have shown ~4.0x. **Nothing to optimise here**; the probe is
kept as the regression test for the property.
