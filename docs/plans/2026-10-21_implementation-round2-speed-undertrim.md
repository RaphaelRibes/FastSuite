# FastSuite round-3 implementation plan — speed levers + undertrimming reduction

> **For Hermes:** execute with `subagent-driven-development` (one subagent per
> task, two-stage review: spec compliance, then code quality). Every task that
> touches `tools/*/src` runs under the `karpathy-guidelines` skill and its
> complete diff goes into the task reply.

Date: 2026-10-21. Tree at plan time: `b8719c0` (round-2 audit landed, no `src/`
changes). Baseline test manifest: 740 (143/134/203/93/96/70 + 1 pathology).

**Goal:** ship the 31 round-2 propositions — 23 speed/UX levers across the six
tools and 8 undertrimming fixes (UT-1..UT-5 plus the quality contract) — so
FastSuite beats fastp on CPU-seconds per pair, memory, *and* residual adapter,
with every claim backed by a two-axis measurement.

**Architecture:** six independent tools, no shared crate (D-1). Each tool's
changes land in its own crate and its own test files; nothing is extracted
"for reuse". The undertrimming work is confined to `tools/fastprep`
(`adapter.rs`, `pe_overlap.rs`, `algorithm.rs`) and the quality contract to
`tools/fasthelper` (`adapter_guard.rs`, already present).

**Tech stack:** Rust (pinned `>=1.99,<2` — always build through `pixi run`,
never ambient cargo), clap 4, rayon, flate2, Python 3 for the bench harnesses,
Typst for the report.

---

## 0. Standing constraints (read before task 1)

| # | Constraint | Where it comes from | How it is enforced in this plan |
|---|---|---|---|
| C-1 | **No shared crate, no workspace.** Anything two tools need exists twice, copied per tool. | `docs/DECISIONS.md` D-1 | No task creates a common library. Task 2.11 (fastcheck) and 2.13 (fastmaps) deliberately duplicate the same idea in two crates. |
| C-2 | **Two-axis discipline.** A trimming change is only measured when it is quoted on BOTH axes: residual adapter (undertrim) and false trims on clean reads (overtrim). | `docs/UNDERTRIM-REDUCTION.md` §2 | Every task in Phase 1 has a gate that fails if either axis regresses. A knob quoted on undertrim alone is not a measurement. |
| C-3 | **Cell A is a hard invariant.** On adapter-free real corpora, `--pe-overlap` must not change any trim decision (`pe_trims ≈ 0`). | `docs/OVERTRIM-PAIRED.md` §Cell A; `docs/plans/2026-10-19_004500-pe-overlap-with-without.md` | Phase 1 Task 1.6 re-runs Cell A before any PE behaviour change is allowed to ship. |
| C-4 | **Never ship data.** No real reads in the repo, no fixtures derived from real corpora. | D-11 | Phase 1 harnesses generate synthetic corpora; Cell A runs against the desktop corpus, out of tree. |
| C-5 | **Quote CPU-seconds and busy cores, never the broken "1-core wall" cell.** The harness's `pinned` rows escaped `taskset` and burned 2.44 cores. | `docs/optim/raw/round2/_salvage-fastprep-speed.md`; `references/optim-dossier.md` item 1 | Every speed task's A/B records wall + CPU + busy cores; the retired cell is named as retired in `docs/FASTP-COMPARISON.md` (Task 3.2). |
| C-6 | **No win-rate until `winners.tsv` and `fastp-comparison.tsv` are reconciled.** They disagree on 93/101 datasets. | round-2 review 2026-10-21 | Task 3.3 reconciles them or the claim is withheld. |
| C-7 | **Test counts are a gate.** `scripts/expected-test-counts.txt` must match reality; a SHORT count means cargo skipped a rebuild. | `scripts/ci.sh` | Every task that adds tests updates the manifest in the same commit. |
| C-8 | **Version bumps need the user's OK on the diff.** Commit locally, never push. | user standing rule | Task 4.3 shows the diff and waits. |
| C-9 | **Build through pixi, not ambient cargo.** The laptop's ambient cargo is 1.95 (1.44x slower binary, silent). | `pixi.toml` comment | All build/test commands in this plan are `pixi run …`. |

### The two honesty defects this plan must not reintroduce

1. `docs/STATUS.md` still says *"679 cargo tests (143/122/187/89/77/60)"* while
   `scripts/expected-test-counts.txt` says 740 (143/134/203/93/96/70). The
   status file is **61 tests stale**. Task 4.2 fixes it from the manifest, never
   from memory.
2. `docs/optim/round2-fastprep.md` cites the three-stage pipeline at
   `pipeline.rs:445-501` — that range is `detect_phred_offset`, not the
   parse/trim stage. The real anchors are `pipeline.rs:602-645` (the
   `raw_tx`/`raw_rx` chunk pump) and `pipeline.rs:253-269` (`process_chunk`).
   Task 4.2 corrects the citation.

---

## 1. Proposition map — every round-2 item and where it lands

| ID | Proposition | Tool | Task | Gate class |
|---|---|---|---|---|
| UT-1 | PE boundary as primary evidence | fastprep | 1.1 | behaviour change (flagged) |
| UT-2 | G-overprint retry | fastprep | 1.2 | behaviour change (flagged) |
| UT-3 | Corroborated short tails (SE, 6–13 bp) | fastprep | 1.3 | behaviour change (flagged, high risk) |
| UT-4 | Quality/N-aware mismatch scoring | fastprep | 1.4 | byte-identical on clean data |
| UT-5 | Adapter-only reads emit empty/drop | fastprep | 1.5 | contract fix |
| — | Residual-adapter quality contract | fasthelper | 1.7 | behaviour change |
| F2-1 | Three-stage pipeline (parse\|trim\|commit) | fastprep | 2.10 | byte-identical intent |
| F2-2 | Fast default gzip level | fastprep + fastdedup + fasthelper | 2.4, 2.1 | behaviour change (.gz bytes) |
| F2-3 | Clamp `-j` past saturation; report utilisation | fastprep | 2.9 | byte-identical |
| F2-4 | PE-path allocation trims | fastprep | 2.7 | byte-identical |
| F2-5 | Metadata writer staging | fastprep | 2.8 | byte-identical |
| FC2-1 | Single-file parallelism | fastcheck | 2.11 | byte-identical |
| FC2-2 | 785 MB peak → compact table + fused scan | fastcheck | 2.12 | byte-identical |
| FC2-3 | Profile fitter width from read length | fastcheck | 2.12b | behaviour change |
| FC2-4 | Fixed startup cost on tiny inputs | fastcheck | 2.12b | byte-identical |
| FD2-1 | Emit-order determinism | fastdedup | 2.5 | documented contract |
| FD2-2 | Default `--gz-level 1` | fastdedup | 2.4 | behaviour change (.gz bytes) |
| FD2-3 | Faster deflate backend | fastdedup | 2.6 | behaviour change (.gz bytes) |
| FD2-4 | Explain the `-j 1..4` flatline | fastdedup | 2.6 | diagnosis |
| FD2-5 | 128-bit hash visibility | fastdedup | 2.6b | documented contract |
| FM2-1 | Cheap hasher / `Vec` for per-file counts | fastmaps | 2.13 | byte-identical |
| FM2-2 | `lto = "fat"`, `codegen-units = 1` | fastmaps | 2.14 | byte-identical |
| FM2-3 | Fused per-read passes | fastmaps | 2.15 | byte-identical |
| FM2-4 | 8 B index record | fastmaps | 2.16 | proven-equivalent |
| FS2-1 | Runner-level `--gz-level 1` for intermediates | fasthelper + fastsim | 2.1 | behaviour change (.gz bytes) |
| FS2-2 | Parallel simulation, per-pair RNG | fastsim | 2.2 | behaviour change (GT version bump) |
| FS2-3 | Fuse per-base generation loops | fastsim | 2.3 | byte-identical |
| FS2-4 | `--stats` counters zero-cost when absent | fastsim | 2.3b | byte-identical |
| FH2-1 | Sample-level parallelism in aggregate/render | fasthelper | 2.17 | byte-identical |
| FH2-2 | Reuse the step's own fastcheck pass in the guard | fasthelper | 2.18 | proven-equivalent |
| FH2-3 | Crash-resume cleanliness (tmp + rename) | fasthelper | 2.19 | behaviour change |
| FH2-4 | Parse each report once | fasthelper | 2.20 | byte-identical |
| FH2-5 | Two reproduced output defects | fasthelper | 2.21 | correctness |

---

## 2. Phase 0 — the harnesses must exist before any code moves

Nothing in Phase 1 may start until its measurement exists. The decisive UT-1
experiment was run ad hoc in a scratch script and is **not reproducible from the
repo**; the paired AB corpus does not exist at all. These three tasks build them.

### Task 0.1 — Commit the random-tail-mate guard probe

**Objective:** make the decisive UT-1 experiment reproducible from the repo.

**Why this design:** R1 = fragment + TruSeq[:t]; R2 = revcomp(fragment) + a
**purely random** t-base tail (never adapter). Any exact cut on R2 can only come
from the PE-overlap path, so R2 is a pure probe of detector + guard.

**Files:**
- Create: `bench/scripts/undertrim_probe_pe_guard.py`
- Reference: `bench/scripts/undertrim_probe_pe.py` (paired scaffolding, `write_pair`, `run_paired`)

**Step 1: write the script**

```python
#!/usr/bin/env python3
"""The decisive UT-1 guard A/B: does the PE detector's own floor discard the
boundary it just pinned?

R1 = fragment + TruSeq[:t]; R2 = revcomp(fragment) + a RANDOM t-base tail
(never adapter). Any exact cut on R2 can only come from the PE-overlap path,
because R2 carries no adapter at all.

Arms (t, --min-adapter-overlap):
    A  t=20, floor 14  -> both mates cut exactly   (baseline: the guard passes)
    B  t=10, floor 14  -> 0/200                    (the floor blocks the guard)
    C  t=10, floor 10  -> both mates cut exactly   (same geometry, lower floor)

Usage: python3 bench/scripts/undertrim_probe_pe_guard.py [--n 200] [--out JSON]
"""
```

Arms A/B/C exactly as above; `--n` records per arm (default 200); write JSON to
`bench/results/undertrim-probe-pe-guard.json` with the schema
`{"arm": {"t":, "floor":, "r1_exact":, "r2_exact":, "n":}}`.

**Step 2: run it against the current binary and confirm it reproduces the
published numbers**

Run: `python3 bench/scripts/undertrim_probe_pe_guard.py --n 200`
Expected (current `main`, floor 14): A → 200/200 and 200/200; B → 0/200 and
0/200; C → 200/200 and 200/200.
**If B is not 0/200, stop** — the mechanism this plan is built on has moved and
the whole Phase 1 design must be re-derived first.

**Step 3: commit**

```bash
git add bench/scripts/undertrim_probe_pe_guard.py bench/results/undertrim-probe-pe-guard.json
git commit -m "bench: the random-tail-mate guard probe, reproducible from the repo"
```

---

### Task 0.2 — Paired mode for the two-axis A/B

**Objective:** give the AB corpus a paired mode, so UT-1's claim can be priced on
real PE geometry instead of argued.

**Files:**
- Modify: `bench/scripts/undertrim_ab.py` (add `--paired`, `build_corpus_paired`)
- Test: none (Python harness; the numbers are the test)

**Step 1: implement `--paired`**

Add a corpus builder that emits mate pairs with a known insert `I < L`, so each
mate carries `L - I` bases of its own adapter readthrough, and scores **both**
mates against the exact truth:

```
exact   cut == t
under   cut < t     (residual adapter left in the output read)
over    cut > t     (real insert removed)
```

The clean axis must be measured on the same run: a population with `I >= L`
(no readthrough at all) whose false-trim rate is the counterweight. Reuse
`TRUSEQ_R1`, `TRUSEQ_R2`, `dna`, `write_fastq`, `parse_pxmt` from
`undertrim_probe.py` (already imported by `undertrim_ab.py`).

**Step 2: verify the harness reproduces the published single-end numbers**

Run: `python3 bench/scripts/undertrim_ab.py --reads 200000 --out bench/results/undertrim-ab.json`
Expected: shipped config (floor 14, eps 0.2) → residual **25 105 reads =
41.95 % of carriers**, overtrim **20 reads = 0.007 %** of clean.
If these numbers move, the corpus builder changed — diff it before continuing.

**Step 3: run the paired mode and record the baseline**

Run: `python3 bench/scripts/undertrim_ab.py --paired --reads 200000 --out bench/results/undertrim-ab-paired.json`
Expected: the paired residual is **worse** than the single-end residual for the
sub-floor bands (the PE path is currently a dead end below floor 14), and the
clean-axis false-trim rate is ≈ the single-end rate. Record both numbers in
`docs/UNDERTRIM-REDUCTION.md` as the "before" row.

**Step 4: commit**

```bash
git add bench/scripts/undertrim_ab.py bench/results/undertrim-ab-paired.json docs/UNDERTRIM-REDUCTION.md
git commit -m "bench: paired mode for the two-axis adapter A/B"
```

---

### Task 0.3 — One driver that emits the comparison table

**Objective:** make "both axes, always" a single command, so no future trimming
change can be quoted on one axis.

**Files:**
- Create: `bench/scripts/undertrim_report.py`
- Read: `bench/scripts/undertrim_{probe,probe_pe,probe_pe_guard,ab}.py`

**Step 1: implement**

Runs, in order: `undertrim_probe.py` (arms A–D), `undertrim_probe_pe.py`,
`undertrim_probe_pe_guard.py`, `undertrim_ab.py` (SE), `undertrim_ab.py --paired`.
Emits a markdown table with one row per metric and a **PASS/FAIL** against
thresholds declared at the top of the file:

```
RESIDUAL_SE_MAX   = 0.25    # residual/carrier must not exceed 25 %
FALSE_CLEAN_MAX   = 0.0002  # false trims on clean reads must stay <= 0.02 %
PE_EXACT_MIN      = 0.90    # tail10/tail6 exact rate on the PE probe
GUARD_B_MUST_BE   = 0       # arm B (floor 14, t=10) must cut nothing
```

Exit non-zero if any threshold fails — this is the command that gates Phase 1.

**Step 2: run it against current `main` and confirm it passes**

Run: `python3 bench/scripts/undertrim_report.py`
Expected: PASS on every row (the current binary is the "before" state), with
residual 41.95 %, false 0.007 %, guard-B 0.
**A FAIL here means a threshold is wrong, not that the binary regressed** — fix
the threshold and say why in the commit message.

**Step 3: commit**

```bash
git add bench/scripts/undertrim_report.py
git commit -m "bench: one driver for the two-axis undertrim report"
```

---

## 3. Phase 1 — undertrimming reduction (the flagship)

Order matters: UT-1 and UT-2 are independent and each is gated on its own; UT-3
is the risky one and lands last of the trimming changes, behind its own
corroboration evidence.

Every task in this phase follows the same shape:

1. implement behind a **new flag**, default = today's behaviour;
2. measure with `bench/scripts/undertrim_report.py` — both axes;
3. only then, in Task 1.6, decide whether the flag's default flips.

That staging is deliberate: it means a failed hypothesis costs a flag, not a
release.

### Task 1.1 — UT-1: the PE boundary as primary evidence

**Objective:** stop the PE-overlap guard from discarding a boundary the detector
already pinned exactly.

**The defect, in one line:** `readthrough_is_adapter` re-verifies the tail by
calling `match_adapters(tail, set)`, which re-imposes `min_overlap_se = 14` on
a tail of length `t` — so a 10-base tail is *provably* unmatchable there, and
the detector's perfect boundary is thrown away.

**Evidence:** `bench/results/undertrim-probe-pe-guard.json` — arm B (t=10,
floor 14) → 0/200; arm C (same geometry, floor 10) → 200/200.

**Files:**
- Modify: `tools/fastprep/src/adapter.rs` (add a boundary-anchored prefix check)
- Modify: `tools/fastprep/src/pe_overlap.rs:275-288` (`readthrough_is_adapter`)
- Modify: `tools/fastprep/src/cli.rs` (new flag, `--adapter-evidence`)
- Test: `tools/fastprep/tests/pe_overlap.rs` (integration) + `tools/fastprep/src/pe_overlap.rs` unit block

**Step 1: write the failing test**

In `tools/fastprep/src/pe_overlap.rs`'s `#[cfg(test)] mod tests`:

```rust
#[test]
fn a_pinned_boundary_survives_a_tail_below_the_overlap_floor() {
    // t = 10: below the default floor of 14, but the mate geometry pinned
    // the boundary exactly. The guard must not re-impose the floor on it.
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let set = AdapterSet::new(
        vec![Adapter { name: "TruSeq/R1".into(), seq: adapter.to_vec() }],
        0.2,
        14, // the floor the guard must NOT re-impose
    );
    let mut seq = b"ACGT".repeat(30);
    seq.extend_from_slice(&adapter[..10]); // a 10-base readthrough
    assert!(
        readthrough_is_adapter(&seq, 10, &set),
        "a boundary the detector pinned must survive the tail check"
    );
}
```

Run: `cd tools/fastprep && cargo test --offline pinned_boundary_survives`
Expected: FAIL — "a pinned boundary must survive the tail check".

**Step 2: implement the prefix check in `adapter.rs`**

Add next to `match_adapters` (which stays untouched — it is the general
3'-anchored matcher and its semantics are load-bearing elsewhere):

```rust
/// Boundary-anchored adapter-prefix check.
///
/// The tail is the last `t` bases of a read and the adapter's own prefix is
/// compared against it: this is exactly the 3'-anchored alignment at
/// `p = len - t` with `k = t`, WITHOUT the `min_overlap_se` floor that
/// `match_adapters` applies.
///
/// The floor exists to keep a short 3' prefix from matching by chance. This
/// function is only called where independent evidence already fixed the
/// boundary (mate geometry), so the floor is the wrong guard here: it
/// discards the very boundary the detector pinned.
///
/// Mismatch budget: `floor(t * error_rate)`, the same rule `match_adapters`
/// uses, so the two paths never disagree about a tail of the same length.
pub fn prefix_matches_adapter(
    tail: &[u8],
    adapter: &[u8],
    error_rate: f64,
) -> bool {
    let t = tail.len();
    if t == 0 || t > adapter.len() {
        return false;
    }
    let max_mm = max_mismatches(t, error_rate);
    mismatches_upto(tail, &adapter[..t], max_mm) <= max_mm
}

/// True if ANY adapter in the set matches the tail as a boundary-anchored
/// prefix.
pub fn any_prefix_matches(tail: &[u8], set: &AdapterSet) -> bool {
    set.adapters
        .iter()
        .any(|a| prefix_matches_adapter(tail, &a.seq, set.error_rate))
}
```

**Step 3: rewire `readthrough_is_adapter` (the guard's two rules, both kept)**

The guard has two independent rules today. Rule 2 (the biology protection) is
sound and **stays exactly as it is**: the adapter must *start at* the boundary
(matched position == 0). Rule 1 (the floor) is the bug and changes.

```rust
pub fn readthrough_is_adapter(
    seq: &[u8],
    readthrough: usize,
    set: &crate::adapter::AdapterSet,
) -> bool {
    if readthrough == 0 || readthrough > seq.len() {
        return false;
    }
    let tail = &seq[seq.len() - readthrough..];
    // PRIMARY: boundary-anchored prefix check. The caller already pinned this
    // boundary from mate geometry; the tail only has to be adapter-LIKE at
    // the boundary, not long enough to pass a general overlap floor.
    if crate::adapter::any_prefix_matches(tail, set) {
        return true;
    }
    // FALLBACK: the general matcher, unchanged, so every pre-existing
    // behaviour stays reachable (it also covers tails longer than any
    // adapter, where the prefix check cannot apply).
    match crate::adapter::match_adapters(tail, set) {
        Some(m) => m.t3 == readthrough,
        None => false,
    }
}
```

Note the FALLBACK keeps `m.t3 == readthrough` — the "starts at the boundary"
rule. The PRIMARY path must apply the same rule: the tail IS the last
`readthrough` bases, so a prefix match at the boundary is the same condition.

**Step 4: add the flag and its dual-mate corroboration**

In `tools/fastprep/src/cli.rs`, add:

```rust
/// How much independent evidence the adapter trim requires.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum AdapterEvidence {
    /// Today's behaviour: a 3'-anchored match of at least
    /// `--min-adapter-overlap` bases.
    Floor,
    /// Accept a shorter tail when the mate geometry pinned the boundary
    /// (paired-end only). Tails of 5 bases or fewer additionally require
    /// both mates to agree.
    Boundary,
}

#[arg(long, value_enum, default_value_t = AdapterEvidence::Floor,
      help_heading = "Adapter Trimming")]
pub adapter_evidence: AdapterEvidence,
```

The `Floor`/`Boundary` split needs to reach `readthrough_is_adapter`, which
today takes only the set. Add the policy as a parameter rather than a global:
`readthrough_is_adapter(seq, readthrough, set, evidence)`. `Floor` keeps today's
code path byte-for-byte.

**Dual-mate corroboration for t <= 5:** the caller in `pe_overlap.rs`
(`detect_pe_overlap_trim`) knows both mates; when the readthrough is ≤ 5
bases, require that the *other* mate's tail also matches the adapter prefix
(TruSeq/R1 for one mate, TruSeq/R2 for the other). At t ≤ 5 a random tail hits
a 5-base prefix by chance with probability ≈ 4⁻⁵ ≈ 1e-3 per mate; requiring
both mates drops that to ≈ 1e-6 per pair.

**Step 5: unit tests**

- `a_pinned_boundary_survives_a_tail_below_the_overlap_floor` (above) — PASS.
- `the_floor_policy_still_refuses_a_short_tail` — with `AdapterEvidence::Floor`,
  the same 10-base tail is refused. This is the byte-identical property.
- `a_tail_that_is_not_the_adapter_prefix_is_refused` — random 10-base tail,
  `Boundary` policy, refused.
- `a_random_tail_below_five_needs_both_mates_to_agree` — one mate adapter-like,
  the other random ⇒ refused.

Run: `cd tools/fastprep && cargo test --offline` — Expected: all pass, count
rises by 4 (203 → 207).

**Step 6: integration test in `tools/fastprep/tests/pe_overlap.rs`**

End-to-end: build a mate pair with insert `I = L - 10`, R1 carrying TruSeq[:10],
R2 carrying a random 10-base tail; run with `--adapter-evidence boundary` and
assert R1 is cut exactly and R2 is cut at the pinned boundary. Then run the same
input with `--adapter-evidence floor` and assert the output is **byte-identical
to today's binary** (the flag defaults to `Floor`).

**Step 7: measure both axes**

Run: `python3 bench/scripts/undertrim_report.py`
Expected: PE-probe exact rate for tail6/tail10 goes **0/500 → ~500/500**;
the SE residual does not move (UT-1 is PE-only); the clean-axis false-trim
rate stays at 0.007 %.
If the clean axis moves, the corroboration rule is not tight enough — tighten
and re-measure. Do not proceed with a regression.

**Step 8: commit**

```bash
git add tools/fastprep/src/adapter.rs tools/fastprep/src/pe_overlap.rs \
        tools/fastprep/src/cli.rs tools/fastprep/tests/pe_overlap.rs \
        scripts/expected-test-counts.txt
git commit -m "feat(fastprep): accept the PE boundary as primary evidence (UT-1)"
```

---

### Task 1.2 — UT-2: G-overprint retry

**Objective:** recover the adapter that a terminal no-signal G-run hides.

**The mechanism:** `match_adapters` is 3'-anchored — the window `k = min(a, l-p)`
must reach the read's terminus — so a terminal `G×6`/`G×8` sits **inside** the
window and consumes the mismatch budget `floor(k·eps) = 4` at k=20. Measured
(`bench/results/undertrim-probe.json`, Arm C): g=0/2/4 → 2 000/2 000 exact;
g=6/8 → **0 exact / 2 000 untrimmed**; g=10/12 → partial.

**Files:**
- Modify: `tools/fastprep/src/algorithm.rs:1304-1320` (the adapter step of `trim`)
- Modify: `tools/fastprep/src/adapter.rs` (add a G-run-aware retry helper)
- Test: `tools/fastprep/tests/undertrim.rs` (new file)

**Step 1: write the failing test**

New file `tools/fastprep/tests/undertrim.rs`:

```rust
//! Undertrimming regressions: reads whose adapter is hidden by a terminal
//! no-signal G-run, and the two-axis expectations that keep the fix honest.

#[test]
fn a_terminal_g_run_does_not_hide_the_adapter() {
    // 20-base TruSeq readthrough whose last 6 bases are G (NovaSeq no-signal).
    // Today: 0/2000 exact. After UT-2: cut exactly at the adapter start.
    // ...
}
```

Use the same scaffolding as `tools/fastprep/tests/regressions.rs` (`tmpdir`,
`run`, `records` helpers) — copy them, do not import across test binaries.

**Step 2: implement**

In `algorithm.rs`, the adapter step becomes:

```rust
    // --- Step 1: adapter detection ---
    let (t5_adapter, t3_adapter, adapter) = match params.pe_adapter {
        Some((t5, t3)) => (t5, t3, Some(crate::adapter::PE_OVERLAP_ADAPTER)),
        None => match params.adapter.as_deref() {
            None => (0, 0, None),
            Some(set) => {
                match crate::adapter::match_adapters(bases, set) {
                    Some(m) => (m.t5, m.t3, Some(m.adapter)),
                    None => match crate::adapter::match_adapters_ignoring_terminal_run(
                        bases, set, crate::adapter::G_RUN_MIN,
                    ) {
                        Some(m) => (m.t5, m.t3, Some(m.adapter)),
                        None => (0, 0, None),
                    },
                }
            }
        },
    };
```

And in `adapter.rs`:

```rust
/// Length of the terminal homopolymer run that a no-signal overprint can
/// produce at a 3' end. Six is the measured boundary: g=4 still matches
/// (the mismatch budget absorbs it), g=6 does not.
pub const G_RUN_MIN: usize = 6;

/// Retry the 3'-anchored match on the slice BEFORE a terminal G-run.
///
/// NovaSeq's two-colour no-signal overprint replaces the last bases of a read
/// with G. Because the matcher is 3'-anchored, those G's sit INSIDE the
/// alignment window and spend the mismatch budget, so an adapter that is
/// really there is not found.
///
/// The retry matches the pre-G slice only, so the full `k >= min_overlap`
/// bar still applies to the adapter bases themselves; the G-run is consumed
/// with the adapter (it is removed by the polyX stage anyway).
pub fn match_adapters_ignoring_terminal_run(
    seq: &[u8],
    set: &AdapterSet,
    run_min: usize,
) -> Option<AdapterMatch> {
    let run = terminal_run_len(seq, b'G');
    if run < run_min {
        return None;
    }
    let cut = seq.len() - run;
    let m = match_adapters(&seq[..cut], set)?;
    // The match must END at the pre-G terminus: the adapter occupies the
    // window right up to where the G-run starts. A match that stops earlier
    // means the G-run hides something else.
    if m.t5 + m.t3 != cut {
        return None;
    }
    Some(AdapterMatch { adapter: m.adapter, t5: m.t5, t3: m.t3 + run })
}
```

`terminal_run_len(seq, b'G')` counts the trailing run of `G` (case-insensitive).

**Step 3: verify the unit test passes and the existing suites do not move**

Run: `cd tools/fastprep && cargo test --offline`
Expected: new test passes; **the 203 existing tests are unchanged in outcome**.

**Step 4: measure both axes**

Run: `python3 bench/scripts/undertrim_report.py`
Expected: Arm C g=6/8 goes **0 → ~100 % exact**; the AB corpus' `g under`
column (3 673 reads, 14.6 % of residual) drops to ≈ 0; the clean-axis
false-trim rate stays at 0.007 % (the G-run is not adapter, and a clean read
with a terminal G-run has no adapter prefix to match).

**Step 5: commit**

```bash
git add tools/fastprep/src/algorithm.rs tools/fastprep/src/adapter.rs \
        tools/fastprep/tests/undertrim.rs scripts/expected-test-counts.txt
git commit -m "feat(fastprep): retry the adapter match before a terminal G-run (UT-2)"
```

---

### Task 1.3 — UT-3: corroborated short tails (single-end, 6–13 bases)

**Objective:** recover the 58.5 % of residual that sits in bands 6–13 on
single-end data, **without** the false-trim explosion of a lower floor.

**Why this is the risky one:** the uncorroborated version of this idea is
exactly the mo=6 arm — 1.29 % false trims on clean reads, a 180× regression
on the overtrim axis. The corroboration is the whole point.

**Corroboration rule (two independent signals, either sufficient):**
1. **Quality drop at the boundary**: the mean quality of the 4 bases before
   the candidate boundary is at least `Q_DROP` (default 8) below the mean
   quality of the 10 bases before that. A real adapter readthrough is where
   the sequencer starts hitting the adapter; biology does not produce a sharp
   quality cliff at an arbitrary position.
2. **Dataset-level 3'-prefix enrichment**: the adapter prefix occurs at the 3'
   end of more reads than chance predicts. The auto-detection sample
   (`cli.rs` `det.sampled`, ~4 000 reads) already exists; compute the prefix
   frequency there and require it to exceed `4^-t` by a margin.

**Files:**
- Modify: `tools/fastprep/src/algorithm.rs` (the adapter step; needs qualities)
- Modify: `tools/fastprep/src/adapter.rs` (prefix candidate enumeration)
- Modify: `tools/fastprep/src/cli.rs` (flags: `--adapter-short-tail`,
  `--adapter-short-tail-min`, `--adapter-q-drop`)
- Test: `tools/fastprep/tests/undertrim.rs`

**Step 1: enumerate candidates without changing behaviour**

Add to `adapter.rs`:

```rust
/// Enumerate 3'-terminal adapter-prefix candidates of length
/// `min_len..=max_len` for every adapter in the set, returning
/// (adapter_index, tail_len) pairs that match EXACTLY (0 mismatches).
///
/// Exactness is the point: a 6-base prefix matches by chance with probability
/// 4^-6 = 2.4e-4, and the corroboration rules decide the rest.
pub fn exact_prefix_candidates(
    seq: &[u8],
    set: &AdapterSet,
    min_len: usize,
    max_len: usize,
) -> Vec<(u8, usize)> {
    let mut out = Vec::new();
    let l = seq.len();
    for (ai, a) in set.adapters.iter().enumerate() {
        let hi = max_len.min(a.seq.len()).min(l);
        for t in min_len..=hi {
            if seq[l - t..] == a.seq[..t] {
                out.push((ai as u8, t));
            }
        }
    }
    out
}
```

**Step 2: add the quality-drop test**

```rust
#[test]
fn a_short_exact_prefix_without_corroboration_is_not_trimmed() {
    // Clean read whose last 6 bases happen to equal an adapter prefix, with
    // FLAT quality. Today: not trimmed (the floor refuses it). After UT-3
    // with no corroboration: still not trimmed.
}
```

and the positive case: same read but with a quality cliff at the boundary →
trimmed exactly.

**Step 3: implement the corroborated path in `algorithm.rs`**

Only when `params.short_tail.is_some()` (the flag defaults to `None` = off):

```rust
// A short exact 3'-terminal adapter prefix, accepted only with independent
// corroboration (a quality cliff at the boundary, or dataset-level prefix
// enrichment). Without corroboration a 6-base prefix is indistinguishable
// from chance: the measured mo=6 arm costs 1.29 % false trims on clean
// reads (docs/UNDERTRIM-REDUCTION.md §2).
```

The quality-drop helper:

```rust
/// True when the quality cliff at `boundary` is large enough to corroborate
/// an adapter start: the mean quality of the 4 bases before the boundary is
/// at least `drop` below the mean of the 10 bases before that.
fn quality_cliff_at(quals: &[u8], boundary: usize, drop: u8) -> bool {
    if boundary < 14 || quals.len() < boundary {
        return false;
    }
    let near = &quals[boundary - 4..boundary];
    let far = &quals[boundary.saturating_sub(14)..boundary - 4];
    if far.is_empty() {
        return false;
    }
    let mean = |v: &[u8]| v.iter().map(|&q| q as u32).sum::<u32>() / v.len() as u32;
    far_mean >= near_mean + drop as u32
}
```

**Step 4: measure both axes, hard**

Run: `python3 bench/scripts/undertrim_report.py`
Expected: SE residual drops by a large fraction of the 6–13 bands; clean-axis
false-trim rate stays **≤ 0.02 %** (the threshold in `undertrim_report.py`).

**If the clean axis exceeds 0.02 %, the fix is rejected.** Do not tune the
threshold to make it pass; tighten the corroboration (raise `Q_DROP`, require
both signals) or drop UT-3 and say so in the round-3 report.

**Step 5: commit**

```bash
git add tools/fastprep/src/algorithm.rs tools/fastprep/src/adapter.rs \
        tools/fastprep/src/cli.rs tools/fastprep/tests/undertrim.rs \
        scripts/expected-test-counts.txt
git commit -m "feat(fastprep): corroborated short adapter tails (UT-3, flag-gated)"
```

---

### Task 1.4 — UT-4: quality- and N-aware mismatch scoring

**Objective:** stop an ambiguous base or a low-quality base from consuming the
whole mismatch budget.

**The defect:** `mismatches_upto` (adapter.rs:700) counts every differing byte
as one mismatch, so an `N` (a wildcard) costs the same as a real substitution.
The AB corpus shows 1 185 n-class + 851 sub-class residual reads (~5–8 % of
residual). Separately, `floor(k·eps)` at eps=0.2 allows only 2 mismatches at
k=14–19, so a tail with 3 sequencing errors is missed (Arm B: t=14, e=3 →
27.2 % exact).

**Files:**
- Modify: `tools/fastprep/src/adapter.rs` (add `scored_mismatches`)
- Modify: `tools/fastprep/src/algorithm.rs` (use the scored variant on the
  anchored path; the general matcher keeps the raw count unless the
  measurement says otherwise)
- Test: `tools/fastprep/src/adapter.rs` unit block + `tools/fastprep/tests/undertrim.rs`

**Step 1: implement the scored comparison**

```rust
/// Mismatch count where an ambiguous base costs HALF a mismatch and a
/// low-quality base costs half too.
///
/// N is not a disagreement, it is an absence of information: charging it a
/// full mismatch spends the budget on a base the sequencer never called.
/// A base whose quality is at or below `LOW_Q` is called, but badly, so it
/// costs half.
///
/// The wildcards are only allowed when the rest of the window agrees: a
/// window of ambiguous bases is not evidence, so an `N` counts as a FULL
/// mismatch unless at least `MIN_UNAMBIGUOUS` bases are real calls.
pub const LOW_Q: u8 = 2;
pub const MIN_UNAMBIGUOUS: usize = 12;

pub fn scored_mismatches(
    a: &[u8],
    b: &[u8],
    quals: Option<&[u8]>,
    low_q: u8,
) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut half = 0usize;
    let mut unambig = 0usize;
    for i in 0..a.len() {
        if a[i] == b[i] {
            unambig += 1;
            continue;
        }
        let ambig = a[i] == b'N' || b[i] == b'N';
        let low = quals
            .map(|q| q.get(i).copied().unwrap_or(u8::MAX) <= low_q)
            .unwrap_or(false);
        if ambig || low {
            half += 1;
        }
    }
    if unambig < MIN_UNAMBIGUOUS {
        // Not enough real calls: an ambiguous window is not evidence.
        return a.len() as f64;
    }
    half as f64 * 0.5
}
```

**Step 2: use it on the anchored path only, first**

Wire it into `prefix_matches_adapter` (Task 1.1's new function) as a second
variant, behind the same `--adapter-evidence` flag, so the general matcher —
whose semantics other code depends on — is untouched until measured.

**Step 3: measure**

Run: `python3 bench/scripts/undertrim_report.py`
Expected: the n-class and sub-class residual columns shrink; the clean-axis
rate does not move (a clean read has no adapter prefix to match).

**Step 4: commit**

```bash
git add tools/fastprep/src/adapter.rs tools/fastprep/src/algorithm.rs \
        tools/fastprep/tests/undertrim.rs scripts/expected-test-counts.txt
git commit -m "feat(fastprep): quality- and N-aware mismatch scoring (UT-4)"
```

---

### Task 1.5 — UT-5: adapter-only reads must not come out partial

**Objective:** fix a visible contract bug: a read that is 100 % adapter is
emitted **partial** (cut somewhere) instead of dropped or emitted empty.

**Evidence:** `undertrim_probe.py` Arm C0 (insert 0, t=150) → 200/200 reads
emitted partial.

**Files:**
- Modify: `tools/fastprep/src/algorithm.rs` (the `TrimResult` decision when
  the insert length reaches 0)
- Test: `tools/fastprep/tests/undertrim.rs`

**Step 1: write the failing test** — a read that is entirely adapter must be
dropped (or emitted with an empty sequence) per the `--min-len` semantics, and
must be **counted** in the drop counters, not emitted as a 1-base record.

**Step 2: implement** — when `insert_len == 0` after the adapter step, return
`TrimResult::Drop { reason: DropReason::Adapter }` (the reason enum already
exists — `stats.rs:170` counts `reads_dropped_adapter`).

**Step 3: verify + commit**

```bash
git add tools/fastprep/src/algorithm.rs tools/fastprep/tests/undertrim.rs \
        scripts/expected-test-counts.txt
git commit -m "fix(fastprep): a fully adapter read is dropped, not emitted partial (UT-5)"
```

---

### Task 1.6 — The measurement wave and the default decision

**Objective:** decide, with numbers, whether UT-1..UT-4 become defaults — and
re-run the hard invariants.

**This is the task that earns the "more appealing than fastp" claim. It may
end with "no default flips" — that is a legitimate outcome.**

**Step 1: re-run Cell A (hard invariant, C-3)**

The desktop holds the corpus (`/mnt/Fast/fastsuite-scratch`, reachable over
Tailscale; **not** mounted on the laptop). On the desktop:

```bash
cd ~/FastSuite && pixi run build
# the Cell A driver from bench/work/cellB/ against the desktop corpus
python3 bench/work/cellB/reconcile_pe_trims.py <gt.tsv> <R1.fastq.gz> <R2.fastq.gz> <outdir>
```

Expected: `pe_trims ≈ 0` on the adapter-free corpus, i.e. **no trim decision
changes** relative to the `fastprep-hd0` baseline. The 102-dataset corpus has
no adapter readthrough at all, so this measures the false-fire rate only.

**A non-zero `pe_trims` on this corpus is a ship blocker.** Tighten the
corroboration (raise the minimum corroborated tail length, require both mates)
and re-measure. Record the final number in `docs/OVERTRIM-PAIRED.md`.

**Step 2: run the full two-axis suite, before and after**

```bash
python3 bench/scripts/undertrim_report.py
```

Fill this table into `docs/UNDERTRIM-REDUCTION.md` §5 (the "after" row):

| metric | before | after | threshold |
|---|---|---|---|
| SE residual / carriers | 41.95 % | ? | ≤ 25 % |
| PE residual / carriers | ? (Task 0.2) | ? | ≤ 25 % |
| false trims / clean reads | 0.007 % | ? | ≤ 0.02 % |
| PE probe tail10 exact | 0/500 | ? | ≥ 90 % |
| Cell A `pe_trims` | 0 | ? | = 0 |

**Step 3: decide the defaults, explicitly**

For each of UT-1, UT-2, UT-3, UT-4 write one of:
- **flip the default** — the measured table justifies it; or
- **keep it opt-in** — with the measured reason (e.g. "UT-3 buys 9 pp of
  residual for 0.015 % of false trims; that is inside the noise band of the
  clean corpus, so it stays opt-in").

A flip is a behaviour change ⇒ a minor version bump (fastprep 0.3.0 → 0.4.0,
D-2) and a `docs/RELEASES.md` entry.

**Step 4: commit the measurement wave**

```bash
git add docs/UNDERTRIM-REDUCTION.md bench/results/undertrim-*.json scripts/expected-test-counts.txt
git commit -m "docs: the round-3 undertrim measurement wave, both axes"
```

---

### Task 1.7 — The quality contract fastp does not ship

**Objective:** make "adapter survived trimming" a **machine-checkable** failure
that any sceptical user can run in one command.

**This already exists** — `fasthelper run --fail-on-adapter-residual`
(`tools/fasthelper/src/cli.rs:104`, wired at `runner.rs:909`, exits 2 at
`runner.rs:1194`). The work here is to make it *usable and advertised*, not to
invent it.

**Files:**
- Modify: `tools/fasthelper/src/adapter_guard.rs` (measure on the step's own
  output only — see Task 2.18; add a standalone `fasthelper guard` verb)
- Modify: `tools/fasthelper/src/cli.rs` (new `guard` subcommand)
- Modify: `tools/fasthelper/src/report.rs` (surface the measured residual in
  the aggregate report)
- Test: `tools/fasthelper/tests/adapter_guard.rs`

**Step 1: add a standalone `fasthelper guard` verb**

Today the guard only exists inside `run`. A standalone verb lets a user check
*their own* trimmed FASTQ without running the whole chain:

```
fasthelper guard -1 trimmed_R1.fastq.gz [-2 trimmed_R2.fastq.gz] \
    --adapter-preset illumina --fail-on-adapter-residual
```

Exit 2 when the residual exceeds the threshold, 0 otherwise, with the measured
percentage on stderr in one line.

**Step 2: measure the guard's own accuracy**

The guard measures adapter content with `fastcheck --adapter-content`. Verify it
detects the residuals the undertrim corpus actually contains (the UT-5 arm: a
read that is 100 % adapter) and that it does not fire on the clean corpus.

**Step 3: document the contract**

In `docs/CLI-SPEC.md` (the fastprep and fasthelper sections) and in the README:
one command, one exit code, one sentence of what it proves.

**Step 4: commit**

```bash
git add tools/fasthelper/src/ tools/fasthelper/tests/adapter_guard.rs docs/CLI-SPEC.md
git commit -m "feat(fasthelper): a standalone residual-adapter guard (exit 2 on trip)"
```

---

## 4. Phase 2 — the speed levers, in leverage order

### Wave A — free wins (no behaviour change, no measurement risk)

### Task 2.14 — fastmaps: `lto` and `codegen-units`

**Objective:** fastmaps' own `Cargo.toml` `[profile.release]` carries **only**
`opt-level = 3` — every other tool in the suite has `lto`/`codegen-units`.

**Files:**
- Modify: `tools/fastmaps/Cargo.toml`

**Step 1: add the two flags**

```toml
[profile.release]
opt-level = 3
lto = "fat"
codegen-units = 1
```

**Step 2: build and verify byte-identity**

Run: `pixi run build`
Run: `cd tools/fastmaps && cargo test --offline` — Expected: all 96 pass.
Run the fastmaps output on a fixture and compare the JSON byte-for-byte with the
pre-change build.

**Step 3: measure**

Run the 5 M-pair profile (`bench/profile/pipeline_5m.csv` protocol):
fastmaps was **3.856 s wall / 7.552 s cpu / 22.7 MB peak**. Record the new
wall and cpu.

**Step 4: commit**

```bash
git add tools/fastmaps/Cargo.toml
git commit -m "build(fastmaps): lto=fat, codegen-units=1"
```

---

### Task 2.13 — fastmaps: cheap hasher and a `Vec` for per-file counts

**Objective:** remove SipHash-1-3 from the per-read hot path.

**Why it is safe:** the keys are `(u16, u16)` pairs derived from read content,
never from untrusted network input — the DoS resistance SipHash provides is
pointless here, and the prototype measured **9.4 ns/read → 5.9 ns/read**
(identity+multiply) and a dense-vec accumulator at **1.0 ns/read**.

**Files:**
- Modify: `tools/fastmaps/src/worker.rs` (the four maps at lines 64, 158, 184, 396)
- Modify: `tools/fastmaps/Cargo.toml` (`fxhash` dependency, or hand-rolled)
- Test: `tools/fastmaps/tests/record_integrity_test.rs` + a new determinism test

**Step 1: add the dependency**

`fxhash = "0.2"` (the suite already uses `rustc-hash` in fastcheck — check
`tools/fastcheck/Cargo.toml` and **copy that dependency and its version pin**;
do not introduce a second hashing crate).

**Step 2: swap the map types**

```rust
pub counts: HashMap<(u16, u16), PixelAcc>,      // -> FxHashMap
pub per_file_counts: HashMap<u16, u64>,         // -> Vec<u64> indexed by file_id
```

`per_file_counts` is keyed by `file_id: u16` — a dense index, not a hash key.
Replace it with a `Vec<u64>` and update the three use sites.

**Step 3: verify byte-identity** — the heatmap JSON must be byte-identical
(the map iteration order changes, so the **serialisation** must sort; check
`tools/fastmaps/src/report.rs` before assuming it already does).

**Step 4: measure** — 5 M pairs, wall and cpu vs `pipeline_5m.csv`.

**Step 5: commit**

```bash
git add tools/fastmaps/src/worker.rs tools/fastmaps/src/report.rs \
        tools/fastmaps/Cargo.toml tools/fastmaps/tests/ \
        scripts/expected-test-counts.txt
git commit -m "perf(fastmaps): cheap hasher for the pixel maps, Vec for per-file counts"
```

---

### Task 2.7 — fastprep: PE-path allocation trims

**Objective:** remove two per-pair allocations.

**Measured:** `pe_overlap.rs:177` allocates a `Vec<u8>` via `revcomp(r2)` per
pair; `pe_overlap.rs:196` allocates `admissible: Vec<usize>` per pair. Small per
pair, whole-chunk in aggregate.

**Files:**
- Modify: `tools/fastprep/src/pe_overlap.rs`

**Step 1: thread-local scratch buffers**

`revcomp` is a public function used by tests; add a private
`revcomp_into(src, dst)` that writes into a caller-provided buffer, and give
`detect_pe_overlap` a reusable `&mut Vec<u8>` scratch owned by the caller
(the pipeline's chunk buffer already recycles arenas — see
`pipeline.rs:618-631`).

**Step 2: reuse the `admissible` vector** — keep one `Vec<usize>` per chunk
arena instead of one per pair.

**Step 3: verify byte-identity** — the detector's output must be identical;
run `tools/fastprep/tests/pe_overlap.rs` (which asserts byte-identical output
when the mode is off).

**Step 4: commit**

```bash
git add tools/fastprep/src/pe_overlap.rs
git commit -m "perf(fastprep): remove the per-pair allocations in the PE path"
```

---

### Task 2.8 — fastprep: metadata writer staging

**Objective:** give the PXMT writer the same 256 KiB staging the FASTQ writer
already has.

**Measured:** metadata costs **+0.11 s wall / +0.20 cpu** at 200k pairs, at
**38.4 bytes per read**.

**Files:**
- Modify: `tools/fastprep/src/metadata.rs` (`write_record`, `write_file_header`)
- Modify: `tools/fastprep/src/pipeline.rs` (the metadata writer construction)
- Test: `tools/fastprep/tests/meta_viewer.rs`

**Step 1: implement the staging buffer** — accumulate records into a 256 KiB
`Vec<u8>` and flush when full; the PXMT format is unchanged (byte-identical
output).

**Step 2: verify byte-identity** — the `.pxmt.gz` file must be byte-identical
to the pre-change output for the same input.

**Step 3: measure** — 200k pairs with `--metadata`; expect ~0.1 s off the
wall time.

**Step 4: commit**

```bash
git add tools/fastprep/src/metadata.rs tools/fastprep/src/pipeline.rs
git commit -m "perf(fastprep): stage the PXMT writer"
```

---

### Task 2.20 — fasthelper: parse each report once

**Objective:** stop `check → aggregate → concatenate → html` from re-reading
every report.

**Files:**
- Modify: `tools/fasthelper/src/report.rs` (+ `concatenate.rs` if it exists)
- Test: `tools/fasthelper/tests/compare_test.rs`

**Step 1: implement a parsed-model cache** — parse each report JSON once into
the in-memory model and render every output from that.

**Step 2: verify byte-identity** — the HTML and TSV outputs must be
byte-identical to today's.

**Step 3: commit** — `git commit -m "perf(fasthelper): parse each report once"`.

---

### Wave B — the gzip lever (the biggest wall-time lever in the suite)

**The measurement that decides this wave:** deflate is **86 %** of the fastsim
run; the fastdedup writer is the memory hog (**4 372 MB at gz6 vs 154 MB at
gz1**, 28×); on 5 M pairs gz1 is **18.6 s vs 22.79 s at gz6** and gz9 is
**121 s**. Level 1 is *faster* than level 6 here, and the output is 88 % larger.

### Task 2.4 — the default gzip level

**Objective:** make the fast default the fast one, and document the trade.

**Three separate defaults, three separate decisions:**

| Where | Today | Proposed | Task |
|---|---|---|---|
| `fastprep --gz-level` (`cli.rs:423`) | 6 | 1 | 2.4a |
| `fastdedup --gz-level` (`cli.rs:83-89`) | 6 | 1 | 2.4b |
| `fasthelper` config `default_gz()` (`config.rs:60-62`) | 6 | 1 | 2.4c |
| `fastsim --gz-level` (`cli.rs:110-114`) | 6 | **unchanged** | — |

**fastsim keeps level 6 deliberately**: its bytes are the ground-truth corpus
for the closed loop (`tools/fastsim/src/fastq.rs:341` asserts
`DEFAULT_GZ_LEVEL == 6` "the default must stay 0.1.0's level"). Changing it
would invalidate every stored campaign corpus for no runtime gain in the
simulator itself (the runner-level change in Task 2.1 is what buys the speed).

**Step 2.4a: fastprep**

Change `default_value_t = 6` → `1` in `tools/fastprep/src/cli.rs:423`; update
the doc-comment to state the trade (".gz bytes grow ~18 %; the decompressed
FASTQ is identical"). **Interaction to verify:** `cli.rs:1050-1053` sizes
`gz_threads` from `--threads`, and `io.rs:52-62` records that single-threaded
deflate was 73 % of wall clock at level 6. At level 1 the same thread count
applies — verify the level-1 default does not lose the parallel-gzip path.

**Step 2.4b: fastdedup**

Change `default_value_t = 6` → `1` in `tools/fastdedup/src/cli.rs:83-89`.
**Critical interaction:** `tools/fastdedup/src/parallel.rs:222` uses gzippy
(parallel gzip) **only when `compression_level >= 6`**. Dropping the default to
1 therefore *disables parallel compression* — the very parallelism that makes
`-j 8` fast (j1 94.84 s, j2 92.13 s, j4 92.2 s, **j8 22.73 s**). The measured
gz1 number (18.6 s) must be re-measured with the level-1 path **and** the
parallel-compression decision resolved, or the default flip is a net loss.

**This interaction is the single most likely way the gzip lever backfires.**
Measure it before flipping anything:
```bash
# 5 M pairs, uncompressed output, to separate parse from compression cost
fastdedup -1 R1 -2 R2 -o /dev/null --gz-level 1 -j 8
fastdedup -1 R1 -2 R2 -o /dev/null --gz-level 6 -j 8
fastdedup -1 R1 -2 R2 -o out.fq.gz --gz-level 1 -j 8
```

**Step 2.4c: fasthelper config**

`default_gz()` returns 6 in `tools/fasthelper/src/config.rs`. The runner
passes `--gz-level` to every tool (`runner.rs:770, 790`). Change the default to
1 and document it in the config doc-comment.

**Step 3 (all three): document the trade in `docs/CLI-SPEC.md`** and in each
tool's `--help` text. A user must be able to see the cost of the default.

**Step 4: version bumps** — a default change is a behaviour change ⇒ minor
bumps (fastprep 0.3.0 → 0.4.0, fastdedup 1.3.0 → 1.4.0, fasthelper 0.1.0 →
0.2.0) with `docs/RELEASES.md` entries. Per C-8, show the diff and get the
user's OK before committing the bumps.

**Step 5: commit** — one commit per tool (D-1: the tools ship independently).

---

**Step 6: commit (one commit per tool — D-1: the tools ship independently)**

```bash
git add tools/fastprep/src/cli.rs tools/fastprep/tests/ scripts/expected-test-counts.txt
git commit -m "perf(fastprep): default the gzip output level to 1"

git add tools/fastdedup/src/cli.rs tools/fastdedup/tests/ scripts/expected-test-counts.txt
git commit -m "perf(fastdedup): default the gzip output level to 1"

git add tools/fasthelper/src/config.rs tools/fasthelper/tests/ scripts/expected-test-counts.txt
git commit -m "perf(fasthelper): default the run gzip level to 1"
```

### Task 2.1 — FS2-1: the runner-level gzip lever for simulated corpora

**Objective:** take the chain's single largest wall-time item — fastsim is
**283.5 s of the ~307 s chain (80 %)** — and cut its deflate cost.

**Why this is separate from Task 2.4:** fastsim's own default stays at 6
(byte-compat with every stored campaign corpus). The *runner* selects a fast
level for the corpora it generates as intermediates.

**Files:**
- Modify: `tools/fasthelper/src/config.rs:60-62` (`default_gz`) + a new `Run.gz_level_intermediate`
- Modify: `tools/fasthelper/src/runner.rs` (`argv_for`, the `simulate` branch)
- Modify: `tools/fasthelper/src/cli.rs` (a `--gz-level-intermediate` override)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: implement**

```toml
# fastsuite.toml
[run]
gz_level = 6            # final outputs (unchanged)
gz_level_intermediate = 1   # simulated corpora: deflate is 86 % of their cost
```

`argv_for` passes `--gz-level 1` for the simulate step and `--gz-level
<gz_level>` for the real steps.

**Step 2: measure the chain end to end**

Run: `bash bench/scripts/chain.sh` (or the demo) and compare against
`bench/profile/pipeline_5m.csv` (chain **283.5 s** → expected ~22.6 s for the
simulate stage at level 1).

**Step 3: commit**

```bash
git add tools/fasthelper/src/ tools/fasthelper/tests/cli.rs
git commit -m "feat(fasthelper): a fast gzip level for simulated intermediates"
```

---

### Wave C — work partitioning (the wall-clock levers)

### Task 2.2 — FS2-2: parallel simulation with per-pair deterministic RNG

**Objective:** make fastsim use the cores. Today `-j` is accepted and ignored
(`cli.rs:294-298` warns "simulation is single-threaded").

**The design constraint that makes this hard:** fastsim's RNG is a single
stream, so pairs cannot be generated independently *and* keep the current
bytes. The fix is to derive each pair's stream from `(seed, pair_index)`:
reproducible AND parallel — but the emitted bytes change even at `-j 1`.

**Files:**
- Modify: `tools/fastsim/src/artefact.rs` (`decide`, `build_read` — take a
  per-pair RNG instead of a shared one)
- Modify: `tools/fastsim/src/main.rs` (the generation loop → a rayon pool over
  pair index ranges, with an ordered emit)
- Modify: `tools/fastsim/src/truth.rs` (**the GT format must carry the RNG
  contract** — a `gt_version` and the seeding rule)
- Modify: `tools/fastsim/src/cli.rs` (the "single-threaded" warning goes away)
- Test: `tools/fastsim/tests/golden_stream.rs`, `tools/fastsim/tests/cli_test.rs`

**Step 1: write the contract test first**

```rust
#[test]
fn the_same_seed_and_pair_count_give_the_same_bytes_at_any_thread_count() {
    // fastsim -j 1 and -j 8 must produce byte-identical FASTQ and an
    // identical GT. This is the property that makes the parallelism safe.
}
```

**Step 2: implement the per-pair stream**

`artefact.rs:49` already notes "needs per-pair seeding and is a separate
change". Derive each pair's RNG as
`ChaCha8Rng::seed_from_u64(seed ^ hash(pair_index))` — the exact function must
be written down in `truth.rs` so the GT is self-describing.

**Step 3: bump the GT version and document the break**

The GT schema changes ⇒ `docs/REPORT-SCHEMA.md` versioning rules (major: a new
schema string, old key aliased for one release) and a `docs/RELEASES.md` entry.
The closed-loop corpus on disk becomes stale — say so explicitly.

**Step 4: measure**

Run the 200k-pair probe (`bench/profile/pipeline_5m.csv` protocol): expect
~8× on 8 cores.

**Step 5: commit** — `git commit -m "feat(fastsim): per-pair deterministic RNG and parallel generation"`.

---

### Task 2.11 — FC2-1: parallelise the single-input path (fastcheck)

**Objective:** fastcheck's most common invocation runs on **one** file and is
therefore serial.

**The defect:** `tools/fastcheck/src/main.rs:240` takes the serial branch when
`files.len() == 1 || threads <= 1`. Measured: a 2×500k pair at `-j 8` runs at
**1.93 busy cores** (8.159 cpu / 4.219 wall).

**Files:**
- Modify: `tools/fastcheck/src/main.rs` (the branch at line 240)
- Modify: `tools/fastcheck/src/reader.rs` (record-boundary chunking)
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: implement chunked single-file parallelism**

Split the input on record boundaries, scan the chunks in parallel, merge in
input order. The merge must be deterministic: the existing multi-file path
already sorts by input position (`main.rs:270`) — reuse that discipline.

**Step 2: verify byte-identity** — the report for a single file at `-j 1` and
`-j 8` must be byte-identical (every report field, not just the metrics).

**Step 3: measure** — same 1 M-read file at `-j 1` vs `-j 8`: cpu/M-read must
stay flat (no duplicated work) and wall must drop.

**Step 4: commit** — `git commit -m "perf(fastcheck): parallelise the single-input path"`.

---

### Task 2.17 — FH2-1: sample-level parallelism in the aggregator

**Objective:** parallelise `fasthelper`'s report aggregation and rendering.

**Measured:** 16 shards × 625 grids — serial sum **8.90 s** → 16 concurrent
processes **1.79 s** = **5.0× on 16 cores**.

**Files:**
- Modify: `tools/fasthelper/src/report.rs`, `tools/fasthelper/src/runner.rs`
- Test: `tools/fasthelper/tests/compare_test.rs`

**Step 1: implement** — reuse the existing `-P/--jobs` knob (currently used
for steps) for the aggregation stage; the outputs must be byte-identical
(aggregate by sample, then merge in a fixed order).

**Step 2: measure + commit** — `git commit -m "perf(fasthelper): sample-level parallelism in aggregation"`.

---

### Task 2.10 — F2-1: the three-stage fastprep pipeline

**Objective:** decouple parsing from the ordered write, so a slow write no
longer stalls the parser.

**This is the largest and riskiest speed change** (it touches the ordering
invariant that the whole tool depends on). It is deliberately last in the wave.

**The current structure** (`tools/fastprep/src/pipeline.rs`):
- `pipeline.rs:602-645` — the chunk pump: `fill_chunk` → `raw_tx.send` →
  `raw_rx.recv` → `pool.install(|| process_chunk(...))`;
- `pipeline.rs:253-269` — `process_chunk`: `chunk.into_par_iter().map(work)`;
- `pipeline.rs:1097+` — `write_pairs`: the **ordered** commit thread that owns
  the dedup decision, the stats update and the emit.

**Proposed structure:**
```
parse thread  ->  ChunkBuf ring (depth 4)  ->  trim pool  ->  ordered commit thread
                     (recycled arenas)        (rayon)        (dedup, stats, emit)
```

**Files:**
- Modify: `tools/fastprep/src/pipeline.rs` (the pump and `process_chunk`)
- Test: `tools/fastprep/tests/regressions.rs` (the byte-identical property)

**Step 1: write the ordering test first**

The output order must be **exactly** the input order regardless of thread
count — this is the invariant the tool's dedup and stats depend on.

```rust
#[test]
fn output_order_is_input_order_at_every_thread_count() {
    // 10k pairs at -j 1, -j 4, -j 8: the output FASTQ must be byte-identical.
}
```

**Step 2: implement the ring** — a bounded channel of `ChunkBuf` (depth 4)
between a dedicated parse thread and the rayon trim pool; the commit thread
consumes results **in chunk order** (the current code already serialises the
commit — keep that).

**Step 3: verify byte-identity** — `-j 1` vs `-j 8` output must be
byte-identical, and identical to the pre-change `-j 8` output.

**Step 4: measure** — 5 M pairs, pinned `-j 8`: busy cores and CPU/M-pair vs
`results_gz6.csv` and `bench/profile/pipeline_5m.csv`.

**Step 5: commit** — `git commit -m "perf(fastprep): parse, trim and commit in three stages"`.

---

### Wave D — memory and honest contracts

### Task 2.12 — FC2-2: shrink fastcheck's 785 MB peak

**Objective:** fastcheck peaks at **785.0 MB** on 5 M pairs, against fastmaps'
22.7 MB on the same input.

**Files:**
- Modify: `tools/fastcheck/src/fingerprint.rs` (the duplication table)
- Modify: `tools/fastcheck/src/accum.rs` (the per-base censuses)
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: measure where the peak is** — the RSS probe showed a single 32 MB
allocation dominating VmHWM (35.3 MB vs 4.5 MB for a 1 MB allocation), so
start with the fingerprint table's pre-allocation (`memory.rs:221`
`MAX_PREALLOC_ENTRIES`, `main.rs:112` `fp_capacity`).

**Step 2: compact the duplication table** — a fixed-width key plus a 32-bit
counter in an open-addressed table, replacing the `FxHashMap` (the same
technique fastdedup uses).

**Step 3: fuse the censuses** — one scan for composition + quality + polyX +
k-mer, following the pattern already used for `--only kmer`.

**Step 4: verify byte-identity** — every report field must be identical.

**Step 5: measure + commit.**

---

**Step 6: commit**

```bash
git add tools/fastcheck/src/fingerprint.rs tools/fastcheck/src/accum.rs \
        tools/fastcheck/tests/cli.rs scripts/expected-test-counts.txt
git commit -m "perf(fastcheck): compact the duplication table and fuse the censuses"
```

### Task 2.5 — FD2-1: fastdedup's emit order is not deterministic

**Objective:** decide, and document, whether fastdedup's output order is part
of the contract.

**The measurement:** the KEEP SET is identical across modes (same sequences,
same counts; sorted-sequence MD5 `f07fb28e2037567da0170f9ecf766620`, 204 675
distinct each), but the **output order** differs (first differing record at
line 826 568 / record 206 641). Anything gating on byte-identity must compare
sorted records — which the probes had to do.

**Two options, and the plan takes the cheap one first:**
1. Serialise the emit order (the writer is already the single ordered stage) —
   a behaviour change with a performance cost.
2. Document that **order is not part of the contract** in
   `docs/REPORT-SCHEMA.md`, and make every gate in this repo compare sorted
   records.

**Decision:** option 2, plus a `fasthelper`/`fastdedup` test that asserts the
*set* is identical. Option 1 becomes a follow-up only if a user needs it.

**Files:**
- Modify: `docs/REPORT-SCHEMA.md` (the contract statement)
- Modify: `tools/fastdedup/README.md`
- Test: `tools/fastdedup/tests/streaming_tests.rs` (set-equality, not order)

**Step 2: fix the repo's own gates** — any harness that compared fastdedup
output byte-for-byte must be corrected to compare sorted records, or it will
fail spuriously once the order changes.

**Step 3: commit** — `git commit -m "docs(fastdedup): output order is not part of the contract"`.

---

### Task 2.9 — F2-3: honest `-j` and visible utilisation

**Objective:** stop users from over-threading, and let them see why.

**Measured:** saturation is at ~8 workers on this box — `-j 8` → `-j 16` costs
**+12.6 % CPU** (5.004 → 5.631 cpu-s) for **+0.2 % wall** (0.649 → 0.651).

**Files:**
- Modify: `tools/fastprep/src/cli.rs` (the `-j` help text and a clamp)
- Modify: `tools/fastprep/src/stats.rs` (a `busy_cores` metric in the JSON report)
- Test: `tools/fastprep/tests/report_v1.rs`

**Step 1: add the metric**

`stats.rs::metrics_flat` gains `busy_cores` (cpu-s / wall-s) and `threads_used`.
REPORT-SCHEMA rule 3 says `metrics` is a flat numeric namespace — a new
numeric metric is a **minor** schema change, so the schema string stays `v1`.

**Step 2: document the plateau** in the `-j` help text: "throughput saturates
around 8 workers; more threads add CPU cost without reducing wall time".

**Step 3: measure + commit.**

---

**Step 4: commit**

```bash
git add tools/fastprep/src/cli.rs tools/fastprep/src/stats.rs \
        tools/fastprep/tests/report_v1.rs scripts/expected-test-counts.txt
git commit -m "feat(fastprep): report busy cores, document the -j plateau"
```

### Task 2.19 — FH2-3: crash-resume cleanliness (fasthelper)

**Objective:** a killed run must leave a step's full output or nothing.

**The defect:** SIGKILL/timeout leaves partial sample dirs (`a_clean_R1.fastq`
etc.) with no terminal manifest record, while `--from` resume exists.

**Files:**
- Modify: `tools/fasthelper/src/runner.rs` (the step output paths)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: write to `*.tmp` and rename on success** — the same discipline
fastprep already uses (`pipeline.rs`'s `part_path` + `TmpFile` guard). **Copy
that pattern** (D-1: no shared crate, copy the code).

**Step 2: verify** — kill a run mid-step; the sample dir must contain either
the complete output or no output.

**Step 3: commit.**

---

**Step 4: commit**

```bash
git add tools/fasthelper/src/runner.rs tools/fasthelper/tests/cli.rs \
        scripts/expected-test-counts.txt
git commit -m "fix(fasthelper): a step writes to a temp path and renames on success"
```

### Task 2.21 — FH2-5: the two reproduced output defects (fasthelper)

**Objective:** verify, then fix, two reported defects.

1. `concatenate` counts the `heatmap.report` wrapper as a phantom sample (stem
   not in the GENERIC name set → "4 report(s) from 3 sample(s)").
2. `check -c` accepts duplicate `[[samples]]` names — which share a work dir, a
   manifest and a resume key.

**The 2026-10-21 wave may have fixed (1) via sample-name propagation — so
verify on the built binary before proposing a patch.**

**Files:**
- Modify: `tools/fasthelper/src/report.rs` (the phantom sample)
- Modify: `tools/fasthelper/src/cli.rs` (duplicate-name refusal)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: reproduce both** on the current build; write a failing test for each.
**Step 2: fix; Step 3: verify; Step 4: commit** (two commits — one per defect).

---

**Step 4: commit (one commit per defect)**

```bash
git add tools/fasthelper/src/report.rs tools/fasthelper/tests/cli.rs
git commit -m "fix(fasthelper): concatenate must not count the report wrapper as a sample"

git add tools/fasthelper/src/cli.rs tools/fasthelper/tests/cli.rs
git commit -m "fix(fasthelper): refuse duplicate sample names in check -c"
```

### Task 2.6 — FD2-4: diagnose the fastdedup `-j 1..4` flatline

**Objective:** explain why 4 threads do no better than 1, and either fix it or
document it.

**The measurement:** `bench/profile/dedup_threads.csv` — j1 94.84 s, j2 92.13 s,
j4 92.2 s (all ≈ 92 s wall, ≈ 100 MB) against j8 **22.73 s / 4 368 MB**. With
uncompressed output the hash pipeline is ~10 s at any `-j`, so the 92 s is
read-side I/O plus writer stall, not hashing.

**Files:**
- Modify: `tools/fastdedup/src/parallel.rs` (the queue bound from `53d6686`)
- Modify: `tools/fastdedup/README.md` (if the flatline is documented instead)
- Test: `tools/fastdedup/tests/streaming_tests.rs`

**Step 1: separate parse from compression**

```bash
# 5 M pairs, plain output: no compression cost at all
fastdedup -1 R1.fq.gz -2 R2.fq.gz -o /tmp/plain.fq --gz-level 1 -j 1
fastdedup -1 R1.fq.gz -2 R2.fq.gz -o /tmp/plain.fq --gz-level 1 -j 4
fastdedup -1 R1.fq.gz -2 R2.fq.gz -o /tmp/plain.fq --gz-level 1 -j 8
```

Expected: all three ≈ the same wall time (the ~10 s hash pipeline). That
isolates the 82 s gap to the compressed-output path.

**Step 2: instrument the writer queue** — measure how long a worker waits for
the compression thread at `-j 1` and `-j 4`. If the wait dominates, the queue
bound is the cause; if the wait is negligible, the cost is the single-threaded
deflate itself (which Task 2.4b's parallel-compression decision also touches).

**Step 3: act**
- If the queue bound is the cause: raise it and re-measure j1/j4.
- If deflate is the cause: this is the same finding as Task 2.4b — keep the
  parallel-compression path for levels ≥ 6 and document the flatline at low `-j`.

**Step 4: commit**

```bash
git add tools/fastdedup/src/parallel.rs tools/fastdedup/README.md
git commit -m "perf(fastdedup): the -j 1..4 flatline, diagnosed and fixed or documented"
```

---

### Task 2.6b — FD2-3: a faster deflate backend

**Objective:** deflate is the suite's largest single wall-time cost; a faster
backend buys it back without changing any semantics.

**Files:**
- Modify: `tools/fastdedup/Cargo.toml` (the `flate2` backend feature)
- Modify: `tools/fastprep/Cargo.toml` (the same, in fastprep's own crate — D-1:
  copy, do not share)
- Test: byte-identity is **not** expected here (the .gz bytes change); the
  decompressed content must be identical.

**Step 1: try `zlib-rs`** behind flate2's `zlib-rs` feature (a drop-in
backend, same API, 2–3× faster deflate at an equal-or-better ratio). Measure
on the `gz_levels.csv` protocol: wall, peak RSS and output size at levels 1, 6
and 9.

**Step 2: verify the content** — `zcat out.fq.gz | md5sum` must equal the
plain-output md5. The .gz bytes differ; that is the declared behaviour change.

**Step 3: commit**

```bash
git add tools/fastdedup/Cargo.toml tools/fastprep/Cargo.toml
git commit -m "perf: zlib-rs deflate backend in fastdedup and fastprep"
```

---

### Task 2.6c — FD2-5: the 128-bit hash path

**Objective:** the 128-bit fingerprint exists but is unreachable without
`--pre-count`; either expose it or document the 64-bit collision bound.

**Files:**
- Modify: `tools/fastdedup/src/cli.rs` (the `--pre-count` gate)
- Modify: `docs/REPORT-SCHEMA.md` (the collision bound, if documented instead)
- Test: `tools/fastdedup/tests/abundance_tests.rs`

**Step 1: find the gate** — `grep -rn 'pre_count' tools/fastdedup/src/` and
read the condition that keeps the 128-bit path unreachable.

**Step 2: choose**
- Expose it: document when the 128-bit path is worth its extra memory.
- Or document: state the 64-bit collision bound (birthday estimate at the
  corpus scale) in `docs/REPORT-SCHEMA.md`, so a reader knows what the
  duplication number means.

**Step 3: commit**

```bash
git add tools/fastdedup/src/cli.rs docs/REPORT-SCHEMA.md
git commit -m "docs(fastdedup): the fingerprint width and its collision bound"
```

---

### Task 2.3 — FS2-3: fuse the per-base generation loops (fastsim)

**Objective:** `quality.rs` and `biology.rs` each walk the read; one walk with
hoisted distribution objects should be faster.

**The constraint:** byte-identical **only if the RNG draw order is kept**. The
closed loop depends on that stream.

**Files:**
- Modify: `tools/fastsim/src/quality.rs`
- Modify: `tools/fastsim/src/biology.rs`
- Test: `tools/fastsim/tests/golden_stream.rs`

**Step 1: pin the draw order with a test** — before touching the loop, add a
test that asserts the exact sequence of RNG draws for one read (the golden
stream test already exists: extend it with a per-read draw trace).

**Step 2: fuse the loops** — one walk over the read, hoisting the distribution
objects out of the per-base call (no per-base RNG construction) and replacing
the quality CDF's per-base search with a LUT.

**Step 3: verify byte-identity** — the output FASTQ and the GT must be
byte-identical to the pre-change build for the same seed.

**Step 4: measure** — the 200k-pair probe.

**Step 5: commit**

```bash
git add tools/fastsim/src/quality.rs tools/fastsim/src/biology.rs         tools/fastsim/tests/golden_stream.rs
git commit -m "perf(fastsim): fuse the per-base generation loops"
```

---

### Task 2.3b — FS2-4: `--stats` counters must be free when absent

**Objective:** verify the counters sit behind the flag in the hot loop, not
collected and discarded.

**Files:**
- Modify: `tools/fastsim/src/stats.rs` (only if the check fails)
- Test: `tools/fastsim/tests/cli_flags_test.rs`

**Step 1: inspect the hot loop** — `grep -n 'stats' tools/fastsim/src/main.rs`
and check each counter write is inside an `if stats { ... }` guard.

**Step 2: measure** — the same 200k-pair run with and without `--stats`; the
wall times must agree within noise.

**Step 3: if the counters are unconditional, move them behind the flag and
re-measure.** Then commit:

```bash
git add tools/fastsim/src/stats.rs tools/fastsim/src/main.rs
git commit -m "perf(fastsim): keep the --stats counters out of the hot loop"
```

---

### Task 2.15 — FM2-3: fused per-read passes (fastmaps)

**Objective:** three accumulator updates in one pass beat three scans.

**Files:**
- Modify: `tools/fastmaps/src/worker.rs`
- Test: `tools/fastmaps/tests/record_integrity_test.rs`

**Step 1: write the byte-identity test** — the heatmap JSON must be identical
before and after the fusion (the accumulator order must not change any sum).

**Step 2: fuse the passes** — one walk over the read updating composition,
quality and the k-mer rolling register together.

**Step 3: measure** — 5 M pairs, wall and cpu vs `pipeline_5m.csv`.

**Step 4: commit**

```bash
git add tools/fastmaps/src/worker.rs tools/fastmaps/tests/record_integrity_test.rs
git commit -m "perf(fastmaps): fuse the per-read accumulator passes"
```

---

### Task 2.16 — FM2-4: pack the index record to 8 bytes

**Objective:** the index record is documented as 8 B but is `(u16,u16,u16,u32)`
= 10 B padded to 12; the spilled record already drops `y`.

**Files:**
- Modify: `tools/fastmaps/src/worker.rs` (the in-memory index record)
- Modify: `tools/fastmaps/src/spiller.rs` (the spilled record, for reference)
- Test: `tools/fastmaps/tests/record_integrity_test.rs`

**Step 1: find the struct** and confirm the 12 B padded size.

**Step 2: pack it** — either drop `y` (as the spilled record already does) or
pack two u16 + a u32 into 8 B. Prove equivalence: the index written to disk
must be byte-identical, and the in-memory size must drop.

**Step 3: measure** — index memory at 5 M pairs.

**Step 4: commit**

```bash
git add tools/fastmaps/src/worker.rs tools/fastmaps/src/spiller.rs
git commit -m "perf(fastmaps): pack the index record to 8 bytes"
```

---

### Task 2.12b — FC2-3: profile the fitter width from the observed read length

**Objective:** 512-wide census rows on short libraries are mostly padding.

**Behaviour change:** fitted values shift within tolerance — flag it loudly.

**Files:**
- Modify: `tools/fastcheck/src/accum.rs` (the `BASE_POSITIONS` width)
- Modify: `tools/fastcheck/src/main.rs` (the `max_pos` derivation)
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: write the tolerance test** — on a short-read corpus, the fitted
profile must stay within the closed loop's tolerances (G fraction within 0.02,
mean Q within 0.5) against the 512-wide baseline.

**Step 2: derive the width** from the observed read length (clamped), keeping
512 where the data reaches it.

**Step 3: measure** — the 5 M-pair peak RSS.

**Step 4: commit**

```bash
git add tools/fastcheck/src/accum.rs tools/fastcheck/src/main.rs
git commit -m "perf(fastcheck): size the census rows from the observed read length"
```

---

### Task 2.12c — FC2-4: the fixed startup cost on tiny inputs

**Objective:** the full pass costs far more than any single census on the same
input — the evidence for the FC2-2 fusion.

**Files:** none (a measurement task; it produces the number Task 2.12's fusion
is judged by).

**Step 1: measure the arms**

```bash
fastcheck --only dup  --json /tmp/dup.json  input.fastq
fastcheck --only comp --json /tmp/comp.json input.fastq
fastcheck             --json /tmp/full.json input.fastq
```

on 1 M reads and on 100 k reads, recording wall and peak RSS for each.

**Step 2: record the table** in `docs/optim/round2-fastcheck.md` (the
`--only` arm timings) and check it against the published numbers (1 M:
`--only dup` 0.437 s / 103.9 MB; 100 k: 0.039 s / 9.9 MB).

**Step 3: commit**

```bash
git add docs/optim/round2-fastcheck.md
git commit -m "docs(fastcheck): the per-census arm timings"
```

### Task 2.18 — FH2-2: reuse the step's own fastcheck pass in the guard

**Objective:** the guard runs an extra `fastcheck --adapter-content` pass on
the reads entering AND leaving every trimming step. Where the chain already
runs fastcheck on that output, reuse it.

**Files:**
- Modify: `tools/fasthelper/src/adapter_guard.rs`
- Modify: `tools/fasthelper/src/runner.rs` (the step wiring)
- Test: `tools/fasthelper/tests/adapter_guard.rs`

**The gate:** "proven-equivalent" — the measured adapter content must be the
same number, obtained from one process instead of two.

**Step 1: find the duplicate pass** — the guard calls `measure()`
(`adapter_guard.rs:280`), which spawns `fastcheck --adapter-content`. Find
where the chain already runs that same check on the same file.

**Step 2: pass the existing measurement in** — change `measure()` to accept an
already-measured value, and have the runner hand it over when the chain
produced one. Keep the spawn as the fallback.

**Step 3: verify equivalence** — the guard's verdict must be identical on a
corpus where the chain ran fastcheck, and the process count must drop by one
per step.

**Step 4: commit**

```bash
git add tools/fasthelper/src/adapter_guard.rs tools/fasthelper/src/runner.rs
git commit -m "perf(fasthelper): reuse the chain's own fastcheck pass in the guard"
```

### Task 3.1 — reconcile `winners.tsv` and `fastp-comparison.tsv`

**Objective:** resolve the 93/101 disagreement so a win-rate can be quoted — or
so the plan can state plainly that none may be.

**Files:**
- Read: `bench/results/phase-d/{winners.tsv,fastp-comparison.tsv}` (inputs)
- Modify: `bench/scripts/aggregate.py` and/or `bench/scripts/fastp_compare.py`
  (the **producer** — never the TSV by hand)
- Create: `bench/results/phase-d/reconciliation.md` (the finding)

**Step 1: find the disagreement** — the two files are produced by different
scripts over the same corpus. Compare them key by key (dataset × variant ×
config) and classify each disagreement into one of:
- a different metric definition (the same quantity measured two ways),
- a different baseline (e.g. `fastprep-hd0` vs the swept winner),
- a different corpus slice (the v1/v2 duplication: `sim.v1.gt.tsv` and
  `sim.v2.gt.tsv` are byte-identical in 100 of 102 datasets, so variant 2 was
  never an independent draw),
- or a genuine bug in one of the two producers.

**Step 2: write the finding down** — `bench/results/phase-d/reconciliation.md`
with the classification and the count per class. This file is the evidence for
whatever the plan is allowed to claim.

**Step 3: fix the producer, not the TSV** — a hand-edited result file is not a
fix, and `scripts/ci.sh`-style staleness would hide it later.

**Step 4: re-run the comparison** over the corpus (desktop: the corpus lives
at `/mnt/Fast/fastsuite-scratch`; the laptop cannot run it).

**Step 5: decide the claim**
- Resolved: the win-rate may be quoted, naming **both** files and the
  reconciliation (C-6).
- Unresolved: the claim is **withheld**, and `docs/FASTP-COMPARISON.md` says
  so explicitly, naming both files.

**Step 6: commit**

```bash
git add bench/scripts/ bench/results/phase-d/reconciliation.md
git commit -m "bench: reconcile winners.tsv against fastp-comparison.tsv"
```

### Task 3.2 — `docs/FASTP-COMPARISON.md`: the honest head-to-head

**Objective:** one page a sceptic can check, with every number traceable to an
artifact and every caveat attached.

**Files:**
- Create: `docs/FASTP-COMPARISON.md`
- Read: `docs/optim/raw/round2/_salvage-fastprep-speed.md`,
  `~/work/FastPrep-Bench/results/results_gz6.csv`,
  `bench/results/phase-d/reconciliation.md` (from Task 3.1)

**Required content:**
1. **CPU-seconds per million pairs**: fastprep **16.36** vs fastp **17.64**
   (7 % less work per pair) — with the note that the wall-time difference is
   utilisation, not work.
2. **Memory**: fastprep **101 MB** peak vs fastp **1 141 MB** (11× lighter).
3. **The retired cell**: the published "1-core 28 s vs 34 s" row is broken
   (the fastprep rows escaped `taskset` and burned 2.44 cores inside a
   supposedly pinned run) and is **retired**, with the reason stated.
4. **The residual-adapter contract**: `fasthelper guard --fail-on-adapter-residual`
   — one command, exit 2, and the measured before/after residual from Phase 1.
5. **Both phase-d files named** if a win-rate is quoted; otherwise none.

**Step 1: write the page** with a table whose every row names its artifact
path. **No number without a file.**

**Step 2: check every number against its artifact** — re-read each source
file and confirm the figure, the units and the conditions. A number copied
from a dossier without re-reading its artifact is exactly the class of error
this plan is built to prevent.

**Step 3: commit**

```bash
git add docs/FASTP-COMPARISON.md
git commit -m "docs: the fastp head-to-head, every number traceable to an artifact"
```

### Task 3.3 — the demo as the proof

**Objective:** make the claim demonstrable in one command.

`bash examples/demo.sh 20000` already runs the whole chain end to end. Extend
it so the run ends with the residual-adapter guard, and the user sees the
contract fire.

**Files:**
- Modify: `examples/demo.sh`
- Modify: `README.md` (the quickstart section)

**Step 1: add the guard to the demo** — after the chain, run the guard on the
trimmed output and print its exit code.

**Step 2: run it** — `bash examples/demo.sh 20000` must end with the guard's
verdict and exit code, and the script's own exit status must reflect a guard
trip.

**Step 3: commit**

```bash
git add examples/demo.sh README.md
git commit -m "feat(demo): the residual-adapter guard closes the demo"
```

### Task 4.1 — the test-count manifest

**Objective:** keep `scripts/ci.sh`'s staleness detector honest.

**Files:**
- Modify: `scripts/expected-test-counts.txt`

**Step 1: run the gate**

```bash
pixi run ci
```

Expected: every per-tool count matches the manifest, GRAND equals the sum of
the tool lines, and the pathology corpus passes.

**Step 2: if a count is SHORT, do not lower the number.** A short count means
cargo skipped a rebuild and the tests silently did not run. Find out why:
`pixi run build` first (the ambient cargo is the wrong compiler), then
`cargo clean` in the affected tool, then re-run.

**Step 3: if a count is LONG, update the manifest** in the same commit — the
file says a long count means tests were added without the manifest being
updated, which is a deliberate, reviewable change.

**Step 4: commit**

```bash
git add scripts/expected-test-counts.txt
git commit -m "test: the expected-count manifest, updated for the round-3 suite"
```

### Task 4.2 — the documentation wave

**Objective:** every document that describes the code must describe the code as
it is now. A stale plan doc is the single biggest source of false claims.

**Files:**
- Modify: `docs/STATUS.md`, `docs/CLI-SPEC.md`, `docs/REPORT-SCHEMA.md`,
  `docs/RELEASES.md`, `docs/UNDERTRIM-REDUCTION.md`, `docs/OVERTRIM-PAIRED.md`,
  `docs/DECISIONS.md`, `docs/optim/round2-*.md`, `README.md`

**Step 1: fix the two known drifts**
1. `docs/STATUS.md` — the test count line is **61 tests stale** (it says 679
   with 143/122/187/89/77/60; the manifest says 740 with 143/134/203/93/96/70).
   Copy the manifest numbers, and add the date of the verification.
2. `docs/optim/round2-fastprep.md` — the F2-1 citation `pipeline.rs:445-501`
   is `detect_phred_offset`, not the parse/trim stage. The real anchors are
   `tools/fastprep/src/pipeline.rs:602-645` (the chunk pump) and
   `tools/fastprep/src/pipeline.rs:253-269` (`process_chunk`).

**Step 2: update the specs for what shipped**
3. `docs/CLI-SPEC.md` — the `--pe-overlap` row gains the new `--adapter-evidence`
   semantics; the `--gz-level` rows state the new default and its cost; the
   new `fasthelper guard` verb is listed.
4. `docs/REPORT-SCHEMA.md` — the new metrics (`busy_cores`, `threads_used`); the
   fastdedup order contract; the GT schema bump if FS2-2 shipped.
5. `docs/RELEASES.md` — one entry per version bump, per D-2.
6. `docs/UNDERTRIM-REDUCTION.md` — the "before"/"after" table from Task 1.6.
7. `docs/OVERTRIM-PAIRED.md` — the final Cell A number.
8. `docs/DECISIONS.md` — a new decision (D-12) recording that the adapter
   instrument is **evidence, not threshold**, and what that forbids.
9. `README.md` — the quickstart shows the guard.

**Step 3: verify every changed claim against the code** — re-read the cited
line for each statement. A doc claim that was not re-read is a claim, not a
fact.

**Step 4: commit**

```bash
git add docs/ README.md
git commit -m "docs: the round-3 documentation wave, verified against the code"
```

### Task 4.3 — version bumps (the user's OK first)

**Objective:** bump the versions of every tool whose behaviour changed, after
the user has seen the diff.

**Files:**
- Modify: `tools/<tool>/Cargo.toml` (the `version = ` line)
- Modify: `docs/RELEASES.md` (the entry)

**Step 1: list the tools that changed behaviour** — any tool whose default,
output bytes or trim decisions changed. Per D-2, each tool carries its own
semver; a behaviour change is at least a minor bump.

**Step 2: show the diff and wait**

```bash
git diff -- tools/*/Cargo.toml docs/RELEASES.md
```

Present it to the user, then wait for the OK. Per C-8 the commit happens
only after that OK, and nothing is ever pushed.

**Step 3: commit**

```bash
git add tools/*/Cargo.toml docs/RELEASES.md
git commit -m "chore: round-3 version bumps"
```

### Task 4.4 — the two-checkout sync

**Objective:** a change made on one machine must exist on the other.

**Files:** none (a procedure).

**Step 1: run the sync**

```bash
bash ~/hermes_scratch/fastsuite-local-artifacts/sync-fastsuite.sh
```

**Step 2: verify both checkouts**

```bash
git rev-parse main                                          # laptop
ssh raphael@100.104.77.4 "git -C /home/raphael/FastSuite rev-parse main"   # desktop
```

Expected: the two SHAs are identical.

**Step 3: no commit** — this task changes no file in the repository; it only
verifies that the two checkouts agree.

**Step 4: check the sync's residual diffs** — the "OUT OF SYNC" line is
expected when the differences are gitignored artefacts (`fastp.html`,
`fastp.json`, `.snakemake/`). Confirm with `git check-ignore -v` on each
residual path; a residual path that is NOT ignored is real drift and must be
resolved before the round is called done.

### Task 4.5 — the round-3 report

**Objective:** one document that carries the whole round's evidence.

**Files:**
- Create: `docs/optim/report/optimisations-round3.typ`
- Create: `docs/optim/report/optimisations-round3.pdf` (built, not committed —
  `docs/optim/report/*.pdf` is gitignored)
- Read: every artifact the numbers come from

**Step 1: draft the source** in the style of
`docs/optim/report/optimisations-round2.typ` (A4, Libertinus Serif 11 pt,
cetz 0.5.2 figures):
1. the undertrim before/after on **both** axes, with the claim tier on every
   number;
2. the speed before/after per tool (wall + CPU + busy cores, never the retired
   1-core cell);
3. the fastp head-to-head, with the retired cell named and the reason stated;
4. the claim tiers and the evidence discipline.

**Step 2: compile it**

```bash
typst compile docs/optim/report/optimisations-round3.typ docs/optim/report/optimisations-round3.pdf
```

**Step 3: render and inspect EVERY page** — `pdftoppm -png -r 110` into a
scratch dir, then look at each page image. The round-2 report had a compile
failure (an unclosed delimiter) and an A4 overflow (long bar labels) that only
a rendered page revealed. Shorten the bar labels to `[283.5 s]`-style tokens
and keep the levers in the propositions table.

**Step 4: deliver the PDF as a MEDIA attachment** and commit the source.

```bash
git add docs/optim/report/optimisations-round3.typ
git commit -m "docs: the round-3 report"
```

### Order of execution

```
Phase 0  (harnesses)              ── blocks everything in Phase 1
  0.1 guard probe
  0.2 paired AB mode
  0.3 one driver

Phase 1  (undertrimming)          ── the flagship; sequential (each task
  1.1 UT-1   PE boundary          depends on 0.1, 0.2, 0.3)
  1.2 UT-2   G-overprint          depends on 0.3
  1.3 UT-3   short tails (risky)  depends on 0.3, after 1.1+1.2
  1.4 UT-4   quality scoring      depends on 1.1 (shares the prefix path)
  1.5 UT-5   adapter-only reads   independent
  1.6 measurement wave + defaults  ← THE GATE for every default flip
  1.7 quality contract            independent

Phase 2  (speed)                  ── parallelisable across tools
  Wave A (free):   2.14, 2.13, 2.7, 2.8, 2.20
  Wave B (gzip):   2.4a/b/c, 2.1
  Wave C (cores):  2.2, 2.11, 2.17, 2.10
  Wave D (memory): 2.12, 2.5, 2.9, 2.19, 2.21, 2.6, 2.3, 2.15, 2.16, 2.18, 2.12b

Phase 3  (fastp case)             ── needs Phase 1's numbers
  3.1 reconcile, 3.2 the doc, 3.3 the demo

Phase 4  (ship)                   ── 4.1 manifest, 4.2 docs, 4.3 versions,
                                    4.4 sync, 4.5 report
```

### Effort estimate

| Phase | Tasks | Implement hours | Measurement hours | Notes |
|---|---|---|---|---|
| 0 | 3 | 3–4 | 1–2 | the harnesses exist as scratch code; they need packaging |
| 1 | 7 | 12–18 | 6–10 | UT-3 may be dropped if the clean axis regresses |
| 2 | ~24 | 20–32 | 12–20 | 2.10 and 2.2 are the two multi-day items |
| 3 | 3 | 3–5 | 4–8 | the corpus re-run dominates |
| 4 | 5 | 4–6 | 3–5 | includes the report |
| **Total** | **42** | **42–65** | **26–45** | |

Wall-clock with one implementer and one measurement pass in parallel: **~5–7
working days**. With the Phase 2 waves fanned out across tools (they are
independent crates, D-1), **~3–4 days**.

**Cost:** the whole plan runs on the session's own model — no paid inference
is needed for the code work. The external costs are compute: the 5 M-pair
profiles (~30 min each on this box), the Cell A corpus run on the desktop
(~1 h), and the phase-d campaign re-run (~2–4 h on 8 cores). If a stronger
model is wanted for the two hard tasks (UT-3's corroboration design and F2-1's
three-stage pipeline), that is a per-task override — say so before dispatch.

### The three tasks most likely to fail

1. **Task 2.4b (fastdedup gzip default)** — the level-1 default disables
   gzippy's parallel compression (`parallel.rs:222`), which is what makes `-j 8`
   fast. Measure before flipping.
2. **Task 1.3 (UT-3)** — the corroborated short tail is the one fix whose
   uncorroborated form costs 180× on the overtrim axis. It may be rejected by
   its own gate; that is an acceptable outcome.
3. **Task 2.2 (fastsim parallel)** — per-pair RNG changes the emitted bytes
   even at `-j 1` and invalidates every stored campaign corpus. The GT schema
   bump and the closed-loop re-run are mandatory, not optional.

---

## 8. Open questions for the user

1. **UT-3:** if the corroborated short tail cannot hold the clean axis at
   ≤ 0.02 %, do we ship UT-1 + UT-2 + UT-4 only (projected residual ≈ 20–25 %
   on PE data), or relax the corroboration and accept a measured false-trim
   rate?
2. **fastdedup's gzip default:** if the level-1 default costs the parallel
   compression path, do we keep level 6 as the default and expose a
   `--gz-level 1` recommendation, or flip the default and pay the single-
   threaded deflate?
3. **The output-order contract (FD2-1):** is "order is not part of the
   contract" acceptable for fastdedup, or must the writer serialise?
4. **fastsim's GT break (FS2-2):** acceptable to invalidate the stored campaign
   corpora for ~8× simulation speed?

## 9. Definition of done

- [ ] `pixi run ci` green, with `scripts/expected-test-counts.txt` matching.
- [ ] Every Phase 1 task measured on **both** axes, with the numbers in
      `docs/UNDERTRIM-REDUCTION.md`.
- [ ] Cell A re-run on the desktop corpus: `pe_trims ≈ 0`.
- [ ] Every Phase 2 task priced with wall + CPU + busy cores; no claim quoted
      from the retired 1-core cell.
- [ ] `docs/FASTP-COMPARISON.md` names an artifact for every number, and
      quotes no win-rate until Task 3.1 reconciles the two phase-d files.
- [ ] `docs/STATUS.md` test counts match the manifest (the 61-test drift fixed).
- [ ] `docs/optim/round2-fastprep.md` F2-1 citation corrected.
- [ ] The round-3 Typst report compiled and every page visually inspected.
- [ ] Version bumps shown as a diff and approved (C-8).
- [ ] Both checkouts synced to the same SHA.
