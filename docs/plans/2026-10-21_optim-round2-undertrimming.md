# Round-2 optimisation audit + undertrimming reduction — technical plan

Date: 2026-10-21. Tree: `30e8fb3` (post fix-wave: `154a39f` landed 8 of round-1's
14 fastprep propositions; `af268d4` made the prefilter sound (`ms_safe`, monotone),
the mismatch tolerance `floor`, the dedup anchor interior-window; `82746b6` tracks
740 tests). Audit only: **no `src/` changes** (`references/optimization-audits.md`),
no commits unless asked.

## Why round 2 is new work

1. Round 1 (`docs/optim/*`, 2026-10-20) is mostly consumed: 8/14 fastprep items
   landed; the deferred ones (FP-1 three-stage pipeline, FP-11 SIMD quality scan,
   FP-12 hash-set sharding) are exactly the "faster without losing efficiency"
   story on the current tree and must be re-priced, not re-stated.
2. The 2026-10-21 fix wave changed the trimming core itself (prefilter soundness,
   tolerance `ceil`→`floor`, anchor interior window). The undertrim/overtrim
   balance moved and has never been re-proposed since.
3. New axis from the user: speed **preserving trimming efficiency**, plus a
   significant **undertrimming** (residual-adapter) reduction — i.e. win the
   head-to-head against fastp on both wall time and residual adapter.

## Steps

### 1. Grounding (lead, ~10 min)

Read `tools/fastprep/src/{adapter.rs,algorithm.rs,pe_overlap.rs,pipeline.rs}` and
the measured basis: `docs/TRIM-ACCURACY*.md`, `docs/OVERTRIM-*.md`,
`bench/scripts/{overtrim_*,closed_loop,fastp_compare}.py`,
`bench/results/phase-d/` + `bench/results/gt2026/`. Output: mechanism notes
appended below — every place residual adapter survives today, and every known
overtrimming ceiling (the counterweight).

### 2. Reviewers in parallel (7 subagents, ~25–35 min wall)

Raw dossiers to `docs/optim/raw/round2/<tool>.md` (on disk = survives timeouts).

- **6 × speed, one per tool** (fastprep, fastcheck, fastdedup, fastmaps, fastsim,
  fasthelper). Brief = `optimization-audits.md`'s "where the per-read cost lives"
  checklist (passes per read, allocations per read, per-base transcendental/RNG
  construction, hot buffer widths, sharded-but-serialised structures, defaults the
  competitor does not run). Every proposition must state: its **output-equivalence
  gate** (byte-identical output / proven-equivalent / deliberate behaviour change),
  its claim tier `[M]`/`[E]`/`[V]`, and the A/B that would price it. Skip
  anything `154a39f`, `1a5380d`, `af268d4` or the fastcheck k-mer fusion already
  landed.
- **1 × undertrimming** on fastprep: `adapter.rs` prefilter (seed width from the
  error bound, `ms_safe` floor at defaults: 4 mismatches at k=21, 5 needs
  `--adapter-error-rate 0.26`), partial overhangs (1–5 bp), `pe_overlap.rs`
  refusal cases and overlap tolerances, quality-cut interplay (where the quality
  trim leaves an adapter fragment or the adapter search skips a low-quality tail),
  poly-G masking order (tail trimming vs adapter search), plus the measurement
  basis: fastsim GT under/over rates, `fasthelper`'s post-trim residual-adapter
  guard (`adapter_guard.rs`) as the validation instrument, fastcheck's adapter
  content module.

Model: MoA think preset (see chat); per-child model pinning is not exposed by
`delegate_task` — a cheaper reviewer tier requires a `delegation.model` config
change.

### 3. Curate (lead, ~20 min)

`docs/optim/<tool>.md` gains a `## Round 2 (2026-10-21)` section; a proposition
survives only with a code line or a committed artifact behind it. Concrete
proposed (unapplied) diffs in `docs/optim/diffs/round2-<tool>.diff` so the
propositions are code-level, per the round-1 shape. Round-1 files stay as the
historical record; nothing is silently deleted.

### 4. Undertrimming strategy (lead, ~15 min)

`docs/UNDERTRIM-REDUCTION.md`:

- Mechanism map: where residual adapter survives by cause (partial adapter
  overhang below the prefilter's effective floor; mismatch-tolerance floor at
  shipped defaults; PE-overlap refusal cases (`--pe-overlap auto` on
  `AdapterSource::Fixed` is off by design since FP-2) and overlap failures on
  short inserts; quality-cut interplay; poly-G tail hiding an upstream adapter;
  adapter-only reads).
- Ranked fixes, each with expected residual-adapter reduction and its
  **overtrimming risk** — the `OVERTRIM-*` series is the counterweight (head
  dropout went 4→0 because it overtrimmed; the anchor interior window; the
  prefilter made *less* aggressive today for soundness). Proposals must be scored
  on both axes (under- AND over-trim rates on fastsim GT), never under alone.
- Validation protocol: fastsim GT over/under rates, `fasthelper run
  --fail-on-adapter-residual` as the shipped gate, fastcheck adapter-content
  curves. **Reconcile `winners.tsv` vs `fastp-comparison.tsv` (93/101 disagree,
  2026-10-21 review) before any win-rate is quoted** — name both, or quote
  neither.

### 5. Report + delivery (lead, ~20–25 min)

One Typst PDF (`docs/optim/report/optimisations-round2.pdf`) with cetz figures
(speed trade-off per tool; undertrim-vs-overtrim trade curve), render-and-inspect
per `typst-cetz-diagrams`, delivered in chat as MEDIA. Then the two-checkout sync
(`fastsuite-local-artifacts/sync-fastsuite.sh`) for every file this run touched —
docs included, even with zero `src/` edits.

## Gates

- Audit only: no `src/` edits, no commits, no version bumps.
- Claim tiers `[M]`/`[E]`/`[V]`; benchmark numbers only from `docs/BENCHMARKS.md`
  or a re-run under `references/bench-rerun.md` — never re-time a committed number.
- The FP-0 harness defect stands: the "1 core" cell in
  `~/work/FastPrep-Bench/results/results_gz6.csv` is not pinned (fastprep
  cpu/wall 2.44). Any fastp speed claim must either fix the harness or quote
  CPU-seconds + busy cores, not that wall column.
- The report never quotes `winners.tsv` or `fastp-comparison.tsv` as canonical
  until reconciled (step 4).

## Effort

Full sweep: 7 reviewer subagents + lead curation, ~60–90 min wall, all on the
same MoA think preset; token spend shows in the monitoring channel. Fastprep-only
variant (2 reviewers, curation, `UNDERTRIM-REDUCTION.md`, PDF): ~35–45 min.


## Step 1 grounding notes (lead, read 2026-10-21 before dossiers)

- `algorithm.rs:1304` `trim()` stage order: (1) adapter detection — PE-overlap-pinned
  boundary wins (`params.pe_adapter`), else `match_adapters` on the FULL read;
  (2) insert_len = len - (t5_adapter + t3_adapter); (3) polyX on the
  adapter-trimmed sub-slice. So the reference-model worry "poly-G runs before
  adapter search and hides it" does NOT apply to fastprep's order — adapter
  search runs FIRST on the full read. The live question for undertrimming is the
  mirror image: a 3' poly-G/quality tail can push an adapter match's right end
  past what `match_adapters` accepts (anchoring semantics — check whether the
  matcher requires the match to reach the read end), and polyX only ever sees
  what the adapter step left.
- `adapter.rs`: `error_rate` (fraction, floor-based `max_mismatches` since
  `af268d4`), `min_overlap_se`, seed table exact-rejection-only (pigeonhole,
  `min_sound_k_for` with the non-monotone-bound fix). Undertrims below the
  effective k-floor fall to the direct scan — cost vs coverage is the trade.
- `pe_overlap.rs`: `MAX_ANCHOR_MISMATCHES = 1` over a 12 bp anchor, verification
  threshold 3 mismatches over the window, refuses ambiguity (more than one
  insert size verifies -> None). Refusal cases return nothing -> read falls back
  to adapter matching; those fallbacks are an undertrim source to quantify.
- FP-2 (landed): `--pe-overlap auto` does NOT run on `AdapterSource::Fixed` —
  PE-overlap trimming is off for users who passed an explicit adapter, which is
  exactly the population whose residual adapter the user wants reduced. Any
  proposal to enable it there must price the overtrim risk (the OVERTRIM-PAIRED
  counterweight).
