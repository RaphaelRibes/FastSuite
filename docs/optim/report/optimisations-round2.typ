// FastSuite, round-2 optimisation audit + undertrimming reduction, compiled report
// Sources: docs/optim/round2-*.md (speed), docs/UNDERTRIM-REDUCTION.md (quality),
//          docs/optim/raw/round2/*.md (evidence), bench/results/undertrim-*.json
// Evidence tiers: [M] measured (artifact named), [E] estimated (pricing A/B named),
//                 [V] verified on today's source or binary.
// Language: English (round 1's compiled report was in French; round 2 matches
// the round-2 artifacts and the repo docs).

#import "@preview/cetz:0.5.2": canvas, draw

#set page(paper: "a4", margin: 2cm)
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
#show figure.caption: set text(size: 9.5pt)
#show raw.where(block: true): it => block(
  width: 100%,
  fill: luma(246),
  inset: 6pt,
  radius: 2pt,
  text(size: 8.5pt, font: "DejaVu Sans Mono", it.text),
)
#set table(stroke: 0.4pt + luma(170), inset: (x: 4pt, y: 2.6pt))
#show table.cell: set text(size: 9.5pt)

#let tag-m = text(fill: rgb("#0b6b3a"), weight: "bold")[M]
#let tag-e = text(fill: rgb("#8a5a00"), weight: "bold")[E]
#let tag-v = text(fill: rgb("#0b4f8a"), weight: "bold")[V]

// ---------------------------------------------------------------------------
// Title page
// ---------------------------------------------------------------------------
#align(center + horizon)[
  #text(size: 21pt, weight: "bold")[FastSuite]
  #v(0.2em)
  #text(size: 15pt)[Round-2 optimisation audit and undertrimming reduction]
  #v(1.2em)
  #text(size: 11pt, style: "italic")[
    One measured quality cliff, the circular guard that hides it,#linebreak()
    and the speed levers ranked by measured leverage
  ]
  #v(2.5em)
  #text(size: 10pt)[
    Read-only audit of `~/hermes_scratch/FastSuite`#linebreak()
    tree `30e8fb3` + `af268d4`, 21 October 2026#linebreak()
    Proposals only --- nothing under `tools/*/src` was changed
  ]
  #v(1em)
  #text(size: 10pt)[Raphaël Ribes]
]

#pagebreak()

// ---------------------------------------------------------------------------
= Verdict

Three families of result, in this order of importance.

*1. The quality contract comes first.* The only claim in this round that a
sceptic can falsify in one command is not a speed number: it is
`fasthelper run --fail-on-adapter-residual`, which exits 2 when adapter
survives trimming (`tools/fasthelper/src/adapter_guard.rs`). fastp ships no
equivalent. Every number below should be read as supporting that contract.

*2. Undertrimming has one dominant cause, and it is not the mismatch
tolerance.* On a 200 000-read probe corpus (29.9 % of reads carrying 1..33
bases of TruSeq readthrough), shipped settings leave *41.95 % of carrier reads*
with residual adapter. The band decomposition is decisive: overhangs of 1--5,
6--9 and 10--13 bases are each $>99.9$ % undertrimmed and together are
*94.44 % of all residual*. The cause is the 14-base overlap floor
(`--min-adapter-overlap` default, `cli.rs:1199`), and the paired-end detector
that should rescue exactly these reads is switched off by a circular guard
(`pe_overlap.rs:275-288` re-imposes the same floor on the tail it is meant to
verify). Lowering the floor globally is the wrong instrument: the two-axis
sweep (figure 2) shows every loosening buys at most 16 points of residual for a
20--50x increase in false trims on clean reads.

*3. The speed story is honest but narrower than the wall column suggests.*
fastprep already does *7 % less CPU work per pair* than fastp
(16.36 vs 17.64 CPU-seconds per million pairs) at *11x lower peak memory*
(101 MB vs 1141 MB) #tag-m. The published "fastp 28 s vs fastprep 34 s at one
core" cell is broken --- the fastprep rows escaped `taskset` (2.44 cores busy
inside a supposedly pinned run) --- so this report quotes CPU-seconds and busy
cores only, never that wall column. The suite's largest wall lever is not in
fastprep at all: `fastsim` is 80 % of the chain's wall time and deflate alone
is ~86 % of `fastsim`'s.

// ---------------------------------------------------------------------------
= Where the 5 M-pair chain spends its wall time

`bench/profile/pipeline_5m.csv`, one machine, one input, `-j 8`, gzip level 6
in and out. Bars are log-scaled. Levers are the measured round-2 propositions;
each is priced in `docs/optim/round2-<tool>.md`.

#figure(
  canvas(length: 1cm, {
    import draw: *
    let c-sim = rgb("#2b6cb0")
    let c-dedup = rgb("#2f855a")
    let c-check = rgb("#b7791f")
    let c-prep = rgb("#c53030")
    let c-maps = rgb("#6b46c1")
    let bar(y, len, fill, label, val) = {
      rect((0, y - 0.26), (len, y + 0.26), fill: fill, stroke: 0.4pt + luma(120))
      content((-0.25, y), text(size: 8.5pt)[#label], anchor: "east")
      content((len + 0.18, y), text(size: 8pt)[#val], anchor: "west")
    }
    bar(4.0, 9.81, c-sim, [fastsim], [283.5 s])
    bar(3.0, 5.41, c-dedup, [fastdedup], [22.6 s])
    bar(2.0, 4.05, c-check, [fastcheck], [10.3 s])
    bar(1.0, 3.47, c-prep, [fastprep], [7.4 s])
    bar(0.0, 2.34, c-maps, [fastmaps], [3.9 s])
    line((-3.4, -0.6), (-3.4, 4.5), stroke: 0.4pt + luma(160))
    line((-3.4, -0.6), (11.6, -0.6), stroke: 0.4pt + luma(160))
    content((4.1, -1.15), text(size: 8pt)[wall seconds, log scale --- 5 M pairs, `-j 8`, gzip 6 in/out])
  }),
  caption: [Chain wall time per stage (`bench/profile/pipeline_5m.csv`), log scale.
  `fastsim` alone is 80 % of the chain. Each stage's measured top lever is in
  the propositions table (figure 4).],
) <fig-chain>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    table.header([*tool*], [*CPU-s / M pairs*], [*peak MB*], [*busy cores*]),
    [fastprep], [16.36], [101], [2.43 #tag-v (harness), 7.71 #tag-m (200 k probe, `-j 8`)],
    [fastp], [17.64], [1141], [3.16 #tag-v (harness)],
    [cutadapt], [27.2], [24], [0.98 #tag-v],
  ),
  caption: [Per-pair position at 5 M pairs, gzip 6 (`~/work/FastPrep-Bench/results/results_gz6.csv`).
  fastprep burns 7 % less CPU per pair and 11x less memory than fastp. The busy-core
  columns: the harness's `pinned` rows are defective (fastprep escaped `taskset`),
  so the honest single-core cell is unmeasured; the 200 k probe shows the current
  build reaching 7.7 busy cores at `-j 8`. Never quote the harness wall column.],
) <fig-vs>

= Speed propositions, ranked by measured leverage

Full propositions with code anchors: `docs/optim/round2-<tool>.md`. Round 1's
deferred items are re-priced there, not repeated.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    table.header([*ID*], [*change*], [*evidence*], [*gate*]),
    [FS2-1], [runner-level `--gz-level 1` for simulated corpora; keep the tool
    default for campaign byte-compat], [#tag-m gz6 13.7 s vs gz1 3.2 s at 200 k;
    deflate is 86 % of the run], [behaviour-change (.gz bytes)],
    [FS2-2], [parallel chunked simulation, per-pair deterministic RNG from
    (seed, pair index)], [#tag-e $tilde$8x on 8 cores; pairs are independent],
    [behaviour-change (versioned GT)],
    [FD2-2], [default `--gz-level` 1 for intermediates], [#tag-m 18.6 s / 154 MB
    vs 22.8 s / 4372 MB at 5 M], [behaviour-change (.gz bytes)],
    [FH2-1], [sample-level parallelism in `concatenate` / render], [#tag-m
    8.90 s serial $arrow$ 1.79 s concurrent = 5.0x], [byte-identical],
    [FC2-1], [parallelise the single-input path (the serial branch), fuse the
    per-read censuses], [#tag-v 1.93 busy cores at `-j 8` on a pair; `--only`
    arms finish 10x faster than the full pass], [byte-identical (deterministic merge)],
    [F2-1], [three-stage pipeline: parse $|$ trim $|$ ordered commit
    (re-priced FP-1)], [#tag-e raise busy cores toward fastp's utilisation at
    unchanged CPU], [behaviour-change (bytes intended identical)],
    [FM2-1], [cheap hasher for the per-pixel maps; `Vec` for per-file counts],
    [#tag-m prototype: 9.4 $arrow$ 5.9 $arrow$ 1.0 ns/read], [byte-identical],
  ),
  caption: [The seven largest speed levers. Everything else (metadata staging,
  index-record packing, deflate backend swap, `--clamp -j` past saturation) is in
  the per-tool files. Every proposition names its A/B; none is implemented.],
) <fig-levers>

// ---------------------------------------------------------------------------
= The undertrimming problem, measured

Probe corpus: `bench/scripts/undertrim_ab.py`, 200 000 synthetic reads, 140 157
adapter-free and 59 843 carrying 1..33 bases of TruSeq readthrough (29.9 %),
with 10 % terminal G-runs, 5 % Ns, 3 % substitutions in the tails. Shipped
settings (`--min-adapter-overlap 14`, `--adapter-error-rate 0.2`). "Residual"
counts a partial cut (0 $<$ cut $<$ t) as undertrim, which is what a user sees:
reads still carrying adapter.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    table.header([*overhang band*], [*reads*], [*exact %*], [*under %*], [*share of residual*]),
    [1--5 bp], [9 028], [0.00], [99.98], [*36.0 %*],
    [6--9 bp], [7 314], [0.00], [99.93], [*29.1 %*],
    [10--13 bp], [7 377], [0.00], [99.96], [*29.4 %*],
    [14--21 bp], [14 428], [94.79], [5.21], [5.2 %],
    [22--33 bp], [21 696], [97.03], [2.97], [2.6 %],
  ),
  caption: [Residual by overhang length (`bench/results/undertrim-bands.json`).
  25 105 residual of 59 843 carriers = 41.95 %; *94.44 % of it is strictly below
  the 14-base floor*. `undertrim_probe.py` Arm A is a clean ladder: perfect tails
  of 1..13 bases are never trimmed (0/2000 at each t), perfect tails of 14..40
  always are (2000/2000). The cliff is the floor and nothing else.],
) <fig-bands>

== Why the naive fix fails

#figure(
  canvas(length: 1cm, {
    import draw: *
    let c-u = rgb("#c53030")
    let c-f = rgb("#2b6cb0")
    // x: configurations, y-left: residual %, y-right: false-trim % (log)
    let xs = (0.0, 1.6, 3.2, 4.8, 6.4, 8.6)
    let under = (4.195, 3.632, 3.074, 2.563, 2.056, 1.5)
    let fals = (0.414, 1.011, 1.975, 2.263, 3.199, 0.6)
    // axes
    line((0, -0.4), (0, 4.6), stroke: 0.5pt + luma(120))
    line((0, -0.4), (9.6, -0.4), stroke: 0.5pt + luma(120))
    line((9.6, -0.4), (9.6, 4.6), stroke: (paint: luma(120), dash: "dashed", thickness: 0.5pt))
    // series
    line((xs.at(0), under.at(0)), (xs.at(1), under.at(1)), (xs.at(2), under.at(2)),
      (xs.at(3), under.at(3)), (xs.at(4), under.at(4)), stroke: 2pt + c-u)
    line((xs.at(0), fals.at(0)), (xs.at(1), fals.at(1)), (xs.at(2), fals.at(2)),
      (xs.at(3), fals.at(3)), (xs.at(4), fals.at(4)), stroke: (paint: c-f, dash: "dashed", thickness: 1.5pt))
    line((xs.at(4), under.at(4)), (xs.at(5), under.at(5)), stroke: (paint: c-u, dash: "dotted", thickness: 1.5pt))
    line((xs.at(4), fals.at(4)), (xs.at(5), fals.at(5)), stroke: (paint: c-f, dash: "dotted", thickness: 1.5pt))
    for i in range(5) {
      circle((xs.at(i), under.at(i)), radius: 0.075, fill: c-u, stroke: none)
      circle((xs.at(i), fals.at(i)), radius: 0.075, fill: c-f, stroke: none)
    }
    circle((xs.at(5), under.at(5)), radius: 0.11, fill: none, stroke: 1.6pt + c-u)
    circle((xs.at(5), fals.at(5)), radius: 0.11, fill: none, stroke: 1.6pt + c-f)
    content((xs.at(5), 1.05), text(size: 7.8pt)[UT-1 + UT-2 \[E\]], anchor: "north")
    // ticks
    for (i, lab) in ("14", "12", "10", "8", "6", "UT").enumerate() {
      content((xs.at(i), -0.75), text(size: 8pt)[#lab], anchor: "north")
    }
    content((4.8, -1.35), text(size: 8pt)[`--min-adapter-overlap` (shipped 14 leftmost)])
    content((-1.15, 2.1), text(size: 8.5pt)[residual % of carriers], angle: -90deg)
    content((10.75, 2.1), text(size: 8.5pt)[false trims on clean, % (log)], angle: -90deg)
    // legend
    line((-0.9, 5.2), (0.1, 5.2), stroke: 2pt + c-u)
    content((0.25, 5.2), text(size: 8pt)[residual (undertrim)], anchor: "west")
    line((4.4, 5.2), (5.4, 5.2), stroke: (paint: c-f, dash: "dashed", thickness: 1.5pt))
    content((5.55, 5.2), text(size: 8pt)[false trims on adapter-free reads], anchor: "west")
    content((8.0, 4.35), text(size: 7.5pt, fill: luma(90))[dotted = proposed fix \[E\]])
  }),
  caption: [The two-axis trade (`bench/results/undertrim-ab.json`,
  `undertrim-probe.json` Arm D). Lowering the floor walks down the red curve and
  up the blue one: 14$arrow$10 buys 11 points of residual for a 20x false-trim
  increase (0.007 % $arrow$ 0.140 %); 14$arrow$6 leaves 20.6 % residual *and*
  multiplies false trims by 200. The proposed UT-1 + UT-2 regime (hollow markers,
  [E]) is not a threshold change: it keeps the per-read bar at 14 and adds
  independent evidence (mate geometry, G-run retry), targeting $~$15 % residual
  at $<=$0.02 % false trims.],
) <fig-trade>

== The circular guard

The paired-end detector is the instrument that should fix sub-floor tails, and
it works: `detect_pe_overlap` (`pe_overlap.rs:156-219`) pins the insert from
mate complementarity --- a 12-base anchor (1 mismatch), a full-window check
($<=$3 mismatches) and a uniqueness refusal (more than one admissible insert
$->$ nothing, refuse rather than guess). Chance of a false pin is
$~10^(-5) times (3/4)^o$.

But `detect_pe_overlap_trim` then gates the *already pinned* boundary on
`readthrough_is_adapter` (`pe_overlap.rs:275-288`), which calls
`match_adapters(tail, set)` on the tail --- and `match_adapters` refuses every
alignment shorter than `min_overlap_se = 14` (`adapter.rs`: `k < min_overlap
$->$ continue`). A 10-base tail is *provably* unmatchable there, so the
detector's exact boundary is discarded. The decisive A/B
(`bench/results/undertrim-probe-pe-guard.json`): R1 carries the adapter tail,
R2 carries a *random* tail, so any exact cut on R2 can only come from the
paired-end boundary path:

#figure(
  table(
    columns: (auto, auto, auto, auto),
    table.header([*arm*], [*tail t / floor*], [*R1 exact*], [*R2 exact (PE-only)*]),
    [A], [20 bp / 14], [200/200], [*200/200* --- PE path fires],
    [B], [10 bp / 14], [0/200], [*0/200* --- guard's floor blocks it],
    [C], [10 bp / 10], [200/200], [*200/200* --- same geometry, floor lowered],
  ),
  caption: [The detector pins the boundary identically in arms B and C; only the
  guard's re-imposed floor differs. Consistent with
  `undertrim-probe-pe.json`: tails of 3, 6 and 10 bases are 0/500 exact in
  *all* of `auto`/`always`/`never`.],
) <fig-guard>

== Mechanisms and ranked fixes

#figure(
  table(
    columns: (auto, auto, auto, auto),
    table.header([*fix*], [*change*], [*undertrim reduction*], [*overtrim risk*]),
    [UT-1], [PE boundary as primary evidence: verify the pinned tail against
    `adapter[..t]` directly (t $>=$ 2 exact; dual-mate corroboration), bypass
    `min_overlap_se`], [41.95 % $arrow$ *12--18 %* on PE libraries \[E\]:
    bands 6--13 (58.5 %) plus most of 1--5], [false pins refused on ambiguity
    and the anchor; gate t $<=$ 5 on both mates; Cell A `pe_trims` $~$ 0],
    [UT-2], [G-overprint retry: on match failure with a terminal G-run $>=$ 4,
    re-match on the pre-G slice], [Arm C g=6..8: 0 $arrow$ $~$100 % exact;
    3 673 reads (14.6 % of residual)], [the G-run is cut by polyX anyway; full
    14-base bar on the pre-G bases],
    [UT-3], [two-tier short tails (SE, t=6..13): exact 3' prefix *with*
    corroboration (3' quality drop or dataset-level prefix enrichment)],
    [SE half of bands 6--13], [uncorroborated form is Arm D mo=6 (1.29 %
    false trims); corroboration must hold $<$ 0.02 %],
    [UT-4], [quality/N-aware scoring in place of raw mismatch counting (N = 0.5
    mismatch when $>=$ 12 unambiguous bases agree)], [1 185 n-class + part of
    851 sub-class; Arm B e=3 gap], [small; no global `eps` change],
    [UT-5], [adapter-only reads come out *partial* (200/200) --- trim-to-empty
    semantics], [correctness wart], [none],
  ),
  caption: [Ranked fixes. Deliberately rejected: global floor or `eps`
  lowering (figure 2), and `--head-dropout` as an undertrim remedy --- it
  over-trims real sequence, and since D-16 (2026-10-22) the shipped default is
  4, a priced 7:1 choice rather than a free lever. Every fix ships with *both*
  axes measured --- a proposition quoted on undertrim alone is not a
  measurement.],
) <fig-fixes>

= Position against fastp

Three claims, in the order a sceptic should attack them:

+ *A machine-checkable quality contract fastp does not ship.*
  `fasthelper run --fail-on-adapter-residual` exits 2 when adapter survives
  trimming. One command; it either fails or it does not.
+ *Residual adapter from 41.95 % to a projected 12--18 % on paired data*
  (UT-1 + UT-2), with clean-read false trims held at the shipped 0.007--0.02 %.
  Projected, not measured: the pricing A/Bs are named in figure 5.
+ *Less work per pair, less memory:* 16.36 vs 17.64 CPU-s per M pairs and
  101 vs 1141 MB peak \[M\], with the broken one-core wall cell retired ---
  and the utilisation gap being closed by the speed propositions in section 2.

// ---------------------------------------------------------------------------
= Evidence discipline

+ Tiers throughout: [M] measured (artifact named), [E] estimated (pricing A/B
  named), [V] verified on today's source or binary.
+ The "1 core" bench cell is defective (fastprep rows escaped `taskset`):
  quote CPU-seconds and busy cores, never that wall column (optim-dossier item 1).
+ `winners.tsv` and `fastp-comparison.tsv` disagree on 93 of 101 datasets
  (2026-10-21 review): name both or neither; neither is canonical.
+ The 41.95 % is conditional on this synthetic corpus (29.9 % carriers, Q37
  tails; partial cuts counted as undertrim). Real-data curves will move; the
  *shape* --- the sub-floor cliff --- is what generalises.
+ Reproduce: `bench/scripts/undertrim_{probe,probe_pe,probe_detect,ab}.py`,
  artifacts `bench/results/undertrim-*.json`; speed evidence in
  `docs/optim/raw/round2/_salvage-*.md`.
