// ─────────────────────────────────────────────────────────────────────────────
// Benchmark results — fastprep.
//
// Provenance. Every number below is read from the Phase D bench tree; none is
// carried over from a prose summary. Cost is the median over datasets; the trim
// decomposition is a mates-weighted share computed over all 20.4 M mates per
// tool; the per-domain retention and pilot-margin bars are the unweighted mean
// of the 16-17 per-dataset values in the domain.
//
//   bench/results/phase-d/fastp-comparison.tsv          speed, peak RSS (101 ds)
//   bench/results/phase-d/trim-report-20261018.tsv      trim decomposition (tol 0)
//   bench/results/phase-d/trim-retention-20261018.tsv   retention per domain (102 ds)
//   bench/results/phase-d/closed_loop.tsv               PASS/FAIL (306 rows = 102 × 3)
//   bench/results/phase-d/pilot.tsv                     sweep margin (101 ds)
//
// The single source of truth for the numbers in the prose is the block of
// #let bindings immediately below; figures and text both read from it so that a
// correction propagates instead of leaving a stale duplicate in a caption.
// ─────────────────────────────────────────────────────────────────────────────

#import "cetz-charts.typ": *

// --- cost: median over 101 datasets, 200,000 read pairs each -----------------
#let sp-fp     = 0.20   // s, fastprep-winner
#let sp-fastp  = 0.50   // s, fastp-default
#let sp-fastpb = 0.50   // s, fastp-best
#let rss-fp    = 34.3   // MB, fastprep-winner
#let rss-fastp = 1134.7 // MB, fastp-default
#let rss-fastpb = 1134.9 // MB, fastp-best

// --- accuracy: shares of all mates, tolerance 0 ------------------------------
// (exact, overtrim, undertrim-recoverable, undertrim-irreducible, dropped)
#let dec-hd0  = (82.73, 0.57, 0.09,  9.48,  7.13)
#let dec-hd4  = (70.03, 20.92, 0.07, 1.76,  7.21)
#let dec-def  = (86.16, 0.74, 0.08, 10.25,  2.77)
#let dec-tun  = (78.20, 9.33, 0.07,  9.70,  2.70)

// --- retention, % of all mates ----------------------------------------------
#let ret-domains = ("bacterial", "generalised", "human_gut", "metagenomics", "rna_seq", "virome")
#let ret-hd0  = (94.67, 92.09, 97.21, 99.55, 89.10, 84.59)
#let ret-fastp = (92.02, 93.82, 99.99, 99.99, 97.80, 99.74)

// --- closed loop -------------------------------------------------------------
#let pass = 142
#let fail = 164

// --- pilot sweep margin ------------------------------------------------------
#let pilot-gap = (0.112, 0.234, 0.006, 0.005, 0.009, 0.061)

// This file is written to be #include-d from paper.typ. The one cross-reference
// to a label defined there is guarded so the fragment also compiles on its own:
// query() returns empty when the label is absent, and the fallback text stands
// in rather than aborting the build with "label does not exist".
#let seclink(l, fallback) = context if query(l).len() > 0 { ref(l) } else { fallback }

#let c-exact = rgb("#2f6f9f")
#let c-over  = rgb("#b23a48")
#let c-und   = rgb("#7d8f9c")
#let c-drop  = rgb("#c9d2d9")
#let c-pass  = rgb("#3f7d5a")

// Fixed-decimal formatting: calc.round drops trailing zeros, and "0.2 s" in a
// value label reads as a different precision from "0.50 s" beside it.
#let fixed(v, n: 2) = {
  let parts = str(calc.round(v, digits: n)).split(".")
  if n == 0 { parts.at(0) }
  else if parts.len() == 1 { parts.at(0) + "." + "0" * n }
  else { parts.at(0) + "." + parts.at(1) + "0" * calc.max(0, n - parts.at(1).len()) }
}

= Benchmark results <sec:bench>

The measurements below come from one harness run over the 102 real Illumina
libraries of the Phase D corpus, drawn from six domains and truncated to 200,000
read pairs each. The cost and sweep tables cover the 101 libraries on which all
three tools produced a complete run; the closed loop covers all 102. Every tool
reads gzip level 6 with four worker threads on the same AMD EPYC 9654 host, and
fastp is passed the same level explicitly so that the timing does not charge one
tool for compressing harder than the other. Time is wall clock for the whole
process, including decompression and output compression; memory is the peak
resident set, sampled from `/proc/<pid>/status` (`VmHWM`) every 50 ms for the
life of the child rather than taken from `getrusage`, which reports a maximum
over every child the process has ever reaped and therefore attributes one large
run to all later ones.

Four configurations appear. `fastprep-winner` is the per-dataset winner of a
55-point sweep over the trimming parameters — the strongest case the method can
make, and the number a user reaches only by calibrating on their own library.
`fastprep-hd0` is the algorithm at `--head-dropout 0` and `fastprep-hd4` is the
`--head-dropout 4` configuration that shipped while the sweep was run. On the
fastp side, `fastp-default` is version 1.4.0 with no flags beyond the input and
output pairs, and `fastp-best` adds a best-effort tuning pass — `--cut_front`,
`--cut_tail`, `--qualified_quality_phred 20`, `--length_required 30`,
`--low_complexity_filter` — the settings its documentation recommends for
two-colour poly-G libraries. `fastp-best` is present so that a loss for fastp
cannot be read as a misconfigured opponent.

The sweep winner and the fixed configuration are both reported because they
differ by more than run-to-run noise, and averaging that difference away would
hide the only question a prospective user has, which is what happens without
calibration.

Two provenance facts bound these numbers, and neither is a defect of the
harness. First, the run used the build predating the 2026-10 scan change, which
charged a run's own matching steps against $Delta_m$; the current scan frees
them (see the
specification's remark on the mercy budget), so the `fastprep-hd4` decomposition
below is a *lower bound* on accuracy. A 15-library re-measure with the same
harness and the same pinned parameters moved the shipped arm from 76.92 % to
79.19 % exact and from 21.70 % to 19.39 % overtrim; the 101-library table is the
pre-change measurement and was not re-run, so the direction, not the magnitude,
is what a reader should take from that comparison. Second, every arm ran at its
swept length floor — a relative 0.9$L$ across most of the grid — not at the
0.67$L$ floor the tool now ships, so the retention and drop columns describe the
sweep's floors. The 0.67$L$ budget has not been re-measured on the full corpus.

== Cost per library

#figure(
  barchart(
    (
      ("fastprep", sp-fp),
      ("fastp-default", sp-fastp),
      ("fastp-best", sp-fastpb),
    ),
    w: 9.0, h: 3.6, ymax: 0.60,
    ylabel: [seconds],
    fmt: v => fixed(v, n: 2),
  ),
  caption: [Wall time per 200,000 paired-end reads, median over 101 libraries,
    gzip level 6, four threads. The two fastp configurations are
    indistinguishable: their trimming flags change where boundaries fall, not
    how much work the process does.],
) <fig:speed>

#figure(
  barchart(
    (
      ("fastprep", rss-fp),
      ("fastp-default", rss-fastp),
      ("fastp-best", rss-fastpb),
    ),
    w: 9.0, h: 3.6, ymax: 1300,
    ylabel: [peak RSS (MB)],
    fmt: v => fixed(v, n: 1),
  ),
  caption: [Peak resident set size per process, same runs as @fig:speed.
    fastprep's footprint is set by its chunk buffers and lookup tables;
    fastp's is set by its read-batch and duplication-detection structures, which
    the trimming flags do not touch, so both fastp bars coincide.],
) <fig:mem>

At this unit size the median library costs #fixed(sp-fp, n: 2) s with fastprep
against #fixed(sp-fastp, n: 2) s with fastp, a factor of #fixed(sp-fastp / sp-fp, n: 1),
and #fixed(rss-fp, n: 1) MB of peak RSS against #fixed(rss-fastp, n: 1) MB, a factor
of #fixed(rss-fastp / rss-fp, n: 0). Absolute times are small because the unit is
small; per pair the medians are 1.0 µs and 2.5 µs, which is the figure that
survives a change of unit. The memory ratio is the more consequential of the
two. A gigabyte per process bounds how many libraries a node can process
concurrently, and it is what makes a per-core deployment of fastp impractical
where #fixed(rss-fp, n: 0) MB is not.

The time ratio is a ratio of constant factors rather than of asymptotics: both
tools are linear in input size, so the ratio is stable across library size,
which is what makes a 200,000-pair unit representative rather than convenient.
The ratio is not attributed here to a single component. The scan loop of
#seclink(<sec:rust-perf>, [the implementation chapter]) is one candidate, but
parse, adapter and output paths have not been profiled separately, and the
component breakdown that would settle it is listed there as outstanding work.

== Where the trim boundary lands

Speed and memory say nothing about whether the boundary is in the right place.
@fig:trim-decomp decomposes every mate into exactly one of five outcomes against
the simulator's ground truth at tolerance 0: written with the exact boundary,
over-trimmed (real bases removed), under-trimmed in a way the measured posterior
says is recoverable, under-trimmed irreducibly, or dropped before scoring. The
shares are of *all* mates, so a tool's five outcomes sum to 100 and the drop
column is comparable with the rest. Recoverable under-trim is under 0.1% for
every tool and is omitted from the bars; it is the smallest column in the
decomposition and the one the project spent the most effort trying to move.

#figure(
  groupedbars(
    ("fastprep-hd0", "fastprep-hd4", "fastp-default", "fastp-tuned"),
    (
      ("exact", c-exact, (dec-hd0.at(0), dec-hd4.at(0), dec-def.at(0), dec-tun.at(0))),
      ("overtrim", c-over, (dec-hd0.at(1), dec-hd4.at(1), dec-def.at(1), dec-tun.at(1))),
      ("underrun", c-und, (dec-hd0.at(3), dec-hd4.at(3), dec-def.at(3), dec-tun.at(3))),
      ("dropped", c-drop, (dec-hd0.at(4), dec-hd4.at(4), dec-def.at(4), dec-tun.at(4))),
    ),
    w: 9.0, h: 4.0, ymax: 100,
    fmt: v => fixed(v, n: 1),
  ),
  caption: [Fate of every mate at tolerance 0, as a share of all mates, so each
    tool's four bars plus the omitted recoverable column sum to 100.
    `underrun` is under-trim whose ground-truth artefact the observable signal
    does not identify: N-runs and adapter-bearing reads outside the posterior.
    Recoverable under-trim, omitted because it never exceeds 0.09%, is 0.087%
    (hd0), 0.071% (hd4), 0.079% (fastp-default) and 0.066% (fastp-tuned).],
) <fig:trim-decomp>

Read across the exact bars, the algorithm at `--head-dropout 0` is at parity
with fastp: 82.7 against 86.2 of all mates, which is 89.1 against 88.6 of kept
mates — the denominator usually quoted, and the one that flatters fastprep
because it drops more. The informative difference is in the columns either
side. fastprep-hd0 leaves 0.8 points less irreducible under-trim than fastp
(9.5 against 10.3) but discards 4.4 points more of the library (7.1 against
2.8), and by the measurement below the extra drops are duplicates that the
dedup default removes and fastp's default retains, not reads the trimmer failed
on — the ledger of @fig:retention rather than a trimming decision. Either way the trade is real, and it is the reason the retention
figure of @fig:retention has to be read beside every exact rate in this section.

Two failure modes are visible in the same figure and neither is the algorithm's
ceiling. `fastprep-hd4`, the configuration that shipped during the sweep, is
worst on the axis that matters most: 20.9% of mates over-trimmed, because a
leading G-run of length 1 to 4 is trimmed unconditionally with no test of
artefact against biology. On a library where a quarter of reads begin with a
genuine G, that rule removes real sequence; the 20.9 points are what
`--head-dropout 0` avoids. It is not the shipped default: the 0.5.0 release
ships 4 again, because an artefact left in the output costs 7x a base trimmed
(D-16), which puts the break-even at 2.6. Re-derived on the two 2026-10-26 fixed
builds the same configuration measures 77.87 / 20.00 / 2.13 and then 77.88 /
20.00 / 2.12 at tolerance 0 (`bench/results/fix2026/` and
`bench/results/wave2/`), so the ordering, and the argument built on it, are
unchanged: the shipped rule is still the worst of the fastprep configurations,
just less so. The second campaign's movement is confined to the under-trim row
(1,690 mates moved to an exact cut, 433 to an over-cut, 20 now dropped instead of
emitted with an untrimmed tail), which is the artefact that used to survive a 3'
cut when a 5' cut had spent the trim budget. `fastp-tuned` shows the mirror image of
the same mistake: its aggressive `--cut_front`/`--cut_tail` raise over-trim from
0.7 to 9.3 points while lowering exactness by 8 points, and unlike boundary
jitter the
over-trim does not decay when the tolerance is relaxed to 2.

== Retention

#figure(
  groupedbars(
    ret-domains,
    (
      ("fastprep-hd0", c-exact, ret-hd0),
      ("fastp-default", c-over, ret-fastp),
    ),
    w: 9.0, h: 4.0, ymax: 110,
    fmt: v => fixed(v, n: 1),
  ),
  caption: [Mates retained as a share of all mates, by domain, at tolerance 0.
    `bacterial` abbreviates `bacterial_isolate`; each domain contributes 17
    libraries. fastprep keeps more only on bacterial isolates; the deficit
    elsewhere is concentrated in virome and RNA-seq.],
) <fig:retention>

fastprep-hd0 retains 92.87% of all mates against fastp's 97.23%, and the deficit
is not spread evenly. It keeps *more* on bacterial isolates (94.7 against 92.0)
and loses 15.2 points on virome (84.6 against 99.7) and 8.7 on RNA-seq. Since
the exact rate of @fig:trim-decomp is computed over kept mates, the virome
comparison is measured on a kept-set roughly a sixth smaller, and the apparent
fastprep advantage there is the least apples-to-apples number in the section.
The bacterial isolate row is the opposite case and the one worth trusting: both
tools keep a comparable set and tie on exactness.

The aggregate deficit is a contract difference, not a trimming one, and it is
attributable to a single stage. With deduplication disabled fastprep retains
98.47% against fastp's 97.23%, and its exact rate over kept mates moves by 0.4
points, from 89.09 to 88.69, which leaves it within 0.1 of fastp's 88.61. The
reads it gives up relative to fastp are therefore duplicates: the harness
charges them to retention because a dropped read leaves the scored set. The shipped default keeps deduplication on, so both
configurations are quoted here rather than the flattering one alone.

== Does the simulator close the loop?

The accuracy figures above are scored against fastsim's ground truth, which is
only meaningful if fastsim reproduces the library it is imitating. The closed
loop tests that directly: each real library is profiled, re-simulated in three
variants, and fastprep's QC report on the simulated copy must match the real
one within per-gate tolerances on GC fraction, duplication ratio, poly-G
fraction, mean quality and per-position base composition. A PASS means the
trimmer's output is indistinguishable, gate by gate, between the real library
and a simulated one built to match it.

#figure(
  barchart(
    (
      ("PASS", pass, c-pass),
      ("FAIL", fail, c-over),
    ),
    w: 9.0, h: 3.6, ymax: 200,
    ylabel: [comparisons],
    fmt: v => str(v),
  ),
  caption: [Closed-loop verdicts over 306 comparisons: 102 libraries times three
    simulated variants, each gated against the real library's QC report.
    #fixed(100 * pass / (pass + fail), n: 1)% of comparisons close. Failures are
    dominated by the per-position base-composition gate, not by trimming.],
) <fig:closed-loop>

#fixed(100 * pass / (pass + fail), n: 1)% of comparisons close. The residue is
concentrated rather than diffuse: 150 rows fail the base-composition gate at a
median excess of 0.086, and those rows are overwhelmingly metagenomics and
RNA-seq, the two domains whose position-dependent composition a single mean
profile cannot represent. A further block fails the duplication gate on
libraries where the simulator's per-pair draw rate and the read-level measured
rate cannot be reconciled by any constant factor. Neither block implicates the
trim boundary, which is why the accuracy decomposition above is scored against
simulated ground truth without a closed-loop caveat: the simulator's deficiency
is in base composition, not in where it places artefacts.

== Margin over the shipped default

The sweep exists to answer whether the parameter surface is worth navigating.
@fig:pilot gives the per-domain mean of `gap_vs_default`, the difference in the
sweep's composite score between the per-dataset winner and the default
configuration.

#figure(
  barchart(
    (
      ("bacterial", pilot-gap.at(0)),
      ("generalised", pilot-gap.at(1)),
      ("human_gut", pilot-gap.at(2)),
      ("metagenomics", pilot-gap.at(3)),
      ("rna_seq", pilot-gap.at(4)),
      ("virome", pilot-gap.at(5)),
    ),
    w: 9.0, h: 3.6, ymax: 0.28,
    ylabel: [score gap],
    fmt: v => fixed(v, n: 3),
  ),
  caption: [Mean composite-score margin of the per-dataset sweep winner over the
    default configuration, by domain, over 101 libraries. The winner is ahead on
    100 of 101 datasets; the overall mean margin is 0.071 and the median 0.005,
    so the mean is carried by a small number of large wins rather than by a
    uniform advantage.],
) <fig:pilot>

The winner beats the default on 100 of 101 libraries, which sounds decisive and
is not. The median margin is 0.005 against a mean of 0.071: most libraries gain
almost nothing from calibration and a handful gain a great deal, with a maximum
of 0.90. The distribution is bimodal in domain terms: generalised (0.23) and bacterial
isolate (0.11) carry almost all of the mean, human gut, metagenomics and RNA-seq
sit at or below 0.01, and virome at 0.06 is the only domain in between. That
matches where the sweep is in fact choosing between two good configurations
rather than between a good and a bad one. A user with a human gut library should expect the default to be within
noise of the best configuration the sweep can find, and should not read the
100-of-101 headline as a reason to calibrate.

The honest summary of this section is therefore narrower than the figures
suggest at a glance. fastprep is #fixed(sp-fastp / sp-fp, n: 1)× faster and
#fixed(rss-fastp / rss-fp, n: 0)× lighter than fastp with the same asymptotic
cost, and at `--head-dropout 0` it matches fastp's exactness at tolerance 0
while dropping more reads. It does not dominate on boundary accuracy, it does
not remove the irreducible under-trim that both tools leave, and its advantage
over its own default is real on two domains out of six.
