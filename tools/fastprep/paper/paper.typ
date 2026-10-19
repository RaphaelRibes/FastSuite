// PolyX trimming — formal specification + Rust implementation guide
// Revised 2026-04-23: asymmetric per-step score, clean indexing,
//                     separated end-scan / interior-scan, tuned naming,
//                     whole-pipeline complexity, captioned figures.
// Extended 2026-05-26: Introduction (§1) and Rust implementation
//                     parts (§8--§10) bridging math to production code.

#set page(paper: "a4", margin: 2cm)
#set text(font: "New Computer Modern", size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.")

#show figure.caption: set text(size: 10pt)
#show figure.where(kind: table): set figure.caption(position: top)

#import "@preview/codly:1.3.0": *
#show: codly-init.with()
#codly(zebra-fill: luma(245), stroke: 0.5pt + gray)

#let Prop = text(weight: "bold", "Proposition")
#let Def  = text(weight: "bold", "Definition")
#let Rem  = text(weight: "bold", "Remark")

#set page(numbering: "1")

#align(center)[
  #v(1.5cm)
  #text(size: 18pt, weight: "bold")[
    Asymmetric quality-weighted polyX trimming
    for short-read sequencing data
  ]

  #v(3mm)
  #text(size: 12pt)[A specification and its Rust implementation]

  #v(9mm)
  #text(size: 12pt)[Raphaël Ribes]

  #v(2mm)
  #text(size: 10pt)[`raphael.ribes@gmail.com`]

  // TODO: add affiliation before circulating.

  #v(3mm)
  #text(size: 10pt)[#datetime.today().display("[day] [month repr:long] [year]")]
]

#v(10mm)

#block(width: 100%, inset: (x: 1.4cm))[
  #set text(size: 10pt)
  #set par(justify: true)
  *Abstract.* Homopolymer runs, adapter readthrough and whole-read
  quality decay are the three artefact classes that short-read
  preprocessors must remove before any downstream analysis. Existing
  tools treat them with independent heuristics, and in particular
  decide where a homopolymer run ends without consulting the basecall
  qualities that determine whether an interrupting base is real. We
  specify a trimming procedure built on a single asymmetric
  quality-weighted per-step score, in which matches are strong
  evidence of a run almost independently of quality while mismatches
  count only when the basecalls are reliable. Threshold crossings of
  the resulting cumulative gap define both the trim length and the
  drop decision, in $O(L + Delta_L)$ time per read with an
  $O(Delta_L)$ early-drop path. Because the score is a path
  functional and therefore not skippable, we pair it with a
  quality-blind composition prefilter that provably over-approximates
  its stopping step, so the expensive decision runs only on reads
  where the boundary is genuinely ambiguous. The document is in two
  parts: a formal specification with proofs, and a description of the
  Rust implementation in which every object of the specification has
  a named counterpart.
]

#v(10mm)

#outline(indent: auto, depth: 2)

#pagebreak()

= Introduction <sec:intro>

== Problem and motivation

Next-generation sequencing reads carry three distinct classes of low-information
artefacts that must be removed before any downstream analysis (variant calling,
assembly, expression quantification):

#set enum(numbering: "(a)")
+ _Homopolymer runs at read ends._ Two-colour Illumina chemistry (NextSeq,
  NovaSeq, MiniSeq) reports a G whenever both colour channels read dark, so
  any signal loss at the 3′ end of a read produces a polyG
  tail@chen2018fastp. PolyA and
  polyT contamination arises from poly(A) tail readthrough in RNA-Seq, from
  oligo-dT priming artefacts, and from low-complexity regions in genomic
  libraries.
+ _Adapter readthrough._ When the cDNA or DNA insert is shorter than the
  read length, the sequencer reads past the insert into the ligated adapter.
  The adapter sequence is a constant per library prep kit, but its prefix
  appears at variable positions in the read depending on insert size.
+ _Whole-read quality decay._ Reads with low overall Phred quality@ewing1998phred
  mislead downstream tools and are best discarded entirely.

Every widely used preprocessor addresses all three, and does so with a
separate heuristic for each: fixed-length end trimming, sliding-window
quality averages@bolger2014trimmomatic, mismatch-rate thresholds on an
adapter alignment@martin2011cutadapt, and counted homopolymer runs with a
fixed error allowance@chen2018fastp. Each works adequately on the case it
was designed for. None of them lets the basecall qualities participate in
the decision of _where a homopolymer run ends_, which is the decision this
document is about.

== Existing tools <sec:related>

Four tools define the current practice, summarised in @tab:tools.

Quality profiling itself is conventionally delegated to a separate tool,
`FastQC`@andrews2010fastqc, and is not discussed further here.

`Trimmomatic`@bolger2014trimmomatic is the oldest of the four still in
routine use. It trims adapters by seed-and-extend with a palindrome mode
for paired-end reads, and controls quality with a sliding-window mean
applied as a separate stage. It has no homopolymer facility, so polyG
tails from two-colour chemistry pass through unless they happen to fail a
quality window.

`cutadapt`@martin2011cutadapt models adapter removal as a semi-global
alignment under unit edit costs@sellers1980, which admits insertions and
deletions and makes it the most general of the four. That generality is
paid for per read: the alignment is a banded dynamic program rather than a
byte comparison. Poly-A trimming is available, but as a separate operation
with its own error tolerance, not as part of the adapter decision.
`Atropos`@didion2017atropos is a fork that retains the model and improves
the implementation.

`fastp`@chen2018fastp@chen2023fastp@chen2025fastp1 is the throughput
reference, and the tool this work must be measured against. It performs
quality control, adapter trimming, quality filtering and polyG/polyX
trimming in a single pass over the data, is written in C++ with genuine
multithreading, and reports itself as two to five times faster than
`Trimmomatic` or `cutadapt`. Its polyX procedure counts run bases from the
read end under a minimum-length threshold and a mismatch allowance, both
of which are quality-independent: a mismatching base interrupts a run by
the same amount whether it was called with high confidence or is
indistinguishable from noise.

#figure(
  table(
    columns: 5,
    stroke: 0.5pt,
    align: (left + horizon,) * 5,
    table.header(
      [*Tool*], [*polyX*], [*Adapter model*],
      [*Quality in trim boundary*], [*Language*],
    ),
    [`Trimmomatic`], [none], [seed and extend, palindrome mode for PE],
    [no], [Java],

    [`cutadapt`], [separate poly-A operation],
    [semi-global alignment, unit edit costs], [no], [Python/Cython],

    [`fastp`], [counted run, fixed allowance],
    [leftmost sequence match; PE overlap], [no], [C++],

    [this work], [asymmetric quality-weighted score],
    [suffix-anchored match; PE overlap], [yes], [Rust],
  ),
  caption: [Homopolymer and adapter handling across current
    preprocessors. The fourth column asks whether basecall quality enters
    the decision of where a trim boundary falls, as opposed to being
    applied as an independent filter before or after trimming. Only the
    procedure specified here answers yes.],
) <tab:tools>

== Why another trimmer <sec:why>

The gap in @tab:tools is narrow but consequential. A polyG tail produced
by signal dropout is, by construction, a region of collapsing quality: the
same physical failure that causes the instrument to call G causes the
basecall to be unreliable. A mismatching base inside such a tail is
therefore much more likely to be an error than a genuine end of the run,
and treating it as a hard interruption truncates the trim early and leaves
artefact in the read. Conversely, a high-quality mismatch in a
high-quality region is real evidence and should stop the scan. A
counted-run procedure cannot express the difference, because it has
discarded the only signal that distinguishes the two cases.

This motivates four design commitments, which the rest of the document
develops:

#set enum(numbering: "(a)")
+ _Quality participates in the boundary._ A single asymmetric per-step
  score weights matches and mismatches differently as a function of
  quality (@sec:intro-approach), and its threshold crossings define both
  the trim length and the drop decision.
+ _The procedure is specified, not merely implemented._ Part I states the
  algorithm with proofs of its complexity, its budget invariant and the
  correctness of its reconciliation step, so that behaviour is a property
  of the specification rather than of a particular revision of the source.
+ _Every read is auditable._ A compact binary record per read
  (@sec:pxmt) attributes every trimmed base to adapter or polyX, at the 5′
  or 3′ end, with the drop reason where applicable. Aggregate counters
  cannot answer why a specific read was lost; this format can.
+ _Throughput is a correctness-preserving concern._ The quality-weighted
  score is a path functional and cannot be skipped, which makes it
  intrinsically more expensive per base than a counted run. The cascade of
  @sec:prefilter resolves this by placing a quality-blind test in front of
  it that provably never rejects a read the full procedure would have
  trimmed, so the cost is paid only on reads where the boundary is
  ambiguous.

+ _Stages are fused, not chained._ Every separate tool in a
  preprocessing chain costs one decompress, parse, write and compress
  round trip, and on compressed input that round trip dominates
  wall-clock time. PCR duplicate removal is therefore specified here as
  a first-stage filter (@sec:dedup) rather than delegated, on the same
  reasoning that led `fastp` to fuse its own stages. What is _not_
  inherited from `fastp` is the probabilistic representation: the
  fingerprint set here is exact, so the false-duplicate rate is bounded
  by a birthday argument on the hash width rather than by filter
  occupancy.

The fourth commitment is the load-bearing one for adoption. A trimmer that
is more principled but slower than `fastp` will not be used on a
production lane, and no argument about evidence models changes that.

== High-level approach <sec:intro-approach>

The procedure specified in this document is built around a single asymmetric
per-step score that distinguishes the evidence asymmetry between matches and
mismatches in a homopolymer run:

- _Matches are strong evidence_ of a polymer extension, almost independently
  of quality. Random substitution errors that coincidentally produce identical
  bases are improbable.
- _Mismatches are informative only when reliable._ A low-quality mismatch may
  correspond to a sequencing error inside a genuine run and should contribute
  little.

This asymmetry is encoded as two parameters $mu, nu in [0, 1]$ scaling the match
and mismatch branches independently. The resulting algorithm produces a single
scalar (the cumulative gap $Gamma$) whose threshold crossings define both the
trim length and the dropping decision. The whole pipeline runs in
$O(L + Delta_L)$ time per read, with an $O(Delta_L)$ early-drop fast path
when the 5′ end is overwhelmed by polyX.

== Pipeline overview

Each read $(R, Q)$ passes through the following stages, in order:

#set enum(numbering: "1.")
+ Deduplication (optional): fingerprint the read, or read pair, over an
  artefact-free anchor window and drop it if that fingerprint has already
  been seen (@sec:dedup).
+ Phred decode: subtract the ASCII offset (33 or 64) from each quality byte.
+ Hard trim (optional): unconditionally remove $T_5^"hard"$ bases from the
  5′ end and $T_3^"hard"$ from the 3′ end before any scanning.
+ Adapter trim: detect adapter contamination via sequence matching
  (single-end) or PE overlap geometry (paired-end), producing
  $(T_5^"adp", T_3^"adp")$.
+ Mean-quality filter: drop the read if $overline(q) < theta_q$, where
  $overline(q)$ is the mean Phred over the full read.
+ PolyX 5′ end-scan: walk from the new 5′ end accumulating the cumulative gap
  $Gamma$ until either a mismatch threshold (mercy clause) or the end of the
  read (boundary clause) is reached. Resolves into $T_5^"px"$.
+ PolyX interior scan: continue accumulating interior homopolymer mass
  $H_("int")$ from individual mercy-terminated qualifying runs in the read body.
+ PolyX 3′ end-scan: a symmetric reverse scan, capped by the remaining trim
  budget, producing $T_3^"px"$.
+ Reconciliation: remove any interior run that overlaps the resolved 3′ trim
  zone (avoids double-counting).
+ Decision rule: drop the read if $T_5 + T_3 + H_("int") > Delta_L$, where
  $Delta_L = L - L_min$ is the trim budget. Otherwise, emit the kept slice
  $R[T_5 : L - T_3]$ with its qualities.

== Document structure

The remainder of this document is structured in two parts. Part I (§2--§7)
gives the formal specification: the asymmetric per-step score, the stopping
rule, the end-scan and interior-scan procedures, the adapter sub-procedures,
and the EPYC-specific optimisations. Part II (§8--§10) bridges the
specification to the Rust implementation: the reproducible build environment,
the type-level encoding of every mathematical object, and the
performance-critical code paths with their micro-optimisations
(quality-gate LUT, fused streaming walk, const-generic monomorphisation).

= PolyX trimming

We specify a procedure for the detection and removal of homopolymer
(polyX) runs at either end of a next-generation sequencing read, with
optional accounting for interior runs. The procedure is built around
an asymmetric, quality-weighted per-step score that distinguishes
evidence asymmetries between matches and mismatches.

== Input

A _read_ is a pair $(R, Q)$ consisting of a base string
$R = (b_0, b_1, dots, b_(L-1)) in {A, C, G, T}^L$ and an aligned
Phred quality vector $Q = (q_0, q_1, dots, q_(L-1)) in NN^L$, where
$L := |R|$ denotes the read length. Indices are $0$-based: $b_(L-1)$
is the terminal base. Bases are assumed to be $2$-bit encoded.

== Parameters

The procedure is configured by the tuple
$(L_min, Delta_m, t, x, mu, nu, m_h)$:

#set enum(numbering: "(a)")
+ $L_min in NN$, $L_min <= L$ — minimum acceptable post-trim length.
  We write $Delta_L := L - L_min$ for the _trim budget_.
+ $Delta_m in RR_(>=0)$ — _mercy budget_ of the stopping rule.
+ $t in NN_(>0)$ — Phred saturation point of the quality gate.
+ $x in RR_(>0)$ — shape parameter of the quality gate.
+ $mu in [0, 1]$ — _match leniency_, controlling the quality
  sensitivity of match contributions.
+ $nu in [0, 1]$ — _mismatch strength_, controlling the magnitude
  of high-quality mismatch contributions.
+ $m_h in NN_(>=2)$ — minimum homopolymer length (in bases) required
  for a run to be recognised.

== Quality gate

#Def (quality gate). The _quality gate_ $g : RR -> [0, 1]$ is
$ g(q) := cases(
  0 &"if" q <= 0,
  (q \/ t)^x quad &"if" 0 < q < t,
  1 &"if" q >= t.
) $

$g$ is non-decreasing, continuous on $RR$, and satisfies $g(t) = 1$.
The shape parameter $x$ tunes severity: $x < 1$ admits sub-threshold
scores with substantial weight; $x > 1$ concentrates weight near
$q = t$.

== Match indicator and joint quality

Fix a scan starting position $s in {0, 1, dots, L - 2}$. For each
step $k in NN$ with $s + k + 1 <= L - 1$, define

$ M_s (k) := cases(
  +1 quad &"if" b_(s+k) = b_(s+k+1),
  -1 &"otherwise",
) quad quad
W_s (k) := (g(q_(s+k)) + g(q_(s+k+1))) / 2. $

Then $M_s (k) in {-1, +1}$ and $W_s (k) in [0, 1]$.

#Rem (anchored variant). When a specific base $X in {A, C, G, T}$ is
targeted (e.g. polyG for two-colour Illumina chemistry), replace $M_s$
by
$ M_s^X (k) := cases(
  +1 quad &"if" b_(s+k+1) = X,
  -1 &"otherwise".
) $
All subsequent definitions and results apply mutatis mutandis.

== Asymmetric per-step score

#Def (per-step score). For each admissible $k$, the _per-step score_
$sigma_s : NN -> RR$ is
$ sigma_s (k) := cases(
  (1 - mu) + mu dot W_s (k) quad &"if" M_s (k) = +1,
  -nu dot W_s (k) &"if" M_s (k) = -1.
) $

The match branch maps to $[1 - mu, 1]$, the mismatch branch to
$[-nu, 0]$.

#Rem (asymmetry rationale). A homopolymer match between adjacent
bases is strong evidence of a polymer run nearly independently of
basecall quality, since random substitution errors coincidentally
producing identical bases are improbable. A homopolymer mismatch, by
contrast, is informative only when the basecalls are reliable:
low-quality mismatches may correspond to sequencing noise within a
genuine run and should contribute negligibly. The parameters
$(mu, nu)$ encode this asymmetry; the default values adopted here are
shown in @tab:asymmetry.

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    align: (left + horizon, center + horizon, center + horizon),
    table.header(
      [], [*$W = 1$* (high quality)], [*$W = 0$* (low quality)],
    ),
    [*$M = +1$* (match)], [$+1.0$], [$+0.9$],
    [*$M = -1$* (mismatch)], [$-0.5$], [$0.0$],
  ),
  caption: [Default per-step score as a function of match indicator
    $M$ and joint quality $W$, corresponding to $mu = 0.1$ and
    $nu = 0.5$. Matches contribute positively with weak quality
    dependence; mismatches contribute negatively only when backed
    by high-quality basecalls.],
) <tab:asymmetry>

== Cumulative score and cumulative gap

#Def (cumulative score). For $K in NN_(>=1)$ with $s + K <= L - 1$,
$ cal(P)_s (K) := sum_(k=0)^(K-1) sigma_s (k). $

#Def (cumulative gap). Equivalently, the _cumulative gap_ is
$ Gamma_s (K) := K - cal(P)_s (K) = sum_(k=0)^(K-1) (1 - sigma_s (k)). $

Per-step gap contributions are
$ 1 - sigma_s (k) = cases(
  mu (1 - W_s (k)) quad &"match",
  1 + nu dot W_s (k) &"mismatch",
) $
hence in $[0, mu]$ for matches and in $[1, 1 + nu]$ for mismatches.

#Prop (bounds). For every admissible $s$ and $K >= 1$,
$ -nu K <= cal(P)_s (K) <= K, quad quad 0 <= Gamma_s (K) <= (1 + nu) K. $
Equality $cal(P)_s (K) = K$ holds iff every pair is a high-quality
match; equality $cal(P)_s (K) = -nu K$ holds iff every pair is a
high-quality mismatch.

_Proof._ Immediate from $sigma_s (k) in [-nu, 1]$, and the equality
conditions follow from the definitions of $M_s$ and $W_s$. $qed$

#Rem. The cumulative gap is monotonically non-decreasing, which is
the property the stopping rule exploits.

== Stopping rule

#Def (stopping step). Given a starting position $s$ and an upper
index $u in {s + 1, dots, L - 1}$, the _stopping step_ $K^*_s (u)$
is the minimal $K >= 1$ satisfying
$ Gamma_s (K) > Delta_m quad quad (text("mercy clause")) $
or
$ s + K = u quad quad (text("boundary clause")). $

We write $K^*_s := K^*_s (L - 1)$ when the upper index is the end of
the read.

_Existence and termination._ The scan runs over a finite window, so
the boundary clause guarantees termination; mercy fires earlier on
sequences dominated by mismatches. Hence $K^*_s (u)$ is well-defined
and satisfies $1 <= K^*_s (u) <= u - s$.

We additionally record the _stopping mode_
$cal(M)_s (u) in {"mercy", "boundary"}$ indicating which clause fired
(mercy, if both are simultaneously satisfied at the same $K$), and
the _match count_
$ m^*_s (u) := |{k in {0, dots, K^*_s (u) - 1} : M_s (k) = +1}|, $
the number of matching steps the scan executed. The corresponding
_run length_ is $m^*_s + 1$, accounting for the anchor base
$b_s$ and one neighbour per matching step. The match count is
necessary because $K^*$ alone conflates run length with the mercy
overshoot: a scan that mercy-terminates after a short run followed by
several mismatches has a $K^*$ much larger than the run itself.

== End-scan procedure

The end-scan produces the trim length at one end of the read.

#Def (5′ trim length). Let $(K^*_0, cal(M)_0, m^*_0)$ be the stopping
step, mode, and match count of the scan with $s = 0$, $u = L - 1$.
Define
$ T_5 := cases(
  0 &"if" m^*_0 + 1 < m_h,
  K^*_0 &"if" cal(M)_0 = "mercy",
  K^*_0 + 1 &"if" cal(M)_0 = "boundary".
) $

#Def (3′ trim length). Apply the same procedure to the reversed read
$tilde(R) := (b_(L-1), b_(L-2), dots, b_0)$ with matching reversed
qualities, and let $(tilde(K)^*_0, tilde(cal(M))_0, tilde(m)^*_0)$ be
its stopping step, mode, and match count. Then
$ T_3 := cases(
  0 &"if" tilde(m)^*_0 + 1 < m_h,
  tilde(K)^*_0 &"if" tilde(cal(M))_0 = "mercy",
  tilde(K)^*_0 + 1 &"if" tilde(cal(M))_0 = "boundary".
) $

#Rem (boundary correction). The mercy clause stops one base _past_
the homopolymer run (so that $b_(K^*)$ is the non-matching base and
$b_0, dots, b_(K^* - 1)$ are trimmed), while the boundary clause
stops within the run and requires the $+1$ correction to include the
final base.

== Interior-scan procedure <sec:interior>

After the end-scans, interior homopolymer mass is accumulated over
the residual window $[T_5, L - 1 - T_3]$. Let
$u_("int") := L - 1 - T_3$.

#Def (interior iteration). Initialise $s_0 := T_5$ and
$H_("int") := 0$. For $i = 0, 1, 2, dots$ while
$s_i <= u_("int") - 1$:

+ Compute $(K^*_(s_i)(u_("int")), cal(M)_(s_i)(u_("int")), m^*_(s_i)(u_("int")))$.
+ If $cal(M)_(s_i)(u_("int")) = "mercy"$ and
  $m^*_(s_i) + 1 >= m_h$:
  $ H_("int") <- H_("int") + K^*_(s_i) . $
+ Advance $s_(i+1) := s_i + K^*_(s_i) + 1$.

#Rem (interpretation). Only mercy-terminated scans contribute to
$H_("int")$: a scan that stops by the interior boundary $u_("int")$
has exhausted the window without a mismatch strong enough to
terminate, and its mass is already accounted for by $T_3$ through
the 3′ end-scan.

== Decision rule

Given $(T_5, T_3, H_("int"))$, define
$T_("tot") := T_5 + T_3 + H_("int")$.

- If $T_("tot") > Delta_L$: _drop_ the read.
- Otherwise: _emit_ $(b_(T_5), dots, b_(L - 1 - T_3))$ with their
  matching qualities.

#block(inset: 0.6em, stroke: 0.5pt + gray, radius: 2pt, [
  _Sanity check._ On a read with no homopolymer at a given end,
  $M = -1$ dominates; each mismatch contributes at least $1$ to
  $Gamma$; the mercy clause fires after at most $floor(Delta_m) + 1$
  steps, with $m^* = 0$ matches accumulated, so $m^* + 1 = 1 < m_h$
  and $T = 0$. On a read that is purely $X$ at that end, every match
  contributes at most $mu$ to $Gamma$; the mercy clause does not
  fire; the boundary clause triggers at $K^* = L - 1$ with
  $m^* = L - 1$ matches; provided $L >= m_h$, the threshold
  $m^* + 1 >= m_h$ is satisfied and $T = L$, yielding near-certain
  drop at the decision rule.
])

== Bidirectional mode and direction assignment

Because polyX runs at the 5′ and 3′ ends correspond to distinct
library artefacts (e.g., polyG from two-colour chemistry at the 3′
end; polyA/T tails at either end), the two end-scans are independent
and may be configured separately. In particular, the direction of
the 5′ scan for one mate ($R_1$ or $R_2$) may be reversed to
prioritise expected artefacts:

- 5′→3′ on $R_2$ prioritises $5'$-end polyG of the Nextera protocol;
- 3′→5′ on $R_1$ prioritises polyA/T of the 3′ end.

When the library is well-characterised, either end-scan (or the
interior scan) may be disabled without affecting the correctness of
the remaining components.

== Complexity <sec:complexity>

We count the number of base-pair evaluations $B$ performed by the
combined 5′ end-scan, interior scan, 3′ end-scan, and decision rule
on a single read of length $L$. Space is measured exclusive of the
read buffer and the precomputed lookup tables.

#Prop (whole-pipeline complexity). $B in O(L + Delta_L)$.

_Proof._ The 5′ end-scan and interior scan together traverse the
prefix $[0, L - 1 - T_3]$ once; the advance rule of @sec:interior
increments $s_i$ strictly, so
each base participates in $O(1)$ pair evaluations. The 3′ end-scan
revisits at most $T_3 <= Delta_L - T_5 <= Delta_L$ bases from the
opposite end. The decision rule is $O(1)$. Summing gives the
bound. $qed$

#Prop (best case: early drop). Suppose the 5′ end-scan encounters a
homopolymer run spanning $[0, Delta_L - 1]$ and the budget
short-circuit (@sec:shortcircuit) is active. Then the drop condition
$T_5 > Delta_L$ is detected at step $K = Delta_L$, the decision rule
fires, and neither the 3′ end-scan nor the interior scan is invoked.
Total cost is $B = Delta_L$, independent of $L$.

_Proof._ Direct from the short-circuit: once $T_5 > Delta_L$, the
drop decision is determined regardless of $T_3$ and $H_("int")$, so
those computations are skipped. $qed$

The two propositions combine into four operational regimes,
summarised in @tab:regimes.

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    align: (left + horizon, center + horizon, left + horizon),
    table.header(
      [*Read content*], [*Cost per mate*], [*Typical trigger*],
    ),
    [5′ polyX run $>= Delta_L$],
    [$O(Delta_L)$],
    [doomed polyA/polyG at 5′ — early drop],

    [No homopolymer at either end],
    [$L + O(m_h + Delta_m)$],
    [clean insert, read kept],

    [Short interior runs only],
    [$O(L)$],
    [typical kept read with low-complexity patches],

    [3′ polyX present],
    [$L + O(Delta_L)$],
    [NextSeq polyG, kept or dropped],
  ),
  caption: [Per-mate cost of the full pipeline across the four
    operationally distinct regimes. Only the first regime decouples
    cost from $L$; all others pay a single traversal of the kept
    window. In the clean-insert regime the 3′ scan exits in
    $O(m_h + Delta_m)$ pairs, because a terminal mismatch triggers
    mercy within $O(Delta_m)$ steps or the candidate run fails the
    $m_h$ threshold after $O(m_h)$ steps.],
) <tab:regimes>

#Rem (space). Auxiliary space is $Theta(1)$. The precomputed
lookup tables for the quality gate, joint weight, and per-step
score are $O(t^2)$ and are shared across all reads in a batch.

#Rem (paired-end). Mates are processed independently; the
worst-case paired-end cost is $2(L + Delta_L)$. Early drop of one
mate does not allow skipping the other, as the decision is per-read.

== Default parameters

For typical Illumina short-read data, the values in @tab:defaults
are reasonable starting points.

#figure(
  table(
    columns: 2,
    stroke: 0.5pt,
    align: (left + horizon, center + horizon),
    table.header(
      [*Parameter*], [*Default value*],
    ),
    [$L_min$ — minimum post-trim length],  [$0.6 L$],
    [$Delta_m$ — mercy budget],            [$3.0$],
    [$t$ — Phred saturation],              [$20$],
    [$x$ — gate shape],                    [$1.0$],
    [$mu$ — match leniency],               [$0.1$],
    [$nu$ — mismatch strength],            [$0.5$],
    [$m_h$ — minimum run length],          [$6$],
  ),
  caption: [Default parameter values for typical Illumina short-read
    input. $nu$ and $Delta_m$ trade off: halving $nu$ roughly doubles
    the number of mismatches absorbed before the mercy clause fires.
    $mu$ has a much weaker effect, since match gaps are bounded by
    $mu$ per step.],
) <tab:defaults>

= Implementation optimisations

This section collects optimisations preserving the specification
above.

== Quality-gate tabulation

Phred scores are bounded non-negative integers, so the gate $g$
admits a complete tabulation over $q in {0, 1, dots, t - 1}$. The
case $q >= t$ returns $1$ and requires no table lookup.
Precomputed cost: $O(t)$ at initialisation. Per-query cost: a
single load.

== Joint-weight tabulation

$W_s (k)$ depends only on the pair of integer quality scores
$(q_(s+k), q_(s+k+1))$. A two-dimensional table of shape $(t, t)$
indexed by these scores yields $W$ in one load, eliminating the
division and power evaluation at scan time.

For $t = 40$ (typical Illumina), the table is
$40 times 40 = 1600$ floating-point entries ($approx 6.4$ kB in
$f 32$), trivially L1-resident.

== Per-step score tabulation <sec:scoretable>

Combining the match indicator and joint weight, the per-step score
$sigma_s (k)$ is a function of the quadruple
$(b_(s+k), b_(s+k+1), q_(s+k), q_(s+k+1))$. With 2-bit bases and
$t$-bounded qualities, a full table has $4 times 4 times t times t$
entries, or $16 t^2$ values. For $t = 40$, this is $approx 102$ kB
in $f 32$: $f 16$ halves it and L2-residency is comfortable. The
scan loop reduces to

$ sigma <- "table"[b_(s+k), b_(s+k+1), q_(s+k), q_(s+k+1)] $

followed by an accumulation into $Gamma$ and one comparison.

#Rem (quality binning). The $t^2$ factor assumes distinct integer
qualities. Recent Illumina basecallers do not emit them: RTA3
(NovaSeq) quantises Phred scores to four values and the HiSeq X
pipeline to eight. Mapping $q$ through a 256-byte binning table
before indexing collapses the per-step score table to
$4 times 4 times 4 times 4 = 256$ entries on NovaSeq input, which is
register-resident rather than merely L2-resident, and reduces the
joint-weight table of @sec:scoretable to a single addition of two
binned half-gates. The number of distinct quality values can be
established from the same sampled prefix used for adapter
auto-detection (@sec:adp-auto), with the unbinned path retained as
a fallback.

== Gap-only accumulation

The stopping rule and trim lengths depend only on $Gamma$, not on
$cal(P)$. Operating directly on $Gamma$ avoids carrying both
quantities and removes the algebraic subtraction $K - cal(P)_s (K)$
at every step. Analogously, a _gap-increment_ table $Delta Gamma$
indexed by the same quadruple as above allows

$ Gamma <- Gamma + Delta Gamma[b_(s+k), b_(s+k+1), q_(s+k), q_(s+k+1)] $

with termination triggered by $Gamma > Delta_m$. Because all gap
increments are non-negative, no cancellation occurs and $f 16$ or
fixed point (e.g. $u 8$ with $1 \/ 128$ resolution) becomes viable,
further shrinking the table and improving vectorisation.

== Two-bit base comparison <sec:xor>

Under $2$-bit base encoding, $b_k = b_(k+1)$ if and only if
$b_k xor b_(k+1) = 0$. The match indicator therefore reduces to one
XOR and one zero-test. When combined with the per-step score
tabulation, the XOR result can be used as a bit of the table index
rather than as a branch.

The match indicator is also the per-step contribution to the match
count $m^*$ used in the $m_h$ threshold (@sec:algorithm). Tracking
$m^*$ during the scan therefore costs one additional dependent
addition per step in the scalar implementation, or — under the SIMD
strategy of @sec:simd — a single popcount of the lane-wise match
mask once mercy fires, with no per-step overhead inside the
vectorised body.

== SIMD vectorisation <sec:simd>

The scan loop has exactly one loop-carried dependency — the running
gap $Gamma$ — and the per-step increment is a gather-load with a
small table. Bases and qualities are narrow (2 bits and $<= 1$ byte
respectively), so a $256$-bit SIMD register processes at least $16$
pairs per batch:

+ Gather the four-channel index $(b_i, b_(i+1), q_i, q_(i+1))$ for
  each lane.
+ Perform a gathered load of the increment table.
+ Compute a prefix sum of increments across lanes.
+ Compare against $Delta_m$ to obtain the first-exit lane as a
  bitmask.
+ Update $Gamma$ from the last unmasked lane.

The mercy-exit condition is evaluated as a lane mask and resolved
by a trailing-zeros count, avoiding data-dependent branches inside
the vectorised body.

#Rem (ordering of effort). Steps (a) and (b) are the expensive part
of this scheme: gather throughput on current x86 cores is roughly
one gather per five to six cycles, so the vectorised body is
gather-bound rather than compute-bound. The composition prefilter
of @sec:prefilter vectorises with plain byte comparisons and no
gather at all, and eliminates the scan on most reads rather than
accelerating it. It is therefore the first vectorisation target;
the gather-based scheme above is worthwhile only for the residual
fraction of reads that survive the prefilter.

== Struct-of-arrays layout

Each scan iteration touches one base and one quality at the same
index. Interleaving them as $(b_i, q_i)$ pairs (AoS) reduces the
two cache lines per iteration to approximately one, at the cost of
slightly less SIMD-friendly loads. For SIMD implementations a
hybrid layout — two aligned parallel arrays with software-prefetch
one cache line ahead — typically performs best.

== Reverse-scan without reverse copy <sec:reverse>

The 3′ end-scan can be implemented as a forward scan over the
original read with decremented indices, using symmetric definitions
$M'_s (k) := "match"(b_(L-1-s-k), b_(L-1-s-k-1))$ and analogous
$W', sigma'$. No allocation or reversal of the input is required.

== Budget short-circuit and early disable <sec:shortcircuit>

Because the decision rule drops a read whenever
$T_("tot") > Delta_L$, any scan that has already committed
$Delta_L$ trim-equivalent bases may terminate immediately, marking
the read for drop. Concretely:

- The 5′ end-scan terminates early once $K = Delta_L$.
- The 3′ end-scan runs with effective budget $Delta_L - T_5$.
- The interior scan runs with effective budget
  $Delta_L - T_5 - T_3$.

#Rem (degenerate budget). When $Delta_L = 0$ (i.e. $L_min = L$, no
trim is permitted), the short-circuit is inapplicable: it would fire
on step $K = 0$ before any pair has been examined and incorrectly
drop reads that have no homopolymer at all. In that regime the scan
must run to its mercy or boundary clause; the decision rule then
drops any read for which the resolved trim is positive. The same
caveat applies to subsequent scans whose effective budget evaluates
to zero. Equivalently: the budget cap is applied only when strictly
positive.

This short-circuit is load-bearing for the $O(Delta_L)$ best-case
bound of @sec:complexity. Additionally, if the library provenance
guarantees artefacts at only one end, the corresponding end-scan is
the dominant cost and skipping its counterpart roughly halves total
work; the interior scan may be skipped entirely when false-negative
drops (low-complexity reads passing through) are acceptable.

== Batch execution

In FASTX processing, many reads are processed sequentially.
Amortisation is aided by:

- Allocating the increment table once per process and sharing it
  across worker threads (read-only, cache-friendly).
- Processing reads in blocks aligned to SIMD width to avoid
  epilogue overhead.
- Prefetching the next read's buffers while processing the current
  one.

== Composition prefilter and the two-tier cascade <sec:prefilter>

The optimisations above reduce the constant factor of the scan but
not its structure: every kept read still pays one per-step score
evaluation per base. On a clean library the overwhelming majority
of reads carry no homopolymer at either end, and for those the
answer $T_5 = T_3 = 0$ can be established without evaluating a
single per-step score.

The observation that makes this possible is that the run-length
threshold of @sec:algorithm depends on the _composition_ of a
window, whereas $Gamma$ is a path functional. Composition supports
prefix sums and constant-time window queries; $Gamma$ does not.
Splitting the procedure along that line yields a two-tier cascade:
a quality-blind first tier that rejects, and the $Gamma$ scan as a
second tier that resolves the boundary of whatever survives.

#Def (anchor count). For a target base $X$ and a window
$W subset.eq [0, L)$, the _anchor count_ is
$ c_X (W) := |{k in W : b_k = X}|, $
and the mismatch count of $W$ is $|W| - c_X (W)$. In the generic
(unanchored) case the same role is played by the adjacency mask
$[b_k = b_(k+1)]$, whose population count over $W$ is the number of
matching steps in $W$.

#Prop (prefilter soundness). Let $K^*_s$ be the stopping step of
$"ScanFrom"$ under mercy budget $Delta_m$ and let $"mm"(K)$ denote
the number of mismatching steps among the first $K$. Then
$ "mm"(K) > floor(Delta_m) ==> Gamma_s (K) > Delta_m, $
and consequently $K^*_s <= l_"max"$, where $l_"max"$ is the largest
$K$ satisfying $"mm"(K) <= floor(Delta_m)$.

_Proof._ By the per-step gap contributions of @sec:algorithm, a
matching step contributes $mu(1 - W_s (k)) >= 0$ and a mismatching
step contributes $1 + nu W_s (k) >= 1$. Summing,
$Gamma_s (K) >= "mm"(K)$ for every admissible $K$. If
$"mm"(K) > floor(Delta_m)$ then, $"mm"(K)$ being an integer,
$"mm"(K) >= floor(Delta_m) + 1 > Delta_m$, hence
$Gamma_s (K) > Delta_m$ and the mercy clause has fired at or before
$K$. The bound on $K^*_s$ follows by contraposition, and $"mm"$ is
non-decreasing in $K$, so ${K : "mm"(K) <= floor(Delta_m)}$ is the
interval $[0, l_"max"]$. $qed$

#Prop (rejection test). Let $w := m_h + floor(Delta_m) - 1$ and let
$c_X$ be the anchor count over the $w$ positions an end-scan visits
first. If $c_X < m_h - 1$, that end-scan resolves to $T = 0$.

_Proof._ By the previous proposition the scan terminates at or
before the $(floor(Delta_m) + 1)$-th mismatching step. Among the
first $w$ steps the number of mismatching steps is at least
$w - c_X > w - (m_h - 1) = floor(Delta_m)$, so termination occurs
within the window. The match count at termination is bounded by the
anchor count of the visited prefix, so $m^* <= c_X < m_h - 1$,
giving $m^* + 1 < m_h$; the run-length threshold of @sec:algorithm
then resolves $T := 0$. $qed$

#Rem (one-sided error). Together the two propositions give exactly
the property a cheap first tier requires: it may pass windows that
the $Gamma$ scan subsequently rejects, but it can never reject a
window that the $Gamma$ scan would have trimmed. A procedure that
computes $l_"max"$ by counting and then runs $"ScanFrom"$
restricted to $[0, l_"max"]$ therefore returns the same $T_5$,
$T_3$ and $H_"int"$ as the unrestricted procedure, base for base.
The cascade is an optimisation, not a change of semantics, and is
testable as such: the conformance property is
$"cascade"(R, Q) = "PolyXTrim"(R, Q)$ over randomised reads and
randomised parameter tuples@claessen2000quickcheck.

#Rem (widening the window). Replacing $w$ by any $w' >= w$ keeps
the rejection test sound, since a wider window can only increase
$c_X$ while the mismatch bound $w' - c_X > floor(Delta_m)$ still
holds. Implementations may therefore round $w$ up to a machine word
or SIMD register width. For the defaults of @tab:defaults
($m_h = 6$, $Delta_m = 3$) the required width is $w = 8$, so a
single 32-byte comparison covers it with margin.

#Rem (SIMD realisation). Unlike the gather-based strategy of
@sec:simd, the first tier vectorises without gathers: a byte-wise
comparison of the terminal window against a broadcast anchor
($"VPCMPEQB"$), a mask extraction ($"VPMOVMSKB"$), and a population
count give $c_X$ over 32 bases in three instructions. The generic
case substitutes a comparison of the window against itself shifted
by one lane. This is the reason to implement the first tier before
attempting to vectorise $Gamma$ itself: it removes most of the work
rather than accelerating it, and it needs no gather.

#Rem (lattice-walk formulation). A cumulative-path encoding@gates1986, in
which each base contributes a fixed unit vector (for instance
$A = (0, 1)$, $C = (1, 0)$, $G = (-1, 0)$, $T = (0, -1)$) and a
homopolymer of length $l$ appears as a net displacement of $l$
along one axis, detects the same event. Over exactly $l$ steps the
condition $Delta y = l$ already forces $n_C + n_G + 2 n_T = 0$ and
hence a pure run, so the criterion is exact in the noiseless case
and the second coordinate is redundant. Under a drift tolerance it
is however strictly weaker than the anchor count: the displacement
$(n_C - n_G, n_A - n_T)$ compresses three composition degrees of
freedom into two, so C and G mismatches cancel and two errors may
be priced as one. It is also strictly more expensive, since window
queries require prefix sums over the whole read, whereas the
rejection test above touches $w$ terminal bases. The lattice
formulation is retained here for exposition and for diagnostic
visualisation, where a polyG tail appears as a straight ray, and
not as a component of the hot path.

The expected per-read cost of the cascade is
$ C = C_1 + Pr[text("tier 1 passes")] dot C_Gamma, $
where $C_1$ is the (constant) rejection test and $C_Gamma$ is now
bounded by $l_"max"$ rather than by $L$. The trade is asymmetric:
on a library where every read carries a polyG tail the cascade pays
$C_1$ on top of the full scan, whereas on a clean library it
removes the scan almost entirely. Since $C_1 lt.double C_Gamma$,
the worst case costs a few percent and the common case saves an
order of magnitude.

#Rem (measurement status). The cost model above is analytic. The
speedup factors of @tab:opt-impact predate the cascade and were
measured against the unfiltered scan; the cascade's realised effect
depends on the polyX-positive fraction of the input and has not yet
been measured on a production lane.

= Full algorithm <sec:algorithm>

This section assembles the primitives of §1 and §2 into a single
executable procedure. The layout — one forward pass accumulating
$T_5$ and $H_("int")$, one budget-capped reverse pass for $T_3$, and
a reconciliation step — achieves the $O(L + Delta_L)$ worst-case
bound of @sec:complexity with a single traversal of the read.

== Scan primitive

All three scans (5′ end, interior, 3′ end) share the same inner
loop, presented in @alg:scan. With the gap-increment table of
@sec:scoretable in $u 8$ fixed-point and the XOR match trick of
@sec:xor, each iteration is one gather-load, one add, one compare.

#figure(
```
procedure ScanFrom(s, u) → (K*, mode, m*)
    Γ ← 0;  K ← 0;  m* ← 0
    while s + K + 1 ≤ u do
        Γ ← Γ + ΔΓ[b[s+K], b[s+K+1], q[s+K], q[s+K+1]]
        if b[s+K] = b[s+K+1] then m* ← m* + 1   // match count
        K ← K + 1
        if Γ > Δ_m then return (K, mercy, m*)
    return (K, boundary, m*)
```,
  caption: [Scan primitive. Walks from position $s$ toward upper
    index $u$, accumulating gap mass until mercy or boundary
    termination. Returns the step count $K$, the stopping mode, and
    the match count $m^*$ (number of matching steps). The match
    count is the load-bearing input to the homopolymer-length
    threshold $m^* + 1 >= m_h$ used downstream; thresholding on
    $K^*$ alone would conflate run length with mercy overshoot.],
) <alg:scan>

The reverse variant, $"ScanFromReverse"(s, ell)$, is structurally
identical but indexes the read from the 3′ end with decremented
indices (@sec:reverse). The parameter $ell$ is a _budget_: the
scan additionally stops as soon as $K = ell$, allowing the 3′ scan
to terminate in $O(Delta_L)$ on a long polyG tail.

== Main pipeline

The full pipeline is given in @alg:main. Pass 1 handles 5′
end-scan and interior scan in a single forward traversal, recording
each mercy-terminated interior run for later reconciliation. Pass 2
runs a budget-capped reverse scan for $T_3$. The reconciliation
step removes any recorded interior run that overlaps the 3′ trim
zone, avoiding double-counting. Three short-circuit points (lines
10, 17, 26) exit the pipeline as soon as the drop decision becomes
inevitable.

#figure(
```
procedure PolyXTrim(R, Q, L)
    s ← 0;  T_5 ← 0;  H_int ← 0;  runs ← [ ]
    (K*, mode, m*) ← ScanFrom(0, L - 1)
    if m* + 1 < m_h then
        T_5 ← 0
    else if mode = mercy then
        T_5 ← K*
    else
        T_5 ← K* + 1
    if T_5 > Δ_L then return DROP              // 5′ short-circuit
    s ← K* + 1
    while s ≤ L - 2 do                         // interior scan
        (K*, mode, m*) ← ScanFrom(s, L - 1)
        if mode = mercy and m* + 1 ≥ m_h then
            H_int ← H_int + K*
            append (s, K*) to runs
            if T_5 + H_int > Δ_L then return DROP
        s ← s + K* + 1
    (K̃*, mode, m̃*) ← ScanFromReverse(0, Δ_L - T_5)   // Pass 2
    if m̃* + 1 < m_h then
        T_3 ← 0
    else if mode = mercy then
        T_3 ← K̃*
    else
        T_3 ← K̃* + 1
    if T_5 + T_3 > Δ_L then return DROP
    for (start, len) in reversed(runs) do      // reconcile
        if start + len - 1 ≥ L - T_3 then
            H_int ← H_int - len
        else break
    if T_5 + T_3 + H_int > Δ_L then return DROP
    return EMIT(R[T_5 : L - T_3], Q[T_5 : L - T_3])
```,
  caption: [Main pipeline. Pass 1 (lines 2–18) fuses 5′ end-scan
    and interior scan in a single forward traversal; Pass 2 (line
    19) is a budget-capped reverse scan; lines 27–30 reconcile
    interior runs that straddle the 3′ boundary. The
    homopolymer-length threshold uses the match count $m^*$ rather
    than $K^*$ to avoid trimming on scans whose distance was
    inflated by mercy overshoot.],
) <alg:main>

#Prop (correctness of reconciliation). The reconciliation step
(lines 27–30) removes exactly those recorded runs whose bases lie
in $[L - T_3, L - 1]$.

_Proof._ Pass 1 visits positions in strictly increasing order, so
the recorded runs are sorted by start position. For a run
$(s_i, K^*_i)$, its last position is $s_i + K^*_i - 1$. The
overlap condition $s_i + K^*_i - 1 >= L - T_3$ is monotone in
$s_i$: if it fails for run $i$, it fails for all $i' < i$ as well.
Iterating in reverse order and breaking at the first non-overlap
therefore examines exactly the straddling and fully-3′ runs. $qed$

#Rem (wholesale vs. partial removal). Lines 28–29 remove a
straddling run's contribution wholesale. The partial alternative
subtracts only the overlap length $"end" - (L - T_3) + 1$. Both
preserve the correctness of the drop decision whenever $T_3$ is
dominant; wholesale removal is slightly conservative (may drop a
read the partial rule would keep) but avoids per-base bookkeeping.
Either choice should be documented and held fixed across a
pipeline.

= Optimisations for AMD EPYC 9654 <sec:epyc>

The EPYC 9654 (Zen 4, 96 cores, 12 CCDs, 384 MB aggregate L3, 12
channels of DDR5-4800) is a useful concrete target: its memory
topology exposes most of the tradeoffs that matter for throughput
on FASTX pipelines. The recommendations below assume a single
socket with NPS4 (NUMA-per-socket = 4) or NPS1 BIOS setting; other
configurations shift the constants but not the advice.

== Memory topology and table placement

Each CCD has a private 32 MB L3 slice; access to a remote CCD's L3
(or DRAM channel in NPS4) traverses Infinity Fabric at roughly
10× the local-L3 latency. The per-step score table (@sec:scoretable)
is $approx 25$ kB in $u 8$ fixed-point, $approx 51$ kB in $f 16$,
$approx 102$ kB in $f 32$. Even the $f 32$ version fits comfortably
in a single CCD's L3, but _not_ in L1D (32 kB), so the table
placement that matters is L3.

Concretely:

- Replicate the increment table _per CCD_, not per process or per
  socket. With 12 CCDs and a 25 kB table, the aggregate cost is
  300 kB — negligible against the 384 MB of L3.
- Pin worker threads to CCDs with `numactl --cpunodebind=N` (in
  NPS4, each CCD is its own NUMA node) or use `hwloc` to enumerate
  CCDs programmatically.
- Allocate per-thread read buffers with first-touch on the worker's
  CCD so the DRAM channels servicing the buffer are local.
- Request 2 MB huge pages for the FASTX input buffer
  (`madvise(MADV_HUGEPAGE)` or `mmap(MAP_HUGETLB)`). For 150-bp
  reads at throughput, the TLB miss rate on 4 kB pages becomes
  measurable; 2 MB pages eliminate this.

Avoid `numactl --interleave=all`: it spreads the increment table
across CCDs, which is the opposite of what you want for a
read-heavy hot table.

== AVX-512 scan loop

Zen 4 supports the full AVX-512 instruction set over a 256-bit
datapath (double-pumped), so 512-bit ops retire at roughly one per
cycle for most lanes. Relevant instructions for the scan loop:

- `VPGATHERDD` / `VPGATHERDQ` for the 16-lane gather into $Delta Gamma$.
  Zen 4 gather throughput is ~1 per 5–6 cycles; above $L approx 20$
  this still beats scalar, below it the scalar loop wins. In
  practice, 150-bp reads always benefit.
- `VPERMB` (from AVX-512 VBMI) unpacks 2-bit base indices into byte
  lanes in ~3 cycles, feeding the gather index register without a
  scalar shuffle.
- `VPCMPUB` + `KMOVW` + `TZCNT` implements the mercy-exit check
  branchlessly: compare $Gamma$ against $Delta_m$ across lanes,
  extract the mask, find the first set bit.
- Unroll the accumulator by 4. The $Gamma$ dependency chain is
  ADD + CMP, roughly 2 cycles; four rotating accumulators hide the
  latency and let the gather port saturate.

The $u 8$ fixed-point increment table (1/128 resolution) is
strongly preferred over $f 16$ or $f 32$: the table is smaller, the
ADD is single-cycle, and saturation arithmetic (`VPADDUSB`) caps
$Gamma$ at 255 for free, so the mercy comparison becomes a fixed
threshold $Delta_m dot 128$ against an 8-bit accumulator.

== I/O and decompression pipeline

At full throughput the scan compute is often not the bottleneck —
gzip decompression is. Recommendations:

- Replace zlib with `libdeflate`@libdeflate (~2× faster scalar decode) or
  `ISA-L`@isal (hardware-accelerated, ~4–5× faster). Both are drop-in.
- For new pipelines, prefer `zstd`@zstd over gzip: comparable ratios,
  faster decode, and the reference implementation already uses
  optimised paths on AMD64.
- Use `io_uring` for async FASTX reads. This overlaps file I/O
  with scan compute and is essential when reading from
  cluster-mounted storage (Lustre, BeeGFS) with non-trivial
  latency.
- Double-buffer input: two 2 MB huge-page buffers rotated between
  reader and scanner so the scan never stalls on I/O.

== Thread topology and SMT

The scan loop is L1D-bandwidth bound once the tables are hot.
Recommendations:

- Start with 1 worker thread per physical core (96 threads), SMT
  disabled or ignored. This is the configuration that maximises
  throughput per thread and avoids L1D contention between sibling
  threads.
- If the pipeline includes decompression or other non-L1-bound
  work, enable SMT and gang-schedule a compute thread with an I/O
  or decompression sibling on the same core. Realised speedups
  from SMT in this configuration are typically 1.3–1.5×.
- Benchmark both. The crossover depends on compression level,
  read length, and how much of the pipeline is fused into a single
  worker.

Pin threads explicitly (`sched_setaffinity` or `taskset`) rather
than relying on the kernel scheduler. On a 96-core socket, thread
migration across CCDs costs more than the scheduler is accounting
for.

== System configuration

A small set of system settings materially affects throughput:

- `cpupower frequency-set -g performance` (or `tuned-adm profile
  throughput-performance`). The default `schedutil` governor
  under-clocks during short bursts, visible as ~5–10% loss on
  mixed-workload FASTX pipelines.
- Disable deep C-states (`C2` and below) on worker cores:
  `cpupower idle-set -D 1`. The wake-up latency from C2 is ~10 µs,
  enough to matter when reads arrive in bursts.
- NPS4 BIOS setting exposes each CCD as a NUMA node, enabling
  fine-grained first-touch allocation. NPS1 is simpler and adequate
  if the workload fits in a single CCD's L3 per worker, which for
  this pipeline it does.
- Transparent huge pages set to `always` (or `madvise` with
  explicit hints) for large buffers.
- Move IRQs off the worker cores: the network and storage IRQs
  should be pinned to a small reservation (e.g., cores 0–3) with
  `irqbalance` disabled.

A well-tuned deployment of this pipeline on a single EPYC 9654 can
sustain several GB/s of compressed FASTQ throughput — generally
I/O-bound rather than compute-bound, which is the intended
operating regime.

= Pre- and post-processing extensions <sec:extras>

The core algorithm of §1 trims homopolymer artefacts. Production
pipelines combine that with three adjacent operations: removal of
PCR duplicates, an unconditional hard trim of fixed-length anchors
(adapter remnants, UMI bases, primer sequence) at either end, and a
per-read mean quality filter. This section specifies all three as
complements to the core procedure rather than modifications of it,
so that the $O(L + Delta_L)$ complexity bound and the correctness
arguments of §§1–3 carry over unchanged. They are presented in the
order in which the pipeline applies them.

== Deduplication <sec:dedup>

PCR amplification produces several sequenced copies of the same
template fragment. Duplicates inflate coverage, distort allele
frequencies, and give false support to variants that arose during
amplification rather than in the sample@ebbert2016. Their removal is
conventionally a separate pipeline stage, which costs one full
decompress, parse, write and compress round trip over the data.
Because that round trip dominates wall-clock time on compressed
input, fusing deduplication into the trimming pass is worth more
than any inner-loop optimisation specified elsewhere in this
document. The procedure below is the one implemented in
`FastDedup`@fastdedup, restricted to the exact-match case and
adapted to run ahead of trimming.

#Def (fingerprint). Let $H : {0,1}^* -> {0,1}^N$, with
$N in {64, 128}$, be a non-cryptographic hash. For a single-end read
the fingerprint is $phi(R) := H(R)$. For a pair, a plain
$H(R_1) xor H(R_2)$ is inadmissible: it vanishes identically
whenever $R_1 = R_2$, so every pair with identical mates would
collide deterministically. The fingerprint is instead
$ phi(R_1, R_2) := H(R_1) xor "rot"(H(R_2), N \/ 2), $
where $"rot"(x, j)$ is a bitwise left rotation by $j$. Rotation is
bijective, so it preserves the entropy of its operand, and the
exclusive-or of two independent uniform variables is uniform; the
combined fingerprint therefore remains uniformly distributed. It is
also order-dependent, which is intended, since $R_1$ and $R_2$ carry
distinct geometric meaning.

#Prop (collision bound). For $x$ distinct fingerprints drawn
uniformly from $N$ bits, the probability that any collision occurs
anywhere in the dataset is approximately
$ p approx x^2 \/ 2^(N+1). $
At $N = 64$ this reaches $10^(-2)$ near $x approx 0.6 times 10^9$
sequences; at $N = 128$ the same bound is reached near
$2.6 times 10^18$. The implementation selects $N$ from an estimate of
$x$ and a user-supplied threshold on $p$.

#Rem (interpretation of the bound). $p$ bounds the existence of at
least one collision over the whole run, not a per-read
false-duplicate rate. A collision removes exactly one read, so at
the threshold the expected number of spuriously dropped reads is
below one, which is negligible against a dataset large enough to
approach it.

=== Anchor window

Fingerprinting the full read before trimming is unsound on
two-colour chemistry. PolyG tails arise from per-cluster signal
dropout, which is stochastic: two amplification copies of one
fragment need not fail on the same cycle, so one may read
$I ‖ G^8$ and the other $I ‖ G^12$ for the same insert $I$. Their
full-read fingerprints differ and both survive. The same divergence
occurs, more weakly, at adapter boundaries. Deduplicating on the
full read before trimming therefore under-detects duplicates in
precisely the artefact class this document is concerned with.

#Def (anchor window). Fix a reference length $L_"ref"$ for the run,
taken from the expected read length or from the first record. Let
$a := T_5^"hard"$ and
$k := min(L_min (L_"ref"), L_"ref" - T_5^"hard" - T_3^"hard")$. The
_anchor window_ is $W_"anc" := [a, a + k)$, and the fingerprint is
computed over $R[W_"anc"]$ rather than over $R$.

#Rem (the width must not depend on $L$). Deriving $k$ from the length
of the individual read defeats the construction. Under a fractional
$L_min$ two copies of one fragment that differ in length receive
different widths, hence different windows, hence different
fingerprints, which is precisely the divergence the anchor exists to
absorb. Fixed-length Illumina output masks the error, because there
$L = L_"ref"$ for every read; already-trimmed input, long-read data and
FASTA do not, and the property should not be conditional on the
platform. Reads shorter than $a + k + T_3^"hard"$ cannot carry the
window and fall back to a full-read fingerprint; these are reads
already below $L_min$, which the decision rule discards in any case.

#Prop (anchor invariance). $phi$ is a function of $R[W_"anc"]$
alone, hence invariant under every modification of $R$ outside
$W_"anc"$. In particular, two reads agreeing on $W_"anc"$ receive
identical fingerprints regardless of the length or content of their
3′ homopolymer tails and adapter remnants.

_Proof._ Immediate from the definition: the hash is evaluated on the
window only. $qed$

#Rem (soundness condition). The proposition is only useful if
$W_"anc"$ is disjoint from the region the configured pipeline
targets. Under $"Direction::Forward"$ with 5′ scanning disabled the
prefix satisfies this. Under $"Direction::Reverse"$, or when the 5′
end-scan is enabled for the Nextera geometry of @sec:complexity, the
anchor must instead be taken from the 3′ end: the rule is to anchor
on the end that is _not_ scanned. When both end-scans are enabled
there is no artefact-free anchor, and the implementation must either
fall back to fingerprinting the trimmed read after the polyX
procedure, at the cost of performing the scan on duplicates, or
accept the reduced sensitivity of a full-read anchor. This choice is
exposed rather than decided here, because it depends on library
provenance.

#Rem (what the anchor discards). Restricting to a prefix discards
fragment-length information. For single-end reads of uniform length
none is available in any case, since two fragments sharing a start
position yield the same read whatever their lengths. For paired-end
input both mates are anchored and hashed together, which pins both
ends of the fragment and recovers the discrimination between two
fragments that share a start but differ in length. Single-end
deduplication is, unavoidably, deduplication by start position.

#Rem (memory). This is the only stage of the pipeline whose space is
not $Theta(1)$ in the number of reads. The fingerprint set grows as
$Theta(U)$ in the number of _unique_ reads, at $N \/ 8$ bytes each
plus open-addressing overhead: approximately $4.6$ GB at
$3 times 10^8$ reads for $N = 64$, and $8.7$ GB for $N = 128$. Every
other stage is chunk-bounded. Deduplication must therefore be
opt-in, and its memory cost must be stated at the point of use
rather than in documentation alone.

#Rem (pair-level semantics). For paired-end input the decision is
taken on the pair. A duplicate pair drops both mates and yields no
singleton. Singletons produced later by the polyX or adapter stages
are not deduplicated against a single-end fingerprint space, since
mixing the two would allow a surviving mate to collide with an
unrelated read.

#figure(
```
procedure Dedup(R1, R2 or ⊥, S, a, k)
    w1 ← R1[a : a + k]
    if R2 = ⊥ then
        φ ← H(w1)
    else
        w2 ← R2[a : a + k]
        φ ← H(w1) xor rot(H(w2), N/2)
    if φ ∈ S then return DROP(Duplicate)
    S ← S ∪ {φ}
    return CONTINUE
```,
  caption: [Deduplication stage. Runs before every other operation on
    the read. $S$ is the fingerprint set, $[a, a+k)$ the anchor window
    of @sec:dedup. The set is the only unbounded state in the
    pipeline; all other stages are chunk-bounded.],
) <alg:dedup>

=== Admission control and saturation <sec:dedup-mem>

Because $|S|$ grows without bound in the number of unique reads, the
stage can exhaust the memory available to the process. Both the
condition under which the run is refused and the behaviour when the
bound is reached mid-run are specification-level properties, not
implementation details, because they change what the tool guarantees
about its output.

#Def (enforced memory bound). Let $M$ denote the quantity of memory
the process may consume before it is terminated. $M$ is *not* the
physical memory of the machine. On a Linux host it is the minimum of
every discoverable constraint: the control-group limit of the process
and of every ancestor group, the allocation granted by the batch
scheduler, the address-space resource limit, and the currently
available physical memory. A job granted $8$ GB on a $1.5$ TB node has
$M = 8$ GB, and a projection made against the node's physical memory
would be wrong by two orders of magnitude.

#Prop (exact projection). The fingerprint set is an open-addressed
table@hashbrown with a maximum load factor of $7 \/ 8$ and one control
byte per bucket. For a requested capacity $c$ and a fingerprint width of $N$
bits, the resident size is
$ B(c, N) = 2^ceil(log_2 (8c \/ 7)) dot (N \/ 8 + 1) "bytes". $
This is an identity, not an estimate: the bucket count is a power of
two fixed by $c$ alone, and no rehashing occurs when the table is
pre-sized.

#Rem (verification of the projection). At $c = 3 times 10^8$ the
formula gives $2^29$ buckets, hence $4608$ MB at $N = 64$ and $8704$ MB
at $N = 128$. Measured peak resident set sizes for the same
configuration are $4616$ MB and $8712$ MB, the $8$ MB difference being
the remainder of the pipeline. Because the projection is exact,
refusing a run on the strength of it is defensible in a way that a
heuristic estimate would not be.

#Rem (the bucket count is a step function). $B$ doubles as $c$ crosses
$c_k := 7 dot 2^k \/ 8$, so a capacity estimate exceeding one of these
points by a single element doubles the allocation. Two consequences
follow. First, the projected bucket count and not merely the projected
byte total should be reported, since a user marginally over a boundary
can act on that information. Second, a clamp on the capacity estimate
should be placed *at* a boundary. Note that the boundaries are not the
powers of two themselves: at $N = 64$ the largest capacity fitting in
$2^29$ buckets is $c_29 = 469 thin 762 thin 048$, projecting $4608$ MB,
whereas a round clamp at $5 times 10^8$ exceeds it by $6%$, rounds to
$2^30$ buckets, and therefore projects $9216$ MB for no benefit.

#Rem (capacity should be sized on unique reads). The table holds
unique fingerprints, whereas an estimate derived from file size counts
all reads. Scaling the estimate by the complement of the expected
duplication rate is a reduction of the order of $10%$, which is
frequently enough to cross a bucket boundary downwards and halve the
allocation. The rate may be taken from the sampled prefix already read
for adapter auto-detection (@sec:adp-auto).

#Def (admission rule). Let $R$ be a reserve covering the chunk arenas,
compression buffers and allocator overhead. The run is admitted when
$B(c, N) + R < M$, and refused otherwise. On refusal, if the width was
$N = 128$ and $B(c, 64) + R < M$, the implementation may instead
proceed at $N = 64$, reporting the collision probability that follows
from the bound of @sec:dedup. This is a quantified degradation and is
preferable to failure.

#Def (saturation). The projection depends on a capacity estimate,
which may be wrong. If resident usage crosses the bound during the
run, the set enters the _saturated_ state: it continues to answer
membership queries but accepts no further insertions. Reads whose
fingerprints were recorded before saturation are still recognised as
duplicates; reads first seen afterwards are never recorded and their
later copies are not detected.

#Prop (saturation is one-sided). Let $S_t$ and $S'_t$ be the
fingerprint sets after $t$ reads in an unsaturated and a saturated run
on the same input. Then $S'_t subset.eq S_t$ for every $t$, and
consequently every read dropped by the saturated run is also dropped
by the unsaturated run. Saturation can therefore fail to remove
duplicates, but can never remove a read that the unsaturated procedure
would have kept.

_Proof._ By induction on $t$. The base case is
$S'_0 = S_0 = nothing$. Assume $S'_(t-1) subset.eq S_(t-1)$ and
consider read $t$ with fingerprint $phi_t$. If
$phi_t in S_(t-1)$ then $S_t = S_(t-1)$, and $S'_t$ is either
$S'_(t-1)$ or $S'_(t-1) union {phi_t}$, both contained in $S_(t-1)$
since $phi_t in S_(t-1)$. If $phi_t in.not S_(t-1)$ then
$S_t = S_(t-1) union {phi_t}$ and
$S'_t subset.eq S'_(t-1) union {phi_t} subset.eq S_t$. In both cases
$S'_t subset.eq S_t$. A read is dropped exactly when its fingerprint
lies in the set at that point, and $S'_(t-1) subset.eq S_(t-1)$, so the
saturated run's drop set is a subset of the unsaturated run's. $qed$

#Rem (saturation must be reported). The proposition makes saturation
safe with respect to data loss, not with respect to interpretation: a
saturated run reports a duplication rate lower than the true one, and a
user who is unaware of this will draw wrong conclusions about library
complexity. The state must therefore appear in the machine-readable
statistics, together with the read index at which it was entered, and
on the diagnostic stream. A silently partial deduplication is the
failure mode this document criticises in probabilistic filters; the
difference here is only that the error is one-sided and declared, and
that difference is lost if it is not declared.

#Rem (source of the bound). Diagnostics should name which constraint
produced $M$. When a job on a machine with a terabyte of memory refuses
a five-gigabyte allocation, the identity of the binding constraint is
the information the user needs, and it is not recoverable from the
numbers alone.

== Hard trim

#Def (hard trim parameters). The procedure accepts two
non-negative integers $T_5^"hard"$, $T_3^"hard" in NN_(>=0)$
configuring fixed cuts at the 5′ and 3′ ends respectively. They are
applied prior to any scan, on the original read $(R, Q)$ of length
$L$.

#Def (hard-trim semantics). Define the working sub-read
$(R', Q')$ by
$ R' := (b_(T_5^"hard"), dots, b_(L - 1 - T_3^"hard")), $
$ Q' := (q_(T_5^"hard"), dots, q_(L - 1 - T_3^"hard")), $
of length $L' := L - T_5^"hard" - T_3^"hard"$. The polyX procedure
of §1 is then run on $(R', Q')$ with trim budget
$Delta_L' := Delta_L - T_5^"hard" - T_3^"hard"$. The total trim
returned to the decision rule is augmented:
$T_5^"out" := T_5^"hard" + T_5$,
$T_3^"out" := T_3^"hard" + T_3$,
with $H_("int")$ unchanged.

#Prop (hard-trim short-circuits). Three degenerate inputs admit
immediate decisions before the scan:

+ If $T_5^"hard" + T_3^"hard" > Delta_L$, the hard trim alone
  exceeds the trim budget; the read is dropped without scanning.
+ If $L' < 2$, the working sub-read is too short to admit even one
  scan step; the read is dropped.
+ If $T_5^"hard" + T_3^"hard" = L$, no read content remains; the
  read is dropped.

In all surviving cases $L' >= 2$ and $Delta_L' >= 0$, satisfying
the preconditions of @alg:main.

#Rem (interaction with adapter trimmers). When $T_5^"hard"$ targets
adapter or primer bases, the polyX scan thereafter starts at the
true insert. Composition is intentional: a homopolymer immediately
following a hard-trimmed adapter is detected and trimmed by the
core algorithm, so the user need not specify the run length in the
hard cut.

== Mean quality filter

#Def (mean quality). For a read $(R, Q)$ of length $L$ the
_mean Phred quality_ is
$ overline(q) := 1/L sum_(i=0)^(L - 1) q_i. $
The mean is computed over raw Phred values; FASTQ inputs encoded
with an ASCII offset (33 for modern Illumina, 64 for legacy) are
decoded once before this step.

#Def (mean-quality filter). Given a threshold
$theta_q in RR_(>=0)$, the read is _filtered_ (dropped) if
$overline(q) < theta_q$. The filter is evaluated on the original
read, before any hard trim or polyX scan: the decision reflects the
quality of the input, not of the kept slice.

#Rem (rationale for filtering before trimming). Computing
$overline(q)$ on the post-trim slice would be circular: the trim
itself is sensitive to quality (through the gate $g$ in §1.3), so a
post-trim mean would describe the result of the algorithm rather
than the read it acts on. Filtering before trimming gives a clean
quality criterion independent of homopolymer content. Pipelines
preferring fastp-style post-trim filtering can disable the
mean-quality filter here and apply it externally.

#Prop (cost). Computing $overline(q)$ is $Theta(L)$ in time and
$Theta(1)$ in space. Combined with the polyX procedure, the
whole-pipeline complexity of @sec:complexity becomes
$Theta(L) + O(L + Delta_L) = O(L + Delta_L)$, unchanged
asymptotically; the constant gains one full pass through $Q$.

#Rem (implementation). The mean-quality pass is the only place
where the full quality vector is unconditionally touched; all other
passes either operate on the working sub-read $(R', Q')$ or are
short-circuited. Implementations on hardware with a streaming
prefetcher (e.g. AMD Zen 4 of @sec:epyc) absorb the cost into the
prefetch shadow of the subsequent scan.

== Composition

The complete pipeline including extensions is:

+ Fingerprint the anchor window per @alg:dedup. If the fingerprint
  has been seen, drop; otherwise record it.
+ Decode Phred ASCII offset to obtain raw $Q$.
+ Compute $overline(q)$. If $overline(q) < theta_q$, drop.
+ Apply hard trim per @alg:hardtrim, possibly dropping early.
+ Run the polyX procedure of @alg:main on $(R', Q')$ with the
  reduced budget $Delta_L'$.
+ Augment $T_5, T_3$ by their hard-trim contributions and apply
  the decision rule.

#figure(
```
procedure HardTrim(R, Q, L, T_5_hard, T_3_hard, Δ_L)
    if T_5_hard + T_3_hard > Δ_L then return DROP
    L' ← L - T_5_hard - T_3_hard
    if L' < 2 then return DROP
    R' ← R[T_5_hard : L - T_3_hard]
    Q' ← Q[T_5_hard : L - T_3_hard]
    Δ_L' ← Δ_L - T_5_hard - T_3_hard
    return (R', Q', Δ_L')
```,
  caption: [Hard-trim wrapper. Returns either a working sub-read
    $(R', Q')$ with reduced budget $Delta_L'$ for the polyX
    procedure, or a drop decision when the hard trim alone makes
    the read non-emittable.],
) <alg:hardtrim>

#Rem (orthogonality). Neither extension affects the asymmetric
score, the stopping rule, or the reconciliation step. They operate
on either side of the polyX procedure and pass through its
correctness guarantees verbatim. A single conformance test
verifying that hard trim with $T_5^"hard" = T_3^"hard" = 0$ and
mean-quality filter with $theta_q = 0$ reproduces the §1 procedure
exactly is sufficient to confirm the extensions do not perturb core
behaviour.

= Adapter trimming <sec:adapter>

Adapter contamination arises when the insert is shorter than the read
length: the sequencer reads through the insert into the adapter
ligated at the far end. Two complementary mechanisms handle this in
production pipelines: _sequence matching_ (single-end, or
paired-end fallback) detects a known adapter sequence by sliding it
over the 3′ tail of the read; _PE overlap detection_ infers the
insert size geometrically by comparing R1 and R2 directly, without
requiring knowledge of the adapter sequence. This section specifies
both, together with a sampling-based _auto-detection_ procedure that
eliminates the need to name the adapter at the command line.

== Input and parameters <sec:adapter-params>

#Def (adapter). An _adapter_ is a pair $(eta, A)$ where $eta$ is a
human-readable name (used in statistics and metadata) and
$A = (a_0, dots, a_(p-1)) in {A, C, G, T}^p$ is the adapter
sequence of length $p$. A finite set $cal(A) = {(eta_j, A_j)}$ of
adapters constitutes an _adapter set_.

#Def (adapter parameters). The adapter sub-procedure is configured
by the tuple $(cal(A), epsilon, o_"SE", o_"PE")$:

#set enum(numbering: "(a)")
+ $cal(A)$ — the adapter set (possibly empty).
+ $epsilon in [0, 1]$ — maximum mismatch rate; a candidate alignment
  with fractional mismatch exceeding $epsilon$ is rejected.
+ $o_"SE" in NN_(>=1)$ — minimum overlap for single-end sequence
  matching (and paired-end fallback).
+ $o_"PE" in NN_(>=1)$ — minimum overlap for paired-end geometric
  detection; should satisfy $o_"PE" <= L$ to be useful.

The adapter trim amounts $T_5^"adp"$ and $T_3^"adp"$ carry the
same positional meaning as the hard-trim parameters of §5.1: they
are applied to the 5′ and 3′ ends of the read respectively, prior
to the polyX scan.

== Sequence matching (SE and PE fallback) <sec:adp-se>

#Def (overlap alignment). For a read $R = (b_0, dots, b_(L-1))$
and an adapter $A = (a_0, dots, a_(p-1))$, a _3′ overlap_ of
length $k$ is the alignment that places $a_0, dots, a_(k-1)$
against the last $k$ bases $b_(L-k), dots, b_(L-1)$ of the read.
The mismatch count of this alignment is
$ "mm"(k) := |{i in {0, dots, k-1} : b_(L-k+i) eq.not a_i}|. $
The alignment is _admissible_ if $k >= o_"SE"$ and
$"mm"(k) / k <= epsilon$.

#Def (best match). For each adapter $(eta_j, A_j) in cal(A)$,
compute the admissible overlap of largest $k$ (greediest trim
first). Among admissible overlaps, select the one with the
largest $k$ across all adapters, breaking ties by adapter order
in $cal(A)$. If no admissible overlap exists for any adapter,
there is no match and $T_3^"adp" := 0$.

#Rem (direction). Only 3′ adapter contamination is modelled here:
$T_5^"adp" := 0$ in single-end mode. In PE mode, the roles of $T_5$
and $T_3$ swap between mates depending on the geometry; see
@sec:adp-pe.

#Rem (complexity). Scanning each adapter against one read is
$O(L + p)$ per adapter in the worst case (decrementing $k$ from
$min(L, p)$ to $o_"SE"$, with an early exit on the first
admissible overlap). Over $|cal(A)|$ adapters the cost per read
is $O((L + p_"max") dot |cal(A)|)$, dominated by $O(L)$ in practice
because $p_"max" lt.double L$ for typical Illumina adapters ($p approx
30$--$60$) and $|cal(A)|$ is small (3–5 in a preset). The
per-base work is one byte comparison and one integer add per step,
with no table lookup.

#Rem (suffix locality). The overlap alignment above is anchored at
the 3′ end, and this is not a simplification but a consequence of
the contamination model of @sec:adapter. If the insert length is
$I < L$, the read is $"insert"[0 : I]$ followed by
$A[0 : L - I]$: the adapter begins at position $I$ and continues to
the end of the read, because there is no further template to
sequence. An adapter occurrence followed by more insert sequence
cannot arise from readthrough. Searching internal start positions
therefore searches for an event the model excludes, at a cost of
$O(L)$ candidate alignments instead of the
$O(min(L, p) - o_"SE")$ suffix alignments of @alg:adp-se, without
any corresponding gain in sensitivity.

#Rem (seeded matching). The $O((L + p_"max") dot |cal(A)|)$ bound
above is pessimistic even for suffix alignments, since every
candidate $k$ re-reads the tail. A seed prefilter removes the
repetition: locate the first $k_"seed" in {8, dots, 12}$ bases of
$A$ in $R$ with a single substring search, then verify the one
suffix alignment implied by the hit. Substring search over a byte
string runs at well under one cycle per byte with vectorised
implementations, so the per-read cost becomes $O(L)$ with a small
constant and one verification of length $<= p$. This is the same
cascade discipline as @sec:prefilter applied one layer out: a cheap
over-approximation in front of an exact resolver.

#Rem (indels). The mismatch model above is Hamming, not edit
distance. For Illumina chemistry this is adequate: the substitution
rate is of order $10^(-3)$ per base while the indel rate is of
order $10^(-6)$, so an indel-tolerant alignment fires on a
negligible fraction of reads at a per-cell cost several times that
of a byte comparison. Pipelines requiring edit-distance semantics
(other chemistries, linked or anchored 5′ adapters) can restore a
banded semi-global alignment as a third tier, invoked only on
candidates where the seed matches but Hamming verification fails.

The sequence-matching procedure is shown in @alg:adp-se.

#figure(
```
procedure MatchAdapters(R, L, 𝒜, ε, o_SE) → (T3_adp, name)
    best_k ← 0;  best_name ← ∅
    for (η, A) with |A| = p in 𝒜 do
        for k ← min(L, p) downto o_SE do          // greediest first
            mm ← count_mismatches(R[L−k : L], A[0 : k])
            if mm / k ≤ ε then
                if k > best_k then
                    best_k ← k;  best_name ← η
                break                               // largest k for this adapter
    return (best_k, best_name)
```,
  caption: [Single-end adapter matching. Tries overlaps from largest
    to smallest for each adapter; accepts the first admissible overlap
    per adapter, then reports the globally largest overlap. Returns
    $T_3^"adp" = 0$ when no adapter matches.],
) <alg:adp-se>

== Paired-end overlap detection <sec:adp-pe>

#block(inset: 0.6em, stroke: 0.5pt + gray, radius: 2pt, [
  *Implemented behind `--pe-overlap`, off by default.* This section specifies the
  paired-end overlap mode, and the current revision carries it as an opt-in flag
  rather than a default. Two findings shaped that choice and are worth stating
  with their numbers.
  #block(inset: 0.6em, [
    *Detection is anchored, not the scan below.* The algorithm as specified accepts
    the first insert size whose mismatch fraction is within $epsilon$; at
    $epsilon = 0.1$ a 150-base window admits 15 mismatches, so it fires on
    unrelated sequence. The implementation uses a 12-base anchor with
    full-window verification and refuses ambiguity instead, following the
    suite's own discipline (`fastcheck`'s insert measurement).
  ])
  #block(inset: 0.6em, [
    *The boundary is applied only when adapter starts exactly at it.* An overlap
    proves where the fragment ends; it does not prove that anything was ligated
    past it. The first implementation missed that and trimmed biology on
    adapter-free data --- the same failure class as the 4-base-overlap false
    trim that `--min-adapter-overlap 8` exists to prevent. The guard is not
    thresholded: any firing on adapter-free data is disqualifying.
  ])
  Measured on the 102-dataset corpus: adapter readthrough is configured in 53 of
  102 profiles but materialises in 11,395 of 10.2M pairs (0.112%), because it
  needs an insert shorter than the read. The corpus therefore cannot answer
  whether the mode improves accuracy, and a synthetic cohort whose ground truth
  carries both G-runs and readthrough is used for that instead
  (`docs/plans/2026-10-19_004500-pe-overlap-with-without.md`).
  #block(inset: 0.6em, [
    *Result: the mode is a measured null.* On that cohort (800,000 mates, 14.1%
    readthrough, 489,546 mates carrying a G-run) exact matches are identical
    (347,171) with and without the flag, and overtrim rises by 106 --- +0.094 pp
    on the readthrough slice, +0.000 pp overall, against a pre-registered 0.5 pp
    threshold. The mode is therefore not restored as a default. The mechanism is
    that the changed decisions are either invisible to the scorer
    (dropped-to-dropped) or net-negative: it trims 1--2 bp *less* on mates the
    ordinary path had already handled, correcting an over-trim rather than
    discovering one, over an addressable set of 1.2% of readthrough mates.
  ])
])

In paired-end sequencing, when the insert length $I$ is shorter
than the read length $L$, R1 reads bases $0, dots, L-1$ of the
fragment while R2 reads the same fragment from the opposite strand.
In standard orientation (R2 reverse-complemented to the same
strand as R1), the first $I$ bases of R1 equal the first $I$ bases
of R2 reverse-complemented:

$ R_1[i] = overline(R_2[I - 1 - i]) quad "for" i = 0, dots, I - 1, $

where $overline(b)$ denotes the complement of base $b$. The
positions $I, dots, L - 1$ in each mate are adapter-derived.

#Def (PE overlap). For a read pair $(R_1, R_2)$ of length $L$, an
_overlap_ of size $I$ is admissible if $I >= o_"PE"$ and
$ "mm"_"PE"(I) :=
  |{i in {0, dots, I-1} : R_1[i] eq.not overline(R_2[I-1-i])}|
  \/ I <= epsilon. $
Define the _PE trim amount_
$ T^"PE" := L - I. $
R1 receives $(T_5^"adp"=0,\ T_3^"adp"=T^"PE")$; R2 receives
$(T_5^"adp"=T^"PE",\ T_3^"adp"=0)$.

#Def (conservative scan). Overlaps are tested in decreasing order
of $I$ (largest insert first, smallest trim first), and the first
admissible overlap is accepted. This minimises trim length: among
all valid insert sizes, the longest is preferred, since a shorter
trim is always at least as correct.

_Termination._ The scan runs over $I in {min(L_1, L_2), dots,
o_"PE"}$ (at most $L - o_"PE"$ iterations) and always terminates.
When no admissible overlap exists (long insert, clean reads), the
procedure returns $T^"PE" = 0$ and sequence matching is applied as
a fallback.

_Rationale._ PE overlap detection requires no knowledge of the
adapter sequence, is robust to adapter sequence variation (e.g.
dual-indexed or custom adapters), and is computationally cheaper
per pair than two SE sequence-matching passes since the per-base
cost is one complement comparison with no table lookup.

The PE overlap procedure is shown in @alg:adp-pe.

#figure(
```
procedure DetectPEOverlap(R1, R2, L1, L2, ε, o_PE) → T_PE
    for I ← min(L1, L2) downto o_PE do
        mm ← 0
        for i ← 0 to I − 1 do
            if R1[i] ≠ complement(R2[I − 1 − i]) then mm ← mm + 1
        if mm / I ≤ ε then
            T_PE ← L1 − I           // adapter length in R1
            if T_PE = 0 then return 0   // no adapter content
            return T_PE
    return 0                         // no overlap found; use SE fallback
```,
  caption: [Paired-end overlap detection. Scans insert sizes in
    decreasing order (conservative: smallest trim first) and
    accepts the first admissible overlap. Returns 0 when no
    overlap is found; the caller applies sequence matching as a
    fallback in that case.],
) <alg:adp-pe>

#Prop (PE complexity). $"DetectPEOverlap"$ costs $O(L^2)$ in the
worst case (no overlap found), $O(L dot I^*)$ in the typical case
where an overlap is found at insert size $I^*$. For a clean read
($I >= L$) the loop body never executes and the cost is $O(L)$
for the outer loop alone.

_Note._ The $O(L^2)$ case is the _common_ case, not a pathological
one, and an earlier draft of this document claimed otherwise. A
clean long-insert pair admits no overlap at any $I$, so the loop
runs over all $min(L_1, L_2) - o_"PE"$ offsets. Each offset exits
only once it has accumulated $ceil(I epsilon) + 1$ mismatches, and
on unrelated sequence mismatches arrive at rate $3 \/ 4$, so the
expected cost per offset is $(ceil(I epsilon) + 1) dot 4 \/ 3$
comparisons rather than $O(1)$. For $L = 150$, $epsilon = 0.1$ and
$o_"PE" = 30$ this is roughly $120 times 21 approx 2500$ byte
comparisons per pair, incurred on every pair that is not adapter
contaminated, and followed by the SE fallback.

#Rem (vectorised mismatch counting). The asymptotic bound is
unchanged, but the constant is not. Evaluating $"mm"_"PE" (I)$ one
byte at a time places a data-dependent branch on the critical path
at every position. Comparing $32$ bases at a time against the
reverse-complemented mate and accumulating population counts, with
the saturation test hoisted out of the inner loop and applied once
per vector, reduces the per-offset cost by more than an order of
magnitude while returning exactly the same $"mm"_"PE" (I)$. This is
the fix to prefer, because it is a pure constant-factor change and
carries no correctness argument.

#Rem (probe gate, heuristic). A sampled probe, in which a fixed
subset of positions is compared before the full evaluation and the
offset is abandoned when too many disagree, removes the quadratic
term outright but is _not_ sound: an admissible offset may carry up
to $ceil(I epsilon)$ mismatches, and nothing prevents them from
falling inside the sample, so no rejection threshold below the
sample size is safe in the worst case. Its false-negative rate is
governed by the true sequencing error rate $p$ rather than by the
tolerance $epsilon$; for an $8$-position probe rejecting on two or
more disagreements and $p = 10^(-2)$, a correct overlap is lost
with probability $1 - (1-p)^8 - 8p(1-p)^7 approx 3 times 10^(-3)$.
That is defensible for a throughput-oriented default and
indefensible as a silent one. If adopted it must be exposed as a
flag, excluded from conformance testing against @alg:adp-pe, and
reported in the statistics block.

== Auto-detection <sec:adp-auto>

When no adapter set is specified, the pipeline samples the first
$N_"peek"$ reads (default $N_"peek" = 10^4$) into an in-memory
buffer before processing begins. A frequency analysis of 12-mers
at the 3′ end of each sampled read identifies candidate adapter
sequences without reference to a known database.

#Def (k-mer frequency analysis). For $k = 12$ and a sample of
$N_"peek"$ reads, extract the terminal 12-mer from each read of
length $>= 12$ and build a frequency table. A 12-mer is a
_candidate_ if it appears in at least $f_"min" := ceil(N_"peek"
dot 0.01)$ reads (i.e. frequency $>= 1%$). The top three
candidates by frequency (ties broken lexicographically) are
returned as adapters named "detected".

#Rem (rationale for 12-mers). Twelve bases provides sufficient
specificity to distinguish adapter signal from random sequence
while remaining short enough to tolerate a few cycles of adapter
degradation at the 3′ end. At 1% frequency in 10 k reads, a
candidate must appear $>= 100$ times, which exceeds the expected
random-sequence frequency by a factor of $4^{12} / 10^4 approx
1.7 times 10^4$ — false positive rate is negligible.

#Rem (PE auto-detection). For paired-end input, 12-mer analysis
is applied to R1 3′ ends. The PE overlap procedure of @sec:adp-pe
is independent of the adapter set and runs unconditionally;
auto-detected adapters serve only as a fallback for pairs where
no overlap is found.

#Rem (peek buffer). The $N_"peek"$ records are held in a
$"VecDeque"$ and processed before any records from the main
reader, so the auto-detected adapter set is available for all
reads in the file and latency is bounded by the cost of reading
$N_"peek"$ compressed records ($O(N_"peek" dot L)$).

== Integration with the polyX pipeline <sec:adp-integration>

#Def (adapter hint). At the call site in the pipeline, an
_adapter hint_ disambiguates the two operating modes:

$ "Hint" := cases(
  "Detect" quad &"(SE sequence matching from" cal(A) ")",
  "Override"(T_5, T_3, eta) &"(pre-computed trim from PE overlap)".
) $

#block(inset: 0.6em, stroke: 0.5pt + gray, radius: 2pt, [
  *Current implementation.* The hint has one inhabitant: every call site passes
  $"Detect"$, and the $"Override"$ branch was never constructed --- it existed
  only to carry a trim from @sec:adp-pe. The implementation therefore takes no
  hint parameter at all, which is why @sec:adp-pe is marked not implemented.
])

#Def (combined trim procedure). The public entry point $"Trim"$
takes a read $(R, Q)$, a parameter set $Theta$ (§1.2), and a hint,
and returns a $"TrimOutput"$ consisting of a $"TrimResult"$ and an
optional adapter name $eta$:

+ _Resolve adapter trim._ If $"Hint" = "Override"(T_5^"adp",
  T_3^"adp", eta)$, use those values directly. Otherwise run
  $"MatchAdapters"$ to obtain $(T_3^"adp", eta)$ with $T_5^"adp"
  = 0$.

+ _Adapter budget check._ Let $T_"adp" := T_5^"adp" +
  T_3^"adp"$. If $T_"adp" > Delta_L$, return
  $"Drop"("Adapter")$ immediately.

+ _PolyX scan on sub-read._ Let $(R', Q')$ be the adapter-trimmed
  sub-read of length $L' = L - T_"adp"$. Run $"PolyXTrim"$ on
  $(R', Q')$ with modified budget $Delta_L' := Delta_L -
  T_"adp"$ and $L_"orig" := L$, so that the min-len invariant is
  enforced relative to the original read length.

+ _Coordinate reconstruction._ Let the polyX result be
  $(T_5^"px", T_3^"px", "runs")$. Return
  $ T_5 = T_5^"adp" + T_5^"px", quad
    T_3 = T_3^"adp" + T_3^"px", $
  with interior runs translated back to original coordinates by
  adding $T_5^"adp"$ to each run's start position.

The procedure is shown in @alg:trim.

#figure(
```
procedure Trim(R, Q, L, Θ, Hint) → TrimOutput
    // Step 1: resolve adapter trim
    if Hint = Override(T5a, T3a, η) then
        T5_adp ← T5a;  T3_adp ← T3a;  name ← η
    else
        (T3_adp, name) ← MatchAdapters(R, L, 𝒜, ε, o_SE)
        T5_adp ← 0
    T_adp ← T5_adp + T3_adp
    // Step 2: adapter budget check
    if T_adp > Δ_L then return Drop(Adapter)
    // Step 3: polyX scan on adapter-trimmed sub-read
    safe5 ← min(T5_adp, L);  safe3 ← min(T3_adp, L − safe5)
    R' ← R[safe5 : L − safe3];  Q' ← Q[safe5 : L − safe3]
    Θ' ← Θ with Δ_L' ← Δ_L − T_adp, L_orig ← L
    res ← PolyXTrim(R', Q', L − T_adp, Θ')
    // Step 4: coordinate reconstruction
    match res with
    | Drop(reason, T5px, T3px) →
        return Drop(reason, T5_adp, T3_adp, T5px, T3px)
    | Emit(T5px, T3px, runs) →
        abs_runs ← [(start + safe5, len) for (start, len) in runs]
        return Emit(T5_adp, T3_adp, T5px, T3px, abs_runs, name)
```,
  caption: [Combined adapter + polyX trim. The adapter trim is
    applied first; any remaining budget is passed to the polyX
    procedure. Interior-run start positions are translated back
    to original-read coordinates by adding the 5′ adapter trim.
    The `Drop(Adapter)` path exits without scanning.],
) <alg:trim>

#Prop (budget invariant). For any input $(R, Q, L)$ and any
admissible hint, the final kept-read length satisfies
$L - T_5 - T_3 >= L_min$.

_Proof._ The adapter budget check ensures $T_"adp" <= Delta_L$.
The polyX procedure receives budget $Delta_L' = Delta_L -
T_"adp" >= 0$ and enforces $T_5^"px" + T_3^"px" + H_"int" <=
Delta_L'$ at the decision rule (§1.8). Adding the adapter
contributions:
$ T_5 + T_3 + H_"int" = T_"adp" + T_5^"px" + T_3^"px" + H_"int"
  <= T_"adp" + Delta_L' = Delta_L = L - L_min. $
Hence $L - T_5 - T_3 >= L_min$. $qed$

#Prop (correctness of coordinate reconstruction). Let $(s_i,
K_i^*)$ be an interior run recorded by $"PolyXTrim"$ in the
adapter-trimmed sub-read frame. Its start in the original-read
frame is $s_i + T_5^"adp"$, and the run occupies original-read
positions $s_i + T_5^"adp", dots, s_i + T_5^"adp" + K_i^* - 1$.
These positions lie strictly inside $[T_5, L - 1 - T_3]$
whenever the reconciliation step of @alg:main has accepted them.

_Proof._ Direct from the definition: the sub-read index
$j = 0$ corresponds to original-read index $T_5^"adp"$, so
position $j$ in the sub-read maps to $j + T_5^"adp"$ in the
original. Reconciliation in the sub-read frame removes runs
straddling the 3′ boundary of the sub-read; those same bases
straddle $L - T_3^"adp" - T_3^"px" - 1$ in the original read,
which is $L - 1 - T_3^"out"$, consistent with the reconciliation
condition of @sec:algorithm. $qed$

== Default adapter presets <sec:adp-presets>

Three named presets are provided:

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    align: (left + horizon, left + horizon, left + horizon),
    table.header(
      [*Preset*], [*Adapters included*], [*Typical use*],
    ),
    [`illumina`],
    [TruSeq R1, TruSeq R2, Nextera],
    [Standard Illumina TruSeq and Nextera/NextEra libraries],

    [`nextera`],
    [Nextera only],
    [Nextera rapid-capture, ATAC-seq, low-input WGS],

    [`small-rna`],
    [Illumina Small RNA 3′ adapter],
    [miRNA, sRNA, CLIP-seq],
  ),
  caption: [Built-in adapter presets. The `illumina` preset covers
    the large majority of standard short-read libraries. Omitting
    a preset activates auto-detection (@sec:adp-auto).],
) <tab:presets>

When a FASTA file is supplied via `--adapters`, its sequences are
loaded as-is; sequence IDs become adapter names. Presets and
auto-detection are mutually exclusive with `--adapters`. Setting
`--no-adapter` disables adapter trimming entirely (equivalent to
$cal(A) = emptyset$), which may be appropriate when upstream
tools (e.g. `cutadapt`, `fastp`) have already performed adapter
removal.

== Updated pipeline composition

Incorporating the adapter step, the complete single-read pipeline
of §5.3 becomes:

+ Fingerprint the anchor window and drop on a repeat fingerprint
  (@alg:dedup).
+ Decode Phred ASCII offset to obtain raw $Q$.
+ Compute $overline(q)$. If $overline(q) < theta_q$, drop.
+ Run $"Trim"(R, Q, L, Theta, "Hint")$ (@alg:trim):
  + Resolve adapter trim (sequence matching or PE override).
  + Adapter budget check; drop on $T_"adp" > Delta_L$.
  + Apply hard trim per @alg:hardtrim on the adapter-trimmed sub-read.
  + Run $"PolyXTrim"$ with reduced budget $Delta_L'$.
  + Reconstruct absolute coordinates.
+ Emit the kept slice $R[T_5 : L - T_3]$.

For paired-end input, step 3 is preceded by
$"DetectPEOverlap"$ on the pair; the result is passed as an
$"Override"$ hint to the $"Trim"$ calls for each mate.


= Rust implementation: environment and architecture <sec:rust-env>

== Reproducible build environment

The implementation is distributed with a `pixi.toml`@pixi that pins the entire
build to byte-level reproducibility on a per-platform basis. The relevant
fragment:

```toml
[project]
name = "fastprep"
version = "0.1.0"
channels = ["conda-forge"]
platforms = ["linux-64", "osx-arm64"]

[dependencies]
rust            = "~=1.83"
pkg-config      = "*"
cmake           = "*"
zlib            = "*"
zstd            = "*"
cargo-nextest   = "*"
cargo-flamegraph = "*"
hyperfine       = "*"

[tasks]
build      = "cargo build --release"
test       = "cargo nextest run --workspace"
bench      = "hyperfine -w 3 'target/release/rpolyx -1 testdata/sample.fq.gz -o /tmp/out.fq.gz -j 8'"
flamegraph = "cargo flamegraph --bin rpolyx --release -- -1 testdata/sample.fq.gz -o /tmp/out.fq.gz -j 8"
```

Bootstrap with `pixi install`, then run tasks through `pixi run build`,
`pixi run test`, etc. The lockfile `pixi.lock` records exact versions of
every transitive dependency, including the Rust compiler and the C ABI of
zlib/zstd.

The Rust-level dependencies are declared in `Cargo.toml`:

```toml
[dependencies]
anyhow            = "1"
clap              = { version = "4", features = ["derive"] }
crossbeam-channel = "0.5"
flate2            = { version = "1", features = ["zlib-ng"], default-features = false }
mimalloc          = { version = "0.1", default-features = false }
needletail        = "0.7.3"
rayon             = "1"
serde             = { version = "1", features = ["derive"] }
serde_json        = "1"
zstd              = "0.13"

[profile.release]
lto           = "fat"
codegen-units = 1
debug         = true
```

Notable choices: `flate2` with `zlib-ng` for 1.5--2× faster gzip decode;
`mimalloc`@leijen2019mimalloc as the global allocator for per-thread bucket
reuse on the many small `Vec<u8>` allocations the pipeline creates;
`needletail`@needletail for FASTX parsing@cock2010fastq with transparent
gzip/zstd/xz/bz2 auto-detection; `rayon`@rayon and
`crossbeam-channel` for parallelism and backpressure; `lto = "fat"` with
single-codegen-unit linking to enable cross-crate inlining of the hot path
into the chunk filler.

== Crate layout

```text
src/
├── lib.rs           # re-exports public surface
├── main.rs          # CLI entry point (rpolyx)
├── bin/pxmt.rs      # metadata viewer binary
├── algorithm.rs     # core polyX + scan primitives
├── adapter.rs       # adapter matching + auto-detect
├── pipeline.rs      # single-end / paired-end orchestration
├── cli.rs           # clap-based argument parsing + validation
├── io.rs            # FASTQ/FASTA output + compression
├── metadata.rs      # binary per-read trim record format (PXMT)
└── stats.rs         # aggregated statistics + reporting
```

The split isolates pure algorithms (`algorithm`, `adapter`) from I/O concerns
(`io`, `pipeline`). The pure modules are tested in isolation; the I/O modules
use the pure modules as opaque oracles.

== Type-level encoding of the mathematical model

Each formal object from §2 has a direct Rust type. The correspondence is
deliberate: every variable name in the math maps to a struct field or enum
variant, and every parameter in the tuple $(L_min, Delta_m, t, x, mu, nu, m_h)$
appears as a named field in `Params`.

=== `MinLen`: encodes $L_min$

The spec allows both an absolute base count and a fraction of read length, so
the type is a sum:

```rust
/// Minimum post-trim length, either absolute or fractional.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MinLen {
    Abs(usize),
    Frac(f64),
}

impl MinLen {
    /// L_min = min(N, L) for absolute, floor(f * L) for fractional.
    pub fn resolve(self, l: usize) -> usize {
        match self {
            MinLen::Abs(n)  => n.min(l),
            MinLen::Frac(f) => ((f * l as f64).floor() as usize).min(l),
        }
    }
}
```

`FromStr` accepts strings like `"100"` (absolute) or `"0.6L"` (fractional).

=== `MatchRule`: encodes the choice between $M_s$ and $M_s^X$

```rust
#[derive(Debug, Clone)]
pub enum MatchRule {
    /// Generic: match iff b[k] == b[k+1].
    Homopolymer,
    /// Anchored: match iff b[k+1] is in the target set.
    /// Box<[bool; 256]> for one-comparison lookup regardless of
    /// how many letters were specified.
    Anchored(Box<[bool; 256]>),
}
```

=== `ScanResult`: output of $"ScanFrom"$

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanResult {
    pub k_star:  usize,
    pub mode:    StopMode,        // Mercy, Boundary, or Budget
    pub matches: usize,           // m*: matching steps; run length = matches + 1
}
```

The `matches` field stores the match count $m^*$. The downstream resolver
thresholds against $m^* + 1 >= m_h$ rather than against $K^*$, matching the
spec §1.7 rationale.

=== `Params`: full parameter tuple plus precomputed gate LUT

```rust
#[derive(Debug, Clone)]
pub struct Params {
    pub min_len:           MinLen,        // L_min
    pub mercy:             f64,           // Δ_m
    pub phred_sat:         u8,            // t
    pub gate_shape:        f64,           // x
    pub match_leniency:    f64,           // μ
    pub mismatch_strength: f64,           // ν
    pub min_run:           usize,         // m_h
    pub rule:              MatchRule,
    pub scan_5:            bool,
    pub scan_3:            bool,
    pub scan_interior:     bool,
    pub direction:         Direction,
    pub mean_q:            Option<f64>,   // θ_q
    pub trim_5:            Option<usize>, // T_5^hard
    pub trim_3:            Option<usize>, // T_3^hard

    /// Arc-wrapped so per-mate clones share the allocation.
    pub adapter:           Option<Arc<AdapterSet>>,
    /// Trim budget already consumed by a prior step (adapter trim).
    pub extra_committed:   usize,
    /// Override length for Δ_L computation: used when polyx_trim runs
    /// on an adapter-sliced sub-read so min-len is enforced relative
    /// to the original read length.
    pub original_len:      Option<usize>,
    /// Precomputed g(q) for q in [0, 255]. Populated by prepare().
    /// Eliminates the powf call from the hot loop.
    #[doc(hidden)]
    pub gate_lut:          [f64; 256],
}
```

=== `TrimResult`: encodes the Drop/Emit decision

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TrimResult {
    Drop {
        reason:     DropReason,
        t5_adapter: usize, t3_adapter: usize,
        t5_poly:    usize, t3_poly:    usize,
        runs:       Vec<(usize, usize)>,
    },
    Emit {
        t5_adapter: usize, t3_adapter: usize,
        t5_poly:    usize, t3_poly:    usize,
        runs:       Vec<(usize, usize)>,
    },
}
```

Both variants carry separate adapter and polyX trim counts, so the pipeline
can attribute every base precisely. `runs` lists interior runs as
`(start, len)` pairs in absolute 0-based coordinates of the original read,
after adapter coordinates have been reconstructed.

=== Adapter detection has a single mode

The current implementation has no dispatcher type: sequence matching against
`params.adapter` is the only adapter-detection path, and `trim` takes no hint
parameter. The `AdapterHint` enum described in earlier revisions of this
document --- whose `Override` variant carried a pre-computed trim from the PE
overlap geometry of @sec:adp-pe --- was removed together with that mode. It was
matched at exactly one call site and never constructed, because every caller
passed `Detect`.

== Concurrency primitives

The pipeline composes three layers of concurrency:

#set enum(numbering: "1.")
+ _Bounded channels_ (`crossbeam_channel::bounded`) connect the reader thread
  to the compute thread, and the compute thread to the writer (main thread).
  The channel capacity of 2 creates backpressure: the reader cannot run more
  than 2 chunks ahead, bounding pipeline memory in the chunk size.
+ _Rayon thread pool_ (`rayon::ThreadPoolBuilder`) parallelises within each
  chunk via `into_par_iter()`. The pool is configured globally in `main.rs`
  or per-test via a scoped pool.
+ _Arc sharing_ (`Arc<AdapterSet>`, `Arc<Params>`) avoids per-thread cloning
  of read-only state. The atomic refcount is incremented once per chunk, not
  once per read.

== RAII guards for crash safety

Every output file is opened to a `.part` temporary path. The `TmpFile` RAII
guard deletes the `.part` file on `Drop` unless `commit()` was called:

```rust
struct TmpFile {
    path:      std::path::PathBuf,
    committed: bool,
}

impl TmpFile {
    fn commit(mut self) { self.committed = true; }
}

impl Drop for TmpFile {
    fn drop(&mut self) {
        if !self.committed {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}
```

On a clean run, after the writer is flushed and the encoder finalised (zstd's
`auto_finish` finalises the frame on drop), the `.part` file is renamed to
the destination via `rename(2)`, which is atomic on the same filesystem. On
any error path (panic, OOM, SIGKILL caught by the runtime, I/O failure
short-circuiting `?`), the guard removes the partial file. A subsequent
invocation starts from a clean slate.

= Rust implementation: hot path and pipeline <sec:rust-impl>

This section presents the performance-critical code paths with explicit
cross-references to the mathematical objects of @sec:intro through
@sec:adapter.

== Quality gate lookup table

The closed-form quality gate of §1.3 admits a complete tabulation over the
range of admissible Phred values, since $q$ is bounded by 256. The LUT is
populated once per `Params`:

```rust
/// Quality gate g(q). Direct correspondence to spec §1.3.
/// Called only at LUT-init; the hot path reads from gate_lut.
#[inline]
fn quality_gate(q: u8, t: u8, x: f64) -> f64 {
    if q == 0       { 0.0 }
    else if q >= t  { 1.0 }
    else            { (q as f64 / t as f64).powf(x) }
}

impl Params {
    /// Populate gate_lut[q] = g(q) for q in [0, 255].
    /// Called by Default::default() and by cli::Cli::finalize().
    pub fn prepare(&mut self) {
        for q in 0u8..=255 {
            self.gate_lut[q as usize] = quality_gate(q, self.phred_sat, self.gate_shape);
        }
    }
}
```

This is the single most important micro-optimisation in the pipeline.
The hot loop calls $g$ twice per base pair; for 150-bp reads at $10^7$ reads/s,
that is $3 dot 10^9$ `powf` evaluations per second, well beyond any CPU's
floating-point throughput. With the LUT, the hot loop performs two
byte-indexed array loads.

== Per-step gap and match indicator

The per-step gap function $1 - sigma_s (k)$ from §1.6 reads as:

- match branch: $1 - sigma = mu (1 - W) in [0, mu]$
- mismatch branch: $1 - sigma = 1 + nu W in [1, 1 + nu]$

Implementation:

```rust
/// Per-step gap + match flag for one (b1, b2, q1, q2) quadruple.
/// Hot path. Per-call cost: 2 LUT loads (q1, q2) + 1 fadd + 1 fmul
/// + 1 branch. Quality bytes are RAW Phred (pipeline pre-decoded
/// the ASCII offset on ingestion).
#[inline]
fn per_step_gap_and_match(
    b1: u8, b2: u8, q1: u8, q2: u8,
    params: &Params,
    rule: &MatchRule,
) -> (f64, bool) {
    // M_s(k): match indicator. Homopolymer = byte equality;
    // Anchored = presence of b2 in target set.
    let matches = match rule {
        MatchRule::Homopolymer    => b1 == b2,
        MatchRule::Anchored(mask) => mask[b2 as usize],
    };
    // W_s(k) = (g(q1) + g(q2)) / 2 via LUT.
    let w = (params.gate_lut[q1 as usize] + params.gate_lut[q2 as usize]) * 0.5;
    let gap = if matches {
        // σ = (1 − μ) + μW  =>  gap = μ(1 − W)
        params.match_leniency * (1.0 - w)
    } else {
        // σ = −νW  =>  gap = 1 + νW
        1.0 + params.mismatch_strength * w
    };
    (gap, matches)
}
```

The two comment lines `σ = ...` link directly to the per-step score definition
of §1.5; the equality `gap = 1 - σ` collapses to the simple expressions
implemented.

#Rem (branchless form). The `if matches` selection above is a
data-dependent branch whose outcome tracks the base composition of the
read. On non-homopolymeric sequence it is close to unpredictable, and a
mispredict costs more than the arithmetic it guards. Both branches are
affine in $W$, so the selection can be replaced by a two-entry table
lookup: with `ADD = [1.0, μ]` and `MUL = [ν, -μ]` indexed by
`(b1 == b2) as usize`, the gap is `ADD[m] + MUL[m] * w`, one fused
multiply-add and no branch.

#Rem (accumulator width). $Gamma$ is the sole loop-carried dependency, so
the loop cannot retire faster than the latency of its accumulation. In
`f64` that is a four-cycle floating-point add; in the $u 16$ fixed point
already recommended in @sec:scoretable it is a single-cycle integer add.
Adopting fixed point here is therefore worth up to a factor of four on the
critical path before any vectorisation, and combines with the binning of
@sec:scoretable to make the whole per-step computation a table load and an
integer add.

== General scan primitive

The general scan walks from `start` toward `upper`, accumulating gap mass
into `gamma` until mercy, boundary, or (optional) budget fires. This is the
direct encoding of @alg:scan, with the optional budget cap added to support
the 5′ short-circuit of @sec:shortcircuit:

```rust
pub fn scan(
    bases: &[u8], quals: &[u8],
    start: usize, upper: usize,
    direction: Direction,
    params: &Params,
    budget: Option<usize>,
) -> ScanResult {
    debug_assert_eq!(bases.len(), quals.len());
    debug_assert!(upper < bases.len());
    debug_assert!(start <= upper);

    let n = bases.len();
    let mut gamma   = 0.0f64;
    let mut k       = 0usize;
    let mut matches = 0usize;

    loop {
        // Boundary clause: next pair would need start + k + 1 > upper.
        if start + k + 1 > upper {
            return ScanResult { k_star: k, mode: StopMode::Boundary, matches };
        }
        // Budget cap (instance of the 5' short-circuit).
        if let Some(b) = budget {
            if k >= b {
                return ScanResult { k_star: k, mode: StopMode::Budget, matches };
            }
        }
        let i1 = idx(start + k,     n, direction);
        let i2 = idx(start + k + 1, n, direction);
        let (gap, is_match) = per_step_gap_and_match(
            bases[i1], bases[i2], quals[i1], quals[i2], params, &params.rule,
        );
        gamma += gap;
        if is_match { matches += 1; }
        k += 1;
        // Mercy clause: Γ_s(K) > Δ_m.
        if gamma > params.mercy {
            return ScanResult { k_star: k, mode: StopMode::Mercy, matches };
        }
    }
}
```

The `idx` helper abstracts forward versus reverse traversal:

```rust
#[inline]
fn idx(i: usize, n: usize, dir: Direction) -> usize {
    match dir {
        Direction::Forward => i,
        Direction::Reverse => n - 1 - i,
    }
}
```

In the main pipeline this branch is removed entirely by const-generic
monomorphisation (next section).

== Fused streaming walk

The naive implementation performs two passes over the read: one to sum
qualities for the mean-quality filter, and one to run the polyX scan. The
production implementation fuses these into a single byte walk, and
additionally monomorphises the direction parameter so the per-iteration
direction check disappears at compile time.

```rust
/// State returned from the streaming forward pass.
struct StreamOut {
    quality_sum:       u64,
    first_trim:        usize,
    h_int:             usize,
    runs:              Vec<(usize, usize)>,
    /// 5' short-circuit fired during the walk.
    dropped_by_budget: bool,
}

/// Inner streaming walk over [0, l). Monomorphised on direction:
/// idx_const::<true> collapses to i, idx_const::<false> to l - 1 - i.
/// The compiler emits two specialised functions, one per dispatch arm
/// in polyx_trim_inner.
#[inline]
fn streaming_walk<const FORWARD: bool>(
    bases:      &[u8],
    quals:      &[u8],
    params:     &Params,
    win_start:  usize,
    win_end:    usize,
    rem_budget: usize,
) -> StreamOut {
    let l = bases.len();
    let mut quality_sum:        u64   = 0;
    let mut gamma:              f64   = 0.0;
    let mut k_curr:             usize = 0;
    let mut m_curr:             usize = 0;
    let mut sub_start_logical:  usize = win_start;
    let mut first_subscan             = true;
    let mut first_trim:         usize = 0;
    let mut h_int:              usize = 0;
    let mut runs:               Vec<(usize, usize)> = Vec::new();

    #[inline(always)]
    fn idx_const<const FWD: bool>(i: usize, n: usize) -> usize {
        if FWD { i } else { n - 1 - i }
    }

    #[inline]
    fn resolve_subscan(k: usize, m: usize, mode: StopMode, min_run: usize) -> usize {
        if m + 1 < min_run              { 0 }
        else if mode == StopMode::Mercy { k }
        else                            { k + 1 }
    }

    for i_logical in 0..l {
        let i_phys = idx_const::<FORWARD>(i_logical, l);

        // Mean-q accumulation: ALWAYS over the whole read.
        quality_sum += quals[i_phys] as u64;

        // PolyX window: pair (b[i-1], b[i]) requires
        // i_logical in [win_start + 1, win_end - 1].
        if i_logical < win_start + 1 || i_logical >= win_end { continue; }

        let prev_phys = idx_const::<FORWARD>(i_logical - 1, l);
        let (gap, is_match) = per_step_gap_and_match(
            bases[prev_phys], bases[i_phys],
            quals[prev_phys], quals[i_phys],
            params, &params.rule,
        );
        gamma += gap;
        if is_match { m_curr += 1; }
        k_curr += 1;

        // 5' budget short-circuit (@sec:shortcircuit).
        if first_subscan && rem_budget > 0 && k_curr >= rem_budget {
            return StreamOut {
                quality_sum,
                first_trim: 0, h_int: 0, runs: Vec::new(),
                dropped_by_budget: true,
            };
        }

        // Mercy clause fired: close out current sub-scan.
        if gamma > params.mercy {
            let sub_trim = resolve_subscan(k_curr, m_curr, StopMode::Mercy, params.min_run);
            if first_subscan {
                // This was the 5' end-scan: record as T_5.
                first_trim = sub_trim;
                first_subscan = false;
            } else if sub_trim > 0 && params.scan_interior {
                // Interior qualifying run.
                let sub_pos = sub_start_logical - win_start;
                h_int += sub_trim;
                runs.push((sub_pos, sub_trim));
                if first_trim + h_int > rem_budget {
                    return StreamOut {
                        quality_sum, first_trim, h_int, runs,
                        dropped_by_budget: true,
                    };
                }
            }
            // Reset state for next sub-scan.
            gamma = 0.0; k_curr = 0; m_curr = 0;
            sub_start_logical = i_logical + 1;
        }
    }

    // End of walk: if the 5' end-scan never mercy-terminated, resolve via boundary.
    if k_curr > 0 && first_subscan {
        first_trim = resolve_subscan(k_curr, m_curr, StopMode::Boundary, params.min_run);
    }

    StreamOut { quality_sum, first_trim, h_int, runs, dropped_by_budget: false }
}
```

#Rem (why monomorphisation matters). The previous version of this code used
a runtime `Direction` enum and an `idx()` function with a `match` inside the
inner loop. Perf-annotated profiles showed that branch alone consuming
roughly 10% of total runtime, not because the branch itself is expensive,
but because it prevented the compiler from hoisting the bounds check on
`bases[i_phys]` and from coalescing loads. With the `const FORWARD: bool`
parameter, the compiler emits two specialised functions; `idx_const`
inlines to a constant offset or `n - 1 - i`, and the bounds checks become
single comparisons that LLVM hoists out of the loop.

#Rem (fusion versus prefiltering). The fusion above was measured against
two _scalar_ passes, and against that baseline it wins. It is a poor
baseline. The mean-quality sum is a reduction over a byte array with no
loop-carried branch, so it vectorises to roughly one vector add per
$32$ bytes, whereas fusing it into the scan forces the scalar per-step
machinery across every base of every read purely to obtain that sum. Once
the prefilter of @sec:prefilter is in place the two components have
opposite cost profiles: the mean-quality pass is unconditional and cheap
when vectorised, while the polyX scan is expensive and rarely needed.
Separating them is then the faster arrangement, and it also frees the scan
order, allowing the 3′ end to be scanned first on two-colour chemistry
where the polyG tail is the expected artefact. The fusion should be
re-measured against a vectorised mean-quality pass before being retained.

#Rem (early exit on mean quality). The filter of @sec:extras is evaluated
over the whole read, but its outcome is often decidable early: if
$"partial" + q_"max" (L - i) < theta_q L$ after $i$ positions, no
continuation can lift the mean above the threshold and the read can be
dropped without reading the remainder. With $q_"max" = 41$ this fires on
reads that are already condemned, which are exactly the reads for which
the remaining work would be wasted.

== Main pipeline

`polyx_trim` is the public entry point for the polyX core; internally it
calls `polyx_trim_inner` with explicit `extra_committed` and `original_len`
overrides, so the adapter wrapper can invoke it without cloning `Params`
(avoiding an `Arc` bump per read, roughly $1.5 dot 10^7$ atomic ops on a
typical Illumina lane).

```rust
pub fn polyx_trim(bases: &[u8], quals: &[u8], params: &Params) -> TrimResult {
    polyx_trim_inner(bases, quals, params, params.extra_committed, params.original_len)
}

fn polyx_trim_inner(
    bases: &[u8], quals: &[u8],
    params: &Params,
    extra_committed: usize,
    original_len: Option<usize>,
) -> TrimResult {
    assert_eq!(bases.len(), quals.len());
    let l = bases.len();

    if l < 2 {
        return TrimResult::Emit {
            t5_adapter: 0, t3_adapter: 0, t5_poly: 0, t3_poly: 0, runs: vec![],
        };
    }

    // Δ_L = L_eff − L_min − extra_committed  (§1.2).
    let delta_l = {
        let l_eff = original_len.unwrap_or(l);
        l_eff
            .saturating_sub(params.min_len.resolve(l_eff))
            .saturating_sub(extra_committed)
    };

    // Hard trims, capped to read length.
    let t_front = params.trim_5.unwrap_or(0).min(l);
    let t_tail  = params.trim_3.unwrap_or(0).min(l - t_front);
    let sub_l   = l - t_front - t_tail;
    if sub_l < 2 || t_front + t_tail > delta_l {
        return TrimResult::Drop {
            reason: DropReason::Hard,
            t5_adapter: 0, t3_adapter: 0,
            t5_poly: t_front, t3_poly: t_tail, runs: vec![],
        };
    }
    let rem_budget = delta_l - (t_front + t_tail);

    // Direction::Forward → streaming pass owns physical 5', reverse pass owns 3'.
    // Direction::Reverse → flipped.
    let dir = params.direction;
    let (stream_enabled, reverse_enabled) = match dir {
        Direction::Forward => (params.scan_5, params.scan_3),
        Direction::Reverse => (params.scan_3, params.scan_5),
    };

    // Dispatch on const generic. Two specialised functions; the branch
    // here is a single resolved jump, no per-iteration cost inside the walk.
    let stream_out = if !stream_enabled {
        StreamOut {
            quality_sum: quals.iter().map(|&q| q as u64).sum(),
            first_trim: 0, h_int: 0, runs: Vec::new(),
            dropped_by_budget: false,
        }
    } else if matches!(dir, Direction::Forward) {
        streaming_walk::<true >(bases, quals, params, t_front, l - t_tail, rem_budget)
    } else {
        streaming_walk::<false>(bases, quals, params, t_front, l - t_tail, rem_budget)
    };

    if stream_out.dropped_by_budget {
        let reason = match dir {
            Direction::Forward => DropReason::Poly5,
            Direction::Reverse => DropReason::Poly3,
        };
        return TrimResult::Drop {
            reason, t5_adapter: 0, t3_adapter: 0,
            t5_poly: t_front, t3_poly: t_tail, runs: vec![],
        };
    }
    let first_trim   = stream_out.first_trim;
    let h_int        = stream_out.h_int;
    let mut runs     = stream_out.runs;
    let quality_sum  = stream_out.quality_sum;

    // Mean-quality filter (§5).
    if let Some(threshold) = params.mean_q {
        let mean = quality_sum as f64 / l as f64;
        if mean < threshold {
            return TrimResult::Drop {
                reason: DropReason::MeanQ,
                t5_adapter: 0, t3_adapter: 0,
                t5_poly: 0, t3_poly: 0, runs: vec![],
            };
        }
    }

    // Reverse pass (budget-capped). 3' end-scan in Forward mode,
    // 5' end-scan in Reverse mode.
    let sub_bases = &bases[t_front..l - t_tail];
    let sub_quals = &quals[t_front..l - t_tail];
    let upper     = sub_l - 1;
    let second_trim = if reverse_enabled {
        let second_dir = match dir {
            Direction::Forward => Direction::Reverse,
            Direction::Reverse => Direction::Forward,
        };
        let budget = rem_budget.saturating_sub(first_trim);
        let budget_opt = if budget > 0 { Some(budget) } else { None };
        let res = scan(sub_bases, sub_quals, 0, upper, second_dir, params, budget_opt);
        if res.mode == StopMode::Budget {
            let reason = match dir {
                Direction::Forward => DropReason::Poly3,
                Direction::Reverse => DropReason::Poly5,
            };
            let (t5p, t3p) = match dir {
                Direction::Forward => (first_trim, 0),
                Direction::Reverse => (0, first_trim),
            };
            return TrimResult::Drop {
                reason, t5_adapter: 0, t3_adapter: 0,
                t5_poly: t_front + t5p, t3_poly: t_tail + t3p, runs: vec![],
            };
        }
        resolve_trim(res, params.min_run)
    } else { 0 };

    // Reconciliation: remove interior runs straddling the 3' trim zone
    // (§3.2, monotone reverse-iterate-with-break).
    let mut h_int_adj = h_int;
    if second_trim > 0 {
        let cut = sub_l.saturating_sub(second_trim);
        let mut drop_count = 0usize;
        for &(start, len) in runs.iter().rev() {
            if start + len > cut {
                h_int_adj = h_int_adj.saturating_sub(len);
                drop_count += 1;
            } else { break; }
        }
        runs.truncate(runs.len() - drop_count);
    }

    // Final decision rule (§1.9).
    if first_trim + second_trim + h_int_adj > rem_budget {
        let (t5p, t3p) = match dir {
            Direction::Forward => (first_trim, second_trim),
            Direction::Reverse => (second_trim, first_trim),
        };
        let hard = t_front + t_tail;
        let max_poly = t5p.max(t3p).max(h_int_adj);
        let reason = if max_poly >= hard {
            if t5p >= t3p && t5p >= h_int_adj { DropReason::Poly5 }
            else if t3p >= h_int_adj          { DropReason::Poly3 }
            else                              { DropReason::PolyInt }
        } else { DropReason::Hard };
        return TrimResult::Drop {
            reason, t5_adapter: 0, t3_adapter: 0,
            t5_poly: t_front + t5p, t3_poly: t_tail + t3p, runs,
        };
    }

    // Translate interior-run positions into absolute coordinates.
    let (t5p, t3p) = match dir {
        Direction::Forward => (first_trim, second_trim),
        Direction::Reverse => (second_trim, first_trim),
    };
    let abs_runs: Vec<(usize, usize)> = runs.into_iter()
        .map(|(sub_pos, len)| match dir {
            Direction::Forward => (t_front + sub_pos, len),
            Direction::Reverse => (l - (t_front + sub_pos) - len, len),
        })
        .collect();

    TrimResult::Emit {
        t5_adapter: 0, t3_adapter: 0,
        t5_poly: t_front + t5p, t3_poly: t_tail + t3p,
        runs: abs_runs,
    }
}
```

`resolve_trim` implements the run-length threshold and boundary correction of
§1.8:

```rust
#[inline]
fn resolve_trim(res: ScanResult, min_run: usize) -> usize {
    debug_assert!(res.mode != StopMode::Budget);
    if res.matches + 1 < min_run         { 0 }
    else if res.mode == StopMode::Mercy  { res.k_star }
    else                                 { res.k_star + 1 } // boundary
}
```

== Combined adapter + polyX wrapper

The public `trim` entry point composes adapter detection with the polyX core,
implementing @alg:trim:

```rust
pub fn trim(bases: &[u8], quals: &[u8], params: &Params) -> TrimOutput {
    let l = bases.len();

    // --- Step 1: adapter detection (sequence matching) ---
    let (t5_adapter, t3_adapter, adapter_name) = match params.adapter.as_deref() {
        None => (0, 0, None),
        Some(set) => match crate::adapter::match_adapters(bases, set) {
            None => (0, 0, None),
            Some(m) => (m.t5, m.t3, Some(m.name)),
        },
    };

    // Step 2: adapter-only budget check.
    let committed = t5_adapter + t3_adapter;
    if committed > 0 {
        let delta_l = params.trim_budget(l);
        if committed > delta_l {
            return TrimOutput {
                result: TrimResult::Drop {
                    reason: DropReason::Adapter,
                    t5_adapter, t3_adapter,
                    t5_poly: 0, t3_poly: 0, runs: vec![],
                },
                adapter_name,
            };
        }
    }

    // Step 3: polyX scan on adapter-trimmed sub-slice.
    let safe_t5 = t5_adapter.min(l);
    let safe_t3 = t3_adapter.min(l.saturating_sub(safe_t5));
    let sub_bases = &bases[safe_t5..l - safe_t3];
    let sub_quals = &quals[safe_t5..l - safe_t3];

    // Pass committed/Some(l) explicitly: budget is enforced relative to
    // original L, not the sub-slice length. No Params clone.
    let poly_result = polyx_trim_inner(sub_bases, sub_quals, params, committed, Some(l));

    // Step 4: coordinate reconstruction.
    let result = match poly_result {
        TrimResult::Drop { reason, t5_poly, t3_poly, runs, .. } => TrimResult::Drop {
            reason, t5_adapter, t3_adapter, t5_poly, t3_poly,
            runs: runs.into_iter().map(|(p, l)| (p + safe_t5, l)).collect(),
        },
        TrimResult::Emit { t5_poly, t3_poly, runs, .. } => {
            let abs_runs = runs.into_iter()
                .map(|(p, l)| (p + safe_t5, l))
                .collect();
            TrimResult::Emit {
                t5_adapter, t3_adapter, t5_poly, t3_poly, runs: abs_runs,
            }
        }
    };
    TrimOutput { result, adapter_name }
}
```

== Adapter sequence matching

The listing below deviates from @alg:adp-se: it generalises to
leftmost-match-wins over internal start positions rather than restricting
to 3′ suffix alignments. By the suffix-locality argument of @sec:adp-se
that generalisation searches for an event the contamination model
excludes, and it raises the per-read cost from
$O(min(L, p) - o_"SE")$ candidate alignments to $O(L)$ of them, each
early-exiting only after $ceil(k epsilon) + 1$ mismatches. For $L = 150$,
$epsilon = 0.1$ and a three-adapter preset this is on the order of
$10^3$ byte comparisons per read, plausibly exceeding the cost of the
polyX scan it precedes. It is recorded here as the current state of the
code, not as the recommended form; the intended replacement is the seeded
suffix matcher of @sec:adp-se, which the listing should be reverted to.

```rust
pub fn match_adapters(seq: &[u8], set: &AdapterSet) -> Option<AdapterMatch> {
    let l = seq.len();
    let mut best: Option<(usize, AdapterMatch)> = None;

    for adapter in &set.adapters {
        let a = adapter.seq.len();
        if a == 0 { continue; }
        let max_p = l.saturating_sub(set.min_overlap_se);
        // Scan adapter start position left to right; leftmost = most trimming.
        for p in 0..=max_p {
            let k = a.min(l - p);
            if k < set.min_overlap_se { break; }
            // Admissibility: mm/k ≤ ε  =>  mm ≤ ceil(k·ε).
            let max_mm = (k as f64 * set.error_rate).ceil() as usize;
            let mut mm = 0usize;
            let mut matched = true;
            for i in 0..k {
                if seq[p + i] != adapter.seq[i] {
                    mm += 1;
                    // Early exit on saturation: avoids O(L) per position,
                    // gives O(L/ε) for non-matching adapters.
                    if mm > max_mm { matched = false; break; }
                }
            }
            if matched {
                let t3 = l - p;
                if best.as_ref().map_or(true, |(prev, _)| t3 > *prev) {
                    best = Some((t3, AdapterMatch {
                        name: adapter.name.clone(), t5: 0, t3,
                    }));
                }
                break;
            }
        }
    }
    best.map(|(_, m)| m)
}
```

== PE overlap detection

Direct encoding of @alg:adp-pe:

#block(inset: 0.6em, stroke: 0.5pt + gray, radius: 2pt, [
  *Specification only --- not present in `src/`.* As @sec:adp-pe states, this mode
  is not implemented in the current revision: it was removed when the six tools
  were unified under one CLI. `adapter.rs` carries no `detect_pe_overlap` and no
  `complement`; sequence matching is the single adapter path. The listing below
  is the direct encoding of the designed procedure, retained as the reference for
  the specification.
])

```rust
pub fn detect_pe_overlap(
    r1: &[u8], r2: &[u8],
    min_overlap: usize, error_rate: f64,
) -> Option<usize> {
    let l1 = r1.len();
    let l2 = r2.len();
    let max_insert = l1.min(l2);
    if max_insert < min_overlap { return None; }

    // Iterate I from largest insert down to o_PE (conservative: minimal trim).
    for insert in (min_overlap..=max_insert).rev() {
        // At least 1 mismatch allowed so very small overlaps with
        // ε = 0 don't auto-reject.
        let max_mm = ((insert as f64 * error_rate).ceil() as usize).max(1);
        let mut mm = 0usize;
        let mut found = true;
        for i in 0..insert {
            if r1[i] != complement(r2[insert - 1 - i]) {
                mm += 1;
                if mm > max_mm { found = false; break; }
            }
        }
        if found {
            let trim = l1.saturating_sub(insert);
            if trim == 0 { return None; } // read entirely insert, no adapter
            return Some(trim);
        }
    }
    None
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'G' | b'g' => b'C',
        b'C' | b'c' => b'G',
        _ => b'N',
    }
}
```

== Pipeline orchestration

The single-end pipeline is a producer/consumer chain with bounded channels
for backpressure:

```text
[reader thread]──raw_tx──>[compute thread]──res_tx──>[main thread, writes]
                              │
                              ▼
                          rayon pool
```

Code skeleton:

```rust
pub fn run_single(
    in1: &Path, out1: &Path,
    params: &Params, rt: &Runtime,
) -> Result<Stats> {
    // ... open input, set up .part outputs with TmpFile guards ...
    let (raw_tx, raw_rx) = bounded::<Vec<ScanInput>>(2);
    let (res_tx, res_rx) = bounded::<Vec<Processed>>(2);

    let params_arc = Arc::new(params);
    let compute = std::thread::spawn(move || {
        while let Ok(chunk) = raw_rx.recv() {
            // Parallel within the chunk via rayon.
            let results: Vec<_> = chunk.into_par_iter()
                .map(|input| {
                    let output = trim(&input.seq, &input.qual, &params_arc);
                    Processed { input, output }
                })
                .collect();
            if res_tx.send(results).is_err() { break; }
        }
    });

    // Main: read next chunk, push, await result of previous chunk, write.
    let mut chunk = Vec::with_capacity(rt.chunk);
    fill_chunk(&mut *reader, &mut chunk, rt.chunk, phred_offset)?;
    while !chunk.is_empty() {
        raw_tx.send(chunk).ok();
        let mut next = Vec::with_capacity(rt.chunk);
        fill_chunk(&mut *reader, &mut next, rt.chunk, phred_offset)?;
        let results = res_rx.recv().context("compute thread died")?;
        for p in &results {
            stats.record(/* ... */);
            p.emit_to(&mut writer, phred_offset)?;
        }
        chunk = next;
    }
    drop(raw_tx);
    compute.join().ok();
    // Flush, atomic rename, commit RAII guards.
    Ok(stats)
}
```

#Rem (pair alignment in PE mode). Rayon's
`into_par_iter().zip(other.into_par_iter())` builds a joined parallel iterator
that splits the combined index range across tasks. Each task receives a
contiguous slice $[i, j)$ and the invariant $R_1[i]$ paired with $R_2[i]$
is preserved regardless of thread count or work distribution. End-to-end,
output records appear in the same order as input records.

== Binary metadata format <sec:pxmt>

Per-read trim decisions are recorded in a compact binary format (PXMT) for
downstream inspection. Header is 4-byte magic `"PXMT"` plus 1-byte version.
Records are variable-length:

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    align: (left + horizon, left + horizon, left + horizon),
    table.header([*Field*], [*Type*], [*Meaning*]),
    [`id_len`],     [`u16`],        [byte length of read ID],
    [`id`],         [`[u8]`],       [raw ID, no terminator],
    [`mate`],       [`u8`],         [0 = R1/SE, 1 = R2],
    [`status`],    [`u8`],         [0 kept, 1 hard, 2 poly5, 3 poly3, 4 poly_int, 5 meanq, 6 adapter, 7 duplicate],
    [`t5_adapter`], [`u32`],        [adapter bases trimmed 5′],
    [`t3_adapter`], [`u32`],        [adapter bases trimmed 3′],
    [`t5_poly`],    [`u32`],        [polyX bases trimmed 5′],
    [`t3_poly`],    [`u32`],        [polyX bases trimmed 3′],
    [`read_len`],   [`u32`],        [original read length],
    [`n_runs`],     [`u8`],         [number of interior run entries],
    [run entries],  [`[(u32, u32)]`], [start, length pairs],
  ),
  caption: [PXMT record layout. Both kept (status 0) and dropped (status
    1--6) reads can carry interior runs; adapter contributions are
    recorded for every touched read, including drops.],
) <tab:pxmt>

The viewer binary `pxmt` decodes this format and supports filtering by ID
substring, drop reason, mate, and adapter-trim presence. Summary mode
holds $O(1)$ memory regardless of file size, suitable for multi-billion-read
metadata files from production lanes.

= Rust implementation: performance and safety <sec:rust-perf>

== Memory model

The pipeline performs three intentional heap allocations per read:

#set enum(numbering: "1.")
+ `seq: Vec<u8>` (bases, copied from needletail's internal buffer)
+ `qual: Vec<u8>` (qualities, copied with Phred ASCII offset decoded)
+ `id: Vec<u8>` (read ID, copied)

These copies are required because needletail invalidates its internal buffer
on each `next()` call, but the pipeline fans out chunks to worker threads
that must outlive any single buffer. The copies are amortised against the
much larger cost of parsing and (when compressed) decompressing the FASTQ.

#Rem (arena alternative). The amortisation argument holds per read but not
per lane: three allocations per read is of order $10^8$ allocations on a
$50$ M-pair run, and `mimalloc` reduces the cost of each rather than the
count. The copies are required, the allocations are not. A single
chunk-scoped byte arena, with records holding offsets and lengths into it,
preserves the outlive-the-parser property that motivates the copies while
reducing the allocation count to one per chunk. This has not been
implemented or measured.

Per-read shared state (`Params`, $approx 2.5$ kB including the gate LUT;
`AdapterSet`, $approx$ kB depending on adapter count) is wrapped in `Arc` and
shared across worker threads. The atomic refcount increments once per chunk,
not per read.

#Rem (allocator). `mimalloc` is enabled as the global allocator via
`#[global_allocator]`. On the hot path of many small `Vec<u8>` allocations,
mimalloc's per-thread bucket reuse outperforms glibc malloc by roughly
1.3--1.5× on typical Illumina workloads.

== Zero-cost abstractions

The pipeline uses three Rust patterns that the compiler eliminates entirely:

- _Const-generic monomorphisation_ (`streaming_walk<const FORWARD: bool>`).
  Produces two specialised functions; the direction branch is resolved at
  compile time, and the inner index expressions collapse to a constant offset
  or `n - 1 - i`.
- _Sealed enum dispatch_ (`TrimResult::{Drop, Emit}`, `StopMode::{Mercy,
  Boundary, Budget}`). Sealed inside the crate, so the compiler proves
  exhaustiveness in `match` arms and eliminates the discriminant check when
  the variant is known statically.
- _Inline LUT lookups_ (`params.gate_lut[q as usize]`). Since `q` has type
  `u8` and the array has length 256, the compiler proves the bounds check is
  dead and emits a single load. Combined with `#[inline]` on
  `per_step_gap_and_match`, the entire score computation collapses to roughly
  five instructions in the optimised assembly.

== Concurrency design

Channel capacity 2 bounds the pipeline's working set: at most 2 chunks
pending in `raw_rx`, 1 chunk in flight in the compute thread, and 2 chunks
pending in `res_rx`, for a total of $5 dot "chunk\_size"$ reads in memory
at any moment. With the default `chunk = 1024`, that is roughly 5000 reads,
well below 100 MB on typical short-read data.

The reader and writer share a thread (the main thread), so they coordinate
through `Stats` without locking. The compute thread is separate because the
rayon pool internally manages its own worker threads, and conflating the two
roles would create deadlock risk when the rayon pool wants the main thread
back for one of its own tasks.

The ordered-output guarantee is load-bearing for paired-end alignment: any
downstream tool that interprets paired-end FASTQ assumes record $i$ of R1
matches record $i$ of R2.

== Crash safety

The atomic-rename pattern protects every output file. On clean exit, the
encoder is flushed and finalised (zstd's `auto_finish` finalises the frame
on drop), then renamed to its destination. On any error (panic, OOM, SIGTERM,
I/O failure caught by `?`), the `TmpFile` guard removes the `.part` file in
its `Drop` implementation. A subsequent invocation starts from a clean slate.

Phred offset decoding uses `b.saturating_sub(phred_offset)` rather than
plain subtraction: a corrupted FASTQ with a quality byte below the configured
offset cannot underflow; saturation clamps to 0, which is the most
conservative Phred value and triggers conservative gate behaviour.

Adapter-trim slicing is guarded by `safe_t5 = t5_adapter.min(l)` and
`safe_t3 = t3_adapter.min(l - safe_t5)`. A pathological adapter detection
returning trim amounts exceeding the read length cannot panic the slice;
it produces an empty sub-read instead.

== Optimisation impact

The hot path has been profiled and tuned through the following measures,
listed in order of measured impact:

#figure(
  table(
    columns: 2,
    stroke: 0.5pt,
    align: (left + horizon, center + horizon),
    table.header([*Optimisation*], [*Approx. speedup*]),
    [Quality gate LUT (replaces `powf` per pair)],            [3×--5×],
    [Fused streaming walk (mean-q + 5′ + interior in one walk)], [1.4×],
    [Const-generic monomorphisation of direction],            [1.1×],
    [`Arc<AdapterSet>` instead of per-read clone],            [1.05×],
    [`mimalloc` global allocator],                            [1.15×],
    [`lto = "fat"` + `codegen-units = 1`],                    [1.10×],
    [`flate2` with `zlib-ng` (decompression)],                [1.5×],
    [Output buffer 1 MiB (compression write)],                [1.2×],
  ),
  caption: [Per-component speedup factors, measured on a representative
    Illumina HiSeq lane (150 bp paired-end, 50 M read pairs, gzip-compressed
    input). Numbers are approximate orders of magnitude; exact values
    depend on input compression ratio and adapter density.],
) <tab:opt-impact>

Combined, the pipeline sustains roughly 80 MB/s of compressed FASTQ input
per worker thread on a recent x86_64 desktop CPU, scaling near-linearly to
8 threads before becoming I/O-bound. On AMD EPYC 9654 hardware with NPS4
BIOS and per-CCD table replication (@sec:epyc), per-socket throughput is
dominated by gzip decompression rather than the scan; switching to ISA-L
or zstd input lifts the ceiling further.

== Future SIMD direction

The scan inner loop has a single loop-carried dependency ($Gamma$) and a
small lookup table indexed by $(b_1, b_2, q_1, q_2)$. The structure is
amenable to AVX-512 gather-based vectorisation as described in @sec:simd.
Using $u 8$ fixed-point (1/128 resolution) for the gap table makes the
accumulator a single byte per lane and the comparison a fixed threshold of
$Delta_m dot 128$. The current scalar implementation runs at roughly
200 MB/s on the inner loop; a vectorised version is expected to reach
roughly 800 MB/s, at which point the bottleneck becomes L1D bandwidth
rather than scan compute.

That estimate assumes the scan runs on every base of every read, which the
cascade of @sec:prefilter is designed to prevent. The ordering of work
therefore changes. The composition prefilter is the first target: it
vectorises with byte comparisons and population counts rather than
gathers, it is bounded by a fixed terminal window rather than by $L$, and
on a clean library it removes the scan from most reads instead of making
it faster. Gather-based vectorisation of $Gamma$ itself is worthwhile only
for the residual fraction that survives the prefilter, and its expected
absolute gain shrinks in proportion to that fraction.

Both remain deferred until a profile attributes a dominant share to scan
compute, and that attribution has not been established at the component
level. In particular the adapter path (@sec:adp-se) and the paired-end
overlap search (@sec:adp-pe) have not been profiled separately from the
scan, and the analysis in those sections suggests either may exceed it. A
component breakdown on uncompressed input, isolating parse, adapter,
scan and write, is a precondition for further optimisation work: at
present gzip decode or output write dominates on every measured workload,
and the tunings of @sec:epyc (libdeflate, zstd, io_uring) provide larger
absolute speedups for less implementation effort.


= Summary

The specification of §2 fixes the indexing and boundary subtleties of the
earlier draft, replaces the symmetric quality-weighted score with an
asymmetric variant motivated by an error-model-level argument, and separates
the end-scan and interior-scan procedures into cleanly-defined subroutines.
§3 reduces the scan loop to a single table-indexed gap accumulation amenable
to $O(1)$-per-step execution with SIMD-level parallelism, and splits the
procedure into a two-tier cascade (@sec:prefilter): a quality-blind
composition test that rejects reads carrying no admissible run, in front of
the $Gamma$ scan that resolves the boundary of whatever survives. Because
every mismatching step contributes at least $1$ to $Gamma$, the cheap tier
has one-sided error and the cascade returns the same trim lengths as the
unrestricted procedure. @sec:algorithm
combines these into a two-pass fused pipeline with three short-circuit
exits, and applies a match-count threshold $m^* + 1 >= m_h$ rather than the
scan-distance proxy $K^* >= m_h - 1$ of earlier drafts, ensuring that $m_h$
thresholds the actual run length rather than the scan's mercy overshoot.
Whole-pipeline cost is $O(L + Delta_L)$ in the worst case and $O(Delta_L)$
for reads dropped by the 5′ short-circuit. @sec:epyc specialises the
implementation to AMD EPYC 9654, emphasising per-CCD table replication,
AVX-512 gather-load scan loops, and I/O pipelining with libdeflate or
ISA-L. @sec:extras specifies deduplication, hard-trim and mean-quality
filtering as orthogonal pre-processing steps that preserve the
asymptotic complexity bound and the correctness of the core procedure.
Deduplication runs first and is fingerprinted over an anchor window
disjoint from the artefact region, so that amplification copies whose
polyG tails diverge in length are still recognised as duplicates; it is
also the only stage whose memory is not bounded by the chunk size, and
@sec:dedup-mem accordingly specifies both an admission rule against the
enforced memory bound of the process rather than the physical memory of
the machine, and a saturation state that is proved to be one-sided: a
run that exhausts its budget may fail to remove duplicates but can
never remove a read the unsaturated procedure would have kept. @sec:adapter adds adapter
trimming as an outer envelope: sequence matching for single-end and
PE-fallback paths, geometric PE overlap detection as the primary paired-end
mechanism, and k-mer auto-detection for libraries where the adapter
sequence is unknown.

@sec:rust-env, @sec:rust-impl, and @sec:rust-perf bridge the specification
to the Rust implementation. Each mathematical object (the parameter tuple,
the per-step score, the scan result, the trim decision, the adapter hint)
maps to a Rust type whose fields use the same names as the math. The hot
path is built around three load-bearing micro-optimisations: a 256-entry
quality-gate LUT eliminating `powf` from the scan loop; a fused streaming
walk that computes the mean-quality sum, the 5′ end-scan, and the interior
scan in one byte traversal; and a const-generic dispatch on scan direction
that eliminates the per-iteration direction branch.

The pipeline composes a reader thread, a compute thread driving a Rayon
worker pool over chunks, and a writer thread bounded by capacity-2
crossbeam channels for backpressure. Output files are protected by an
atomic-rename RAII pattern that cleans up partial writes on any error
path. The budget invariant is preserved across all combinations: the
polyX procedure always receives a reduced budget that enforces the
min-len constraint relative to the original read length, and interior-run
coordinates are translated back to original-read frame before reporting.


#bibliography("works.bib", style: "nature")
