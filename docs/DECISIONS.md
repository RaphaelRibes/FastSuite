# FastSuite decision log

One entry per decision that shapes more than one tool. Status: **decided** /
**open**. Keep this short — it exists so decisions are not re-litigated.

---

## D-1 — FastSuite is a naming scheme, not an architecture. ✅ decided

FastSuite is the *name* for a group of independent preprocessing tools. What
they share is naming, the common option set and the `fastsuite.report.v1`
payload — nothing else.

* No shared crate, no workspace, no `fastsuite-core`. Tools are separate
  binaries that happen to look alike and talk the same report format.
* Anything two tools need (BGZF reader, FASTQ scanner, fingerprint table,
  report writer) exists **twice**, copied per tool, and may drift. Fix per tool
  and say which in the changelog.
* Consistency is enforced by **tests**, not by the compiler:
  `fasthelper check` validates each tool's `--json` against
  `REPORT-SCHEMA.md`; any two tools computing the same quantity (fastdedup vs
  fastcheck duplication) get a fixture-agreement test on the same input file.

---

## D-2 — Versioning ✅ decided (forced by bioconda)

**One semver per tool, plus a suite release train.** bioconda recipes are
per-package, so each tool carries its own version and its own `meta.yaml`.
A `docs/RELEASES.md` lists which tool versions made up each "FastSuite N"
release, together with each tool's maturity (`WIP`).

---

## D-3 — Distribution ✅ decided

Everything ships *in a pixi environment*, published to **bioconda**, and
containerised with **pixitainer**.

Consequences to implement:

* every tool gets a `pixi.toml` (fastcheck is missing one) and the suite root
  `pixi.toml` is the dev environment;
* one bioconda recipe per tool (`meta.yaml`, license file, `run`/`test`
  sections invoking the binary's `--version` and a tiny fixture run);
* one pixitainer recipe per tool for the container images;
* version pinning matters from now on (D-2), and the recipes must not download
  anything at build time.

---

## D-4 — Standard streams and FASTA ✅ decided: both

Tools work **in a shell pipe and on files**: `-` is accepted for FASTQ input
and for FASTQ output, alongside ordinary paths. A chain like
`zcat in.fq.gz | fastdedup -1 - -o - | fastprep -1 - -o clean.fq.gz` must run
with **zero intermediate files** on disk.

Consequences:

* gzip/BGZF detection must come from **magic bytes**, not from the file
  extension (a stream has no name).
* A stream cannot be seeked, so the reader peeks the header bytes and chains
  them back in front of the stream.
* Settled 2026-10-12 (fastdedup, the first tool with output streaming): output
  to a stream carries an explicit decision — `--out-format fastq|fasta` and
  `--out-compression none|gz` (defaults: fastq, none). `-o -` does **not**
  inherit the input's compression: a pipe whose bytes silently depend on the
  producer's format is a trap for the next command in the chain. Named paths
  keep inferring both from their extension.
* FASTA stays tool-dependent and explicit: fastprep and fastdedup accept it;
  fastcheck and fastmaps need qualities and refuse it with a clear message.
* **`fastmaps` is file-only.** It exists to take *many* files and classify
  R1/R2 by filename; a stream removes that information, so a pipe there would
  cost capability instead of saving I/O. The other four tools take `-`.

Reference implementation: `fastcheck`. fastprep and fastdedup follow.

---

## D-5 — fastcheck ↔ fasthelper boundary ✅ decided

* **fastcheck = the FastQC/Falco slot**: one pass over FASTQ files, per-file
  QC report, per-position matrices. It never reads a run directory and never
  renders HTML.
* **fasthelper = the MultiQC slot**: aggregates every tool's reports across
  samples, computes cross-sample statistics, renders the HTML/TSV/MD report,
  runs the chain from a config.

*Correction (2026-10-12, verified 2026-10-18):* fastcheck **does** render its own
**per-file** HTML report (`--html`, MultiQC-style, inline SVG, no JS). The
boundary therefore reads: fastcheck renders the **per-file** report, fasthelper
renders the **cross-sample** report. fastcheck still never reads a run directory,
which is the half of the original statement that held.

---

## D-6 — Behaviour when a detection returns nothing ✅ decided: B

**Take the safe action, label it in the report, exit 0. `--strict` restores
refusal.**

Concretely: fastprep's adapter auto-detection finding nothing now means "trim
no adapter" (the safe action), with the assumption printed on stderr and
recorded in the report, instead of exit 1. Same for any future detection: the
safe default plus an `assumptions` entry
(`{"what":"adapter","value":"none detected","consequence":"no adapter trimming"}`).

`--strict` makes any assumption fatal, for anyone who prefers the old refusal.
Nothing is ever guessed silently: every assumption appears in the JSON and in
the summary line.

Applied so far: fastprep (`--strict` added, safe continuation by default).
fastcheck already behaved this way for the Phred offset
(`phred_offset_certain: 0`), and now also emits an `assumptions` array.

---

## D-7 — Short reads only, for now. ✅ decided

`--max-read-len` cap with exclusion reporting stays. Binned positions and
long-read modes are out of scope until the short-read suite is finished and
published; `--help` says "short-read QC".

---

## D-8 — Compatibility horizon ⏳ open

Old **flags** forever (aliases cost nothing). Old **binary names**
(`rpolyx`, `fdedup`, `novaseq-sim`) as symlinks installed by
`fasthelper install`, not by cargo/bioconda. Legacy **JSON** shapes accepted by
fasthelper for one release, then dropped with a changelog entry.

---

## D-9 — Licensing and attribution ✅ decided in scope, ⏳ license text open

* **Only FastPreprocessing and FastDedup belong to the Mandrillus Project** —
  they carry the project, ISDM/ISEM and Céline Mandier / Alice Baniel credits.
* **Everything else (fastsim, fastmaps, fastcheck, fasthelper) is the author's
  own work**, credited to Raphaël Ribes alone.
* Still to fix: the license text per tool (MIT proposed, fastdedup already is)
  and a `CITATION.cff` per tool so each is independently citable.

---

## D-10 — Repository layout ⏳ open

Monorepo now. The suite is a naming scheme (D-1), so the repo does not need to
mirror it: keep one repo while D-3 (bioconda recipes + pixitainer) is built
out, then split per tool if release automation or citation demands it.

---

## D-11 — NEVER SHIP ANY DATA. ✅ decided (hard rule)

No sequencing data, no subsets, no derived tables of real reads, in any repo,
container, tarball, test fixture or example. Everything demonstrated publicly
uses **generated** data (fastsim) or hand-written literal fixtures.

Files that were in the tree and are now **out** of it:

| file | what it is | disposition |
|---|---|---|
| `tools/fastsim/profiles/mandril_E.json` | statistical profile measured from the mandrill E-dataset (CDFs, quality curves, GC) | moved to `~/hermes_scratch/fastsuite-local-artifacts/profiles/`; also still in the untouched `~/hermes_scratch/novaseq-sim/` tree. Only the **synthetic** `test_small.json` ships |
| `tools/fastsim/validation/closed_loop_validation.png` | validation plot from that dataset | moved to `~/hermes_scratch/fastsuite-local-artifacts/validation/` |

Enforcement: `.gitignore` blocks `tools/fastsim/profiles/*E.json` and
`tools/fastsim/validation/*.png`, so a profile cannot be re-committed by
accident. Real profiles are built on the compute node
(`scripts/build_profile.py STATS_DIR out.json`) and stay there.

Not data, so they stay: hand-written test fixtures (inline literals such as
`"@r1\nACGTACGT\n+\nIIIIIIII\n"`) and `profiles/test_small.json`, which is
synthetic (gc 0.5, no artefacts, `full_length_rate: 0`).


---

## D-12 — What the suite refuses to do ⏳ open

Planned list, one line each in the README with the reason: PE overlap/merge
(removed from the polyx lineage deliberately — the paper's §`adp-pe` now says so),
UMI handling, alignment-based QC, workflow-engine semantics (arrays, resume — that
is Snakemake/Nextflow's job), per-tile modules (planned, plan Task 11), long reads
(D-7).

---

## D-13 — Deduplication ships ON ✅ decided (2026-10-18)

`fastprep --dedup` is on by default. This is a product decision, not a
measurement outcome: exact-fingerprint dedup is one of the things that separates
the suite from fastp, which does not dedup at all.

The consequence to keep in view, measured over 102 datasets and 20.4 M mates: the
trim-accuracy harness scores only the reads a tool **kept**, so a deduplicating
tool is charged for every duplicate it removes. On the worst-retention datasets
**100% of the drops are `Duplicate`** (SRR20754397: 167,162 of 167,162;
DRR015139: 32,138 of 32,138), and `--min-len` does not move that histogram at all
— the earlier "fastprep drops more reads than fastp" reading blamed the length
floor for the dedup default. With `--no-dedup` the aggregate retention is
**98.47% against fastp's 97.23%** (virome 99.74% vs 99.7%, rna_seq 98.90% vs
97.8%): the suite does not keep fewer reads, it deduplicates.

Rules that follow:

* Quote trim accuracy **with the shipped default** and, whenever fastp is in the
  table, **also with `--no-dedup`** — that is the only apples-to-apples pair.
* Read **retention and exact rate together**, never the rate alone. Dropping a
  read removes it from the scored set, and whether that flatters or depresses the
  rate depends on the class dropped: `--min-len 0.9L` shows the *highest* fastprep
  exact rate (89.470% vs 89.088%) with **441k fewer** exact mates, while
  `--no-dedup` raises both.
* Evidence: `bench/results/phase-d/minlen-ladder-20261018.tsv` (+
  `-perdomain` and `-dropreason`), `docs/TRIM-ACCURACY.md` §Read retention.

**Corollary (decision taken with it): `--min-len` keeps its shipped `0.9L`
default.** The ladder showed no rung buying anything (15/0.5L: +0.15 pp
retention, +0.07 pp exact, i.e. noise; 30 neutral; 50 and uniform 0.9L cost
2.7–2.9 pp of both), so the retention question is answered by D-13, not by the
length floor.

---

## D-14 — `--pe-overlap` needs an adapter set to be sound ✅ decided (2026-10-19)

The paired-end overlap mode pins the adapter boundary from the mate pair. Its
safety guard — the boundary applies only when adapter *starts exactly at* it —
needs an adapter set to check against. **With no set, the guard is inert and the
mode trims biology.**

Measured failure: with adapter auto-detection returning nothing (the default
configuration, and a common one), `Params.adapter` is `None`, the guard is
skipped, and the mode trims a short-insert mate to its insert even though nothing
was ligated past the boundary. 20,000 pairs, insert 100 of 150, random tails →
reads trimmed to 100 bp with `--pe-overlap` and auto-detection.

The rule, and why it is about the *source* rather than the flag:

| adapter source | mode runs? | why |
|---|---|---|
| `--adapter-preset` / `--adapters` (explicit) | **yes** (default) | the set is known, so the guard can actually run — the mode is then *safer* than the ordinary path, since it refuses a boundary the matcher would not confirm |
| auto-detected (the default) | **no** | auto-detect may return an empty set, and with no set the guard is inert |
| `--no-adapter` | **no** | nothing to find; the mode would be inventing a boundary |

So the flag is a **tri-state** — `--pe-overlap=never|auto|always`, `auto` being
the default and resolving to "on iff the adapter source is explicit". `always` is
the explicit opt-in for someone who knows the guard is inert and wants the
boundary anyway.

**Why the default tracks the source rather than being a constant:** the mode's
soundness is a property of whether a set exists to check against, not of the
user's intent. A constant default either exposes the bug (always on) or hides a
useful instrument (always off).

**The measurement behind keeping it at all:** on a synthetic cohort built to be
scorable (800,000 mates, 14.1% readthrough, 489,546 carrying a G-run), exact
matches are **identical** (347,171) with and without the mode, and overtrim rises
by 106 — **+0.000 pp overall** against a pre-registered 0.5 pp threshold
(`docs/OVERTRIM-PAIRED.md`, Cell B). The mode is kept as an opt-in instrument
because its guard is real and because it is the only way to re-run the question
on a different corpus. Throughput cost measured: **+3.7% single-threaded, +7.8%
at 4 threads** (best-of, 400k pairs).

**Amendment (2026-10-21, exactness work Tier 3c):** the new opt-in
`--pe-ambiguous-smallest` narrows the *refusal*, it does not remove the source
rule above. With more than one verifying insert it takes the smallest one
(shortest insert = largest readthrough trim, i.e. the overtrim side) **only
when `readthrough_is_adapter` confirms adapter behind that boundary on the
pinned mate**. With no explicit adapter set the flag still cannot fire, so the
tri-state safety property is unchanged.

---

## D-15 — Trim cost model: undertrim priced at 7× overtrim ✅ decided (2026-10-21)

The exactness work needs one number the scorer cannot invent: the exchange rate
between the two ways to miss. **`w_over = 1`, `w_under = 7`** — artefact left in
the output (wrong data) costs seven times real sequence deleted (less data).
This is the standing preference ("we prefer less data than wrong data") made
numeric, fixed *before* the sweep so the margin cannot be tuned to its result.

* Weights live in `bench/weights/trim-weights-v1.json`; every result row
  records `weights_version = v1`. The ladder `r ∈ {1, 3, 5, 7, 10}` is reported
  alongside, plus a per-rule break-even `r* = −Δover/Δunder` (in reads), so any
  verdict can be re-read at another ratio from the same counts — the per-axis
  raw numbers never depend on the weights.
* Scoring is per-read (primary, comparable to the existing trim-accuracy
  numbers); a per-base variant is a secondary column only. A deliberate 2 bp
  margin and a whole-read wrong drop both cost `w_over` in the primary score.
* The formula is homogeneous: only the ratio is a decision, and 7 is the
  ship/no-ship anchor.
* **Provenance of the trade quoted below** (verified 2026-10-22):
  `bench/results/phase-d/fastprep-flags-analysis-20261018.tsv`, commit
  `c46cbd3` (2026-10-18) — measured *after* the sim regeneration (2026-10-15)
  and after the discrepancy fix `73c2385`, over 102 datasets / 20.4 M mates,
  baseline `fastprep-hd0`. It is generation-consistent. The retired
  66.5 / 29.5 / 4.0 round is **not** generation-consistent and is not an
  acceptance value (`docs/TRIM-ACCURACY-DISCREPANCY.md`).

**What 7 does to the one closed band:** the measured 3′ terminal-run trade
(`--min-run-3p 6`: −35,196 undertrim bought at +236,779 overtrim) breaks even
at `r* = 236779/35196 ≈ 6.73`, so at r = 7 it is net-positive by **+9,593
weighted reads**. The band's closure verdict was computed on the symmetric
exact rate and does not carry over: it is re-scored (plan Step 5) and the
result is a decision for the user, not a rule change — the generic polyX form
stays a non-goal until he reopens it (the A/C/T census is a correctness
argument, the weight is a cost argument; both must agree before it ships).

**Measured at the anchor (2026-10-22, desktop, 1 224 rows):** the campaign's
per-dataset winner configuration (`--head-dropout 4`) is the optimum of that knob
at 7:1 — hd0 reaches 82.80% strict
exactness but buys it with 9.58% undertrim and scores **−0.4947** against base,
hd1 −0.1338, hd2 −0.0255. `--boundary-anchored` is **+0.0121** and
`--overtrim-margin` is **0.0000** on this corpus (its target class is 0.003% of
mates). Full table: plan §7.

**Baseline correction (2026-10-22):** the tolerance-0 yardstick for r = 7 is
**75.59% exact / 23.4% overtrim / 1.0% undertrim** (`overtrim-vs-fastp.tsv`;
75.35% after the `_minlen` repair), against fastp-default's 92.19%. The
proposal's original 66.5 / 29.5 / 4.0 figures came from a defective pairing
(generation-1 trims scored against generation-2 ground truth) and must not be
quoted as a baseline.

**Deliberate margins are attributed from the config, never from PXMT:** the
PXMT v3 record stays byte-stable (25+ bench parsers unpack it at fixed
offsets), and the ≤2 bp class is split jitter-vs-margin from the margin-flag
on/off delta per dataset — tolerance decay alone cannot tell them apart,
because both decay at ±2.

## D-16: `--head-dropout` ships at 4, chosen by the 7:1 price (2026-10-22)

**Decision.** The default is `4` again, from `0`. The knob stays exposed; `0`
disables the rule and `1`-`3` shorten it.

**Why the default moved twice.** The rule trims a 1..=4 bp leading G-run
exactly. It is a prior, not a detector, so it is wrong on roughly three of every
four reads it touches. On 2026-10-16 (`3d76a31`) the default went to `0` on a
*total-error* argument: at tolerance 0 of kept mates the rule cost 24.6 % error
against 7.4 % without it, so switching it off looked like a 3.3x improvement.
Re-measured on the current corpus with the shipped binary it is 22.2 % against
10.5 %: a 2.1x improvement, not 3.3x, and still an improvement at 1:1. That argument holds only if a deleted base and an artefact left in
the output are priced equally.

**The price settles it, and the break-even is measurable.** With undertrim
costed at *r* times overtrim, the rule pays whenever `r > P(biology)/P(artefact)`
for the affected class. Measured three ways, that ratio is 2.6:

| route | value |
|---|---|
| quality posterior for a 1..4 bp leading run (0.27-0.29 artefact) | 2.5-2.7 |
| the 102-library four-tool table, Δover / Δunder (19.97 / 8.31) | 2.40 |
| the 101-library head-dropout axis, same ratio | 3.70 |
| the 7:1 campaign (1 224 rows, desktop), same ratio | 2.39 |

At 1:1 (below 2.6) `0` is right; at 7:1 (above it) `4` is right, and the
measurements agree on the margin: `+43.6` of mates with the rule against `+17.1`
without on the four-tool table, i.e. **+26.5 points** for `4` (`+39.9` against
`+15.2` on the twelve-arm campaign, a +24.7-point margin). The standing project preference is
"prefer less data than wrong data", fixed at 7:1 in D-15 — so the shipped
default must be the 7:1 optimum. Shipping `0` made a bare invocation run the
configuration this project's own cost model rejects.

**What does not change.** The measurement itself was already in D-15 and the
plan's §7; only the default followed it. The three Tier-2/3 flags stay off
(`0.0000` on this corpus). `--boundary-anchored` stays on. The 1:1 reading stays
reachable (`--head-dropout 0`) and is documented as the choice for a pipeline
that prices a lost base and a kept artefact equally.

**Evidence.** `bench/results/exactness/` (`sweep-raw.jsonl`, `sweep.tol0.tsv`,
`sweep.tol2.tsv`, `side.tsv`, `qgate.tsv`), `bench/results/phase-d/head-dropout-axis.tsv`,
`bench/results/phase-d/overtrim-vs-fastp.tsv`, and the rendered tables in
`tools/fastprep/CHANGELOG.md` (0.4.1 entry) and `docs/plans/2026-10-21_fastprep-exactness.md` §7.

## D-17: the 5' dropout is priced, and the price is a library property (2026-10-22)

`--head-dropout` was measured as two global policies (trim every 1-4 bp leading
G-run, or none). Neither is the end of the question: the injected share is not
the same in every library, and on a 21-library probe it runs from 0% to 63%.
A per-library decision is therefore available, and it is estimable from the reads
alone -- the excess of the short leading-G rate over the library's own G
frequency (`p - p^5`), divided by the ambiguous rate, thresholded at the cost
model's own break-even `w_over/(w_over+w_under) = 1/8`.

Measured end to end through the real binary over 102 libraries and 20.4 M mates
(`bench/scripts/auto_hd_sweep.py`, raw in `bench/results/exactness/auto-hd-*`):

| policy | exact | overtrim | undertrim | score@7 | total error (1:1) |
|---|---:|---:|---:|---:|---:|
| hd 4, anchored (shipped) | 72.35% | 18.71% | 1.88% | +0.6749 | 20.59% |
| per-library, class threshold | 83.06% | 7.50% | 2.38% | +0.5218 | 9.88% |
| per-library, all-mates threshold | 84.38% | 5.32% | 3.24% | +0.4961 | 8.56% |
| hd 0, anchored | 83.19% | 0.15% | 9.60% | +0.1820 | 9.75% |

Decision: the default does not change. The estimate is worth +10.7 pp of exactness
and dominates `hd 0` outright (same exactness, a quarter of the undertrim), but it
still misses libraries -- an undertrim mistake costs 7x, so 0.5 pp of extra
undertrim outweighs 3.5 pp of overtrim -- and at 7:1 the global policy wins
(+0.6749 against +0.5218). It is not shipped: it needs a prescan and a per-library
decision the binary does not take, and there is nothing for it to win at the
shipped price. `--head-dropout 0` is no longer the recommended way to buy quality;
the estimate is, as soon as someone writes the flag.

Also negative, and recorded because it closes the search: no per-read feature
separates the classes (1 171 418 ambiguous mates, six features, dataset-level
train/test split; the best rule returns 6 594 mate-units on the held-out half, 0.8%
of what that class is worth, and every other threshold is negative), and lowering
`--min-run` to reach the recoverable 3' class loses at the anchor (mr8 72.17%
exact / +0.6720; mr6 70.65% / +0.6228).

## D-18: the floor is a trim budget, and the mercy clause bounds overshoot only (2026-10-23)

Two coupled defaults changed after the 2026-10-23 fastprep review. Both are
recorded here because they move every number in `docs/TRIM-ACCURACY.md`.

**`--min-len 0.67L` (was `0.9L`).** The floor and the trim budget are one
number: a read that needs more trim than the floor allows is DROPPED rather
than trimmed, so `0.9L` permitted a read to lose only 10 % of its length.
Measured on a 150 bp read with a 30 bp 5′ polyG artefact: dropped at `0.9L`,
trimmed to 117 bp at `0.67L`; in paired mode the same library gave **0 pairs
out of 2000** at `0.9L` and 500/500 at `0.67L`. `0.67L` keeps 100 bp of a
150 bp read — a third of the read trimmable — which is the pipeline owner's
stated requirement. It is NOT a sweep result: the 90-dataset campaign and
Phase D both ran at `0.9L`, and their exact-rate metric excludes dropped
reads, so a floor that discards the hard reads was never charged for it. The
resolved floor is printed in bases at the end of every run ("100 bp of a
150 bp mean read; 50 bp of trim allowed"). `--min-len 0.9L` reproduces the old
behaviour; every table measured before this entry must state that floor.

**The mercy clause charges a run's own bases nothing.** Δ_m is meant to bound
how far the cut walks *past* a run's end — the overshoot window, where a chance
match in the flanking biology can extend the run span. Charging matches turned
it into a length cap, because a match costs `μ(1−W)`: at Q2 that is 0.09, so
roughly 33 bases exhausted Δ_m = 3. Measured cap on a single read: 50, 70 and
100 bp Q2 tails all cut at exactly 37 bases; Q5 at 43, Q10 at 64, Q20 uncapped.
For 4-colour chemistry a low-quality polyG tail is the normal artefact, so the
model trimmed *least* on the case it exists for. After the change a match costs
nothing until the first mismatch; from there the walk is charged exactly as
before, so the overshoot window is unchanged (the ~2 flanking bases). Measured
after: tail + 2 at Q0, Q2, Q5, Q10, Q20, Q37 for tails of 10–100 bp.

What this does NOT settle: the corpus cannot measure either change. Its longest
injected 3′ tail is 20 bp (below the old cap), its 5′ artefact is a single G
(`HEAD2 = {0: 8154, 1: 1846}`), and `bench/results/exactness/undertrim-decomp.tsv`
attributes 0.95 % of undertrim to 3′ runs ≥ 6. A bench arm with long low-quality
tails is required before any number is quoted about either.

### Measured on the desktop, 2026-10-23 (15-dataset subset, `--stride 7`)

`overtrim_vs_fastp.py --tools fastprep` on the GT corpus, each dataset at its own
swept winner and with the swept `--min-len` PINNED, so the numbers isolate the
scan change from the floor change. Binary md5 `7b27b2af204181ce10982deb5fb7bb88`
(desktop release build of `6763a4f`), rows in `/tmp/trim-regress-20261023.tsv`,
compared against the committed `bench/results/phase-d/overtrim-vs-fastp.tsv` for
the same 15 datasets.

| arm | exact (was → now) | overtrim | dropped |
|---|---|---|---|
| fastprep-hd4 (shipped) | 76.92 % → **79.19 %** (+2.27 pp) | 21.70 % → **19.39 %** | 9.08 % → 8.86 % |
| fastprep-hd0 | 89.86 % → 89.81 % (−0.05 pp) | 0.02 % → 0.03 % | 8.98 % → 8.86 % |

So the change is a large win on the shipped arm and a wash on `hd0` (two datasets
move 0.2–0.3 pp down, `DRR015180` and `DRR015431`). Every dataset's `scored`
denominator is within 0.1 % of the old one, so this is not a denominator effect.
The floor change (0.67L) is NOT in these numbers — the harness pins the swept
`--min-len` — and still needs its own arm.

## D-19: `--auto-params` ships a decision list, not a ranker (2026-10-23)

**Decision.** The recommendation is a depth-3 decision list over fastcheck's
input-side metrics, exported to JSON, walked in-tree, and **embedded in the
binary** with the same bytes committed at `tools/fastprep/models/`.
`--print-args` prints the argv (stdout) and the rule that fired with its
training coverage and mean regret (stderr). Explicit flags beat the model.

**Why a decision list.** The tool's failure mode is silently wrong reads, so a
recommendation has to be auditable and overridable. A boosted ranker would score
better and be unexplainable; the 8 KB JSON can be read in a minute, and
`--print-args` names the evidence. Depth 3 with ≥ 10 datasets per leaf is the
shape that a 101-dataset corpus can support without producing one-dataset rules.

**Measured acceptance** (group 5-fold by dataset, stratified by domain; the
script refuses to write a model that loses): mean regret vs the per-dataset
oracle **0.0331** against the default's **0.0721**; delta vs the default
**+0.0389**; beats the default on **70.9 %** of held-out datasets; top-5
**0.533** against the default's **0.020**. Leave-one-domain-out is positive on
all six domains. **One of five folds loses to the default** (regret 0.1180,
delta −0.0462) — an average improvement, not a per-library guarantee, and the
reason the printed explanation exists.

**The action space is 39 of 52 configs.** The 13 excluded set `r2_reverse`, a
bench-side revcomp of R2 that is not a fastprep flag; recommending one would
recommend an effect the tool cannot produce. 23 datasets' oracle winner needed
one of those, which is a permanent ceiling on the model's regret.

**Train/serve aliasing is part of the contract.** The training matrix's
`read_len` and `bases_non_acgt_frac` are DERIVED from `mean_len` and
`bases_non_acgt / bases`; the derivation lives in `auto_params::Profile`. Two of
the seven rules test `read_len`, so without it they would parse, never fire, and
fall through to the default without saying so. Verified by walking all 102 real
QC reports in Rust and in the trainer's Python reference: **102/102
byte-identical argv**.

**Known-bad regimes** (documented in `docs/MODEL-CARD-auto-params.md`): 101 bp
reads only (no 150 bp library in the corpus), no 2-colour chemistry, no 3′ tail
longer than 20 bp, six study types. The model's `min_len` values come from the
Phase D grid and therefore REPLACE D-18's 0.67L floor when applied; pass
`--min-len` explicitly to keep it.
