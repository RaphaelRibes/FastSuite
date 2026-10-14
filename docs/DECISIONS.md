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
(removed from the polyx lineage deliberately), UMI handling, alignment-based QC,
workflow-engine semantics (arrays, resume — that is Snakemake/Nextflow's job),
per-tile modules (implemented, off by default), long reads (D-7).
