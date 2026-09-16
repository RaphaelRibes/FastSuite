# FastSuite decision log

One entry per decision that shapes more than one tool. Status: **decided** /
**open**. Keep this short — it exists so decisions are not re-litigated.

---

## D-1 — Shared code: none. Tools are independent. ✅ decided

Each tool is its own crate and its own binary. What they share is **naming and
the common option set**, for continuity: `-1/-2`, `-o/-p`, `-j/--threads`,
`--gz-level`, `--json`, `--force`, `--dry-run`, `-q/-v`, `-V/-h`, and the
`fastsuite.report.v1` payload.

Consequences, accepted:

* No `fastsuite-core`. Anything two tools need (BGZF reader, FASTQ scanner,
  fingerprint table, report writer) exists **twice**, copied per tool.
* The report schema and the CLI spec are the only interfaces. They cannot be
  enforced by the compiler, so they are enforced by **tests**:
  * `fasthelper check -c fastsuite.toml` validates that each tool's `--json`
    matches `REPORT-SCHEMA.md` (phase 2 task).
  * Any two tools that compute the same quantity (fastdedup's duplicate rate,
    fastcheck's duplication module) get a **fixture-agreement test**: same
    input file, identical numbers. That test is what replaces a shared library.
* A bug fixed in one copy does **not** propagate. Fix or accept per tool, and
  say which in the changelog.

---

## D-2 — Versioning and releases: one suite version or one per tool? ⏳ open

**Why it is big:** it decides what "FastSuite 1.0" means, how compatibility is
communicated, and whether the replacement claim ("replaces FastQC/fastp") is
made about a suite or about a binary.

* Option A — **one suite version** (all tools released together, same number).
  Simple to explain and to cite; forces six changelogs into lockstep.
* Option B — **independent versions** (fastdedup 1.3.0, fastcheck 0.1.0, …).
  Honest about maturity differences; needs a compatibility note instead of a
  version number.

Recommendation: **A, with a maturity column per tool** (`stable` / `beta`), so
`FastSuite 1.0` can ship with fastcheck still marked beta.

---

## D-3 — Distribution and install ⏳ open

**Why it is big:** it decides whether anyone but you can use it. Your real
deployment target is a headless HPC cluster with no root and jobs submitted
under SLURM — `cargo install` is not available there.

* Options: static release binaries (musl + gnu) on GitHub releases ·
  bioconda recipe · pixi/conda channel (you already prefer pixi) · Apptainer
  definitions (already used for the mandrill pipeline) · crates.io.
* Recommendation: **static binaries + Apptainer defs** first (works on Io
  today, no root), **bioconda second** (that is where FastQC/fastp users live),
  crates.io only if the names stay free.

---

## D-4 — Standard streams and FASTA: is the fastp-replacement claim honest? ⏳ open

**Why it is big:** `fastp` is used in pipes (`zcat x | fastp | ...`) and on
FASTA. fastcheck currently refuses FASTA and nothing supports `-`. Without
this, "we replace fastp" is a claim about flags, not about fit.

* Options: implement `-` in/out for fastprep/fastmaps/fastcheck (blocked on a
  non-seeking reader), or declare pipes out of scope and say so.
* Recommendation: implement `-` for **input** only (the common case, cheap
  with the existing sequential reader), keep output-to-file, and document the
  limitation rather than half-supporting it.

---

## D-5 — fastcheck ↔ fasthelper boundary ⏳ open

**Why it is big:** it decides where aggregation lives, and MultiQC's job must
not end up implemented twice.

* Proposal: **fastcheck analyses files (one row per file in `per_file`),
  fasthelper aggregates samples/runs.** fastcheck never reads a run directory,
  never renders HTML, and has no cross-sample logic.
* Cross-sample questions ("worst retention", "highest polyG") are fasthelper's,
  driven by its own config, not by fastcheck flags.

---

## D-6 — "Never guess" as a suite-wide stance ⏳ open

**Why it is big:** it is a repeated UX decision (fastprep refuses to guess an
adapter and exits 1; fastcheck reports an uncertain Phred offset instead of
asserting one; fastsim hardcodes the adapter sequence). Users migrating from
FastQC/fastp expect a tool that always produces output.

* Option A — strict by default, one `--lenient` flag that guesses **and says
  it guessed in the report** (`"assumed": true`).
* Option B — guess by default, warn.
* Recommendation: **A**. Your existing tools already behave this way and it is
  the strongest differentiator against FastQC/Fastp.

---

## D-7 — Long reads (ONT/PacBio) ⏳ open

**Why it is big:** per-position matrices are O(read length); the heatmap grid
is O(len²). Short-read-only is a defensible scope, but it must be stated.

* Options: hard cap with exclusion reporting (current: `--max-read-len 1000`) ·
  binned positions as the default beyond a threshold · refuse long reads
  explicitly.
* Recommendation: cap + exclusion now, **binned mode as an explicit opt-in**,
  and a one-line "short-read QC" statement in the tool's `--help`.

---

## D-8 — Compatibility horizon ⏳ open

**Why it is big:** it decides how long `rpolyx`/`fdedup`/`novaseq-sim` names,
`--compression`, `--input/--output` and the legacy JSON shapes must keep
working — and therefore how long fasthelper keeps its adapter layer.

* Recommendation: old **flags** forever (aliases are free), old **binary
  names** as symlinks shipped by `fasthelper install` (not by cargo), legacy
  **JSON** accepted for one release then dropped with a changelog entry.

---

## D-9 — Licensing and attribution ⏳ open

**Why it is big:** it gates publication, packaging and citation, and it is
awkward to change later. Who is credited matters here: the tools were built in
the Mandrillus Project context (ISDM/ISEM, Céline Mandier, Alice Baniel).

* Recommendation: MIT for all tools (fastdedup already is), one `AUTHORS`/
  citation block listing every contributor, and a `CITATION.cff` so the suite is
  citable per tool.

---

## D-10 — Repository layout ⏳ open

**Why it is big:** six independent tools in one repo is convenient now and
painful at release time (six tags, six issue queues, one CI). The earlier
rpolyx split (code / bench / paper) already showed both sides of this.

* Options: keep the monorepo · one repo per tool plus a meta-repo for docs and
  releases · monorepo now, split when the first tool reaches 1.0.
* Recommendation: **monorepo until D-3 (distribution) is settled**, then split
  per tool if bioconda/release automation demands it.

---

## D-11 — Canonical validation data ⏳ open

**Why it is big:** it decides whether anyone else can reproduce the claims.
Candidates: the mandrill E dataset (read-only on Io), fastsim profiles +
ground truth, one small public dataset (a single SRA/ENA run) for
external readers.

* Recommendation: fastsim ground-truth tests in CI (fast, self-contained) +
  one published real run as the acceptance fixture, referenced by accession in
  the docs. Never ship the mandrill reads themselves.

---

## D-12 — What the suite refuses to do ⏳ open

**Why it is big:** it is the defence against feature creep, and it makes the
scope honest in writing. Currently out: PE overlap/merge (deliberately removed
from the polyx lineage), UMI handling, alignment-dependent QC, workflow-engine
semantics (SLURM arrays/resume), per-tile modules (implemented but off).

* Recommendation: keep this list in the README, one line each with the reason.
