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
release, together with each tool's maturity (`stable` / `beta`).

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

## D-4 — Standard streams and FASTA ⏳ OPEN (explaining)

The question: should a tool be usable **in the middle of a shell pipe**, and
should it accept FASTA?

Today:

* `fastcheck reads.fastq.gz` works. `zcat reads.fastq.gz | fastcheck -` does
  **not** — `-` is not accepted, and the reader infers compression from the
  file extension, which a stream does not have.
* `fastcheck ref.fasta` exits with an error (no quality string to check).
  fastprep and fastdedup already accept FASTA; fastcheck and fastmaps are
  FASTQ-only by nature.

Why it could matter: `fastp` is routinely used inside pipes
(`zcat in.fq.gz | fastp ... | ...`), and pipe-based pipelines (Snakemake,
Nextflow) exist partly to avoid writing intermediates to disk — which is
exactly the I/O argument this suite is built on. Composing
`fastdedup | fastprep | fastmaps` without intermediates needs **stdin on the
consumer and stdout on the producer**.

Options: (a) stdin input only, outputs stay files; (b) stdin **and** stdout, so
a full chain runs with zero intermediates; (c) neither, and say so in the docs.
Cost: (a) is small (read sequentially, detect gzip from magic bytes instead of
extension); (b) adds an explicit output-format flag, since a stream has no
filename to infer from.

---

## D-5 — fastcheck ↔ fasthelper boundary ✅ decided

* **fastcheck = the FastQC/Falco slot**: one pass over FASTQ files, per-file
  QC report, per-position matrices. It never reads a run directory and never
  renders HTML.
* **fasthelper = the MultiQC slot**: aggregates every tool's reports across
  samples, computes cross-sample statistics, renders the HTML/TSV/MD report,
  runs the chain from a config.

---

## D-6 — Behaviour when a detection returns nothing ⏳ OPEN (explaining)

The question is not about heuristics; it is about what a tool does when it
cannot determine something from the data:

* **fastprep, adapter detection** — it fingerprints the first ~2 000 reads to
  identify the adapter. On genuinely adapter-free data it finds nothing and
  today **exits 1** with instructions. In a 700-sample array job that turns a
  *good* sample into a failed task.
* **fastcheck, Phred offset** — it can prove `+33` when a quality byte below 64
  exists; otherwise it assumes `+33`, records `phred_offset_certain: 0` and
  keeps going.
* **fastsim, adapter sequence** — a hardcoded Nextera constant in the profile.

So the suite is currently inconsistent, and both extremes are defensible.
The question is which one becomes the rule:

* Option A — **refuse** (fastprep's current behaviour): never produce a result
  that rests on an unproven assumption.
* Option B — **take the safe action, label it, continue**: no adapter
  detected ⇒ trim no adapter (the safe action), record the assumption in the
  report, exit 0. `--strict` restores Option A for anyone who wants it.
* Option C — guess and warn (what FastQC/fastp effectively do).

Recommendation now: **B**, because with `"assumptions": [...]` in the JSON and
in the summary line, nothing is hidden — the report says "adapter: none
detected, trimming skipped" — and batch runs do not fail on good samples.
Option A stays available per run via `--strict`.

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

Files currently in the tree that fall under this rule and need a ruling:

| file | what it is | concern |
|---|---|---|
| `tools/fastsim/profiles/mandril_E.json` | statistical profile measured from the mandrill E dataset (CDFs, quality curves, GC) | derived from Mandrillus Project data, ~96 KB, currently committed |
| `tools/fastsim/validation/*.png` | closed-loop validation plots from that dataset | same |
| `tests/*` inline FASTQ literals | hand-written 6-read strings ("ACGTACGT…") | not data — but confirm the policy is "no real data", not "no FASTQ bytes at all" |

---

## D-12 — What the suite refuses to do ⏳ open

Planned list, one line each in the README with the reason: PE overlap/merge
(removed from the polyx lineage deliberately), UMI handling, alignment-based QC,
workflow-engine semantics (arrays, resume — that is Snakemake/Nextflow's job),
per-tile modules (implemented, off by default), long reads (D-7).
