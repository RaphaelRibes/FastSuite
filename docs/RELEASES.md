# Releases

Per **D-2** (versioning, forced by bioconda): one semver per tool, plus this
suite release train. bioconda recipes are per-package, so each tool carries
its own version and its own `meta.yaml`; this file records which tool
versions made up each "FastSuite N" release, with each tool's maturity
(`WIP`).

Attribution (**D-9**): **fastprep and fastdedup belong to the Mandrillus
Project** — they carry the project, ISDM/ISEM and Céline Mandier /
Alice Baniel credits. Everything else (fastsim, fastmaps, fastcheck,
fasthelper) is Raphaël Ribes' own work. See each tool's `CITATION.cff`.

---

## FastSuite 1 (2026-10-12)

| tool | version | maturity | notes |
|---|---|---|---|
| fastdedup | 1.3.0 | **WIP** | exact PCR-dedup by 64/128-bit fingerprints; published (bioRxiv 2026.04.29.721745) |
| fastprep | 0.1.0 | WIP | swept-winner default configuration shipped (validated on 90/90 datasets) |
| fastsim | 0.1.0 | WIP | measured insert CDF from the profile; single-end mode; ground-truth TSV |
| fastcheck | 0.1.0 | WIP | exact one-pass QC; the Python differential harness; `--profile` fitter |
| fastmaps | 0.1.0 | WIP | composition plane; `-j/--threads` added in this release |
| fasthelper | 0.1.0 | WIP | run / concatenate / compare / doctor / schema |

Suite-level work in this release: the six tools converge on the canonical
CLI option set (gaps tracked in `docs/CLI-SPEC.md` §4, closed by the plans
under `docs/plans/2026-10-12_*.md`), `fastcheck` is wired into
`fasthelper doctor` and the run manifest, and `fasthelper compare` no longer
fails on a missing output directory.

---

## Before the first bioconda publication (checklist)

* [ ] Every tool has a `meta.yaml` recipe (`recipes/<tool>/meta.yaml`) —
      `run:` and `test:` sections invoke `--version` and a **synthetic**
      fixture run only (**D-11**: never ship any data).
* [ ] Every tool has a `LICENSE` (MIT; fastprep's carries the Mandrillus
      Project copyright line).
* [ ] Every tool has a `CITATION.cff`.
* [ ] Recipes download nothing at build time.
* [ ] The suite is only claimed to replace what `docs/REPLACEMENT-MAP.md`
      grants: "FastQC's modules, replaced and extended, minus per-tile and
      exact over-representation" and "MultiQC's aggregation, replaced".
