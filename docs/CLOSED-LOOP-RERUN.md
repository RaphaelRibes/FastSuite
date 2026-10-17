# The closed loop, re-measured on a regenerated corpus

**Result: 142/306 = 46.4%**, replacing the 167/306 = 54.6% that three papers
cite. The old figure was inflated by a corpus artefact; the new one is lower
for a reason that is itself a finding.

## 1. Why the old number was wrong

The Phase D corpus stores three variants per dataset, meant to be independent
draws (`chain.sh` uses `SEED = 42 + VARIANT*1000`, so 42/1042/2042). They were
not independent:

```
$ md5sum work/<domain>/<acc>/v1/sim.v1.gt.tsv  work/<domain>/<acc>/v2/sim.v2.gt.tsv
```

Byte-identical in **100 of 102 datasets**, written in the same second. Two
seeds cannot produce that. The evidence is preserved in
`results/STALE-BEFORE-RERUN.txt` (written before the tree was rebuilt).

**Cause.** `run_campaign.sh` skips any `(dataset, variant)` with a status file,
and `chain.sh` guards every artefact with `if [[ ! -s file ]]`. When v0 was
regenerated after an earlier fix, v1 and v2 were silently reused. The seed
handling is correct in the current code — three seeds produce three different
outputs, verified — so this was a stale-artefact problem, not an RNG bug.

**Effect on the number.** All 18 duplicate-rate failures came from the
duplicated v1/v2 pair, and v0 had **zero**. So 18 of the 306 rows were the same
9 datasets counted twice, and the pass rate was never measured over 306
independent comparisons.

## 2. The re-run

Three independent changes invalidated every derived artefact, so a patched
re-run was not possible:

1. the polyG fitter fix (`a933980`) changes what `fastcheck --profile` emits;
2. the N-emission work (`fc1bc78`) adds RNG draws when a profile carries
   `n_run`, shifting the stream;
3. the v1/v2 duplication above.

`bench/scripts/rerun-campaign.sh` deletes each variant's work dir before
rebuilding, which `chain.sh` cannot do for itself. 102 datasets × 3 variants,
57-config grid, ~40 min at jobs=8.

### Acceptance

| check | before | after |
|---|---|---|
| `v1==v2` on every metric | 100/102 | **0/102** |
| duplicate-rate failures | 18 | **0** |
| closed loop | 167/306 = 54.6% | **142/306 = 46.4%** |

Per domain: human_gut 94.1%, bacterial_isolate 64.7%, generalised 58.8%,
virome 35.3%, rna_seq 19.6%, metagenomics 5.9%.

`g_frac` now appears in **all 164 failing rows** — it is the sole blocker.

## 3. Why it went down: the fitter trades one metric for another

Removing the contamination did not raise the number. On v0 (never duplicated)
there are **14 pass→fail and 0 fail→pass** — a one-way regression, with the
g_frac mismatch 3–6× worse (e.g. `DRR463958` 0.0033 → 0.0213).

**Isolated by stripping `n_run` from the profile and re-running, same seed:**

| effect of removing the N emission | datasets |
|---|---|
| both fail — emission irrelevant | **12 / 14** |
| tips a marginal case to FAIL | 2 / 14 |

Every one of the 12 still fails at 0.0205–0.0395 with the emission removed,
and the values barely move (in three cases the emission makes it *better*).
**So the N emission is not the cause; the polyG fitter fix is.**

**Mechanism.** The fitter's boundary deconvolution assumes terminal bases are
i.i.d. at the measured background G rate. That is *exact* for the simulator's
random DNA — which is why fitter recovery went **8/10 → 10/10** — and an
*approximation* on real libraries with strong 3′ composition bias. The closed
loop fits from **real** data, and the failure lands on the read's **last base**,
where polyG lives. `DRR001194` moves 0.0069 → 0.0198 from the fitter alone.

**The honest reading:** the fitter became better at the metric it was tuned
against and worse at the metric the closed loop scores. Both are true of the
same change, and this campaign is the first time they have been separated.

## 4. Where the fix goes

Not a constant, and not a wider tolerance. The deconvolution needs a
**per-position boundary model**: instead of one background G rate over the read
interior, model the terminal composition directly. The affected domains are the
3′-biased ones (rna_seq, virome), which is consistent with the mechanism.

## 5. What may not be cited

Anything generated from a profile carrying `n_run` before commit `fc1bc78` is a
different realisation. The staleness guard and the input-identity resume key
added in `73c2385` turn that into a refusal rather than a silently wrong
number, but any sim FASTQ, `sweep.c*.pxmt`, or trim-accuracy row produced
before this re-run must be regenerated, not reused.

## 6. The tooling lesson: three silencers, one shape

Every wrong number in this re-run came from a tool that **quietly reused or
mixed state instead of refusing**. That is a pattern, not bad luck, and each
case has now been closed:

| tool | what it did | guard |
|---|---|---|
| `chain.sh` | `if [[ ! -s file ]]` — reused any existing artefact, which is how v1/v2 survived a fix | `rerun-campaign.sh` deletes each work dir before rebuilding |
| `closed_loop.py` | refused to overwrite `--row`, so a re-score that never ran printed the *previous* run's numbers as if fresh | already had `--force`; the re-run must pass it |
| `overtrim_vs_fastp.py` | `emit` opens `--out` in mode `"a"`, so a second run **appended** to the first — a file where every key appeared twice, once per corpus | now refuses unless `--force`, checked *before* any work |

The general rule these three arrived at: **a script must fail loudly when its
output already exists.** An artefact that silently outlives its inputs is
indistinguishable from a correct one, and the cost of finding out is a full
campaign rebuild — as this one required.

Anything else added under `bench/scripts/` that writes a result should refuse
by default and take an explicit overwrite flag.

