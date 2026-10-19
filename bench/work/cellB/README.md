# Cell B — the synthetic PE-overlap value-test cohort

A `fastsim` cohort whose ground truth carries **both** G-runs and adapter
readthrough, built to answer whether fastprep's `--pe-overlap` improves trimming.
It is the *value* cell of
`docs/plans/2026-10-19_004500-pe-overlap-with-without.md`; the real corpus is the
*guard* cell and cannot answer this question (readthrough in 11,395 of 10.2M
pairs, 0.11%).

---

## 1. The confound this cohort exists to avoid

`bench/scripts/trim_report.py` classifies every mate by comparing the tool's
trim positions against the **ground truth's G-run lengths**:

```python
if abs(t5-h) <= tol and abs(t3-t) <= tol:  'exact'
elif t5 > h or t3 > t:                    'overtrim'
```

An **adapter** trim is not a G-run. On a cohort whose GT has no G-runs, every
adapter trim the tool makes is labelled `overtrim` — a false verdict. That is
what Cell A measured (8,533 of 11,299 readthrough pairs labelled overtrim) and
why Cell A is unreadable.

So the cohort must carry both halves. It does — verified on the built cohort
(`python3 verify_cellB.py …`, sha256 of the re-simulated GT identical):

| | |
|---|---|
| mates carrying an injected G-run | **489,546 of 800,000** |
| mates with a 3′ run the polyX scan can reach (`t >= 12`) | **25,526** |
| pairs carrying adapter readthrough | **112,658 of 400,000 (14.08%)** |
| readthrough mates that ALSO carry a G-run | 74,916 |
| readthrough pairs inside the overlap window | **0 — structural, see §5** |

**The confound is avoided.** An adapter trim on this cohort is no longer forced
into `overtrim`: the G-run boundary is a real, scorable reference.

---

## 2. What was built

### 2.1 The profile

`tools/fastsim/profiles/cellB-pe-overlap.json`, authored by
`bench/work/cellB/make_profile.py` (re-run it to regenerate the JSON).

| field | value | why |
|---|---|---|
| `read_len` | 150 | matches the corpus |
| `adapter.mode` | `readthrough` | the mode's trigger |
| `adapter.seq` | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` | the **real Illumina TruSeq R1** sequence — the same one fastprep's built-in catalogue carries, so the mode's *"the readthrough past the boundary really is adapter"* guard has something to match |
| `head` | 0:1000, 1:170, 2:170, 12–17:20/20/20/20/10/10 | mass at 0 **and** at ≥1. The 12–17 bp class exists so the plain polyX scan (`min_run 12`) has a 5′ reference; the 1–2 bp class is what `--head-dropout` exists for |
| `tail` | 0:720, 4–17 spanning `min_run` | mass at 0 **and** at ≥1, on both mates, 72 % clean / 28 % with a run |
| `joint_tail` | the diagonal of `tail` | correlated mates, as real NovaSeq tails are |
| `full_length_rate` | 0.0 | an all-G read is a drop class, not a trim class — it would only add a large `dropped` bucket with no bearing on the PE question |
| `n_run` | absent | N-runs are forced into `undertrim_irreducible` by construction and would only blur the readout |
| `dup_rate` | 0.0 | a duplicate replays a stored pair verbatim, which would distort the exact-rate denominators |
| `qual_noise` | 0.0 (default) | emitted qualities are the profile's own integers |

### 2.2 The insert distribution

Realised as fastsim's `{length: count}` histogram (the profile builder's own
convention: a saturating bin for lengths past the measurement ceiling
`2L-12 = 288`):

| component | share | range | role |
|---|---|---|---|
| overlap | 49.97 % | 150..288 | mates overlap; no readthrough |
| short | 15.03 % | 112..149 | readthrough 1..38 — **the mode's target** |
| saturated | 35.00 % | 301 | one bin, as faithful as a fitted tail |

The `short` component is what makes readthrough occur, and it is a *fixed 15 %*
of pairs rather than a rare tail — that is the whole point of the cohort.

### 2.3 Determinism

```
fastsim --pairs 400000 \
        --profile tools/fastsim/profiles/cellB-pe-overlap.json \
        --out <cohort>/work/synthetic/cellB-pe-overlap/v0/sim.v0 \
        --seed 4242 --force --drop-rate 0 --dup-rate 0
```

`--seed 4242` is fixed. Re-simulating produced a **byte-identical GT**
(sha256 `3b854613b7f7bba00744b1ace3f1c77edba76ba683f957d816f55f173221ed6f`),
so the ground truth every number below is scored against is reproducible.

### 2.4 Cohort layout

The layout `trim_report.py` discovers — it reads
`<work>/<domain>/<acc>/v0/{sim.v0.gt.tsv,sim.v0_R1.fastq.gz,sim.v0_R2.fastq.gz}`
plus `<results>/winners.tsv`:

```
<root>/cohort/work/synthetic/cellB-pe-overlap/v0/sim.v0_R1.fastq.gz
<root>/cohort/work/synthetic/cellB-pe-overlap/v0/sim.v0_R2.fastq.gz
<root>/cohort/work/synthetic/cellB-pe-overlap/v0/sim.v0.gt.tsv
<root>/cohort/results/winners.tsv          one row: the published winner config
```

`winners.tsv` carries the suite's shipped default winner (min_run 12, mercy 3,
phred_sat 20, gate_shape 1.0, min_len 0.9L, mean_q 25, adapter_err 0.2,
min_overlap 14) — the config the corpus sweep picked on 100 of 101 datasets, so
the baselines reproduce the real harness rather than an arbitrary argv.

---

## 3. The measurement

Run by `bench/work/cellB/build_and_run_cellB.py` (or `run_cellB_remote.sh` on
the desktop). The harness's own mechanisms, no bespoke scorer:

```bash
# WITH the mode
FASTSCRATCH=<scratch>/with-pe \
FASTPREP_VARIANTS='[["fp-pe",["--pe-overlap"]]]' \
python3 bench/scripts/trim_report.py --mode decomp \
    --work <root>/cohort/work --results <root>/cohort/results \
    --out <out>/cellB-with-pe.tsv --tol 0

# WITHOUT (the baseline; FASTPREP_VARIANTS unset)
FASTSCRATCH=<scratch>/without-pe \
python3 bench/scripts/trim_report.py --mode decomp \
    --work <root>/cohort/work --results <root>/cohort/results \
    --out <out>/cellB-without-pe.tsv --tol 0
```

The report refuses to overwrite an existing TSV, so each run gets its own
`--out`. The instrument checks pass in both runs
(`per-mate sums == pooled: PASS`, `slice sums == pooled: PASS`).

Outputs: `<out>/cellB-{with-pe,without-pe}.tsv`, `-permate.tsv`, `-pe.tsv`,
`-permate-pe.tsv` (8 files, under
`/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap/results/`).

### 3.1 Pooled decomposition, tolerance 0

`cellB-with-pe.tsv` (the without-pe run is identical for hd0 / hd4 / fastp, as
expected — those rows carry no flag):

| tool | exact | overtrim | undertrim_rec | undertrim_irr | dropped | scored |
|---|---|---|---|---|---|---|
| fastprep-hd0 | 347,171 | 18,021 | 228,567 | 74,968 | 131,273 | 668,727 |
| fastprep-hd4 | 381,196 | 131,411 | 141,135 | 9,809 | 136,449 | 663,551 |
| **fp-pe** | **347,171** | **18,127** | **228,567** | **74,809** | **131,326** | **668,674** |
| fastp-default | 364,139 | 129,497 | 239,595 | 66,769 | 0 | 800,000 |
| fastp-tuned | 364,139 | 129,497 | 239,595 | 66,769 | 0 | 800,000 |

`fastprep-hd0` is bit-identical between the two runs — the baseline reproduces,
as it must.

### 3.2 The readthrough slice (the primary metric)

`cellB-with-pe-pe.tsv`, `n = 112,658` readthrough mates:

| tool | exact_rt | overtrim_rt | exact rate |
|---|---|---|---|
| fastprep-hd0 | 19,892 | 3,713 | 17.657 % |
| fastprep-hd4 | 22,363 | 9,965 | 19.850 % |
| **fp-pe** | **19,892** | **3,819** | **17.657 %** |
| fastp-default | 3,677 | 100,804 | 3.264 % |

**`fp-pe` == `fastprep-hd0` on the readthrough slice to the mate.** exact is
identical (19,892); overtrim moves by +106 of 112,658 (+0.094 pp). The mode
fires (`pe_trims = 5,328`) but its effect on the scored counters is ~0.

### 3.3 Per mate

`cellB-with-pe-permate.tsv` (tolerance 0):

| tool | exact_R1 | exact_R2 | overt_R1 | overt_R2 | rec_R1 | rec_R2 | irr_R1 | irr_R2 |
|---|---|---|---|---|---|---|---|---|
| fastprep-hd0 | 173,661 | 173,510 | 9,096 | 8,925 | 114,074 | 114,493 | 37,602 | 37,366 |
| fastprep-hd4 | 190,893 | 190,303 | 65,686 | 65,725 | 70,369 | 70,766 | 4,887 | 4,922 |
| **fp-pe** | **173,661** | **173,510** | **9,096** | **9,031** | **114,074** | **114,493** | **37,602** | **37,207** |

The only mate-level movement is `overt_R2` (+106) against `irr_R2` (−159).

---

## 4. Why the effect is ~0 — measured, not asserted

Three measurements explain it. All are scripts in this directory.

### 4.1 The mode does fire, and on the right mates

`explain_pe_effect.py` re-runs both configs and joins the PXMTs to the GT. All
5,328 changed decisions are on **readthrough mates** — the mode does not fire on
clean reads (the guard's question).

`verify_deltas.py` recomputes the whole decomposition from the PXMTs and the GT
with no reference to the report, and **reconciles exactly** against the report's
own deltas (`RECONCILED: YES`):

| class | baseline | fp-pe | delta | report delta |
|---|---|---|---|---|
| exact | 347,171 | 347,171 | +0 | +0 |
| overtrim | 18,021 | 18,127 | **+106** | +106 |
| undertrim_recoverable | 285,003 | 285,003 | +0 | +0 |
| undertrim_irreducible | 18,532 | 18,373 | **−159** | −159 |
| dropped | 131,273 | 131,326 | **+53** | +53 |

The class movement behind those deltas:

```
dropped -> dropped                    4866   (invisible to the scorer)
dropped -> overtrim                    198
overtrim -> dropped                    105
undertrim_irreducible -> dropped       146
undertrim_irreducible -> overtrim       13
```

**4,866 of the 5,328 changed decisions are `dropped -> dropped`** — the tool
drops the mate in both arms, but the drop record's trim positions differ (all
4,866 are a 1–2 bp shorter `t3` on the drop record). That is a real change in
what the tool *decided* that the scorer cannot see at all, because a dropped
mate is scored by nobody. `invisible_decisions.py` isolates it.

The 462 that do move a scored class are net-negative: +211 into `overtrim`,
−251 into `dropped`, i.e. the mode makes the scored decisions slightly worse,
not better.

The direction is the important part. The mode trims the 3′ end by 1–2 bases
*less* than the baseline did, on mates that already carried an adapter the
baseline had already found. That is a correction of an over-trim, not the
discovery of an adapter — and it is worth 0.000 pp of `exact`.

### 4.2 The addressable set is tiny — and mostly unreachable

`why_small_effect.py` buckets the 109,013 readthrough mates on two mechanisms:

| bucket | mates | what it means |
|---|---|---|
| `t>=12 / a<14` | **1,328** | scan active, baseline matcher blind — **the addressable set** |
| `t>=12 / a>=14` | 1,076 | scan active, baseline matcher already sees it |
| `t<12 / a<14` | 39,820 | scan blind — no run for a PE boundary to interact with |
| `t<12 / a>=14` | 66,789 | scan blind, matcher sees it |

The mode can only help where the baseline's matcher is blind (`a < 14`) *and*
the polyX scan is active (`t >= 12`) — **1,328 of 109,013 readthrough mates
(1.2 %)**. And of those, the ones whose readthrough is long enough for the
mode's own detector to fire (`readthrough_is_adapter` needs `t3 == readthrough`
with `>= 14` matched bases) are a subset again.

The reason the rest is unreachable is structural, and worth stating plainly:
the readthrough's length is `L - h - t - min(I, L-h-t)`, i.e. it **shrinks as
the tail run grows**. A mate with a long 3′ polyG tail has very little adapter
left, so both the baseline matcher and the mode's own adapter guard look at a
short tail. The cohort's design puts a 3′ run on 28 % of mates, which is what
makes the G-run boundary a scorable reference — and that same mass is what
starves the adapter.

### 4.3 `readthrough` and `in_window` cannot overlap

`window_readthrough_check.py`, on the built cohort:

```
  readthrough & in_window            0
  readthrough & below_window         109013
    matches I < L - t                109013
    CONTRADICTION                    0
```

`readthrough` needs `I < L - t`; `in_window` needs `I >= L`. The two are
**disjoint by definition** — the `in_window=1, rt=1` cell is empty because the
definitions cannot both hold, not because of a bug in the cohort or the report.
The plan's Cell B wording ("inserts concentrated in `L <= I <= 2L-12` *and*
mates carrying adapter readthrough") describes a set that does not exist; the
PE mode's real operating region is the **below-window** readthrough mates.

### 4.4 The `dropped` bucket is large and mostly correct

`why_dropped.py` on the baseline PXMT:

```
  emitted                                              668727
  dropped                                              131273
  dropped: keep window < min_len (correct drop)         93687
  dropped: GT says KEEP, tool dropped anyway            37586
  dropped: GT says full (correct drop)                      0
```

71 % of drops are required by `--min-len 0.9L` (keep window `< 135`), 29 % are
drops the GT would not require. Both runs' drop rates are within 0.04 pp of
each other, so the kept-set the exact rate is earned on is effectively the same
in both arms.

---

## 5. Verdict on the mode

**Null result.** On a cohort built specifically to give `--pe-overlap` something
to find — 14 % of pairs carrying real adapter readthrough, 50 % of pairs with
overlapping mates — the mode changes the scored counters by ~0:

* pooled `exact` **347,171 = 347,171** (identical);
* readthrough-slice `exact` **19,892 = 19,892** (identical);
* readthrough-slice `overtrim` 3,713 → 3,819 (**+106 of 112,658, +0.094 pp**);
* `pe_trims` = 5,328, all on readthrough mates — **4,866 of them `dropped ->
  dropped`**, which the scorer cannot see at all; the 462 that do move a class
  are net-negative (+211 overtrim, −251 dropped).

Against the plan's proposed threshold (exact-rate gain on readthrough pairs
≥ 0.5 pp), the measured gain is **+0.000 pp**. The decision rule in the plan
says: *do not restore; record both cells.*

Two limits on the claim, stated rather than buried:

1. **Synthetic, so conditional on the profile.** The measured value is "on a
   library with readthrough in the overlap window", not "on real Illumina
   data". The corpus cannot answer the latter because it has almost no
   readthrough.
2. **The addressable set is 1.2 % of readthrough mates on this profile**, and
   that is driven by the profile's own 3′-run mass, which is what the confound
   avoidance requires. A profile with readthrough but no 3′ runs would be
   unreadable for the opposite reason. The two requirements pull against each
   other; this cohort resolves the tension in favour of readability, and the
   price is a small addressable set.

---

## 6. Files

| file | role |
|---|---|
| `tools/fastsim/profiles/cellB-pe-overlap.json` | the profile (regenerate with `make_profile.py`) |
| `bench/work/cellB/make_profile.py` | authors the profile |
| `bench/work/cellB/build_and_run_cellB.py` | builds the cohort and runs both reports (local or desktop) |
| `bench/work/cellB/run_cellB_remote.sh` | the desktop invocation |
| `bench/work/cellB/audit_cohort.py` | the confound/value check (G-runs? readthrough?) |
| `bench/work/cellB/verify_cellB.py` | determinism (same seed → same GT) + the confound, on the built cohort |
| `bench/work/cellB/explain_pe_effect.py` | `pe_trims` by mate and direction, joined to the GT |
| `bench/work/cellB/verify_deltas.py` | recomputes the decomposition from PXMT + GT and asserts it matches the report's deltas |
| `bench/work/cellB/reconcile_pe_trims.py` | changed decisions vs scored-class movement, on the full cohort |
| `bench/work/cellB/invisible_decisions.py` | the `dropped -> dropped` decisions the scorer cannot see |
| `bench/work/cellB/why_small_effect.py` | the addressable-set buckets |
| `bench/work/cellB/window_readthrough_check.py` | proves `readthrough ∩ in_window = ∅` |
| `bench/work/cellB/why_dropped.py` | decomposes the `dropped` bucket against the GT |
| `bench/work/cellB/probe_pxmt.py` | PXMT reader for the raw per-mate decisions |

Cohort + results roots (desktop), both written by the scripts in this dir:

```
/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap/cellB/cohort/work/synthetic/cellB-pe-overlap/v0/
/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap/results/cellB-*.tsv              (the measurement)
/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap/results-fresh/cellB-*.tsv        (the re-run check)
```
