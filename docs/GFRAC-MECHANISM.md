# g_frac: which per-position mechanism fails the metagenomics closed loop

**Result: the metagenomics `g_frac` blocker is not the fitter's polyG boundary
deconvolution. It is three separate defects in the simulator's read assembly
(two of them in the artefact/base-plane geometry of `fastsim`'s profile
consumer), amplified or masked by a fourth, harness-side, effect: the report
pools both mates while the simulator always emits `read_len` bases for both.**

Everything below is measured on the Phase D corpus
(`/mnt/Fast/fastsuite-scratch`, 16 failing + 1 passing metagenomics datasets
plus five cross-domain controls). Result artifacts (timestamp 20261018-1815):

| file | what it is |
|---|---|
| `bench/results/phase-d/gfrac-attr-20261018-1815.tsv` | per (dataset, variant) attribution: curves, argmax, zone means, per-mate terminal G, length-matched counterfactual |
| `bench/results/phase-d/gfrac-mech-20261018-1815.tsv` | the same rows joined with the fitted run-length rates, run-structure agreement and the mechanism columns |
| `bench/results/phase-d/gfrac-arms-20261018-1815.log` | the profile-mutation arms (re-simulated at the same seed) |
| `bench/results/phase-d/gfrac-frameshift-20261018-1815.log` | the frame-shift model checked against the measured simulation, position by position |

Scripts: `bench/scripts/gfrac_attribution.py` (replaces
`verify_gfrac_deconv.py`), `gfrac_mechanism_summary.py`, `probe_plane_arms.py`,
`probe_frame_shift.py`, `probe_head_gfrac.py`, `probe_mate_lengths.py`.

**No Rust file was changed by the diagnosis itself.** The reasons are in §5. Two
follow-up commits then acted on it: `7631292` implemented the fix order's steps
1(a) and 1(b) — the harness's length-matched/per-mate gate, and all four
(plane × artefact) deconvolutions in fastsim — and everything in §1–§4 describes
the state measured **before** that. Where §3 says a pair is "not deconvolved",
read it as "not deconvolved at diagnosis time"; the update is inline below.

## 1. Where the mismatch sits

The closed loop scores `g_frac_max_diff` = max over 1-based read positions of
|G_real − G_sim| on the report's `base_position` plane, which is **pooled over
both mates**. Per failing dataset, the argmax splits the cluster cleanly:

| group | datasets | argmax | real/sim mate lengths | scored | length-matched |
|---|---|---|---|---|---|
| A last base | `DRR083063 065 212 291 326 337 353 390 411` | pos 101 (terminal) | 101 / **99** | 0.0275–0.0411 | **0.098–0.112** |
| B first base | `DRR070915 083082 083283 083295 083328 083438` | pos 1 | 150/150 or 101/101 | 0.0657–0.0774 | same |
| C pure mix | `DRR083073` | pos 97 | 101 / **91** | 0.1029 | **0.0358** |

Controls (all symmetric lengths, no plane defect in play): passing
metagenomics `DRR049542` 0.0044; human_gut `DRR425490` 0.0051 / `DRR425489`
0.0054; bacterial `DRR000978` 0.0126; generalised `DRR000426` 0.0139; virome
`DRR173086` 0.0162; rna_seq `DRR002131` 0.0184.

## 2. Mechanism A — sample mix: the simulator pads the second mate

`profile.read_len` is one number, taken from the first mate
(`fastcheck/src/profile.rs::write_profile`), and `fastsim` emits exactly
`read_len` bases for both mates (`main.rs`: "Sequencer output is always exactly
read_len bases", a `debug_assert`). Real libraries need not be symmetric:
group A is R1 101 / R2 **99**, group C is R1 101 / R2 **91**. At the positions
only R1 reaches in the real library, the *pooled* simulated value is a two-mate
average, so the comparison is diluted and the gate reads a mix effect, not a
composition error.

Measured, DRR083073 (101/91): scored **0.1029** at position 97 → length-matched
(re-weight the sim's per-mate curves by the *real* per-mate depths) **0.0358**.
65% of that row's gate is the pooling artefact. Group A shows the opposite
sign: scored 0.028–0.041 but length-matched 0.098–0.112 — the pooling *hides* a
2.5–3× larger per-mate error (group A also has mechanism B).

## 3. Mechanism B — the injected run is counted twice, on both planes

`profile.base[m]` is the **measured** per-position composition, so at the
positions an injected G-run reaches it already contains those Gs. `fastsim`
deconvolves exactly one of the four (plane × artefact) pairs: `base.R1` against
the **head** run (`profile.rs`: `p_mass = P(h > p)`). Not deconvolved:

* `base.R1` against the **3′** run → R1's terminal base is drawn from a plane
  value that already contains the tail artefact, and the tail run is injected
  on top (group A).
* `base.R2` (consumed revcomp'd as the fragment's 3′ end) against R2's **5′**
  run → same double count at R2's first base (group B). R2's 3′ run is a third,
  smaller case.

**Update (`7631292`): all four pairs are now deconvolved** — `profile.rs` builds
`head_cover`/`tail_cover` and applies them to both mates, revcomp'ing R2's tail
plane, which is why the code cites this section. Measured on the corpus, old →
new: DRR083326's R1 terminal **0.5270 → 0.4166** against a real 0.4245 (and the
old hand-built `fixed_planes` arm read 0.4140), its length-matched gate
**0.1026 → 0.0268**; group A 0.1026–0.1092 → 0.0268–0.0306; group B 0.0655–0.0756
→ 0.0387–0.0409. The **pooled** gate for group A moves *up* (0.0387 → 0.0960)
because those rows sit on two cancelling errors — which is why step 1(a) had to
land first, and why the metagenomics pass rate is unchanged at 3/51: what remains
is mechanism C (the head-run frame, §4) and the per-mate read length (§5).

Measured per mate (`gfrac_attribution.py --per-mate` for the pooled plane,
`probe_plane_arms.py` for the per-mate positions; "model" is the simulator's own
arithmetic `E[terminal G] = Σ_h P(h)·[P(t≥1) + P(t<1)·plane[L−1−h]]`, which
includes the head-length shift):

| dataset | real R1 terminal G | sim R1 terminal G | model | real R2 pos-1 G | sim R2 pos-1 G |
|---|---|---|---|---|---|
| DRR083326 (A) | 0.4245 | **0.5270** | 0.5285 | 0.1804 | 0.1811 |
| DRR083063 (A) | 0.4166 | **0.5243** | 0.5528† | 0.1788 | 0.1804 |
| DRR083082 (B) | 0.1779 | 0.1778 | 0.1779 | 0.3946 | **0.5461** |
| DRR083295 (B) | 0.1817 | 0.1913 | 0.1893 | 0.3918 | **0.5389** |
| DRR083438 (B) | 0.2104 | 0.2177 | 0.2164 | 0.4138 | **0.5518** |

† the plain product `P(t≥1) + P(t<1)·plane[L−1]` without the head mixture; with
the mixture the model lands on the measured value (DRR083326: 0.5285 model
against 0.5270 measured).

The signature is visible without any per-mate work: the share of reads whose
**first** base is G (`polyx 5p G` total / reads) is inflated by ~0.07 on every
group-B dataset, and by ≤0.006 on the controls:

| dataset | real 5′G share | sim 5′G share |
|---|---|---|
| DRR083082 | 0.3925 | **0.4681** |
| DRR083328 | 0.4244 | **0.4957** |
| DRR083438 | 0.4088 | **0.4773** |
| DRR083295 | 0.3888 | **0.4632** |
| DRR083283 | 0.4178 | **0.4884** |
| DRR070915 | 0.4188 | **0.4847** |
| DRR049542 (passes) | 0.2491 | 0.2534 |

Independent cross-domain confirmation of the same double count on the 3′ side:
`bacterial_isolate/DRR000466` fails at 0.0254 (position 110, terminal) with a
fitted tail rate of 0.0319 and `plane_G[last] = 0.2441` →
P(t≥1)·(1 − plane_G) = **0.0241** against the measured **0.0254**.

## 4. Mechanism C — the head run is prepended, so the composition is smeared

`build_read` composes `[G; h] + insert[...] + adapter + [G; t]`: a read whose
head run is `h` shows the fragment's base `p − h` at its own position `p`. The
`base` plane is a **read-position** measurement, so the realised read-position
curve is the measured curve convolved with the head-run length distribution —
which is why the library's own periodic G structure (period ~5–6, e.g.
DRR083082 real G: 0.1044 @5, 0.2547 @6, 0.1355 @10, 0.3085 @11, 0.1756 @12) is
attenuated in the simulation.

Tested as a model, not asserted: `sim[p] = Σ_h P(h)·real[p−h]` reproduces the
measured simulation curve to **≤0.0037** (DRR083082; worst raw difference
0.0405 @11 → residual 0.0014) and ≤0.0068 (DRR083295; raw 0.0392 @11). Full
per-position table in `gfrac-frameshift-*.log`. Removing the head run entirely
(`head = {0}`, a control arm, not a fix) reproduces the whole real curve to
**0.0066** max on DRR083082 (0.0054–0.0288 across the four group-B rows
measured, §6.4) — the plane itself is right, and everywhere outside the
defective regions the sim matches the real curve to ≤0.008.

## 5. What the evidence rules out, and why no code was changed

* **The fitter's polyG boundary deconvolution — ruled out.** (a) On all 16
  failing metagenomics datasets the round trip of the run-length distribution
  (measured {@len: counts} → fitted artefact distribution → injected → measured)
  agrees with the real library to TV ≤ 0.055 (`tail_run_TV`, `head_run_TV` in
  `gfrac-mech`), mostly ≤ 0.03. (b) Group B has `tail = {"0": ~10000}` — no 3′
  artefact exists in those libraries, so a 3′-deconvolution defect cannot
  produce a position-1 error. (c) The only visibly mis-fitted run distribution
  in the corpus is virome `DRR173086` (tail TV **0.1805**) and it **passes**
  (0.0162) — a 0.18 TV on the *run-length* distribution does not by itself break
  the gate.
* **The N overlay — ruled out.** `n_overlay_max_absd` = 0.0000 on every row
  (the N-aware denominator moves no position's G fraction).
* **GC under-model — ruled out.** body-zone mean |Δ| ≤ 0.008 on every failing
  row (median 0.0015); `gc_diff` ≤ 0.0013.
* **Why not the fitter's `profile.rs`:** the fitter is emitting what it should —
  the measured plane and the measured run-length distributions. The defect is in
  how the *consumer* uses them, and the precedent for the fix is the 5′
  deconvolution, which lives in `tools/fastsim/src/profile.rs`, not in the
  fitter. Pre-deconvolving the plane in the fitter would be a regression: an
  arm that deconvolved `base.R1` for the head run (mirroring what fastsim
  already does) drove the sim's R1 position 1 to 0.246 against a real 0.3913 —
  a double subtraction.

## 6. The arms: what a fix would do, measured end to end

Profile-mutation arms, re-simulated at the same seed with the same binary, so
the only difference is the mutated block (`probe_plane_arms.py`; the pooled
gate is the closed loop's, the length-matched gate is §2):

| dataset | arm | pooled gate | length-matched | R1 pos-1 (real) | R2 pos-1 (real) | R1 terminal (real) |
|---|---|---|---|---|---|---|
| DRR083326 | as fitted | **0.0383** @101 | 0.1047 @101 | 0.4114 (0.4099) | 0.1811 (0.1804) | 0.5292 (0.4245) |
| DRR083326 | `base.R1` deconvolved for the 3′ run | **0.0969** @101 | 0.0288 @100 | 0.4124 | 0.1779 | 0.4140 (0.4245) |
| DRR083082 | as fitted | **0.0749** @1 | 0.0750 | 0.3896 (0.3913) | 0.5461 (0.3946) | 0.1785 (0.1779) |
| DRR083082 | `base.R2` deconvolved for the R2 5′ run | **0.0401** @11 | 0.0401 | 0.3931 | 0.3934 | 0.1779 |
| DRR083438 | as fitted | **0.0677** @1 | 0.0677 | 0.4016 (0.4044) | 0.5518 (0.4138) | 0.2193 (0.2104) |
| DRR083438 | same | **0.0403** @11 | 0.0403 | 0.4077 | 0.4115 | 0.2191 |
| DRR083073 | as fitted | **0.1030** @97 | 0.0372 @98 | 0.3462 (0.3496) | 0.3654 (0.3068) | 0.1582 (0.1580) |
| DRR083073 | every plane arm | 0.1024–0.1030 | 0.0341–0.0370 | unchanged | unchanged | unchanged |

Three findings that decide what the fix must be:

1. **The 3′ deconvolution is correct and its obvious application makes the
   judged metric worse.** On DRR083326 it repairs the per-mate terminal
   (0.5292 → 0.4140 against a real 0.4245) and improves the length-matched gate
   (0.1047 → 0.0288), but the *pooled* gate the closed loop scores goes
   0.0383 → **0.0969**, because the pooled row currently sits on the
   cancellation of two errors (double count up, second-mate dilution down).
   Group A is 9 datasets of exactly this shape.
2. **Fixing one plane is not enough.** The R2-head deconvolution removes the
   position-1 error on group B (0.075 → 0.040) but the residual then sits at
   position 11 — mechanism C, the frame shift — still 2× the 0.02 bound. No
   single-plane patch passes a group-B dataset.
3. **Group C is untouched by any plane arm** (0.1030 → 0.1024): it is the
   sample mix alone.
4. **Group B needs the head defects *and* (on some rows) the R2 3′ one.** The
   `head_none` control arm removes the head run everywhere — the limit in which
   both the double count *and* the frame shift are gone, but so is the artefact
   — and it lands at pooled max **0.0054** (DRR070915), **0.0066**
   (DRR083082), **0.0151** (DRR083295) but **0.0288** (DRR083438, still over
   the bound). DRR083438's residual is on the other side: its real R2 terminal
   G is 0.2993 against a simulated 0.3503, and its fitted R2 3′ run rate is
   0.0893; the `base.R2`-for-3′-run arm takes the pooled terminal from 0.2862
   to 0.2635 (real 0.2548). So a correct head fix sits between `head_none` and
   the as-fitted curve, and three of four group-B rows need the R2 3′
   deconvolution on top of it.

## 7. Where the fix goes, in order

1. **Harness first (no simulator change).** The scored metric pools two mates;
   the simulator can only ever emit two equal-length mates today. Until the
   profile carries per-mate read lengths, `g_frac` for asymmetric libraries
   measures the pooling convention. Either compare per mate, or state the
   length-matched gate as the primary statistic — on the current corpus the two
   differ by up to 0.08 in both directions.
2. **Double counts (small, consumer side).** `fastsim/src/profile.rs`, next to
   the existing 5′ deconvolution: deconvolve `base.R1` for the 3′ run and
   `base.R2` for the R2 5′ and 3′ runs (`probe_plane_arms.py::deconvolve` is a
   working reference implementation, and the arms above are its measured
   effect). It must be paired with (1), and with (3), or it lowers the scored
   number while raising the true accuracy.
3. **Frame shift (geometry, not small).** `artefact.rs::build_read` prepends the
   head run, shifting the biology by `h`; the measured plane is a read-position
   composition, so the shift smears it. Fixing it means placing the head run as
   an overwrite of the read's first `h` bases, which moves `bio_len`, the trim
   window in the ground truth and the insert window — a model change, not a
   patch.
4. **Per-mate read length (profile schema + geometry).** Needed for group A and
   C honestly; `read_len` is one scalar today and `main.rs` asserts both mates
   equal it.

## 8. Open / not determined

* Whether the real libraries' first-cycle Gs replace the read's first bases or
  shift the frame is a chemistry question this corpus can only bound: the
  measured read-position profile is the thing the loop compares, and the
  simulator does not reproduce it (mechanism C) — but which convention is
  physically right is not settled here.
* Group A's residual after both (2) and (3) is not measured (that needs a
  fastsim geometry patch, out of scope for a profile mutation).
* Only 17 metagenomics datasets and 5 controls were run per-mate; the virome
  and rna_seq domains (also 3′-biased, 35% and 20% pass) were sampled, not
  swept, and their tail-run TV (virome 0.18) suggests mechanism B/C reach there
  too.
* The `sample_mix_component` column is a counterfactual pooling, not a
  re-simulation: it re-weights the *simulated* per-mate curves, which is exact
  for the pooling effect but does not model what the simulator would emit if it
  could emit a 99 bp R2.

---

## 9. Mechanism C: two fixes attempted, and why the naive one fails (2026-10-19)

Mechanism C was implemented twice, independently, and both attempts failed the
same way. Recording that here so the next attempt starts from the obstacle rather
than rediscovering it.

### The model is confirmed

`build_read` composes `[G*h] + insert + adapter + [G*t]` — the head run is
**prepended**, so a read whose head run is `h` shows the fragment's base `p − h`
at its own position `p`. The `base` plane is a read-position measurement, so what
fastsim consumes is the real curve convolved with the head-run distribution:

```text
sim[p] = Σ_h P(h) · real[p − h]
```

Verified analytically and against the diagnosis's own numbers: a fixed-`h`
subtraction leaves the smearing in place (a delta convolved with `P = {0: 0.5,
1: 0.5}` stays at 0.5 at `p = 1`), while the weighted deconvolution recovers it
exactly. So the inverse is a **deconvolution of the whole distribution P(h)**,
not the subtraction of a single `h`.

### Attempt 1 — fixed-`h` subtraction (subagent, 2026-10-19)

Subtracted one `h`'s worth of mass from each position. Measured end to end:

| dataset | arm | pooled | position |
|---|---|---|---|
| DRR083082 | asis | 0.0401 | 11 |
| DRR083082 | fix | **0.1298** | 1 |

**3× worse.** The subagent hit a model repetition loop before it could iterate,
so no further arms were run. This is the same failure the diagnosis predicted in
§5 for the fitter-side variant: subtracting where the artefact already was.

### Attempt 2 — distributional deconvolution (parent, 2026-10-19)

Implemented `deconvolve_frame` in `tools/fastsim/src/profile.rs`, mirroring the
structure of the existing `deconvolve_row`: a sequential (Z-domain) inversion
`real[p] = (sim[p] − Σ_{h≥1} P(h)·real[p−h]) / P(0)`, applied after the row
deconvolutions in `biased_plane`.

Correct in isolation — four new unit tests pass, including one that pins the
fixed-`h` arm as *wrong* (residue left at position 1) and one that checks the
round trip on a known analytic smearing.

**But it does not compose with the row deconvolution.** The pre-existing test
`profile_deconvolves_both_planes_for_both_injected_runs` fails:

```
R1 terminal biology G should be 0.2857, got 1
```

Position 3's G is *amplified* to 1.0 rather than inverted. The cause is the
composition, not the arithmetic: `deconvolve_row` has already removed the
artefact's mass **at** each position, and `deconvolve_frame` then divides by
`P(0)` **across** positions, so the correction is applied twice. This is the same
double-subtraction failure class as attempt 1, reached from the other side.

The work was reverted rather than committed: `tools/fastsim` is untouched (93
tests today) and the gate is unchanged at 3/51.

### What the obstacle actually is

The two operations are not independent, and neither order is obviously right:

* **frame-then-row**: the frame step needs an un-smeared plane to work on, but
  the row step's coverage inputs (`head_cover`, `tail_cover`) are defined on the
  *read-position* plane — i.e. the smeared one.
* **row-then-frame** (attempt 2): double-corrects, because the row step has
  already removed mass the frame step then re-divides.

A correct joint inverse needs a decision about what the row deconvolution's
coverage should be computed on — the smeared plane (as now) or the un-smeared
one — and that is a modelling question about what `head_cover` *means*, not a
bug to be patched. §5's warning that a fitter-side deconvolution is a regression
is the same shape of problem.

### Recommendation

Do not attempt a third variant without first settling that question. The
diagnosis already notes that a correct head fix "sits between `head_none` and the
as-fitted curve" and that three of four group-B rows also need an R2 3′ fix — so
even a correct frame fix is unlikely to close the gate alone. The gate has been
stable at 3/51 through mechanisms A and B; a partial fix that moves it while
risking a double-subtraction regression is not obviously worth it.

What is worth keeping: the model is confirmed and quantified, and both failed
arms have their numbers recorded. Any future attempt should start by deciding
what `head_cover` is measured on.
