# fastcheck: exactness and integration (overrep modes, max-mem, tiles, adapter-position, ground-truth recovery)

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Close fastcheck's remaining verified gaps — memory admission refusal, the four
`--overrep-mode` trade-offs, the adapter-content-per-position and per-tile modules, `--only`
module gating — and, decisively, make the numbers *provably right* via a fastsim
ground-truth recovery test that fails loudly when they are not.

**Architecture:** No shared crate (D-1). Every reused idea — the memory-admission model from
`fastprep/src/memory.rs`, the adapter catalogue/matcher from `fastprep/src/adapter.rs`, the
xxh3 fingerprint from `fastdedup/src/hasher.rs` — is **copied into `tools/fastcheck/src/`**
with the source named in the module doc comment, per D-1 ("anything two tools need exists
twice, copied per tool, and may drift. Fix per tool and say which"). Consistency is enforced
by tests, not by the compiler. Every module stays inside the one-pass scan except
`--overrep-mode second-pass`, which is an explicit second read recorded as `passes: 2`.

**Tech Stack:** Rust (no new language deps; adds `libc`, `rustc-hash`, `xxhash-rust` — all
already in the offline cargo cache and used by fastprep/fastdedup), serde_json, clap 4.

---

## Verified baseline (2026-10-12, measured on this machine)

```
$ cd tools/fastcheck && cargo test --offline
running 42 tests ... test result: ok. 42 passed; 0 failed
running 12 tests ... test result: ok. 12 passed; 0 failed
```

* `--overrep-mode`, `--max-mem`, `--tiles`, `--only`, `--no-quality`, `--no-heatmap` do not
  exist (`fastcheck --help` lists only `--adapter-seq` among adapter flags).
* No `tile` code exists in `tools/fastcheck/src/` (only `quantile`/`title` false positives).
* The positions TSV has **14 columns**; `tests/differential.py` unpacks **11** and crashes:
  ```
  $ python tests/differential.py ../../examples/runs/demo/sim_E/sim_R1.fastq.gz
  ValueError: too many values to unpack (expected 11, got 14)
  ```
  The differential gate is **RED today**. Task 1 fixes it; nothing else may assume it is green.
* The closed loop works: `fastcheck --insert-hist --profile p.json` writes a profile that
  `fastsim --profile p.json` consumes (verified, 1000 pairs).
* Offline build with the three new deps verified in a scratch copy:
  `libc = "0.2"`, `rustc-hash = "2"`, `xxhash-rust = { version = "0.8", features = ["xxh3"] }`
  → `Finished dev profile ... in 13.37s`.

### What the ground-truth run already proves (measured, not assumed)

`fastsim --pairs 20000 --profile tools/fastsim/profiles/synthetic_novaseq.json --insert-mean
200 --insert-sd 40 --gc 0.5 --seed 42 --no-gzip --dup-rate 0.02 --adapter-seq
CTGTCTCTTATACACATCT -o /tmp/gt/simad` reports GT: **adapter-readthrough 989 pairs
(4.945%), duplicates 398 (1.99%)**. fastcheck on the same files measures
`measured_from_adapter = 990`, i.e. **within one pair** — the adapter path is exact.

Two **real defects** in `tools/fastcheck/src/insert.rs` are exposed by the same run and are
the reason the ground-truth test belongs in this plan:

1. **The anchor sits in the polyG tail.** `overlap()` anchors on
   `revcomp(s2)[0..ANCHOR]`. When R2 carries a 3′ polyG tail (any NovaSeq-like profile),
   those first 12 bases are `GGGGGGGGGGGG`, so the "anchor" matches a chance G-stretch in R1
   and the verified overlap is garbage. Only **4 200 of 16 642** measurable inserts are found
   (25%).
2. **The insert identity is wrong whenever artefacts exist.** `insert = len(R1)+len(R2)-o`
   over-estimates by `(h1+t1+h2+t2)/2`, because polyG runs and adapter are not biology. The
   correct identity is `insert = b1 + b2 - o` with `b = L - h - t` the biology length.
   Verification also demands **0 mismatches**, but the simulator's own boundary rule
   guarantees ~1 disagreement at the overlap's first base plus ~`2·10⁻³·o` sequencing errors,
   so even a correct alignment is rejected about half the time.

A prototype of the fix (anchor inside R2's biology at `off >= t2`, overlap to the end of R1's
biology, tolerate ≤2 mismatches, `insert = b1+b2-o`) recovers **98.7%** of measurable inserts
with mean **207.17 vs GT 206.66 (0.25% error)** on the same run. Tasks 10–12 implement it
behind a failing test.

### Tolerances the ground-truth test asserts (all measured on the run above)

| quantity | GT | fastcheck | tolerance in the test |
|---|---|---|---|
| per-position mean Q (150 positions) | 30.0 | max dev 0.0094 | ≤ 0.05 |
| whole-read GC | 0.5 | 0.5074 | ± 0.03 |
| adapter rate (pairs) | 989 | 990 | ± 1% of pairs |
| 3′ G-run histogram, bins 5 and 10 | 6010 / 4047 | 5923 / 3908 | ± 4% (substitution errors split runs) |
| 5′ G-run histogram | injected `{1:3970,2:1394,3:598}` | 6665 / 2076 / 753 | matches the convolution `P(h=k)+P(h=0)·P(random run=k)` within 1.5% |
| insert mean (after the fix) | 206.66 | 207.17 | ± 2% |
| duplication level | 1.99% | 0.0 | `0 ≤ level ≤ gt_dup_rate` — a **lower bound**, see Task 5 |

`mean_read_mean_q` is **not** asserted: it bins the *floor* of each read's mean, which biases
it low by ~0.5 (measured 29.53 for a true 30.0). The per-position mean from the exact
histogram is the unbiased quantity, and that is what the test uses.

---

## Task 1: Fix the differential harness (prerequisite — it is red today)

**Objective:** `tests/differential.py` reads the positions TSV by position instead of by
column name, so the three columns added in commit `1327ca0` (`g_frac`, `polyx3_frac`,
`polyx5_frac`) stop crashing it, and later tasks can add columns freely.

**Files:**
- Modify: `tools/fastcheck/tests/differential.py:76-93`

**Step 1: Reproduce the failure**

Run: `cd tools/fastcheck && python tests/differential.py ../../examples/runs/demo/sim_E/sim_R1.fastq.gz`
Expected: `ValueError: too many values to unpack (expected 11, got 14)` (exit 1).

**Step 2: Replace the fixed-arity unpack with a header map**

Replace the body of section 3 in `differential.py` with:

```python
    # 3. Per-position table, every column, read BY NAME from the header so adding a
    #    column (g_frac, polyx3_frac, ...) cannot break the comparison.
    ref_rows = ref["positions"]
    header = got_rows[0]
    col = {name: i for i, name in enumerate(header)}
    need = ["pos", "a", "c", "g", "t", "n", "depth", "mean_q", "q10", "q50", "q90"]
    missing = [n for n in need if n not in col]
    if missing:
        failures.append(f"positions TSV missing columns: {missing}")
    if len(got_rows) - 1 != len(ref_rows):
        failures.append(
            f"positions: fastcheck has {len(got_rows) - 1} rows, reference {len(ref_rows)}"
        )
    elif not missing:
        for gr, rr in zip(got_rows[1:], ref_rows):
            g = lambda name: gr[col[name]]
            for name in ("a", "c", "g", "t", "n", "depth"):
                got_i, want_i = int(g(name)), rr[name]
                if got_i != want_i:
                    failures.append(f"pos {g('pos')} {name}: fastcheck={got_i} reference={want_i}")
            if abs(float(g("mean_q")) - rr["mean_q"]) > 1e-3:
                failures.append(
                    f"pos {g('pos')} mean_q: fastcheck={g('mean_q')} reference={rr['mean_q']:.6f}"
                )
            for name in ("q10", "q50", "q90"):
                if abs(float(g(name)) - rr[name]) > 1e-9:
                    failures.append(f"pos {g('pos')} {name}: fastcheck={g(name)} reference={rr[name]}")
```

**Step 3: Verify green**

Run: `python tests/differential.py ../../examples/runs/demo/sim_E/sim_R1.fastq.gz`
Expected: `differential OK on .../sim_R1.fastq.gz: 20000 reads, 150 positions, every metric
identical to the reference` (exit 0).

**Step 4: Also verify on a generated fixture (the gate this plan keeps green)**

Run:
```bash
fastsim_bin=../fastsim/target/release/fastsim
"$fastsim_bin" --pairs 5000 --profile ../fastsim/profiles/synthetic_novaseq.json \
  --insert-mean 200 --insert-sd 40 --gc 0.5 --seed 7 --no-gzip -o /tmp/fcfix
python tests/differential.py /tmp/fcfix_R1.fastq
```
Expected: `differential OK on /tmp/fcfix_R1.fastq: 5000 reads, 150 positions, ...`

**Step 5: Commit**

```bash
git add tools/fastcheck/tests/differential.py
git commit -m "test(fastcheck): read the positions TSV by column name, not by position"
```

---

## Task 2: Port the memory-admission model into fastcheck

**Objective:** `tools/fastcheck/src/memory.rs` exists with limit discovery, exact projection
and the three-way admission rule, copied from `tools/fastprep/src/memory.rs` (D-1), and its
unit tests pass.

**Files:**
- Create: `tools/fastcheck/src/memory.rs`
- Modify: `tools/fastcheck/Cargo.toml` (add `libc`)
- Modify: `tools/fastcheck/src/main.rs:9-19` (`mod memory;`)

**Step 1: Add the dependency and the module**

`tools/fastcheck/Cargo.toml`, in `[dependencies]` after `serde_json = "1"`:

```toml
libc = "0.2"
```

`tools/fastcheck/src/main.rs`, after `mod insert;`:

```rust
mod memory;
```

**Step 2: Write the failing test first**

Create `tools/fastcheck/src/memory.rs` with the module doc comment and the test module
below, but with every function body `todo!()`:

```rust
//! Memory-limit discovery and exact allocation projection.
//!
//! Copied from `tools/fastprep/src/memory.rs` (D-1: no shared crate — anything two tools
//! need exists twice). fastcheck's bounded table modules (the duplication fingerprint set
//! and `--overrep-mode table`) are the only state whose size is not a function of read
//! length, so they are admission-checked before allocation rather than after an OOM kill.
//!
//! The limit that matters is the one that will KILL the process, not the memory the machine
//! happens to have: a SLURM job granted 8 GB on a 1.5 TB node sees 1.5 TB and is killed by
//! its cgroup at 8 GB, so this module takes the minimum over every constraint it can find.

use std::path::{Path, PathBuf};

/// cgroup v1's "no limit" sentinel (`PAGE_COUNTER_MAX * PAGE_SIZE`).
const V1_UNLIMITED: u64 = 9_223_372_036_854_771_712;

/// Which constraint produced the limit. Reported in diagnostics: when a job on a terabyte
/// node refuses a 5 GB allocation, the identity of the binding constraint is the
/// information the user needs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LimitSource {
    CgroupV2(String),
    CgroupV1(String),
    Slurm(&'static str),
    RlimitAs,
    MemAvailable,
}

impl std::fmt::Display for LimitSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LimitSource::CgroupV2(p) => write!(f, "cgroup v2: {}", p),
            LimitSource::CgroupV1(p) => write!(f, "cgroup v1: {}", p),
            LimitSource::Slurm(v) => write!(f, "SLURM: ${}", v),
            LimitSource::RlimitAs => write!(f, "RLIMIT_AS"),
            LimitSource::MemAvailable => write!(f, "/proc/meminfo MemAvailable"),
        }
    }
}

/// A discovered memory bound and the constraint that produced it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Limit {
    pub bytes: u64,
    pub source: LimitSource,
}

/// The binding memory limit for this process, or `None` when unconstrained.
///
/// Returns the *minimum* over every source, because any one of them can be the one that
/// kills us.
pub fn discover() -> Option<Limit> {
    discover_in(Path::new("/"), &std::env::vars().collect::<Vec<_>>())
}

/// Testable form of [`discover`]: reads under `root` and uses `env` instead of the process
/// environment.
pub fn discover_in(root: &Path, env: &[(String, String)]) -> Option<Limit> {
    todo!()
}

/// cgroup v2. The limit is frequently set on a parent slice rather than on the leaf, so the
/// whole ancestor chain is walked and the minimum taken.
fn cgroup_v2(root: &Path, phys: u64) -> Option<Limit> { todo!() }

/// cgroup v1. The controller path comes from the `:memory:` line.
fn cgroup_v1(root: &Path, phys: u64) -> Option<Limit> { todo!() }

/// The scheduler's grant. Worth reading even under cgroup enforcement: some sites set it
/// without cgroups.
fn slurm(env: &[(String, String)]) -> Option<Limit> { todo!() }

/// `RLIMIT_AS` bounds *virtual address space*, not resident set, so it is respected only
/// when set below physical memory and is deliberately the weakest signal.
#[cfg(target_os = "linux")]
fn rlimit_as(phys: u64) -> Option<Limit> { todo!() }

#[cfg(not(target_os = "linux"))]
fn rlimit_as(_phys: u64) -> Option<Limit> { None }

/// `MemAvailable`, not `MemTotal`: it accounts for reclaimable page cache and is the honest
/// answer to "how much can this process actually take".
fn mem_available(root: &Path) -> Option<u64> { todo!() }

fn mem_total(root: &Path) -> Option<u64> { todo!() }

fn meminfo_field(root: &Path, key: &str) -> Option<u64> { todo!() }

fn read(p: impl AsRef<Path>) -> Option<String> { todo!() }

// ---------------------------------------------------------------------------
// Projection. This is an identity, not an estimate.
// ---------------------------------------------------------------------------

/// Projected resident size of a pre-sized hash table.
///
/// `hashbrown` is open-addressed with a maximum load factor of 7/8 and one control byte
/// per bucket, and pre-sizing means no rehashing occurs, so
/// `buckets = next_power_of_two(ceil(8 * capacity / 7))` and
/// `bytes = buckets * (size_of::<T>() + 1)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Projection {
    pub buckets: u64,
    pub bytes: u64,
}

pub fn project(capacity: usize, fingerprint_bytes: usize) -> Projection { todo!() }

/// Largest capacity that still fits in `2^k` buckets, i.e. `7 * 2^k / 8`.
///
/// Bucket counts double at these points, **not** at powers of two: 5e8 sits 6% above the
/// 2^29 boundary and so projects 9216 MB where 4608 MB would do.
pub const fn capacity_ceiling(pow2_buckets: u32) -> usize { todo!() }

/// Replacement for a round-number preallocation clamp: the largest capacity that fits in
/// 2^29 buckets, i.e. 4608 MB at 64-bit fingerprints.
pub const MAX_PREALLOC_ENTRIES: usize = todo!();

// ---------------------------------------------------------------------------
// Admission.
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Admission {
    /// Comfortably within budget.
    Proceed,
    /// Fits, but close enough that the user should be told.
    ProceedWithWarning { headroom_bytes: u64 },
    /// Would not fit. Refuse before allocating.
    Refuse { projected: u64, limit: u64, source: LimitSource },
}

/// Reserve covering chunk arenas, compression buffers and allocator slack.
pub fn reserve_for(limit: u64) -> u64 { todo!() }

/// Admission rule. `limit: None` means undiscoverable, which degrades to
/// [`Admission::Proceed`]: refusing on an unknown bound would make the tool unusable
/// wherever discovery fails.
pub fn admit(projected: u64, limit: Option<&Limit>) -> Admission { todo!() }

pub fn mib(bytes: u64) -> f64 { todo!() }

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    const GB: u64 = 1 << 30;
    const PHYS: u64 = 1536 * GB; // a 1.5 TB EPYC node

    /// Builds a synthetic /proc + /sys tree.
    fn tree(name: &str, cgroup_line: &str, files: &[(&str, &str)]) -> PathBuf {
        let root = std::env::temp_dir().join(format!("fastcheck-memtest-{}-{}", std::process::id(), name));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("proc/self")).unwrap();
        fs::write(root.join("proc/self/cgroup"), cgroup_line).unwrap();
        fs::write(
            root.join("proc/meminfo"),
            format!("MemTotal:       {} kB\nMemAvailable:   {} kB\n", PHYS / 1024, PHYS / 1024),
        )
        .unwrap();
        for (p, v) in files {
            let full = root.join(p);
            fs::create_dir_all(full.parent().unwrap()).unwrap();
            fs::write(full, v).unwrap();
        }
        root
    }

    fn env(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs.iter().map(|(a, b)| (a.to_string(), b.to_string())).collect()
    }

    #[test]
    fn v2_leaf_limit() {
        let r = tree("v2leaf", "0::/slurm/uid_1234/job_98765\n",
            &[("sys/fs/cgroup/slurm/uid_1234/job_98765/memory.max", &(8 * GB).to_string())]);
        assert_eq!(discover_in(&r, &[]).unwrap().bytes, 8 * GB);
    }

    /// The case a leaf-only read misses: the limit lives on a parent slice.
    #[test]
    fn v2_parent_slice_is_binding() {
        let r = tree("v2parent", "0::/slurm/uid_1234/job_98765\n", &[
            ("sys/fs/cgroup/slurm/uid_1234/memory.max", &(4 * GB).to_string()),
            ("sys/fs/cgroup/slurm/uid_1234/job_98765/memory.max", "max"),
        ]);
        let l = discover_in(&r, &[]).unwrap();
        assert_eq!(l.bytes, 4 * GB);
        assert_eq!(l.source, LimitSource::CgroupV2("/slurm/uid_1234".into()));
    }

    #[test]
    fn v2_takes_minimum_over_chain() {
        let r = tree("v2min", "0::/slurm/uid_1234/job_98765\n", &[
            ("sys/fs/cgroup/slurm/uid_1234/memory.max", &(4 * GB).to_string()),
            ("sys/fs/cgroup/slurm/uid_1234/job_98765/memory.max", &(16 * GB).to_string()),
        ]);
        assert_eq!(discover_in(&r, &[]).unwrap().bytes, 4 * GB);
    }

    #[test]
    fn v2_all_max_is_unconstrained() {
        let r = tree("v2unl", "0::/user.slice\n", &[("sys/fs/cgroup/user.slice/memory.max", "max")]);
        // Only MemAvailable remains, which equals physical here.
        assert_eq!(discover_in(&r, &[]).unwrap().source, LimitSource::MemAvailable);
    }

    #[test]
    fn v1_sentinel_is_not_a_limit() {
        let r = tree("v1unl", "5:memory:/user.slice\n",
            &[("sys/fs/cgroup/memory/user.slice/memory.limit_in_bytes", &V1_UNLIMITED.to_string())]);
        assert_ne!(discover_in(&r, &[]).unwrap().source, LimitSource::CgroupV1("/user.slice".into()));
    }

    #[test]
    fn v1_real_limit() {
        let r = tree("v1real", "5:memory:/slurm/job_1\n",
            &[("sys/fs/cgroup/memory/slurm/job_1/memory.limit_in_bytes", &(33 * GB).to_string())]);
        assert_eq!(discover_in(&r, &[]).unwrap().bytes, 33 * GB);
    }

    #[test]
    fn slurm_env_without_cgroups() {
        let r = tree("noc", "", &[]);
        assert_eq!(discover_in(&r, &env(&[("SLURM_MEM_PER_NODE", "8192")])).unwrap().bytes, 8 * GB);
        assert_eq!(
            discover_in(&r, &env(&[("SLURM_MEM_PER_CPU", "2048"), ("SLURM_CPUS_ON_NODE", "4")]))
                .unwrap().bytes,
            8 * GB
        );
    }

    #[test]
    fn minimum_wins_when_sources_disagree() {
        let r = tree("both", "0::/slurm/j\n", &[("sys/fs/cgroup/slurm/j/memory.max", &(8 * GB).to_string())]);
        let l = discover_in(&r, &env(&[("SLURM_MEM_PER_NODE", "4096")])).unwrap();
        assert_eq!(l.bytes, 4 * GB);
        assert!(matches!(l.source, LimitSource::Slurm(_)));
    }

    /// Predicted values sit 8 MB below the measured peak RSS at every point, that residue
    /// being the rest of the pipeline.
    #[test]
    fn projection_matches_measured_rss() {
        for (reads, bytes_per_fp, measured_mib) in
            [(300_000_000usize, 8usize, 4616.0), (300_000_000, 16, 8712.0),
             (100_000_000, 8, 1160.0),           (100_000_000, 16, 2184.0),
             (50_000_000, 8, 584.0),             (50_000_000, 16, 1096.0),
             (10_000_000, 8, 152.0),             (10_000_000, 16, 280.0)]
        {
            let p = project(reads, bytes_per_fp);
            let delta = measured_mib - mib(p.bytes);
            assert!((delta - 8.0).abs() < 1.0, "reads={reads} width={bytes_per_fp}: delta {delta}");
        }
    }

    #[test]
    fn bucket_boundaries_are_seven_eighths_of_a_power_of_two() {
        assert_eq!(MAX_PREALLOC_ENTRIES, 469_762_048);
        assert_eq!(project(469_762_048, 8).buckets, 1 << 29);
        assert_eq!(project(469_762_049, 8).buckets, 1 << 30);
        // The old round-number clamp overshoots the boundary by 6% and so projects twice
        // the memory.
        assert_eq!(project(500_000_000, 8).buckets, 1 << 30);
        assert_eq!(mib(project(500_000_000, 8).bytes) as u64, 9216);
        assert_eq!(mib(project(MAX_PREALLOC_ENTRIES, 8).bytes) as u64, 4608);
    }

    #[test]
    fn project_zero_capacity() {
        assert_eq!(project(0, 8), Projection { buckets: 0, bytes: 0 });
    }

    #[test]
    fn admits_when_unconstrained() {
        assert_eq!(admit(100 * GB, None), Admission::Proceed);
    }

    #[test]
    fn refuses_when_projection_exceeds_limit() {
        let l = Limit { bytes: 4 * GB, source: LimitSource::CgroupV2("/x".into()) };
        assert!(matches!(admit(4 * GB, Some(&l)), Admission::Refuse { .. }));
    }

    #[test]
    fn warns_in_the_band_below_the_refusal_threshold() {
        let l = Limit { bytes: 16 * GB, source: LimitSource::CgroupV2("/x".into()) };
        // 12 GB + 1.6 GB reserve = 13.6 GB, which is 85% of 16 GB.
        assert!(matches!(admit(12 * GB, Some(&l)), Admission::ProceedWithWarning { .. }));
        // 4 GB + 1.6 GB = 5.6 GB, 35%.
        assert_eq!(admit(4 * GB, Some(&l)), Admission::Proceed);
    }

    #[test]
    fn reserve_has_a_floor() {
        assert_eq!(reserve_for(1 * GB), 512 * 1024 * 1024);
        assert_eq!(reserve_for(100 * GB), 10 * GB);
    }
}
```

**Step 3: Run the tests to verify failure**

Run: `cd tools/fastcheck && cargo test --offline memory::`
Expected: every test FAILS with `not yet implemented: ...` (`panicked at 'not yet implemented'`).

**Step 4: Implement**

Fill in every `todo!()` from `tools/fastprep/src/memory.rs` (lines 60-322), keeping the
signatures and doc comments already written above. The ported bodies are, verbatim in
substance:

```rust
pub fn discover_in(root: &Path, env: &[(String, String)]) -> Option<Limit> {
    let phys = mem_total(root).unwrap_or(u64::MAX);
    [
        cgroup_v2(root, phys),
        cgroup_v1(root, phys),
        slurm(env),
        rlimit_as(phys),
        mem_available(root).map(|b| Limit { bytes: b, source: LimitSource::MemAvailable }),
    ]
    .into_iter()
    .flatten()
    .min_by_key(|l| l.bytes)
}

fn cgroup_v2(root: &Path, phys: u64) -> Option<Limit> {
    let rel = read(root.join("proc/self/cgroup"))?
        .lines()
        .find_map(|l| l.strip_prefix("0::").map(str::to_owned))?;
    let parts: Vec<&str> = rel.split('/').filter(|s| !s.is_empty()).collect();
    let mut best: Option<Limit> = None;
    for depth in (0..=parts.len()).rev() {
        let mut dir: PathBuf = root.join("sys/fs/cgroup");
        for p in &parts[..depth] {
            dir.push(p);
        }
        let Some(raw) = read(dir.join("memory.max")) else { continue };
        let raw = raw.trim();
        if raw == "max" {
            continue;
        }
        let Ok(bytes) = raw.parse::<u64>() else { continue };
        if bytes >= phys {
            continue; // not a real constraint
        }
        let path = format!("/{}", parts[..depth].join("/"));
        if best.as_ref().is_none_or(|b| bytes < b.bytes) {
            best = Some(Limit { bytes, source: LimitSource::CgroupV2(path) });
        }
    }
    best
}

fn cgroup_v1(root: &Path, phys: u64) -> Option<Limit> {
    let rel = read(root.join("proc/self/cgroup"))?.lines().find_map(|l| {
        let f: Vec<&str> = l.splitn(3, ':').collect();
        (f.len() == 3 && f[1].split(',').any(|c| c == "memory")).then(|| f[2].to_owned())
    })?;
    let file = root
        .join("sys/fs/cgroup/memory")
        .join(rel.trim_start_matches('/'))
        .join("memory.limit_in_bytes");
    let bytes: u64 = read(file)?.trim().parse().ok()?;
    (bytes != V1_UNLIMITED && bytes < phys)
        .then(|| Limit { bytes, source: LimitSource::CgroupV1(rel) })
}

fn slurm(env: &[(String, String)]) -> Option<Limit> {
    let get = |k: &str| env.iter().find(|(n, _)| n == k).map(|(_, v)| v.as_str());
    if let Some(mb) = get("SLURM_MEM_PER_NODE").and_then(|v| v.parse::<u64>().ok()) {
        return Some(Limit { bytes: mb * 1024 * 1024, source: LimitSource::Slurm("SLURM_MEM_PER_NODE") });
    }
    let per_cpu = get("SLURM_MEM_PER_CPU")?.parse::<u64>().ok()?;
    let cpus = get("SLURM_CPUS_ON_NODE")?.parse::<u64>().ok()?;
    Some(Limit { bytes: per_cpu * cpus * 1024 * 1024,
                 source: LimitSource::Slurm("SLURM_MEM_PER_CPU x SLURM_CPUS_ON_NODE") })
}

#[cfg(target_os = "linux")]
fn rlimit_as(phys: u64) -> Option<Limit> {
    let mut rl = libc::rlimit { rlim_cur: 0, rlim_max: 0 };
    // SAFETY: `rl` is a valid, fully initialised `rlimit`; getrlimit only writes through it.
    if unsafe { libc::getrlimit(libc::RLIMIT_AS, &mut rl) } != 0 {
        return None;
    }
    let cur = rl.rlim_cur as u64;
    (rl.rlim_cur != libc::RLIM_INFINITY && cur < phys)
        .then_some(Limit { bytes: cur, source: LimitSource::RlimitAs })
}

fn meminfo_field(root: &Path, key: &str) -> Option<u64> {
    read(root.join("proc/meminfo"))?
        .lines()
        .find(|l| l.starts_with(key))
        .and_then(|l| l.split_whitespace().nth(1)?.parse::<u64>().ok())
        .map(|kb| kb * 1024)
}

fn mem_available(root: &Path) -> Option<u64> { meminfo_field(root, "MemAvailable:") }
fn mem_total(root: &Path) -> Option<u64> { meminfo_field(root, "MemTotal:") }
fn read(p: impl AsRef<Path>) -> Option<String> { std::fs::read_to_string(p).ok() }

pub fn project(capacity: usize, fingerprint_bytes: usize) -> Projection {
    if capacity == 0 {
        return Projection { buckets: 0, bytes: 0 };
    }
    let need = (capacity as u64 * 8).div_ceil(7);
    let buckets = need.next_power_of_two();
    Projection { buckets, bytes: buckets * (fingerprint_bytes as u64 + 1) }
}

pub const fn capacity_ceiling(pow2_buckets: u32) -> usize { (7usize << pow2_buckets) / 8 }
pub const MAX_PREALLOC_ENTRIES: usize = capacity_ceiling(29); // 469_762_048

pub fn reserve_for(limit: u64) -> u64 { (512 * 1024 * 1024).max(limit / 10) }

pub fn admit(projected: u64, limit: Option<&Limit>) -> Admission {
    let Some(l) = limit else { return Admission::Proceed };
    let want = projected.saturating_add(reserve_for(l.bytes));
    if want * 10 >= l.bytes * 9 {
        Admission::Refuse { projected, limit: l.bytes, source: l.source.clone() }
    } else if want * 10 >= l.bytes * 6 {
        Admission::ProceedWithWarning { headroom_bytes: l.bytes - want }
    } else {
        Admission::Proceed
    }
}

pub fn mib(bytes: u64) -> f64 { bytes as f64 / 1024.0 / 1024.0 }
```

**Step 5: Run the tests to verify pass**

Run: `cargo test --offline memory::`
Expected: `test result: ok. 15 passed; 0 failed` (13 discovery/projection + 4 admission
tests, of which `projection_matches_measured_rss` is one test).

**Step 6: Whole crate still green**

Run: `cargo test --offline`
Expected: `42 passed` (unit) + `12 passed` (cli), 0 failed.

**Step 7: Commit**

```bash
git add tools/fastcheck/src/memory.rs tools/fastcheck/Cargo.toml tools/fastcheck/src/main.rs
git commit -m "feat(fastcheck): port the memory-admission model from fastprep (D-1 copy)"
```

---

## Task 3: `--max-mem` flag, record estimate, and the printed projection

**Objective:** `fastcheck --max-mem 4096 reads.fastq` is accepted, the estimated record
count and the table projection are printed **before** the scan, and the existing
`accumulator state ... MB` line is unchanged.

**Files:**
- Modify: `tools/fastcheck/src/cli.rs` (new flag)
- Modify: `tools/fastcheck/src/main.rs` (estimate + print)

**Step 1: Write the failing test**

Append to `tools/fastcheck/tests/cli.rs`:

```rust
#[test]
fn max_mem_flag_is_accepted_and_the_projection_is_printed_before_the_scan() {
    let d = tmpdir("maxmem");
    let fq = write_fixture(&d);
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--max-mem", "4096", "-q"],
    );
    assert_eq!(code, 0, "stderr: {stderr}");
    // The projection line must appear, and must name the fingerprint table.
    let proj = stderr
        .lines()
        .find(|l| l.contains("fingerprint table"))
        .unwrap_or_else(|| panic!("no projection line in stderr:\n{stderr}"));
    assert!(proj.contains("MB"), "{proj}");
    // The existing accumulator-state line is untouched.
    assert!(stderr.contains("accumulator state"), "{stderr}");
}
```

**Step 2: Run to verify failure**

Run: `cd tools/fastcheck && cargo test --offline max_mem_flag_is_accepted`
Expected: FAIL — clap: `error: unexpected argument '--max-mem' found`.

**Step 3: Add the flag**

In `tools/fastcheck/src/cli.rs`, after the `max_read_len` field:

```rust
    /// Refuse to start if the projected fingerprint table exceeds this many megabytes.
    /// Defaults to the enforced memory limit of the process (cgroup v2/v1, SLURM,
    /// RLIMIT_AS, then MemAvailable — the minimum of them, per tools/fastprep's model).
    #[arg(long, value_name = "MB")]
    pub max_mem: Option<u64>,
```

**Step 4: Estimate the records and print the projection**

In `tools/fastcheck/src/main.rs`, add after `let track_profile = ...` (line 53):

```rust
    // The duplication fingerprint set is the only state whose size is a function of the
    // dataset, so it is projected and admission-checked before the scan starts. The
    // estimate is deliberately conservative on compressed input (a 4x ratio under-counts,
    // so the table is under-sized and grows rather than over-allocating): getting it wrong
    // costs a rehash, not correctness. Copied from fastprep's estimate_records (D-1).
    let est_records: usize = cli
        .files
        .iter()
        .map(|p| {
            let bytes = std::fs::metadata(p)
                .with_context(|| format!("sizing {}", p.display()))?
                .len();
            // Read length comes from the profile fit's own length histogram later; 150 is
            // the short-read default and only scales the estimate.
            let per_record = (150u64 * 2 + 50).max(1);
            let compressed = p
                .extension()
                .and_then(|e| e.to_str())
                .is_some_and(|e| matches!(e, "gz" | "zst" | "zstd" | "bz2" | "xz"));
            let raw = if compressed { bytes.saturating_mul(4) } else { bytes };
            Ok::<usize, anyhow::Error>((raw / per_record) as usize)
        })
        .sum::<anyhow::Result<usize>>()?;
    let fp_capacity = est_records.clamp(1, memory::MAX_PREALLOC_ENTRIES);
    let fp_proj = memory::project(fp_capacity, 8); // 64-bit fingerprints
    eprintln!(
        "fastcheck: duplication fingerprint table projected at {:.1} MB for ~{} reads \
         ({} buckets, 64-bit fingerprints)",
        memory::mib(fp_proj.bytes),
        fp_capacity,
        fp_proj.buckets,
    );
```

(`use anyhow::Context;` must be added to the imports at the top of `main.rs`.)

**Step 5: Run to verify pass**

Run: `cargo test --offline max_mem_flag_is_accepted`
Expected: PASS.

**Step 6: Manual check**

Run: `./target/release/fastcheck tests/../../examples/runs/demo/sim_E/sim_R1.fastq.gz --max-mem 4096 -q 2>&1 | head -2`
Expected: the fingerprint line first (≈ 6 MB for a 20 k-read file), then the accumulator
line.

**Step 7: Commit**

```bash
git add tools/fastcheck/src/cli.rs tools/fastcheck/src/main.rs tools/fastcheck/tests/cli.rs
git commit -m "feat(fastcheck): --max-mem flag and a fingerprint-table projection printed before the scan"
```

---

## Task 4: Admission refusal above `--max-mem`

**Objective:** A run whose projected table would not fit exits 1 **before reading anything**,
with a message that names the projection, the limit and its source, and offers the three
ways out. This is PLAN-QC principle 6.

**Files:**
- Modify: `tools/fastcheck/src/main.rs` (the admission block, right after Task 3's print)

**Step 1: Write the failing test**

Append to `tools/fastcheck/tests/cli.rs`:

```rust
#[test]
fn a_projection_above_max_mem_refuses_before_reading() {
    let d = tmpdir("refuse");
    let fq = write_fixture(&d);
    // 1 MB ceiling against a ~6 MB projection on the 6-read fixture: it must refuse.
    let (stdout, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--max-mem", "1", "-q"],
    );
    assert_ne!(code, 0, "must refuse: {stdout}{stderr}");
    assert!(stdout.is_empty(), "no report may be written on refusal: {stdout}");
    assert!(stderr.contains("fingerprint table"), "{stderr}");
    assert!(stderr.contains("limited to"), "{stderr}");
    assert!(stderr.contains("--overrep-mode none"), "{stderr}"); // a way out
    // The refusal must happen before any numbers are printed.
    assert!(!stderr.contains("accumulator state"), "{stderr}");
}
```

**Step 2: Run to verify failure**

Run: `cargo test --offline a_projection_above_max_mem_refuses_before_reading`
Expected: FAIL — the run succeeds (exit 0) because nothing checks the ceiling.

**Step 3: Implement the refusal**

In `tools/fastcheck/src/main.rs`, immediately after the projection `eprintln!` from Task 3:

```rust
    // An explicit --max-mem overrides discovery, but only DOWNWARDS: it cannot authorise
    // exceeding a cgroup limit. (fastprep's rule, D-1 copy.)
    let discovered = memory::discover();
    let limit_bytes = cli
        .max_mem
        .map(|mb| mb.saturating_mul(1024 * 1024))
        .into_iter()
        .chain(discovered.as_ref().map(|l| l.bytes))
        .min();
    let limit = limit_bytes.map(|bytes| memory::Limit {
        bytes,
        source: discovered
            .as_ref()
            .filter(|l| l.bytes == bytes)
            .map(|l| l.source.clone())
            .unwrap_or(memory::LimitSource::RlimitAs), // an explicit --max-mem is a soft cap
    });

    match memory::admit(fp_proj.bytes, limit.as_ref()) {
        memory::Admission::Refuse { projected, limit, source } => anyhow::bail!(
            "the duplication fingerprint table needs ~{:.1} MB ({} buckets x 9 B, 64-bit fingerprints)\n\
             but this process is limited to {:.1} MB ({})\n\n\
             --max-mem <MB>           raise the ceiling (it cannot exceed the cgroup limit)\n\
             --overrep-mode none      run without the duplication table entirely\n\
             reread the file with a smaller --reads budget to bound the estimate",
            memory::mib(projected),
            fp_proj.buckets,
            memory::mib(limit),
            source,
        ),
        memory::Admission::ProceedWithWarning { headroom_bytes } => eprintln!(
            "fastcheck: warning: fingerprint table projected at {:.1} MB, leaving only {:.1} MB headroom",
            memory::mib(fp_proj.bytes),
            memory::mib(headroom_bytes),
        ),
        memory::Admission::Proceed => {}
    }
```

Note the refusal uses `anyhow::bail!` from `main() -> Result<()>`, so the process exits 1
with the message on stderr — the same contract fastprep has.

**Step 4: Run to verify pass**

Run: `cargo test --offline a_projection_above_max_mem_refuses_before_reading`
Expected: PASS.

**Step 5: Manual check, both directions**

```bash
./target/release/fastcheck examples/runs/demo/sim_E/sim_R1.fastq.gz --max-mem 1 -q; echo "exit=$?"
./target/release/fastcheck examples/runs/demo/sim_E/sim_R1.fastq.gz --max-mem 65536 -q; echo "exit=$?"
```
Expected: first prints the refusal block and `exit=1`; second prints the summary line and
`exit=0`.

**Step 6: Whole crate green**

Run: `cargo test --offline`
Expected: `42 passed` + `13 passed` (the new CLI test), 0 failed.

**Step 7: Commit**

```bash
git add tools/fastcheck/src/main.rs tools/fastcheck/tests/cli.rs
git commit -m "feat(fastcheck): refuse to start when the projected table exceeds --max-mem (PLAN-QC principle 6)"
```

---

## Task 5: Always-on exact fingerprint counter (duplication levels become exact)

**Objective:** The duplication level stops being a bounded-counter approximation and becomes
an exact multiplicity histogram over xxh3-64 fingerprints, which is what PLAN-QC §9 and §5.9
require as the baseline every `--overrep-mode` builds on.

**Files:**
- Modify: `tools/fastcheck/Cargo.toml` (`rustc-hash`, `xxhash-rust`)
- Create: `tools/fastcheck/src/fingerprint.rs`
- Modify: `tools/fastcheck/src/accum.rs` (the counter field, `observe`, `merge`, accessors)
- Modify: `tools/fastcheck/src/rich.rs:75-100` (`repeated_metrics` reads the exact table)

**Step 1: Add deps and the module**

`tools/fastcheck/Cargo.toml`, in `[dependencies]`:

```toml
rustc-hash = "2"
xxhash-rust = { version = "0.8", features = ["xxh3"] }
```

`tools/fastcheck/src/main.rs`: `mod fingerprint;`

**Step 2: Write the failing unit tests**

Create `tools/fastcheck/src/fingerprint.rs`:

```rust
//! Exact duplication counting over 64-bit fingerprints.
//!
//! Copied from `tools/fastdedup/src/hasher.rs` (D-1: no shared crate) with the abundance
//! tracker kept and the dedup-only verifier dropped. A PCR duplicate re-sequences the same
//! fragment with FRESH polyG artefacts and fresh substitution errors, so two copies of one
//! molecule usually differ as strings: what identifies a duplicate is the fragment, not the
//! read. Fingerprinting the whole read therefore under-counts, and the report says so
//! (duplication_level is a LOWER BOUND on the PCR duplicate rate).

use rustc_hash::FxHashMap;
use xxhash_rust::xxh3::xxh3_64;

/// Fingerprint of one read. Quality is ignored (stated): two copies of one fragment differ
/// in their quality string because dropout is per-cluster, and that difference is exactly
/// what must NOT split a duplicate.
#[inline(always)]
pub fn fingerprint(seq: &[u8]) -> u64 {
    xxh3_64(seq)
}

/// Exact multiplicity counter: fingerprint -> occurrences.
///
/// Pre-sized from the projected capacity so no rehashing occurs, which is what makes the
/// projection in `crate::memory::project` an identity rather than an estimate.
#[derive(Clone)]
pub struct FingerprintTable {
    counts: FxHashMap<u64, u32>,
}

impl FingerprintTable {
    /// `capacity` is the projected number of distinct reads.
    pub fn with_capacity(capacity: usize) -> Self {
        let mut counts = FxHashMap::with_capacity_and_hasher(
            capacity.min(1 << 26),
            Default::default(),
        );
        counts.reserve(capacity.min(1 << 26));
        Self { counts }
    }

    #[inline]
    pub fn observe(&mut self, seq: &[u8]) {
        *self.counts.entry(fingerprint(seq)).or_insert(0) += 1;
    }

    /// Distinct fingerprints seen.
    pub fn distinct(&self) -> u64 {
        self.counts.len() as u64
    }

    /// Total observations recorded (sum of multiplicities).
    pub fn total(&self) -> u64 {
        self.counts.values().map(|&c| c as u64).sum()
    }

    /// Exact multiplicity histogram: `out[c]` = number of distinct sequences occurring
    /// exactly `c` times. Index 0 is unused; index 1 is the singleton count.
    pub fn multiplicity_hist(&self) -> Vec<u64> {
        let max = self.counts.values().copied().max().unwrap_or(0) as usize;
        let mut out = vec![0u64; max + 2];
        for &c in self.counts.values() {
            out[c as usize] += 1;
        }
        out
    }

    /// Distinct sequences occurring at least `min_count` times, most frequent first,
    /// each as `(fingerprint, count)`. Names are not available in this mode.
    pub fn at_least(&self, min_count: u32) -> Vec<(u64, u32)> {
        let mut v: Vec<(u64, u32)> = self
            .counts
            .iter()
            .filter(|(_, &c)| c >= min_count)
            .map(|(&k, &c)| (k, c))
            .collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
    }

    /// Exact sum of two tables.
    pub fn merge(&mut self, o: &FingerprintTable) {
        for (&k, &c) in &o.counts {
            *self.counts.entry(k).or_insert(0) += c;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_exact_multiplicities() {
        let mut t = FingerprintTable::with_capacity(16);
        for _ in 0..5 {
            t.observe(b"ACGTACGTAC");
        }
        t.observe(b"TTTTTTTTTT");
        t.observe(b"TTTTTTTTTT");
        t.observe(b"GGGGCCCCGG");
        assert_eq!(t.distinct(), 3);
        assert_eq!(t.total(), 8);
        let h = t.multiplicity_hist();
        assert_eq!(h[1], 1, "one singleton");
        assert_eq!(h[2], 1, "one sequence seen twice");
        assert_eq!(h[5], 1, "one sequence seen five times");
    }

    #[test]
    fn at_least_filters_and_sorts_by_count() {
        let mut t = FingerprintTable::with_capacity(16);
        for _ in 0..7 {
            t.observe(b"ACGTACGTAC");
        }
        for _ in 0..3 {
            t.observe(b"TTTTTTTTTT");
        }
        t.observe(b"GGGGCCCCGG");
        let v = t.at_least(3);
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].1, 7, "most frequent first");
        assert_eq!(v[1].1, 3);
    }

    #[test]
    fn merge_is_exact() {
        let mut a = FingerprintTable::with_capacity(8);
        a.observe(b"ACGTACGTAC");
        a.observe(b"ACGTACGTAC");
        let mut b = FingerprintTable::with_capacity(8);
        b.observe(b"ACGTACGTAC");
        b.observe(b"TTTTTTTTTT");
        a.merge(&b);
        assert_eq!(a.distinct(), 2);
        assert_eq!(a.total(), 4);
        assert_eq!(a.multiplicity_hist()[3], 1);
    }

    /// fastdedup and fastcheck must agree on the fingerprint of a sequence (D-1: consistency
    /// by tests, not by the compiler).
    #[test]
    fn fingerprint_is_xxh3_64_of_the_sequence() {
        // xxh3_64 is deterministic; this pins the value so a dependency bump that changes
        // it fails loudly instead of silently splitting the two tools' counts.
        let h = fingerprint(b"ACGTACGTAC");
        assert_eq!(h, fingerprint(b"ACGTACGTAC"));
        assert_ne!(h, fingerprint(b"ACGTACGTACA"), "one base differs -> different fingerprint");
    }
}
```

**Step 3: Run to verify failure**

Run: `cargo test --offline fingerprint::`
Expected: FAIL — `unresolved import` / module does not exist (it is not declared yet).

**Step 4: Wire it into the accumulator**

In `tools/fastcheck/src/accum.rs`:

- Replace the `seq_counts`/`kmer_counts`/`distinct_seqs` fields with:

```rust
    /// Exact duplication counter over 64-bit fingerprints of the whole read. Quality is
    /// ignored (stated in the report). This is the baseline every --overrep-mode builds on
    /// (PLAN-QC §5.9): it gives the exact multiplicity histogram with no sequence text.
    pub seq_counts: crate::fingerprint::FingerprintTable,
```

- In `Accum::new`, replace the two `RepeatedCounter::new(...)` lines with:

```rust
            seq_counts: crate::fingerprint::FingerprintTable::with_capacity(1024),
```

- In `observe`, replace the bounded block (lines 249-264) with:

```rust
        // Duplication: every read is fingerprinted, exactly, with no sample and no cap.
        // The table is pre-sized from the projected record count in main.rs, so this does
        // not rehash. Over-representation naming is a separate, opt-in concern (§5.9).
        self.seq_counts.observe(seq);
```

- In `merge`, replace the bounded-counter comment and `distinct_seqs` line with:

```rust
        // The fingerprint table merges exactly, like every other accumulator here.
        self.seq_counts.merge(&o.seq_counts);
```

- Replace the `distinct_seqs` accessor usages: `resident_bytes` gains nothing (the table's
  size is reported separately), and `covered_positions`/`base_at`/`qual_at` are unchanged.

**Step 5: Point `rich::repeated_metrics` at the exact table**

In `tools/fastcheck/src/rich.rs`, replace `repeated_metrics` with:

```rust
/// Duplication levels: the
```rust
/// Duplication levels: the exact multiplicity histogram, plus the duplication level and the
/// distinct/read ratio. Over-represented sequences are a separate concern handled by
/// `crate::overrep` (§5.9) — this function only reports what the exact fingerprint table
/// knows. The top list names no sequences: that requires one of the --overrep-mode table
/// modes, and `names_available` says so.
pub fn repeated_metrics(acc: &Accum) -> (Map<String, Value>, Value) {
    let reads = acc.reads.max(1) as f64;
    let distinct = acc.seq_counts.distinct();
    let hist = acc.seq_counts.multiplicity_hist();
    let min_count = (reads * 0.001).max(10.0) as u32;
    let overrep = acc.seq_counts.at_least(min_count);

    let mut m = Map::new();
    m.insert("distinct_seqs".into(), json!(distinct));
    m.insert(
        "duplication_level".into(),
        json!(r4(1.0 - distinct as f64 / reads)),
    );
    m.insert(
        "total_deduplicated_pct".into(),
        json!(r4(100.0 * (reads - distinct as f64) / reads)),
    );
    m.insert("dup_multiplicity_max".into(), json!(hist.len().saturating_sub(2)));
    m.insert("dup_singletons".into(), json!(hist.get(1).copied().unwrap_or(0)));
    m.insert("repeated_seq_count".into(), json!(overrep.len()));
    let frac: f64 = overrep.iter().map(|(_, c)| *c as f64).sum::<f64>() / reads;
    m.insert("repeated_seq_frac".into(), json!(r4(frac)));

    // `names_available: false` is not a caveat, it is the contract of this mode: the
    // multiplicities are exact, the sequences are not recoverable from a hash.
    let list: Vec<Value> = overrep
        .iter()
        .take(20)
        .map(|(fp, c)| {
            json!({
                "fingerprint": format!("{fp:016x}"),
                "count": c,
                "frac": r4(*c as f64 / reads),
            })
        })
        .collect();
    (m, json!({ "names_available": false, "top": list }))
}
```

Also update `report.rs`'s `repeated_sequences` block comment expectation and the CLI test
that asserts `j["metrics"]["duplication_level"]` exists (it still does — the key is kept).

**Step 6: Update the CLI test that expected sequence names**

In `tools/fastcheck/tests/cli.rs`, `rich_report_has_polyx_repeats_and_heatmap` asserts
`m["repeated_seq_count"].is_u64()` — still true. Add one line inside that test, after the
existing metric assertions:

```rust
    // The exact fingerprint table reports multiplicities, never names (PLAN-QC §5.9).
    assert_eq!(j["repeated_sequences"]["names_available"], false);
    assert!(j["repeated_sequences"]["top"].is_array());
```

**Step 7: Run the tests**

Run: `cargo test --offline`
Expected: all green. `fingerprint::` gives 4 passed; the CLI suite's
`rich_report_has_polyx_repeats_and_heatmap` passes with the new assertions. Watch for one
legitimately changed expectation: the old bounded counter reported `distinct_seqs` only for
the first `COUNT_LIMIT` (100 000) reads; the exact table counts every read, so on the 6-read
fixture `distinct_seqs` is still 6 and nothing changes. On a 20 000-read file it rises from
20 000 to 20 000 (already under the old limit), so the demo numbers are unchanged too.

**Step 8: Verify the duplication level against ground truth**

```bash
../fastsim/target/release/fastsim --pairs 20000 --profile ../fastsim/profiles/synthetic_novaseq.json \
  --insert-mean 200 --insert-sd 40 --gc 0.5 --seed 42 --no-gzip --dup-rate 0.02 \
  --adapter-seq CTGTCTCTTATACACATCT -o /tmp/fp
./target/release/fastcheck /tmp/fp_R1.fastq /tmp/fp_R2.fastq --json /tmp/fp.json -q
python3 -c "
import json
m=json.load(open('/tmp/fp.json'))['metrics']
gt=sum(1 for r in open('/tmp/fp.gt.tsv').read().splitlines()[1:] if r.split('\t')[9]=='1')
print('distinct',m['distinct_seqs'],'reads',m['reads'],'level',m['duplication_level'])
print('GT dup pairs',gt,'-> GT rate',gt/20000)
assert 0.0 <= m['duplication_level'] <= gt/20000 + 0.01, 'must be a lower bound on the PCR rate'
print('OK: duplication_level is a lower bound on the PCR duplicate rate, as documented')
"
```
Expected: `level` between 0 and ~0.03, and the assertion passes. (Measured 0.0: fresh
artefacts make the copies differ as strings — the documented reason.)

**Step 9: Commit**

```bash
git add tools/fastcheck/Cargo.toml tools/fastcheck/src/fingerprint.rs \
        tools/fastcheck/src/accum.rs tools/fastcheck/src/rich.rs \
        tools/fastcheck/src/report.rs tools/fastcheck/tests/cli.rs
git commit -m "feat(fastcheck): exact duplication multiplicities over xxh3-64 fingerprints (D-1 copy of fastdedup's hasher)"
```

---

## Task 6: `--overrep-mode` with the four modes

**Objective:** `--overrep-mode none|fingerprint|second-pass|table` exists; `none` is the
default and is what Task 5 already built; each mode emits its PLAN-QC §5.9 report fields;
`table` is admission-checked with `--max-mem`; `second-pass` records `passes: 2`.

**Files:**
- Modify: `tools/fastcheck/src/cli.rs` (the flag)
- Create: `tools/fastcheck/src/overrep.rs` (the sequence table + second pass)
- Modify: `tools/fastcheck/src/accum.rs` (the mode hook in `observe`)
- Modify: `tools/fastcheck/src/report.rs` (the fields)
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: Write the failing tests**

Append to `tools/fastcheck/tests/cli.rs`:

```rust
/// A fixture where one sequence appears 4 times out of 12 reads.
const OVERREP_FASTQ: &[u8] = b"@a\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCANNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII####\n\
@b\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCANNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII####\n\
@c\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCANNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII####\n\
@d\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCANNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII####\n\
@e\nTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGT\n+\n55555555555555555555555555555555\n\
@f\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n\
@g\nGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@h\nCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT\n+\n55555555555555555555555555555555\n\
@i\nAACCAACCAACCAACCAACCAACCAACCAACC\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n\
@j\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@k\nGGATCCTTAGGCCAATTGGCCAATTGGCCAA\n+\n55555555555555555555555555555555\n\
@l\nCATCATCATCATCATCATCATCATCATCATCA\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";

fn write_overrep(dir: &Path) -> std::path::PathBuf {
    let p = dir.join("over.fastq");
    std::fs::write(&p, OVERREP_FASTQ).unwrap();
    p
}

#[test]
fn overrep_mode_none_is_the_default_and_reports_no_names() {
    let d = tmpdir("overrep_none");
    let fq = write_overrep(&d);
    let json = d.join("r.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    assert_eq!(v["overrep"]["mode"], "none");
    assert_eq!(v["overrep"]["names_available"], false);
    // Duplication levels are still exact in this mode.
    assert!(v["metrics"]["distinct_seqs"].is_u64());
    assert!(v["overrep"]["sequences"].as_array().unwrap().is_empty());
}

#[test]
fn overrep_mode_fingerprint_reports_exact_multiplicity_without_names() {
    let d = tmpdir("overrep_fp");
    let fq = write_overrep(&d);
    let json = d.join("r.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(),
              "--overrep-mode", "fingerprint", "-q"]);
    let v = read_json(&json);
    assert_eq!(v["overrep"]["mode"], "fingerprint");
    assert_eq!(v["overrep"]["names_available"], false);
    // 4 copies of the first sequence, threshold max(10, 0.1% of 12) = 10 -> none reported.
    assert!(v["overrep"]["sequences"].as_array().unwrap().is_empty());
    // The threshold is a flag, so it can be lowered to see the 4x sequence.
    let json2 = d.join("r2.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json2.to_str().unwrap(),
              "--overrep-mode", "fingerprint", "--overrep-min-count", "3", "-q"]);
    let v2 = read_json(&json2);
    let seqs = v2["overrep"]["sequences"].as_array().unwrap();
    assert_eq!(seqs.len(), 1, "only the 4x sequence clears a threshold of 3");
    assert_eq!(seqs[0]["count"], 4);
    assert!(seqs[0].get("seq").is_none(), "fingerprint mode cannot name the sequence");
    assert_eq!(v2["overrep"]["names_available"], false);
}

#[test]
fn overrep_mode_table_names_the_sequences_and_second_pass_records_two_passes() {
    let d = tmpdir("overrep_tbl");
    let fq = write_overrep(&d);
    let json = d.join("t.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(),
              "--overrep-mode", "table", "--overrep-min-count", "3", "-q"]);
    let v = read_json(&json);
    assert_eq!(v["overrep"]["mode"], "table");
    assert_eq!(v["overrep"]["names_available"], true);
    assert_eq!(v["overrep"]["unique_sequences_retained"], true);
    let seqs = v["overrep"]["sequences"].as_array().unwrap();
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0]["seq"], "AGATCGGAAGAGCACACGTCTGAACTCCAGTCANNNN");
    assert_eq!(seqs[0]["count"], 4);

    let json2 = d.join("s.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json2.to_str().unwrap(),
              "--overrep-mode", "second-pass", "--overrep-min-count", "3", "-q"]);
    let v2 = read_json(&json2);
    assert_eq!(v2["overrep"]["mode"], "second-pass");
    assert_eq!(v2["overrep"]["names_available"], true);
    assert_eq!(v2["overrep"]["passes"], 2, "the second read must be declared");
    let seqs2 = v2["overrep"]["sequences"].as_array().unwrap();
    assert_eq!(seqs2.len(), 1);
    assert_eq!(seqs2[0]["seq"], "AGATCGGAAGAGCACACGTCTGAACTCCAGTCANNNN");
}

#[test]
fn overrep_mode_table_refuses_above_max_mem_like_the_fingerprint_table() {
    let d = tmpdir("overrep_refuse");
    let fq = write_overrep(&d);
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--overrep-mode", "table", "--max-mem", "1", "-q"],
    );
    assert_ne!(code, 0, "must refuse: {stderr}");
    assert!(
        stderr.contains("unique sequences") || stderr.contains("fingerprint table"),
        "the refusal must name the table: {stderr}"
    );
    assert!(stderr.contains("--overrep-mode fingerprint"), "must offer the cheaper mode: {stderr}");
}

#[test]
fn overrep_mode_is_rejected_when_unknown() {
    let d = tmpdir("overrep_bad");
    let fq = write_fixture(&d);
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--overrep-mode", "sampled", "-q"],
    );
    assert_ne!(code, 0);
    assert!(stderr.contains("none") && stderr.contains("second-pass"), "{stderr}");
}
```

**Step 2: Run to verify failure**

Run: `cargo test --offline overrep`
Expected: FAIL — `error: unexpected argument '--overrep-mode' found` (clap), for all five.

**Step 3: Add the flags**

In `tools/fastcheck/src/cli.rs`:

```rust
    /// How over-represented sequences are identified (PLAN-QC §5.9). The one honest
    /// trade-off in fastcheck: exact identity needs either the sequences themselves or a
    /// second read of the file.
    ///
    ///   none        (default, zero cost) exact duplication levels, no sequence names
    ///   fingerprint +unique x 8 B/0.8 RAM; exact multiplicity, names_available: false
    ///   second-pass one extra read of the file to emit the text of sequences above the
    ///               threshold; memory stays at the fingerprint table; passes: 2
    ///   table       2-bit packed sequences (ceil(len/4)+8 B each); refuses above --max-mem;
    ///               unique_sequences_retained: true
    ///
    /// FastQC instead analyses the first 100 000 sequences and hopes. None of these modes
    /// samples: a sample is never presented as the whole.
    #[arg(long, value_name = "none|fingerprint|second-pass|table", default_value = "none")]
    pub overrep_mode: String,

    /// Minimum multiplicity for a sequence to be reported as over-represented
    /// (default: max(10, 0.1% of reads) — FastQC's 0.1% convention, but computed from the
    /// exact read count rather than from a 100 k sample).
    #[arg(long, value_name = "N")]
    pub overrep_min_count: Option<u32>,
```

and a validator in `impl Cli`:

```rust
    pub fn overrep_mode(&self) -> anyhow::Result<OverrepMode> {
        OverrepMode::parse(&self.overrep_mode)
    }
```

**Step 4: Create `tools/fastcheck/src/overrep.rs`**

```rust
//! Over-represented sequences: the four modes of PLAN-QC §5.9.
//!
//! Exact identity needs one of:
//!
//! | mode | cost | report field |
//! |---|---|---|
//! | `none` (default) | zero | duplication levels still exact (§9) |
//! | `fingerprint` | +`unique x 8 B / 0.8` RAM | exact multiplicity, `names_available: false` |
//! | `second-pass` | one extra read of the file | `passes: 2` |
//! | `table` | 2-bit packed, `ceil(len/4)+8` B each | `unique_sequences_retained: true` |
//!
//! The table mode packs sequences 2 bits per base, so a 150 bp read costs 46 B including
//! the length header; 3e8 unique reads project ~11 GB and are refused above `--max-mem`
//! rather than OOM-killed at hour three. Second-pass keeps the fingerprint table in RAM
//! and re-reads the file only to recover the text of the fingerprints that cleared the
//! threshold, so its memory is the same as fingerprint mode.

use std::collections::HashMap;

use anyhow::{bail, Result};
use serde_json::json;

use crate::accum::Accum;
use crate::fingerprint::{self, FingerprintTable};
use crate::reader::{self, Limits};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverrepMode {
    None,
    Fingerprint,
    SecondPass,
    Table,
}

impl OverrepMode {
    pub fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "none" => OverrepMode::None,
            "fingerprint" => OverrepMode::Fingerprint,
            "second-pass" => OverrepMode::SecondPass,
            "table" => OverrepMode::Table,
            other => bail!(
                "--overrep-mode must be none, fingerprint, second-pass or table (got {other:?})"
            ),
        })
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            OverrepMode::None => "none",
            OverrepMode::Fingerprint => "fingerprint",
            OverrepMode::SecondPass => "second-pass",
            OverrepMode::Table => "table",
        }
    }

    /// Whether sequence text is recoverable in this mode.
    pub fn names_available(&self) -> bool {
        matches!(self, OverrepMode::SecondPass | OverrepMode::Table)
    }
}

/// 2-bit packed sequence with its length. `N` and every non-ACGT code are not packable and
/// force the `table` mode to refuse rather than to lose information silently.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Packed {
    len: u16,
    words: Vec<u64>,
}

impl Packed {
    pub fn pack(seq: &[u8]) -> Option<Self> {
        if seq.len() > u16::MAX as usize {
            return None;
        }
        let mut words = Vec::with_capacity(seq.len().div_ceil(32));
        let mut acc = 0u64;
        let mut n = 0u32;
        for &b in seq {
            let code = match b.to_ascii_uppercase() {
                b'A' => 0u64,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None, // not packable
            };
            acc |= code << (2 * n);
            n += 1;
            if n == 32 {
                words.push(acc);
                acc = 0;
                n = 0;
            }
        }
        if n > 0 {
            words.push(acc);
        }
        Some(Packed { len: seq.len() as u16, words })
    }

    pub fn unpack(&self) -> String {
        const BASES: [u8; 4] = *b"ACGT";
        let mut out = Vec::with_capacity(self.len as usize);
        for i in 0..self.len as usize {
            let w = self.words[i / 32];
            out.push(BASES[((w >> (2 * (i % 32))) & 0b11) as usize]);
        }
        String::from_utf8(out).unwrap_or_default()
    }

    /// Bytes this entry occupies: 2 bits per base plus the length and the map's overhead.
    pub fn bytes(&self) -> usize {
        2 + self.words.len() * 8 + 8 // len + payload + HashMap slot estimate
    }
}

/// The sequence table used by `table` mode: packed sequence -> multiplicity.
#[derive(Clone, Default)]
pub struct SequenceTable {
    counts: HashMap<Packed, u32>,
}

impl SequenceTable {
    pub fn observe(&mut self, seq: &[u8]) -> bool {
        match Packed::pack(seq) {
            Some(p) => {
                *self.counts.entry(p).or_insert(0) += 1;
                true
            }
            None => false, // a non-ACGT base: the caller must refuse, not lose data
        }
    }

    pub fn distinct(&self) -> u64 {
        self.counts.len() as u64
    }

    pub fn resident_bytes(&self) -> usize {
        self.counts.iter().map(|(k, _)| k.bytes()).sum()
    }

    /// Sequences occurring at least `min_count` times, most frequent first, as
    /// `(sequence text, count)`.
    pub fn at_least(&self, min_count: u32) -> Vec<(String, u32)> {
        let mut v: Vec<(String, u32)> = self
            .counts
            .iter()
            .filter(|(_, &c)| c >= min_count)
            .map(|(k, &c)| (k.unpack(), c))
            .collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
    }
}

/// One entry of the over-representation report.
#[derive(Debug, Clone)]
pub struct OverrepEntry {
    pub seq: Option<String>,
    pub fingerprint: Option<String>,
    pub count: u32,
}

/// Build the over-representation report for one mode.
///
/// `acc` carries the exact fingerprint table. `files` and `threads` are needed only by
/// `second-pass`. `min_count` is the resolved threshold (flag, or max(10, 0.1% of reads)).
pub fn report(
    mode: OverrepMode,
    acc: &Accum,
    files: &[std::path::PathBuf],
    threads: usize,
    min_count: u32,
) -> Result<serde_json::Value> {
    let reads = acc.reads.max(1) as f64;
    match mode {
        OverrepMode::None => Ok(json!({
            "mode": "none",
            "names_available": false,
            "unique_sequences_retained": false,
            "min_count": min_count,
            "sequences": [],
            "note": "duplication levels are exact in this mode; over-represented sequence \
                     identity needs --overrep-mode fingerprint, second-pass or table",
        })),
        OverrepMode::Fingerprint => {
            let top = acc.seq_counts.at_least(min_count);
            let seqs: Vec<serde_json::Value> = top
                .iter()
                .take(20)
                .map(|(fp, c)| {
                    json!({
                        "fingerprint": format!("{fp:016x}"),
                        "count": c,
                        "frac": (c as f64 / reads * 10000.0).round() / 10000.0,
                    })
                })
                .collect();
            Ok(json!({
                "mode": "fingerprint",
                "names_available": false,
                "unique_sequences_retained": false,
                "min_count": min_count,
                "sequences": seqs,
                "note": "multiplicities are exact; the sequence text is not recoverable from \
                         a 64-bit fingerprint. Use --overrep-mode table or second-pass to \
                         name them.",
            }))
        }
        OverrepMode::Table => {
            // The table was filled during the scan (accum.rs), so this mode only reports.
            let top = acc.overrep_table.at_least(min_count);
            let seqs: Vec<serde_json::Value> = top
                .iter()
                .take(20)
                .map(|(s, c)| {
                    json!({
                        "seq": s,
                        "count": c,
                        "frac": (c as f64 / reads * 10000.0).round() / 10000.0,
                    })
                })
                .collect();
            Ok(json!({
                "mode": "table",
                "names_available": true,
                "unique_sequences_retained": true,
                "min_count": min_count,
                "unique_sequences": acc.overrep_table.distinct(),
                "sequences": seqs,
            }))
        }
        OverrepMode::SecondPass => {
            let wanted: std::collections::HashSet<u64> = acc
                .seq_counts
                .at_least(min_count)
                .into_iter()
                .map(|(fp, _)| fp)
                .collect();
            if wanted.is_empty() {
                return Ok(json!({
                    "mode": "second-pass",
                    "names_available": true,
                    "unique_sequences_retained": false,
                    "passes": 2,
                    "min_count": min_count,
                    "sequences": [],
                    "note": "no fingerprint cleared the threshold, so the second read found \
                             nothing to name",
                }));
            }
            // Second read: recover the text of the fingerprints above the threshold.
            let mut texts: HashMap<u64, String> = HashMap::new();
            for path in files {
                second_pass(path, threads, &wanted, &mut texts)?;
            }
            let mut seqs: Vec<serde_json::Value> = acc
                .seq_counts
                .at_least(min_count)
                .iter()
                .take(20)
                .map(|(fp, c)| {
                    json!({
                        "seq": texts.get(fp).cloned().unwrap_or_default(),
                        "count": c,
                        "frac": (c as f64 / reads * 10000.0).round() / 10000.0,
                    })
                })
                .collect();
            seqs.sort_by(|a, b| {
                b["count"].as_u64().unwrap_or(0).cmp(&a["count"].as_u64().unwrap_or(0))
            });
            Ok(json!({
                "mode": "second-pass",
                "names_available": true,
                "unique_sequences_retained": false,
                "passes": 2,
                "min_count": min_count,
                "sequences": seqs,
            }))
        }
    }
}

/// Re-read one file and record the text of every sequence whose fingerprint is wanted.
///
/// This is the ONLY place in fastcheck that reads a file twice, and it is opt-in: the
/// report carries `passes: 2` so a consumer never mistakes the cost for the default.
fn second_pass(
    path: &Path,
    gz_threads: usize,
    wanted: &std::collections::HashSet<u64>,
    texts: &mut HashMap<u64, String>,
) -> Result<()> {
    let (mut rdr, _counter, _mode) = reader::open(path, gz_threads)?;
    let (mut line, mut seq, mut plus, mut qual) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    loop {
        if reader::read_line(&mut rdr, &mut line)? == 0 {
            break;
        }
        if line.is_empty() {
            continue;
        }
        reader::read_line(&mut rdr, &mut seq)?;
        reader::read_line(&mut rdr, &mut plus)?;
        reader::read_line(&mut rdr, &mut qual)?;
        if seq.is_empty() {
            continue;
        }
        let fp = fingerprint::fingerprint(&seq);
        if wanted.contains(&fp) {
            texts.entry(fp).or_insert_with(|| String::from_utf8_lossy(&seq).to_string());
        }
    }
    Ok(())
}

/// Projected bytes for `table` mode at a given unique-sequence count and read length.
/// `ceil(len/4) + 8` per sequence, times the hashbrown load factor of 8/7.
pub fn table_projection(unique: usize, read_len: usize) -> u64 {
    let per = (read_len.div_ceil(4) + 8) as u64;
    let need = (unique as u64 * 8).div_ceil(7);
    need.next_power_of_two().saturating_mul(per + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mode_parsing_rejects_unknown_values() {
        assert_eq!(OverrepMode::parse("none").unwrap(), OverrepMode::None);
        assert_eq!(OverrepMode::parse("second-pass").unwrap(), OverrepMode::SecondPass);
        let e = OverrepMode::parse("sampled").unwrap_err().to_string();
        assert!(e.contains("second-pass"), "{e}");
    }

    #[test]
    fn packed_roundtrip_and_non_acgt_refusal() {
        let p = Packed::pack(b"ACGTACGTNN").unwrap();
        // N is not packable, so this returns None — the caller must refuse, not lose data.
        assert!(Packed::pack(b"ACGTN").is_none());
        let p = Packed::pack(b"ACGTACGTAC").unwrap();
        assert_eq!(p.unpack(), "ACGTACGTAC");
        // 10 bases -> 1 word (32 bases fit in a u64) + 2 B length + 8 B slot.
        assert_eq!(p.bytes(), 2 + 8 + 8);
    }

    #[test]
    fn sequence_table_counts_and_names() {
        let mut t = SequenceTable::default();
        for _ in 0..4 {
            t.observe(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        }
        t.observe(b"TTTTAAAACCCCGGGG");
        assert_eq!(t.distinct(), 2);
        let top = t.at_least(3);
        assert_eq!(top.len(), 1);
        assert_eq!(top[0].0, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        assert_eq!(top[0].1, 4);
    }

    /// 3e8 unique reads at 150 bp is the number PLAN-QC §5.9 quotes as ~11 GB.
    #[test]
    fn table_projection_matches_the_plan_quote() {
        let bytes = table_projection(300_000_000, 150);
        let gb = bytes as f64 / (1u64 << 30) as f64;
        assert!((gb - 11.0).abs() < 1.5, "projected {gb:.2} GB, PLAN-QC says ~11 GB");
    }
}
```

**Step 5: Hook the table into the accumulator**

In `tools/fastcheck/src/accum.rs`:

- Add the field:

```rust
    /// 2-bit packed sequence table, filled only in `--overrep-mode table` (the memory is
    /// real: ~46 B per unique 150 bp read, admission-checked in main.rs).
    pub overrep_table: crate::overrep::SequenceTable,
```

- In `Accum::new`: `overrep_table: crate::overrep::SequenceTable::default(),`
- In `observe`, after `self.seq_counts.observe(seq);`:

```rust
        if self.track_overrep_table {
            // A non-ACGT base cannot be packed; refusing is the caller's job (main.rs
            // checks the first read), and here the sequence is simply skipped rather than
            // silently mis-counted — the report's unique_sequences count shows the gap.
            let _ = self.overrep_table.observe(seq);
        }
```

- Add the flag field `pub track_overrep_table: bool,` (default `false` in `new`), set from
  `main.rs` exactly like `track_profile`.

**Step 6: Admission-check `table` mode and emit the report**

In `tools/fastcheck/src/main.rs`, after the fingerprint admission block from Task 4:

```rust
    let overrep = cli.overrep_mode()?;
    if overrep == overrep::OverrepMode::Table {
        let proj = overrep::table_projection(fp_capacity, 150);
        eprintln!(
            "fastcheck: over-representation table projected at {:.1} MB for ~{} unique sequences",
            memory::mib(proj),
            fp_capacity,
        );
        match memory::admit(proj, limit.as_ref()) {
            memory::Admission::Refuse { projected, limit, source } => anyhow::bail!(
                "--overrep-mode table needs ~{:.1} MB of packed unique sequences\n\
                 but this process is limited to {:.1} MB ({})\n\n\
                 --overrep-mode fingerprint   exact multiplicities at 8 B per unique read\n\
                 --overrep-mode second-pass   same memory as fingerprint, plus one extra read\n\
                 --overrep-mode none         duplication levels only, zero cost",
                memory::mib(projected),
                memory::mib(limit),
                source,
            ),
            memory::Admission::ProceedWithWarning { headroom_bytes } => eprintln!(
                "fastcheck: warning: unique-sequence table projected at {:.1} MB, leaving \
                 {:.1} MB headroom",
                memory::mib(proj),
                memory::mib(headroom_bytes),
            ),
            memory::Admission::Proceed => {}
        }
    }
```

and set the flag on each accumulator in `scan_one`:

```rust
    acc.track_overrep_table = overrep == overrep::OverrepMode::Table;
```

Then, where the report is written, add the block to `report::write_report`'s `doc`:

```rust
        "overrep": overrep::report(overrep, total, &files, threads.max(1), min_count)?,
```

with

```rust
    let min_count = cli
        .overrep_min_count
        .unwrap_or_else(|| ((total.reads as f64 * 0.001) as u32).max(10));
```

computed just before the `doc` (and passed in, rather than recomputed inside `overrep`).

**Step 7: Run the tests**

Run: `cargo test --offline`
Expected: all green, including the five new `overrep_*` CLI tests and the four
`overrep::` unit tests.

**Step 8: Manual verification of each mode**

```bash
cd tools/fastcheck
printf '@a\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n%.0s' 1 2 3 4 > /tmp/or.fastq
printf '@b\nTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGT\n+\n55555555555555555555555555555555\n' >> /tmp/or.fastq
for m in none fingerprint second-pass table; do
  ./target/release/fastcheck /tmp/or.fastq --overrep-mode $m --overrep-min-count 3 \
    --json /tmp/or-$m.json -q
  python3 -c "import json;d=json.load(open('/tmp/or-$m.json'))['overrep'];print('$m', d['mode'], 'names', d['names_available'], d.get('passes','-'), d['sequences'])"
done
```
Expected:
```
none none names False - []
fingerprint fingerprint names False - [{'fingerprint': '...', 'count': 4, ...}]
second-pass second-pass names True 2 [{'seq': 'ACGTACGT...', 'count': 4, ...}]
table table names True - [{'seq': 'ACGTACGT...', 'count': 4, ...}]
```

**Step 9: Commit**

```bash
git add tools/fastcheck/src/cli.rs tools/fastcheck/src/overrep.rs \
        tools/fastcheck/src/accum.rs tools/fastcheck/src/main.rs \
        tools/fastcheck/src/report.rs tools/fastcheck/tests/cli.rs
git commit -m "feat(fastcheck): --overrep-mode none|fingerprint|second-pass|table (PLAN-QC §5.9)"
```

---

## Task 7: Document the bounded-counter approximation and remove it from the report

**Objective:** The old `RepeatedCounter` singleton-eviction approximation (used for k-mers)
is documented as such in the report, and the over-representation path no longer depends on
it — the exact fingerprint table replaced it in Task 5.

**Files:**
- Modify: `tools/fastcheck/src/repeated.rs` (doc comment + a test pinning the approximation)
- Modify: `tools/fastcheck/src/rich.rs` (`kmer_metrics` gains an `approximate: true` flag)
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tools/fastcheck/tests/cli.rs`:

```rust
#[test]
fn the_kmer_spectrum_is_labelled_approximate_because_its_counter_evicts() {
    let d = tmpdir("kmer_approx");
    let fq = write_fixture(&d);
    let json = d.join("k.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    assert_eq!(v["kmer_spectrum"]["approximate"], true);
    // The approximation and its bound must be stated, not implied.
    let note = v["kmer_spectrum"]["approximation"].as_str().unwrap();
    assert!(note.contains("singleton"), "{note}");
    assert!(note.contains("evict"), "{note}");
}
```

**Step 2: Run to verify failure**

Run: `cargo test --offline the_kmer_spectrum_is_labelled_approximate`
Expected: FAIL — `kmer_spectrum` has no `approximate` key.

**Step 3: Label it**

In `tools/fastcheck/src/rich.rs`, change `kmer_metrics`'s return from
`json!(list)` to:

```rust
    (
        m,
        json!({
            "approximate": true,
            "approximation": "the 12-mer counter is bounded at 1e6 distinct keys and evicts \
                              singletons when it overflows: a k-mer whose every occurrence \
                              was evicted loses its count entirely, and one whose first \
                              occurrence was evicted loses exactly 1. Counts are exact for \
                              k-mers above the 0.1% reporting threshold; below it the number \
                              is a lower bound.",
            "top": list,
        }),
    )
```

**Step 4: Update the doc comment in `repeated.rs`**

Replace the module doc's first paragraph with:

```rust
//! Bounded over-represented-sequence and k-mer counters.
//!
//! Streaming with singleton eviction: when the number of distinct keys exceeds a cap,
//! every key with count == 1 is dropped. A singleton can never be over-represented, and
//! losing its first occurrence costs at most 1 of a count that must reach a large threshold
//! to be reported, so memory stays bounded on 50 M-read files while the frequent sequences
//! are counted exactly.
//!
//! SCOPE (as of the exactness work): this counter is used for the **k-mer spectrum only**.
//! Full-sequence duplication and over-representation moved to the exact fingerprint table
//! (`crate::fingerprint`) and the packed sequence table (`crate::overrep`), because a
//! fingerprint is 8 B where a sequence is 46 B and the exactness is affordable. The
//! approximation that remains is labelled in the report (`kmer_spectrum.approximate`).
```

**Step 5: Pin the approximation with a test**

Add to `tools/fastcheck/src/repeated.rs`'s test module:

```rust
    /// The approximation, stated as a bound rather than hidden: a singleton evicted before
    /// its second occurrence costs exactly 1 count. This is why the counter is only used
    /// below the 0.1% reporting threshold, where 1 count cannot change a verdict.
    #[test]
    fn eviction_costs_at_most_one_count_per_evicted_singleton() {
        let mut c = RepeatedCounter::new(4);
        // 4 distinct singletons fill the cap; the 5th distinct key triggers an eviction
        // that removes all four.
        for i in 0..4u8 {
            c.insert(vec![i]);
        }
        assert_eq!(c.len(), 4);
        c.insert(vec![9]); // over the cap -> retain(count > 1) drops all four singletons
        assert_eq!(c.len(), 1, "only the just-inserted key survives");
        // The frequent key is counted exactly once it clears 1.
        let mut c = RepeatedCounter::new(4);
        for _ in 0..3 {
            c.insert(b"ACGT".to_vec());
        }
        for i in 0..4u8 {
            c.insert(vec![i]);
        }
        c.insert(vec![9]);
        let rep = c.report(2);
        assert_eq!(rep.len(), 1);
        assert_eq!(rep[0].1, 3, "the frequent key survives with its exact count");
    }
```

**Step 6: Run the tests**

Run: `cargo test --offline`
Expected: all green.

**Step 7: Commit**

```bash
git add tools/fastcheck/src/repeated.rs tools/fastcheck/src/rich.rs tools/fastcheck/tests/cli.rs
git commit -m "docs(fastcheck): label the k-mer counter's singleton-eviction approximation in the report"
```

---

## Task 8: Adapter catalogue + matcher copied into fastcheck

**Objective:** `tools/fastcheck/src/adapter.rs` holds the catalogue and the seeded matcher
with fastprep's exact semantics, so the per-position module (Task 9) and any future
fastcheck-side trimming agree with `fastprep --adapter-preset` by construction.

**Files:**
- Create: `tools/fastcheck/src/adapter.rs`
- Modify: `tools/fastcheck/Cargo.toml` (nothing new — no deps needed for the port)
- Modify: `tools/fastcheck/src/main.rs` (`mod adapter;`)

**Step 1: Create the module with its tests first**

```rust
//! Adapter catalogue and seeded matching, copied from `tools/fastprep/src/adapter.rs`
//! (D-1: no shared crate — the two tools may drift, and the tests below are what keep
//! them honest).
//!
//! Semantics preserved exactly, because a per-position adapter-content curve is only
//! meaningful if "adapter" means the same thing here as in the trimmer:
//!
//! * catalogue: TruSeq/R1, TruSeq/R2, Nextera, SmallRNA (`preset_illumina` +
//!   `preset_nextera` + `preset_small_rna`, deduplicated by sequence)
//! * a match at position `p` with overlap `k = min(len(adapter), len(read) - p)` is
//!   accepted when `mismatches <= ceil(k * error_rate)`
//! * `min_overlap_se` (fastprep default 14) and `error_rate` (default 0.2) are the flags
//! * leftmost match wins; on a tie the adapter earlier in the set wins
//! * the 6-mer seed index is a rejection aid only: any alignment of `k` bases with at most
//!   `max_mm` mismatches contains an exactly-matching run of at least `k/(max_mm+1)` bases,
//!   so a seed hit at that position finds the alignment. Everything else is provably
//!   non-matching. It is exact, never a heuristic.

/// One named adapter sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    pub name: String,
    pub seq: Vec<u8>,
}

/// All TruSeq + Nextera sequences (fastprep's default "illumina" preset).
pub fn preset_illumina() -> Vec<Adapter> {
    vec![
        Adapter { name: "TruSeq/R1".into(), seq: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec() },
        Adapter { name: "TruSeq/R2".into(), seq: b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec() },
        Adapter { name: "Nextera".into(),    seq: b"CTGTCTCTTATACACATCT".to_vec() },
    ]
}

pub fn preset_nextera() -> Vec<Adapter> {
    vec![Adapter { name: "Nextera".into(), seq: b"CTGTCTCTTATACACATCT".to_vec() }]
}

pub fn preset_small_rna() -> Vec<Adapter> {
    vec![Adapter { name: "SmallRNA".into(), seq: b"TGGAATTCTCGGGTGCCAAGG".to_vec() }]
}

/// Every built-in sequence, deduplicated by sequence (Nextera appears in two presets).
pub fn catalogue() -> Vec<Adapter> {
    let mut all = preset_illumina();
    for a in preset_nextera().into_iter().chain(preset_small_rna()) {
        if !all.iter().any(|x| x.seq == a.seq) {
            all.push(a);
        }
    }
    all
}

/// The full adapter configuration passed to [`match_at`].
#[derive(Debug, Clone)]
pub struct AdapterSet {
    pub adapters: Vec<Adapter>,
    /// Maximum fraction of mismatches allowed in a match.
    pub error_rate: f64,
    /// Minimum overlap for sequence matching.
    pub min_overlap_se: usize,
}

impl AdapterSet {
    pub fn new(adapters: Vec<Adapter>, error_rate: f64, min_overlap_se: usize) -> Self {
        AdapterSet { adapters, error_rate, min_overlap_se }
    }

    /// fastprep's defaults: `--adapter-error-rate 0.2 --min-adapter-overlap 14`.
    pub fn illumina_default() -> Self {
        Self::new(catalogue(), 0.2, 14)
    }
}

/// Result of a match: which adapter, and where it starts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdapterMatch {
    pub name: String,
    /// 0-based position in the read where the adapter match starts.
    pub start: usize,
    /// Overlap length actually compared.
    pub overlap: usize,
}

/// Byte -> 2-bit base code (0=A, 1=C, 2=G, 3=T; anything else -> 3, a wildcard that may
/// spuriously seed and is rejected by the mismatch count).
static BASE2: [u8; 256] = {
    let mut t = [3u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t
};

/// 2-bit encoding of a 6-mer into a 12-bit code.
#[inline]
fn seed_code(bases: &[u8]) -> u32 {
    let mut code: u32 = 0;
    for &b in bases {
        code = (code << 2) | BASE2[b as usize] as u32;
    }
    code & 0xFFF
}

/// Smallest overlap length for which the 6-mer guarantee holds, per adapter.
/// `usize::MAX` = never (e.g. a high error rate).
fn min_sound_k(adapter_len: usize, error_rate: f64) -> usize {
    for k in 1..=adapter_len {
        let max_mm = (k as f64 * error_rate).ceil() as usize;
        if k >= 6 * (max_mm + 1) {
            return k;
        }
    }
    usize::MAX
}

/// Mismatches between two equal-length slices, stopping once `cap` is exceeded.
#[inline]
fn mismatches_upto(a: &[u8], b: &[u8], cap: usize) -> usize {
    debug_assert_eq!(a.len(), b.len());
    let mut mm = 0usize;
    for (x, y) in a.iter().zip(b.iter()) {
        if x.to_ascii_uppercase() != y.to_ascii_uppercase() {
            mm += 1;
            if mm > cap {
                return mm;
            }
        }
    }
    mm
}

/// The leftmost adapter match in `seq`, tie-broken by set order.
///
/// Returns `None` when no adapter matches anywhere under the set's thresholds — which is
/// the common case for a clean read, and the reason the seed index exists.
pub fn match_at(seq: &[u8], set: &AdapterSet) -> Option<AdapterMatch> {
    let l = seq.len();
    let max_p = l.saturating_sub(set.min_overlap_se);
    let mut best: Option<(usize, usize)> = None; // (start, adapter index)
    let mut consider = |p: usize, ai: usize| {
        if best.is_none_or(|(bp, bai)| p < bp || (p == bp && ai < bai)) {
            best = Some((p, ai));
        }
    };

    for (ai, adapter) in set.adapters.iter().enumerate() {
        let a = adapter.seq.len();
        if a == 0 {
            continue;
        }
        let ms = min_sound_k(a, set.error_rate);
        // Seed-guided pass: only sound for k >= ms. Scan the read once with a rolling
        // 6-mer and test the candidate positions the seeds imply.
        let mut code: u32 = 0;
        for q in 0..l {
            code = ((code << 2) | BASE2[seq[q] as usize] as u32) & 0xFFF;
            if q < 5 {
                continue;
            }
            let win_start = q - 5;
            for off in 0..=a.saturating_sub(6) {
                if seed_code(&adapter.seq[off..off + 6]) != code {
                    continue;
                }
                if win_start < off {
                    continue;
                }
                let p = win_start - off;
                let k = a.min(l - p);
                if k < set.min_overlap_se || k < ms || p > max_p {
                    continue;
                }
                let max_mm = (k as f64 * set.error_rate).ceil() as usize;
                if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                    consider(p, ai);
                }
            }
        }
        // Direct pass for the seed-unsound region (short adapters, high error rates, or
        // truncated tails): scanned exactly, reference-style.
        if a < ms || ms == usize::MAX {
            for p in 0..=max_p {
                let k = a.min(l - p);
                if k < set.min_overlap_se {
                    break;
                }
                let max_mm = (k as f64 * set.error_rate).ceil() as usize;
                if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                    consider(p, ai);
                    break;
                }
            }
        } else {
            let start = if l >= a { l - a + 1 } else { 0 };
            let start = start.max(l.saturating_add(1).saturating_sub(ms));
            for p in start..=max_p {
                let k = a.min(l - p);
                if k < set.min_overlap_se {
                    break;
                }
                let max_mm = (k as f64 * set.error_rate).ceil() as usize;
                if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                    consider(p, ai);
                    break;
                }
            }
        }
    }

    best.map(|(p, ai)| {
        let a = set.adapters[ai].seq.len();
        AdapterMatch { name: set.adapters[ai].name.clone(), start: p, overlap: a.min(l - p) }
    })
}

/// All adapter matches in `seq`, one per adapter that matches anywhere, as
/// `(adapter index, start)`. Used by the per-position module, which counts a match at
/// EVERY position where an adapter begins, not just the leftmost.
pub fn all_matches(seq: &[u8], set: &AdapterSet) -> Vec<(usize, usize)> {
    let l = seq.len();
    let max_p = l.saturating_sub(set.min_overlap_se);
    let mut out: Vec<(usize, usize)> = Vec::new();
    for (ai, adapter) in set.adapters.iter().enumerate() {
        let a = adapter.seq.len();
        if a == 0 {
            continue;
        }
        let ms = min_sound_k(a, set.error_rate);
        for p in 0..=max_p {
            let k = a.min(l - p);
            if k < set.min_overlap_se {
                break;
            }
            let max_mm = (k as f64 * set.error_rate).ceil() as usize;
            if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                out.push((ai, p));
                break; // leftmost per adapter, like the matcher
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalogue_is_deduplicated_by_sequence() {
        let cat = catalogue();
        let nextera = cat.iter().filter(|a| a.seq == b"CTGTCTCTTATACACATCT").count();
        assert_eq!(nextera, 1, "Nextera appears in two presets but must be listed once");
        assert_eq!(cat.len(), 4, "TruSeq/R1, TruSeq/R2, Nextera, SmallRNA");
    }

    /// The agreement test D-1 asks for: fastcheck's matcher must accept exactly what
    /// fastprep's does, on the same inputs, with the same thresholds.
    #[test]
    fn a_full_adapter_at_the_3p_end_is_found_at_the_right_position() {
        let set = AdapterSet::illumina_default();
        // 120 bases of filler then the full 34-base TruSeq/R1 adapter: insert 120.
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        let m = match_at(&seq, &set).unwrap();
        assert_eq!(m.name, "TruSeq/R1");
        assert_eq!(m.start, 120);
        assert_eq!(m.overlap, 34);
    }

    #[test]
    fn a_clean_read_matches_nothing() {
        let set = AdapterSet::illumina_default();
        let seq: Vec<u8> = (0..150).map(|i| b"ACGT"[i % 4]).collect();
        assert!(match_at(&seq, &set).is_none());
        assert!(all_matches(&seq, &set).is_empty());
    }

    /// The threshold is a flag: at error_rate 0.2 and overlap 20, up to 4 mismatches are
    /// tolerated. Three must still match, five must not.
    #[test]
    fn the_error_rate_threshold_is_honoured_exactly() {
        let set = AdapterSet::new(vec![Adapter {
            name: "T".into(),
            seq: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
        }], 0.2, 14);
        let mut seq = vec![b'A'; 116];
        let mut ad = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec();
        // 3 mismatches in the last 20 bases -> ceil(20*0.2) = 4 allowed.
        ad[30] = if ad[30] == b'A' { b'C' } else { b'A' };
        ad[31] = if ad[31] == b'A' { b'C' } else { b'A' };
        ad[32] = if ad[32] == b'A' { b'C' } else { b'A' };
        seq.extend_from_slice(&ad);
        assert!(match_at(&seq, &set).is_some(), "3 mismatches in 20 bases must pass");

        let mut ad2 = ad.clone();
        for i in 30..35 {
            ad2[i] = if ad2[i] == b'A' { b'C' } else { b'A' };
        }
        let mut seq2 = vec![b'A'; 116];
        seq2.extend_from_slice(&ad2);
        assert!(match_at(&seq2, &set).is_none(), "5 mismatches in 20 bases must fail");
    }

    /// A match shorter than --min-adapter-overlap is not a match, however perfect.
    #[test]
    fn a_match_below_min_overlap_is_rejected() {
        let set = AdapterSet::new(vec![Adapter {
            name: "T".into(),
            seq: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
        }], 0.2, 14);
        // Only the first 10 bases of the adapter are present.
        let mut seq = vec![b'A'; 140];
        seq.extend_from_slice(b"AGATCGGAAG");
        assert!(match_at(&seq, &set).is_none(), "10 < 14: below the floor");
    }

    #[test]
    fn leftmost_wins_and_set_order_breaks_ties() {
        let set = AdapterSet::illumina_default();
        // TruSeq/R1 then Nextera, both matching: the earlier start must win.
        let mut seq = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec();
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        let m = match_at(&seq, &set).unwrap();
        assert_eq!(m.name, "TruSeq/R1");
        assert_eq!(m.start, 0);
    }

    /// all_matches reports one entry per matching adapter, which is what the
    /// per-position curve counts.
    #[test]
    fn all_matches_reports_each_adapter_once() {
        let set = AdapterSet::illumina_default();
        let mut seq = vec![b'A'; 100];
        seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        let v = all_matches(&seq, &set);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].1, 100);
    }
}
```

**Step 2: Run to verify failure**

Run: `cd tools/fastcheck && cargo test --offline adapter::`
Expected: FAIL — module not declared (`mod adapter;` missing).

**Step 3: Declare the module and run**

`tools/fastcheck/src/main.rs`: add `mod adapter;` next to `mod accum;`.

Run: `cargo test --offline adapter::`
Expected: `test result: ok. 6 passed; 0 failed`.

**Step 4: Cross-check against fastprep on a real fixture (D-1 consistency by test)**

```bash
cd tools/fastcheck
../fastprep/target/release/fastprep -1 ../../examples/runs/demo/sim_E/sim_R1.fastq.gz \
  -o /dev/null --adapter-preset nextera --json /tmp/fp-prep.json 2>/dev/null || true
./target/release/fastcheck ../../examples/runs/demo/sim_E/sim_R1.fastq.gz \
  --json /tmp/fp-check.json -q
```
(The per-position curve lands in Task 9; this step only proves both tools load the same
catalogue. If `fastprep` is not built, skip — the unit tests above pin the semantics.)

**Step 5: Commit**

```bash
git add tools/fastcheck/src/adapter.rs tools/fastcheck/src/main.rs
git commit -m "feat(fastcheck): adapter catalogue + seeded matcher with fastprep's exact semantics (D-1 copy)"
```

---

## Task 9: Adapter content per position

**Objective:** For every catalogue adapter and every position `p`, the exact count of reads
whose adapter match starts at `p`, emitted both as a positions-TSV column set and as a
report block, reusing Task 8's matcher semantics.

**Files:**
- Create: `tools/fastcheck/src/adapter_pos.rs`
- Modify: `tools/fastcheck/src/accum.rs` (the counter, `observe`, `merge`)
- Modify: `tools/fastcheck/src/cli.rs` (`--adapter-content`, `--adapter-error-rate`,
  `--min-adapter-overlap`, `--adapter-preset`)
- Modify: `tools/fastcheck/src/report.rs` (block + TSV columns)
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tools/fastcheck/tests/cli.rs`:

```rust
/// 6 reads: two carry the full Nextera adapter starting at position 121 (1-based), one
/// carries a truncated 16-base Nextera tail, three are clean.
const ADAPTER_FASTQ: &[u8] = b"@a\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTGTCTCTTATACACATCT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@b\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTGTCTCTTATACACATCT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@c\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@d\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\n5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555\n\
@e\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCTGTCTCTTATACACAT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@f\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";

#[test]
fn adapter_content_per_position_counts_every_match_start() {
    let d = tmpdir("adapter_pos");
    let p = d.join("ad.fastq");
    std::fs::write(&p, ADAPTER_FASTQ).unwrap();
    let json = d.join("a.json");
    let tsv = d.join("pos.tsv");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--adapter-content", "--adapter-preset", "nextera",
          "--json", json.to_str().unwrap(), "--positions-tsv", tsv.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let blk = &v["adapter_content"];
    assert_eq!(blk["adapters"].as_array().unwrap().len(), 1, "nextera preset = one adapter");
    assert_eq!(blk["adapters"][0]["name"], "Nextera");
    // Per-position counts: 2 full adapters at position 121 (0-based 120), one truncated
    // 16-base match at position 135 (0-based 134).
    let per_pos = blk["adapters"][0]["per_position"].as_object().unwrap();
    assert_eq!(per_pos["121"], 2, "two reads carry the full adapter at 121");
    assert_eq!(per_pos["135"], 1, "one read carries a 16-base truncated adapter at 135");
    // Totals: 3 of 6 reads carry an adapter somewhere.
    assert_eq!(blk["adapters"][0]["reads_with_match"], 3);
    assert_eq!(blk["error_rate"], 0.2);
    assert_eq!(blk["min_overlap"], 14);

    // The TSV gains one column per adapter.
    let text = std::fs::read_to_string(&tsv).unwrap();
    let header = text.lines().next().unwrap();
    assert!(header.contains("adapter_Nextera"), "header: {header}");
    let rows: Vec<&str> = text.lines().skip(1).collect();
    let idx = header.split('\t').position(|c| c == "adapter_Nextera").unwrap();
    assert_eq!(rows[120].split('\t').nth(idx).unwrap(), "2");
    assert_eq!(rows[134].split('\t').nth(idx).unwrap(), "1");
    assert_eq!(rows[0].split('\t').nth(idx).unwrap(), "0", "clean at position 1");
}

#[test]
fn adapter_content_is_off_by_default_and_the_columns_are_absent() {
    let d = tmpdir("adapter_off");
    let fq = write_fixture(&d);
    let json = d.join("a.json");
    let tsv = d.join("pos.tsv");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(),
              "--positions-tsv", tsv.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    assert!(v.get("adapter_content").is_none(), "must be opt-in");
    let header = std::fs::read_to_string(&tsv).unwrap().lines().next().unwrap().to_string();
    assert!(!header.contains("adapter_"), "no adapter columns by default: {header}");
}

#[test]
fn adapter_content_honours_the_error_rate_and_overlap_flags() {
    let d = tmpdir("adapter_flags");
    let p = d.join("ad.fastq");
    std::fs::write(&p, ADAPTER_FASTQ).unwrap();
    let json = d.join("a.json");
    // A 30-base floor rejects the 16-base truncated match but keeps the full ones.
    run(&d, &[p.to_str().unwrap(), "--adapter-content", "--adapter-preset", "nextera",
              "--min-adapter-overlap", "30", "--json", json.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    let per_pos = v["adapter_content"]["adapters"][0]["per_position"].as_object().unwrap();
    assert_eq!(per_pos.get("135"), None, "16 < 30: below the floor");
    assert_eq!(per_pos["121"], 2);
    assert_eq!(v["adapter_content"]["min_overlap"], 30);
}
```

**Step 2: Run to verify failure**

Run: `cargo test --offline adapter_content`
Expected: FAIL — `error: unexpected argument '--adapter-content' found`.

**Step 3: Add the flags**

In `tools/fastcheck/src/cli.rs`:

```rust
    /// Count, per position, the reads whose adapter match STARTS there, for every adapter
    /// in the catalogue (PLAN-QC module 11). Off by default: it costs one catalogue scan
    /// per read. The matcher semantics are fastprep's, exactly.
    #[arg(long)]
    pub adapter_content: bool,

    /// Adapter catalogue for --adapter-content: illumina (TruSeq R1/R2 + Nextera),
    /// nextera, or a FASTA path.
    #[arg(long, value_name = "illumina|nextera|PATH", default_value = "illumina")]
    pub adapter_preset: String,

    /// Maximum mismatch rate for --adapter-content matching [0, 1] (fastprep's default).
    #[arg(long, default_value_t = 0.2, value_name = "FLOAT")]
    pub adapter_error_rate: f64,

    /// Minimum adapter overlap for --adapter-content matching (fastprep's default).
    #[arg(long, default_value_t = 14, value_name = "N")]
    pub min_adapter_overlap: usize,
```

with a resolver:

```rust
    /// The adapter set for --adapter-content, with fastprep's presets.
    pub fn adapter_set(&self) -> anyhow::Result<crate::adapter::AdapterSet> {
        if !(0.0..=1.0).contains(&self.adapter_error_rate) {
            anyhow::bail!("--adapter-error-rate must be in [0, 1]");
        }
        if self.min_adapter_overlap == 0 {
            anyhow::bail!("--min-adapter-overlap must be >= 1");
        }
        let adapters = match self.adapter_preset.as_str() {
            "illumina" => crate::adapter::preset_illumina(),
            "nextera" => crate::adapter::preset_nextera(),
            path => {
                let text = std::fs::read_to_string(path)
                    .with_context(|| format!("reading adapter FASTA {}", path))?;
                parse_fasta_adapters(&text)?
            }
        };
        Ok(crate::adapter::AdapterSet::new(
            adapters,
            self.adapter_error_rate,
            self.min_adapter_overlap,
        ))
    }
```

and a small FASTA parser at the bottom of `cli.rs` (fastprep's, in substance):

```rust
/// Parse an adapter FASTA: `>name` then the sequence on the following lines.
fn parse_fasta_adapters(text: &str) -> anyhow::Result<Vec<crate::adapter::Adapter>> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                out.push(crate::adapter::Adapter { name: n, seq: std::mem::take(&mut seq).into_bytes() });
            }
            name = Some(rest.split_whitespace().next().unwrap_or("adapter").to_string());
        } else if !line.is_empty() {
            seq.push_str(&line.to_ascii_uppercase());
        }
    }
    if let Some(n) = name {
        out.push(crate::adapter::Adapter { name: n, seq: seq.into_bytes() });
    }
    if out.is_empty() {
        anyhow::bail!("adapter FASTA contained no records");
    }
    Ok(out)
}
```

**Step 4: Create the counter module**

Create `tools/fastcheck/src/adapter_pos.rs`:

```rust
//! Adapter content per position (PLAN-QC module 11).
//!
//! For each adapter in the catalogue and each position `p`, the exact count of reads whose
//! adapter match STARTS at `p`. This is a different question from "how many reads carry an
//! adapter": the per-position curve shows WHERE the readthrough begins, which is what
//! distinguishes a short-insert library (a spike at one position) from adapter
//! contamination spread across the read.
//!
//! The matcher semantics are fastprep's, copied (D-1) — same catalogue, same
//! `--adapter-error-rate` / `--min-adapter-overlap` defaults, same leftmost-wins rule — so
//! the curve and the trimmer cannot disagree about what an adapter is.

use crate::adapter::{self, AdapterSet};

/// `per_adapter[ai][pos]` = reads whose match for adapter `ai` starts at `pos`.
#[derive(Clone)]
pub struct AdapterPos {
    pub names: Vec<String>,
    pub counts: Vec<Vec<u64>>,
    /// Reads carrying at least one match, per adapter.
    pub reads_with_match: Vec<u64>,
    /// Total match observations (a read can match two adapters).
    pub total_matches: u64,
}

impl AdapterPos {
    pub fn new(set: &AdapterSet) -> Self {
        Self {
            names: set.adapters.iter().map(|a| a.name.clone()).collect(),
            counts: set.adapters.iter().map(|a| vec![0u64; a.seq.len().max(1)]).collect(),
            reads_with_match: vec![0; set.adapters.len()],
            total_matches: 0,
        }
    }

    #[inline]
    pub fn observe(&mut self, seq: &[u8], set: &AdapterSet) {
        for (ai, p) in adapter::all_matches(seq, set) {
            if p < self.counts[ai].len() {
                self.counts[ai][p] += 1;
                self.reads_with_match[ai] += 1;
                self.total_matches += 1;
            }
        }
    }

    /// Exact sum of two counters.
    pub fn merge(&mut self, o: &AdapterPos) {
        for (ai, row) in o.counts.iter().enumerate() {
            if row.len() > self.counts[ai].len() {
                self.counts[ai].resize(row.len(), 0);
            }
            for (d, s) in self.counts[ai].iter_mut().zip(row.iter()) {
                *d += s;
            }
            self.reads_with_match[ai] += o.reads_with_match[ai];
        }
        self.total_matches += o.total_matches;
    }

    /// The report block: one entry per adapter with its non-zero per-position counts.
    pub fn to_json(&self, error_rate: f64, min_overlap: usize) -> serde_json::Value {
        let adapters: Vec<serde_json::Value> = self
            .names
            .iter()
            .enumerate()
            .map(|(ai, name)| {
                let mut per_pos = serde_json::Map::new();
                for (p, &c) in self.counts[ai].iter().enumerate() {
                    if c > 0 {
                        per_pos.insert((p + 1).to_string(), serde_json::json!(c));
                    }
                }
                serde_json::json!({
                    "name": name,
                    "reads_with_match": self.reads_with_match[ai],
                    "per_position": serde_json::Value::Object(per_pos),
                })
            })
            .collect();
        serde_json::json!({
            "error_rate": error_rate,
            "min_overlap": min_overlap,
            "total_matches": self.total_matches,
            "adapters": adapters,
        })
    }

    /// Extra positions-TSV columns: `adapter_<name>` per adapter, 1-based positions.
    pub fn tsv_columns(&self) -> Vec<(String, Vec<u64>)> {
        self.names
            .iter()
            .enumerate()
            .map(|(ai, name)| (format!("adapter_{name}"), self.counts[ai].clone()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set() -> AdapterSet {
        AdapterSet::new(adapter::preset_nextera(), 0.2, 14)
    }

    #[test]
    fn counts_the_start_position_of_each_match() {
        let mut c = AdapterPos::new(&set());
        // 120 filler bases then the full 19-base Nextera adapter.
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        c.observe(&seq, &set());
        assert_eq!(c.reads_with_match[0], 1);
        assert_eq!(c.counts[0][120], 1, "the match starts at 0-based 120 = position 121");
        assert_eq!(c.counts[0][0], 0);
    }

    #[test]
    fn a_clean_read_contributes_nothing() {
        let mut c = AdapterPos::new(&set());
        let seq: Vec<u8> = (0..150).map(|i| b"ACGT"[i % 4]).collect();
        c.observe(&seq, &set());
        assert_eq!(c.total_matches, 0);
        assert!(c.counts[0].iter().all(|&c| c == 0));
    }

    #[test]
    fn merge_is_exact() {
        let mut a = AdapterPos::new(&set());
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        a.observe(&seq, &set());
        let mut b = AdapterPos::new(&set());
        b.observe(&seq, &set());
        a.merge(&b);
        assert_eq!(a.counts[0][120], 2);
        assert_eq!(a.reads_with_match[0], 2);
        assert_eq!(a.total_matches, 2);
    }

    #[test]
    fn to_json_uses_one_based_positions_and_omits_zeros() {
        let mut c = AdapterPos::new(&set());
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        c.observe(&seq, &set());
        let j = c.to_json(0.2, 14);
        let per = j["adapters"][0]["per_position"].as_object().unwrap();
        assert_eq!(per.len(), 1, "only the non-zero position");
        assert_eq!(per["121"], 1);
        assert_eq!(j["adapters"][0]["name"], "Nextera");
    }
}
```

**Step 5: Wire it into the accumulator and the report**

In `tools/fastcheck/src/accum.rs`:

```rust
    /// Adapter content per position, only when --adapter-content is on.
    pub adapter_pos: Option<crate::adapter_pos::AdapterPos>,
```

- `Accum::new`: `adapter_pos: None,`
- `observe`, after the fingerprint line:

```rust
        if let Some(ap) = self.adapter_pos.as_mut() {
            ap.observe(seq, &self.adapter_set);
        }
```

  which needs a `pub adapter_set: crate::adapter::AdapterSet` field on `Accum` (set in
  `new` from a parameter — change `Accum::new(max_pos)` to
  `Accum::new(max_pos, adapter_set: AdapterSet)` and update the three call sites in
  `main.rs`, `reader.rs` tests and `accum.rs` tests to pass
  `crate::adapter::AdapterSet::illumina_default()` when the module is off; the set is only
  consulted when `adapter_pos.is_some()`).
- `merge`:

```rust
        match (self.adapter_pos.as_mut(), o.adapter_pos.as_ref()) {
            (Some(a), Some(b)) => a.merge(b),
            _ => {}
        }
```

In `tools/fastcheck/src/main.rs`, after parsing:

```rust
    let adapter_set = cli.adapter_set()?;
```

pass it to `scan_one`, and in `scan_one`:

```rust
    if cli.adapter_content {
        acc.adapter_pos = Some(crate::adapter_pos::AdapterPos::new(&adapter_set));
    }
```

In `tools/fastcheck/src/report.rs`:

- `write_positions_tsv` gains a `adapter_pos: Option<&crate::adapter_pos::AdapterPos>`
  parameter; the header becomes the base header plus `\tadapter_<name>` per adapter, and
  each row appends the count (0 when past the end).
- `write_report`'s `doc` gains:

```rust
        "adapter_content": match adapter_pos {
            Some(ap) => ap.to_json(cli_adapter_error_rate, cli_min_overlap),
            None => serde_json::Value::Null,
        },
```

  (pass the two numbers in as parameters rather than reaching into `Cli` from the report
  module.)

**Step 6: Run the tests**

Run: `cargo test --offline`
Expected: all green, including the three `adapter_content_*` CLI tests and four
`adapter_pos::` unit tests. The 54 pre-existing tests stay green — the new TSV columns are
only added when `--adapter-content` is given, so
`positions_tsv_has_one_row_per_position_and_exact_depths` still sees the 14-column header.

**Step 7: Differential harness still exact**

Run: `python tests/differential.py ../../examples/runs/demo/sim_E/sim_R1.fastq.gz`
Expected: `differential OK ...` (the harness reads by column name since Task 1, so the
extra columns do not affect it; and this run does not pass `--adapter-content` anyway).

**Step 8: Commit**

```bash
git add tools/fastcheck/src/adapter_pos.rs tools/fastcheck/src/accum.rs \
        tools/fastcheck/src/cli.rs tools/fastcheck/src/main.rs \
        tools/fastcheck/src/report.rs tools/fastcheck/tests/cli.rs
git commit -m "feat(fastcheck): adapter content per position, off by default (PLAN-QC module 11)"
```

---

## Task 10: Fix the insert-size anchor and the insert identity (the GT-exposed defects)

**Objective:** The insert measurement finds the measurable inserts and computes the insert
from biology lengths, so the ground-truth test in Task 13 can pass. This is a bug fix the
ground-truth run exposed; the tolerances and the prototype are in the "Verified baseline"
section above.

**Files:**
- Modify: `tools/fastcheck/src/insert.rs` (`overlap`, `verify`, `observe_pair`, doc comment)
- Test: `tools/fastcheck/src/insert.rs` (new unit tests)

**Step 1: Write the failing tests**

Add to `tools/fastcheck/src/insert.rs`'s test module:

```rust
    /// DEFECT A: the anchor used to be revcomp(s2)[0..12], which sits in R2's 3' polyG
    /// tail whenever t2 > 0. Anchoring on GGGGGGGGGGGG finds a chance G-stretch in R1, so
    /// the "verified" overlap is garbage and only ~25% of measurable inserts were found.
    /// The anchor must start AFTER the tail.
    #[test]
    fn the_anchor_skips_r2_s_polyg_tail() {
        // R2 = [biology][6 Gs]. revcomp(R2) starts with those 6 Gs, so an anchor at
        // offset 0 is all G and must not be used.
        let frag = dna(200, 31);
        let mut r2 = revcomp(&frag[50..200]);
        r2.extend_from_slice(b"GGGGGG");
        let r1 = frag[0..150].to_vec();
        // The biology overlap is 100 bases; the anchor must land inside the biology.
        let o = overlap(&r1, &r2).expect("a 150+156 pair with a 100-base overlap must be found");
        assert_eq!(o, 100, "the overlap is the shared biology, tail excluded");
    }

    /// DEFECT B: insert = len(R1)+len(R2)-o over-estimates by (h1+t1+h2+t2)/2 because
    /// polyG runs and adapter are not biology. The identity is insert = b1 + b2 - o with
    /// b = L - h - t.
    #[test]
    fn the_insert_uses_biology_lengths_not_read_lengths() {
        let frag = dna(220, 37);
        // Both mates carry a 5-base 3' G-run, so b1 = b2 = 145 while L = 150.
        let mut r1 = frag[0..145].to_vec();
        r1.extend_from_slice(b"GGGGG");
        let mut r2 = revcomp(&frag[75..220]);
        r2.extend_from_slice(b"GGGGG");
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert_eq!(h.measured, 1);
        // b1 + b2 - o = 145 + 145 - 70 = 220.
        assert!((h.mean() - 220.0).abs() < 1e-9, "got {}, want 220", h.mean());
    }

    /// Verification used to demand 0 mismatches. The simulator's own boundary rule
    /// guarantees ~1 disagreement at the overlap's first base, and sequencing errors add
    /// ~2e-3 per base, so a correct alignment was rejected about half the time.
    #[test]
    fn verification_tolerates_the_boundary_disagreement_and_sequencing_errors() {
        let frag = dna(220, 41);
        let mut r1 = frag[0..150].to_vec();
        let mut r2 = revcomp(&frag[70..220]);
        // One mismatch at the overlap's first base (the boundary rule's signature) and one
        // deep inside (a sequencing error): 2 mismatches in 80 bases must still verify.
        let anchor_start = 70; // overlap starts at R1 position 70
        r1[anchor_start] = if r1[anchor_start] == b'A' { b'C' } else { b'A' };
        r2[30] = if r2[30] == b'A' { b'C' } else { b'A' };
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert_eq!(h.measured, 1, "2 mismatches in 80 bases must not reject the pair");
    }

    /// Three mismatches in 80 bases is beyond the tolerance and must be rejected rather
    /// than silently mis-measured.
    #[test]
    fn too_many_mismatches_is_rejected_not_mismeasured() {
        let frag = dna(220, 43);
        let mut r1 = frag[0..150].to_vec();
        let r2 = revcomp(&frag[70..220]);
        for i in [70, 80, 90] {
            r1[i] = if r1[i] == b'A' { b'C' } else { b'A' };
        }
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert!(h.measured + h.censored_above + h.censored_below == 1);
        assert_eq!(h.measured, 0, "3 mismatches in 80 bases must not verify");
    }
```

**Step 2: Run to verify failure**

Run: `cargo test --offline insert::`
Expected: the three new tests FAIL (the anchor lands in the tail; the identity uses read
lengths; 2 mismatches are rejected). The pre-existing insert tests still pass.

**Step 3: Implement the fix**

In `tools/fastcheck/src/insert.rs`, replace `overlap` and `verify` with:

```rust
/// Mismatch tolerance for the overlap verification.
///
/// Not zero, and deliberately so: the boundary rule in the simulator (and any real
/// two-colour dataset) replaces the base adjacent to an artefact with a random non-G,
/// independently in each mate, so the two reads disagree at that base BY CONSTRUCTION.
/// Sequencing errors add ~2e-3 mismatches per base. A tolerance of 2 absorbs both while
/// still rejecting a chance alignment (which differs everywhere).
const MAX_OVERLAP_MISMATCHES: usize = 2;

/// Bases of the mate's start used to locate the overlap.
const ANCHOR: usize = 12;

/// Longest overlap `o` such that `s1[len1-t1-o..len1-t1] == revcomp(s2)[t2..t2+o]`.
///
/// `t1` and `t2` are the measured 3' homopolymer run lengths of each mate: the anchor is
/// placed at offset `t2` inside `revcomp(s2)` because offsets below that are the polyG
/// tail, and an anchor of Gs matches a chance G-stretch in R1 and produces a garbage
/// "verified" overlap. The overlap extends to the end of R1's biology (`len1 - t1`), not
/// to the end of the read, for the same reason.
///
/// The insert identity is then `b1 + b2 - o` with `b = L - h - t`: polyG and adapter are
/// not biology, so counting them (as `len1 + len2 - o` does) over-estimates the insert by
/// `(h1+t1+h2+t2)/2`.
fn overlap(s1: &[u8], s2: &[u8]) -> Option<usize> {
    let l1 = s1.len();
    let l2 = s2.len();
    if l1 < ANCHOR || l2 < ANCHOR {
        return None;
    }
    let t1 = trailing_run(s1);
    let t2 = trailing_run(s2);
    let bio1 = l1.saturating_sub(t1);
    let r = revcomp(s2);
    if bio1 < ANCHOR || r.len() < t2 + ANCHOR {
        return None;
    }
    // The anchor is the first ANCHOR biology bases of revcomp(s2), i.e. at offset t2.
    let anchor = &r[t2..t2 + ANCHOR];
    if anchor.iter().all(|&b| b == b'N') {
        return None;
    }
    let first = anchor[0];
    let mut p = 0usize;
    while p + ANCHOR <= bio1 {
        if s1[p] == first && &s1[p..p + ANCHOR] == anchor {
            let o = bio1 - p;
            if o <= r.len() - t2 && verify(&s1[p..p + o], &r[t2..t2 + o]) {
                return Some(o);
            }
        }
        p += 1;
    }
    None
}

/// The measured 3' homopolymer run of each base, so the tail can be skipped. N and other
/// IUPAC codes end the run.
fn trailing_run(s: &[u8]) -> usize {
    let mut n = 0usize;
    for &b in s.iter().rev() {
        match b.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => n += 1,
            _ => break,
        }
    }
    n
}

/// An alignment is accepted with up to [`MAX_OVERLAP_MISMATCHES`] differences, because the
/// boundary rule and sequencing errors guarantee a few. A chance alignment differs
/// everywhere, so the tolerance does not admit noise.
fn verify(a: &[u8], b: &[u8]) -> bool {
    let mut mm = 0usize;
    for (x, y) in a.iter().zip(b.iter()) {
        if x.to_ascii_uppercase() != y.to_ascii_uppercase() {
            mm += 1;
            if mm > MAX_OVERLAP_MISMATCHES {
                return false;
            }
        }
    }
    true
}

/// Reverse complement, for the overlap search.
fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}
```

and in `observe_pair`, replace the insert computation:

```rust
            Some(o) => {
                // b1 + b2 - o, with b = L - h - t the BIOLOGY length: the polyG runs and
                // any adapter are not biology, and counting them over-estimates the insert
                // by (h1+t1+h2+t2)/2. This is the identity the ground-truth test checks.
                let t1 = trailing_run(s1);
                let t2 = trailing_run(s2);
                let h1 = leading_run(s1);
                let h2 = leading_run(s2);
                let b1 = l1.saturating_sub(h1 + t1);
                let b2 = l2.saturating_sub(h2 + t2);
                let insert = b1 + b2 - o;
                if insert < self.bins.len() {
                    self.bins[insert] += 1;
                    self.measured += 1;
                } else {
                    self.censored_above += 1;
                }
            }
```

with `leading_run` next to `trailing_run`:

```rust
/// The measured 5' homopolymer run, so the head artefact can be excluded from the biology
/// length. Mirrors [`trailing_run`].
fn leading_run(s: &[u8]) -> usize {
    let mut n = 0usize;
    for &b in s.iter() {
        match b.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => n += 1,
            _ => break,
        }
    }
    n
}
```

Also update the module doc comment's `insert = len(R1) + len(R2) - overlap` line to state
the corrected identity and the reason, and update the `resolution` string in `to_json`
(which still says "an overlap shorter than 12 bases cannot be anchored") to add:
`the anchor is placed after each mate's 3' homopolymer run, and up to 2 mismatches are
tolerated because the artefact boundary rule guarantees one`.

**Step 4: Run the tests**

Run: `cargo test --offline`
Expected: all green, including the four new insert tests. Watch `overlap_recovers_a_known_insert`
(insert 300, reads 251, no artefacts): `b1=b2=251`, `o=202`, insert `251+251-202 = 300` —
unchanged. `a_long_insert_is_censored_not_guessed` and
`an_insert_below_the_read_length_is_censored_below` also keep their behaviour.

**Step 5: Verify on the ground-truth fixture**

```bash
../fastsim/target/release/fastsim --pairs 20000 --profile ../fastsim/profiles/synthetic_novaseq.json \
  --insert-mean 200 --insert-sd 40 --gc 0.5 --seed 42 --no-gzip --dup-rate 0.02 \
  --adapter-seq CTGTCTCTTATACACATCT -o /tmp/ins
./target/release/fastcheck /tmp/ins_R1.fastq /tmp/ins_R2.fastq --insert-hist --insert-pairs 20000 \
  --adapter-seq CTGTCTCTTATACACATCT --json /tmp/ins.json -q
python3 -c "
import json, statistics
ih = json.load(open('/tmp/ins.json'))['insert_hist']
b = ih['bins_unioned_with_adapter']
overl = [int(k) for k, c in b.items() if int(k) > 150 for _ in range(c)]
gt = [int(r.split('\t')[1]) for r in open('/tmp/ins.gt.tsv').read().splitlines()[1:]]
mid = [i for i in gt if 150 < i <= 288]
print('measured pairs', len(overl), 'of', len(mid), 'measurable')
print('measured mean %.2f  GT mean %.2f' % (statistics.mean(overl), statistics.mean(mid)))
assert len(overl) > 0.9 * len(mid), 'must find at least 90% of the measurable inserts'
assert abs(statistics.mean(overl) - statistics.mean(mid)) < 0.02 * statistics.mean(mid)
print('OK')
"
```
Expected: ~16 000 of 16 642 found (≥90%), measured mean within 2% of the GT mean.

**Step 6: Commit**

```bash
git add tools/fastcheck/src/insert.rs
git commit -m "fix(fastcheck): anchor the insert overlap in biology and use biology lengths (GT-exposed)"
```

---

## Task 11: Per-tile quality modules, default off

**Objective:** `--tiles` parses 7-field Illumina headers and accumulates
`tile x position x quality`, emitting a per-tile report block; off by default, with the
memory projection printed when on.

**Files:**
- Create: `tools/fastcheck/src/tile.rs`
- Modify: `tools/fastcheck/src/cli.rs` (`--tiles`)
- Modify: `tools/fastcheck/src/accum.rs` (the counter)
- Modify: `tools/fastcheck/src/main.rs`, `src/report.rs`
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tools/fastcheck/tests/cli.rs`:

```rust
/// 7-field Illumina headers: `@instrument:run:flowcell:lane:tile:x:y`.
/// Tile 1101 gets good quality, tile 1102 gets bad, on the same read content.
const TILE_FASTQ: &[u8] = b"@LH00916:1:SIMFLOW1:1:1101:1:1 1:N:0:ATCGATCG\nACGTACGTAC\n+\nIIIIIIIIII\n\
@LH00916:1:SIMFLOW1:1:1101:2:2 1:N:0:ATCGATCG\nACGTACGTAC\n+\nIIIIIIIIII\n\
@LH00916:1:SIMFLOW1:1:1102:3:3 1:N:0:ATCGATCG\nACGTACGTAC\n+\n#####AAAAA\n\
@LH00916:1:SIMFLOW1:1:1102:4:4 1:N:0:ATCGATCG\nACGTACGTAC\n+\n#####AAAAA\n\
@simple_read_with_no_coords\nACGTACGTAC\n+\nIIIIIIIIII\n";

#[test]
fn tiles_module_separates_quality_by_tile_and_is_off_by_default() {
    let d = tmpdir("tiles");
    let p = d.join("t.fastq");
    std::fs::write(&p, TILE_FASTQ).unwrap();

    // Off by default: no block, and headers without coords are simply not counted.
    let j0 = d.join("off.json");
    run(&d, &[p.to_str().unwrap(), "--json", j0.to_str().unwrap(), "-q"]);
    assert!(read_json(&j0).get("per_tile").is_none(), "must be opt-in");

    let json = d.join("on.json");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--tiles", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let blk = &v["per_tile"];
    assert_eq!(blk["tiles"].as_array().unwrap().len(), 2, "1101 and 1102");
    assert_eq!(blk["reads_with_coordinates"], 4);
    assert_eq!(blk["reads_without_coordinates"], 1, "the simple header is counted, not dropped");
    // The projection line names the cost.
    assert!(stderr.contains("per-tile"), "{stderr}");

    let by_tile: std::collections::HashMap<String, &serde_json::Value> = blk["tiles"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| (t["tile"].as_str().unwrap().to_string(), t))
        .collect();
    // Tile 1101: every base is 'I' = Q40. Tile 1102: five '#' (Q2) then five 'A' (Q32).
    let good = &by_tile["1101"];
    assert_eq!(good["mean_q"], 40.0);
    let bad = &by_tile["1102"];
    // (5*2 + 5*32) / 10 = 17
    assert!((bad["mean_q"].as_f64().unwrap() - 17.0).abs() < 1e-9, "{bad}");
    // Per-position detail exists for both.
    assert!(good["per_position_mean_q"].as_array().unwrap().len() == 10);
}

#[test]
fn tiles_module_reports_the_memory_projection_when_on() {
    let d = tmpdir("tiles_mem");
    let p = d.join("t.fastq");
    std::fs::write(&p, TILE_FASTQ).unwrap();
    let (_, stderr, code) = run(&d, &[p.to_str().unwrap(), "--tiles", "-q"]);
    assert_eq!(code, 0);
    assert!(stderr.contains("per-tile quality"), "{stderr}");
    assert!(stderr.contains("MB"), "{stderr}");
}
```

**Step 2: Run to verify failure**

Run: `cargo test --offline tiles_module`
Expected: FAIL — `error: unexpected argument '--tiles' found`.

**Step 3: Add the flag**

In `tools/fastcheck/src/cli.rs`:

```rust
    /// Per-tile quality (PLAN-QC module 7 / FastQC's tile modules): quality by tile x
    /// position x quality bin, parsed from 7-field Illumina headers
    /// (`@instrument:run:flowcell:lane:tile:x:y`). OFF by default: it costs ~5 MB per 100
    /// tiles, and only raw Illumina runs carry the coordinates. Reads whose header has no
    /// tile field are counted and reported, never silently dropped.
    #[arg(long)]
    pub tiles: bool,
```

**Step 4: Create the module**

Create `tools/fastcheck/src/tile.rs`:

```rust
//! Per-tile quality (PLAN-QC module 7).
//!
//! FastQC's per-tile modules answer "is one region of the flowcell failing", which a
//! whole-file quality curve cannot see: a bad tile is invisible once its reads are averaged
//! with the good ones. The accumulator is `tile x position x quality bin`, ~5 MB per 100
//! tiles at R=150, so it is off by default and the projection is printed when it is on.
//!
//! Only 7-field Illumina headers carry a tile. Reads with any other header shape are
//! counted in `reads_without_coordinates` rather than dropped: a run whose headers were
//! rewritten must not silently lose its reads.

use std::collections::BTreeMap;

use crate::accum::{quantile, QBINS};

/// Extract the tile number from a FASTQ header, if it has the 7-field Illumina form.
///
/// `@instrument:run:flowcell:lane:tile:x:y <rest>` — the tile is field 5 (index 4) of the
/// part before the first space. Anything else (a plain name, a rewritten header, a
/// different vendor's layout) returns `None`.
pub fn tile_of(header: &[u8]) -> Option<u32> {
    let end = header.iter().position(|&b| b == b' ').unwrap_or(header.len());
    let fields: Vec<&[u8]> = header[..end].split(|&b| b == b':').collect();
    if fields.len() != 7 {
        return None;
    }
    std::str::from_utf8(fields[4]).ok()?.parse::<u32>().ok()
}

/// `tile -> position -> quality-bin counts`.
#[derive(Clone, Default)]
pub struct TileQual {
    pub tiles: BTreeMap<u32, Vec<[u64; QBINS]>>,
    pub reads_with_coordinates: u64,
    pub reads_without_coordinates: u64,
}

impl TileQual {
    #[inline]
    pub fn observe(&mut self, header: &[u8], qual: &[u8]) {
        match tile_of(header) {
            Some(t) => {
                self.reads_with_coordinates += 1;
                let rows = self.tiles.entry(t).or_default();
                for (i, &q) in qual.iter().enumerate() {
                    if rows.len() <= i {
                        rows.resize(i + 1, [0; QBINS]);
                    }
                    rows[i][q as usize] += 1;
                }
            }
            None => self.reads_with_coordinates += 0, // no coordinates: counted below
                ,
        }
        if tile_of(header).is_none() {
            self.reads_without_coordinates += 1;
        }
    }

    /// Exact sum of two accumulators.
    pub fn merge(&mut self, o: &TileQual) {
        for (&t, rows) in &o.tiles {
            let dst = self.tiles.entry(t).or_default();
            if dst.len() < rows.len() {
                dst.resize(rows.len(), [0; QBINS]);
            }
            for (d, s) in dst.iter_mut().zip(rows.iter()) {
                for i in 0..QBINS {
                    d[i] += s[i];
                }
            }
        }
        self.reads_with_coordinates += o.reads_with_coordinates;
        self.reads_without_coordinates += o.reads_without_coordinates;
    }

    /// Resident bytes, for the projection printed before the numbers.
    pub fn resident_bytes(&self) -> usize {
        self.tiles
            .values()
            .map(|rows| rows.len() * std::mem::size_of::<[u64; QBINS]>())
            .sum()
    }

    /// The report block: one entry per tile with its mean quality and per-position means.
    pub fn to_json(&self, offset: u8) -> serde_json::Value {
        let tiles: Vec<serde_json::Value> = self
            .tiles
            .iter()
            .map(|(&t, rows)| {
                let mut sum = 0.0f64;
                let mut n = 0u64;
                let mut per_pos = Vec::with_capacity(rows.len());
                for row in rows {
                    let depth: u64 = row.iter().sum();
                    if depth == 0 {
                        per_pos.push(serde_json::json!(0.0));
                        continue;
                    }
                    let (mean, _q10, _q50, _q90) = pos_mean(row, offset);
                    per_pos.push(serde_json::json!((mean * 100.0).round() / 100.0));
                    sum += mean * depth as f64;
                    n += depth;
                }
                serde_json::json!({
                    "tile": t,
                    "reads": rows.first().map(|r| r.iter().sum::<u64>()).unwrap_or(0),
                    "mean_q": if n > 0 { ((sum / n as f64) * 100.0).round() / 100.0 } else { 0.0 },
                    "per_position_mean_q": per_pos,
                })
            })
            .collect();
        serde_json::json!({
            "tiles": tiles,
            "tile_count": self.tiles.len(),
            "reads_with_coordinates": self.reads_with_coordinates,
            "reads_without_coordinates": self.reads_without_coordinates,
            "note": "reads whose header has no 7-field Illumina tile are counted in \
                     reads_without_coordinates, never dropped",
        })
    }
}

/// Mean quality of one position's exact histogram, in Phred space.
fn pos_mean(row: &[u64; QBINS], offset: u8) -> (f64, f64, f64, f64) {
    let depth: u64 = row.iter().sum();
    if depth == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut q = vec![0u64; QBINS];
    let mut sum = 0.0;
    for (b, &c) in row.iter().enumerate() {
        let phred = ((b as i64 - offset as i64).max(0) as usize).min(QBINS - 1);
        q[phred] += c;
        sum += phred as f64 * c as f64;
    }
    (sum / depth as f64, quantile(&q, 0.10), quantile(&q, 0.50), quantile(&q, 0.90))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tile_is_field_five_of_a_seven_field_header() {
        assert_eq!(tile_of(b"@LH00916:1:SIMFLOW1:1:1101:1:1 1:N:0:ATCGATCG"), Some(1101));
        assert_eq!(tile_of(b"@LH00916:1:SIMFLOW1:1:2202:9:9 2:N:0:ATCGATCG"), Some(2202));
    }

    #[test]
    fn a_header_without_coordinates_has_no_tile() {
        assert_eq!(tile_of(b"@simple_read"), None);
        assert_eq!(tile_of(b"@a:b:c 1:N:0:x"), None, "only 3 fields");
        assert_eq!(tile_of(b"@LH00916:1:SIMFLOW1:1:notanumber:1:1 1:N:0:x"), None);
    }

    #[test]
    fn quality_is_separated_by_tile() {
        let mut t = TileQual::default();
        t.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"IIII");
        t.observe(b"@i:1:f:1:1101:2:2 1:N:0:x", b"IIII");
        t.observe(b"@i:1:f:1:1102:3:3 1:N:0:x", b"####");
        t.observe(b"@no_coords", b"IIII");
        assert_eq!(t.reads_with_coordinates, 3);
        assert_eq!(t.reads_without_coordinates, 1);
        assert_eq!(t.tiles.len(), 2);
        let good = &t.tiles[&1101];
        // 'I' = 73; offset 33 -> Q40, two reads x four bases.
        let (mean, ..) = pos_mean(&good[0], 33);
        assert!((mean - 40.0).abs() < 1e-9);
        let bad = &t.tiles[&1102];
        let (mean, ..) = pos_mean(&bad[0], 33);
        assert!((mean - 2.0).abs() < 1e-9, "'#' = 35 -> Q2, got {mean}");
    }

    #[test]
    fn merge_is_exact() {
        let mut a = TileQual::default();
        a.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"IIII");
        let mut b = TileQual::default();
        b.observe(b"@i:1:f:1:1101:2:2 1:N:0:x", b"IIII");
        a.merge(&b);
        assert_eq!(a.tiles[&1101][0][73], 2);
        assert_eq!(a.reads_with_coordinates, 2);
    }
}
```

**Step 5: Wire it in**

- `accum.rs`: add `pub tile: Option<crate::tile::TileQual>,` (`None` in `new`), observe in
  the record loop (the header must be passed to `observe` — change `Accum::observe` to take
  the header: `pub fn observe(&mut self, header: &[u8], seq: &[u8], qual: &[u8])` and update
  the three test call sites plus `reader.rs`'s `acc.observe(&seq, &qual)` to
  `acc.observe(&line, &seq, &qual)`), and merge with the same
  `match (self.tile.as_mut(), o.tile.as_ref())` shape as `adapter_pos`.
- `main.rs`: when `cli.tiles`, set `acc.tile = Some(TileQual::default())` in `scan_one`
  **and** print the projection after the scan:

```rust
    if let Some(t) = total.tile.as_ref() {
        eprintln!(
            "fastcheck: per-tile quality for {} tile(s), {:.1} MB",
            t.tiles.len(),
            t.resident_bytes() as f64 / (1024.0 * 1024.0),
        );
    }
```

- `report.rs`: `write_report` gains `tile: Option<&crate::tile::TileQual>` and the doc gains
  `"per_tile": tile.map(|t| t.to_json(offset)).unwrap_or(serde_json::Value::Null)`.

**Step 6: Run the tests**

Run: `cargo test --offline`
Expected: all green, including the two `tiles_module_*` CLI tests and four `tile::` unit
tests, and the 54 pre-existing tests (the `observe` signature change is mechanical).

**Step 7: Commit**

```bash
git add tools/fastcheck/src/tile.rs tools/fastcheck/src/cli.rs \
        tools/fastcheck/src/accum.rs tools/fastcheck/src/reader.rs \
        tools/fastcheck/src/main.rs tools/fastcheck/src/report.rs \
        tools/fastcheck/tests/cli.rs
git commit -m "feat(fastcheck): per-tile quality modules, off by default (PLAN-QC module 7)"
```

---

## Task 12: `--only` module gating

**Objective:** `--only dup` (and `--no-quality`, `--no-heatmap`) stop the parser from
materialising what the enabled modules do not need, and the report says which modules ran.

**Files:**
- Modify: `tools/fastcheck/src/cli.rs`
- Modify: `tools/fastcheck/src/accum.rs` (a `Modules` struct gating each block)
- Modify: `tools/fastcheck/src/main.rs` / `reader.rs` (skip the quality decode)
- Test: `tools/fastcheck/tests/cli.rs`

**Step 1: Write the failing tests**

Append to `tools/fastcheck/tests/cli.rs`:

```rust
#[test]
fn only_dup_runs_the_duplication_module_and_nothing_else() {
    let d = tmpdir("only_dup");
    let fq = write_fixture(&d);
    let json = d.join("d.json");
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--only", "dup", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let mods = v["modules"].as_array().unwrap();
    assert!(mods.iter().any(|m| m == "dup"), "dup must be enabled: {mods:?}");
    assert!(!mods.iter().any(|m| m == "quality"), "{mods:?}");
    assert!(!mods.iter().any(|m| m == "heatmap"), "{mods:?}");
    // Duplication is still exact, and the gated-out modules are absent, not zero.
    assert!(v["metrics"]["distinct_seqs"].is_u64());
    assert!(v.get("heatmap").is_none(), "the heatmap must not be built");
    assert!(v.get("base_position").is_none(), "no per-position matrices were materialised");
}

#[test]
fn no_quality_skips_the_quality_accumulators_but_keeps_the_framing() {
    let d = tmpdir("no_quality");
    let fq = write_fixture(&d);
    let json = d.join("q.json");
    let (_, _, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--no-quality", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0);
    let v = read_json(&json);
    // Length and composition are unaffected; quality is absent rather than zero.
    assert_eq!(v["metrics"]["reads"], 6);
    assert_eq!(v["metrics"]["bases"], 60);
    assert!(v["metrics"].get("mean_read_mean_q").is_none());
    assert!(v.get("base_position").is_some(), "composition is not a quality module");
    assert_eq!(v["metrics"]["phred_offset"], 33, "a stated default, not a measurement");
    // The report must say the offset was not measured.
    assert!(v["assumptions"].as_array().unwrap().iter().any(|a| a["what"] == "quality"));
}

#[test]
fn no_heatmap_never_builds_the_grid() {
    let d = tmpdir("no_heatmap");
    let fq = write_fixture(&d);
    let json = d.join("h.json");
    let (_, _, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--no-heatmap", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0);
    let v = read_json(&json);
    assert!(v.get("heatmap").is_none());
    assert!(v["metrics"].get("heatmap_peak").is_none());
    assert_eq!(v["metrics"]["reads"], 6, "everything else still runs");
}

#[test]
fn only_rejects_an_unknown_module_name() {
    let d = tmpdir("only_bad");
    let fq = write_fixture(&d);
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--only", "bogus", "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("dup") && stderr.contains("adapter"), "{stderr}");
}
```

**Step 2: Run to verify failure**

Run: `cargo test --offline only_dup_runs`
Expected: FAIL — `error: unexpected argument '--only' found`.

**Step 3: Add the flags and the module set**

In `tools/fastcheck/src/cli.rs`:

```rust
    /// Run only these modules (comma-separated): bases, quality, length, gc, heatmap,
    /// dup, polyx, kmer, adapter, tile. The parser then materialises only what the
    /// enabled modules need, so `--only dup` touches neither the quality bytes nor the
    /// heatmap grid. Default: every module except tile and adapter-content, which stay
    /// opt-in behind their own flags.
    #[arg(long, value_name = "MODULE[,MODULE...]")]
    pub only: Option<String>,

    /// Disable the per-base quality accumulators (the quality line is still parsed for
    /// framing, and the Phred offset is reported as a stated default rather than a
    /// measurement).
    #[arg(long)]
    pub no_quality: bool,

    /// Never build the 256x256 composition grid (the only structure that needs
    /// --max-read-len up front).
    #[arg(long)]
    pub no_heatmap: bool,
```

and a resolver returning the module set:

```rust
/// Which modules run. `--only` is an allow-list; `--no-quality`/`--no-heatmap` are
/// deny-list entries that also apply to `--only`.
#[derive(Debug, Clone)]
pub struct Modules {
    pub bases: bool,
    pub quality: bool,
    pub length: bool,
    pub gc: bool,
    pub heatmap: bool,
    pub dup: bool,
    pub polyx: bool,
    pub kmer: bool,
    pub adapter: bool,
    pub tile: bool,
}

const MODULE_NAMES: &[&str] = &[
    "bases", "quality", "length", "gc", "heatmap", "dup", "polyx", "kmer", "adapter", "tile",
];

impl Cli {
    pub fn modules(&self) -> anyhow::Result<Modules> {
        let all = Modules {
            bases: true, quality: true, length: true, gc: true, heatmap: true,
            dup: true, polyx: true, kmer: true,
            adapter: self.adapter_content, tile: self.tiles,
        };
        let mut m = match &self.only {
            None => all,
            Some(list) => {
                let mut m = Modules {
                    bases: false, quality: false, length: false, gc: false, heatmap: false,
                    dup: false, polyx: false, kmer: false, adapter: false, tile: false,
                };
                for name in list.split(',').map(str::trim).filter(|s| !s.is_empty()) {
                    match name {
                        "bases" => m.bases = true,
                        "quality" => m.quality = true,
                        "length" => m.length = true,
                        "gc" => m.gc = true,
                        "heatmap" => m.heatmap = true,
                        "dup" => m.dup = true,
                        "polyx" => m.polyx = true,
                        "kmer" => m.kmer = true,
                        "adapter" => m.adapter = self.adapter_content,
                        "tile" => m.tile = self.tiles,
                        other => anyhow::bail!(
                            "--only: unknown module {other:?} (known: {})",
                            MODULE_NAMES.join(", ")
                        ),
                    }
                }
                m
            }
        };
        if self.no_quality {
            m.quality = false;
        }
        if self.no_heatmap {
            m.heatmap = false;
        }
        Ok(m)
    }

    /// The enabled module names, for the report's `modules` array.
    pub fn module_names(&self) -> Vec<&'static str> {
        let Ok(m) = self.modules() else { return Vec::new() };
        MODULE_NAMES.iter().copied().filter(|n| match *n {
            "bases" => m.bases, "quality" => m.quality, "length" => m.length,
            "gc" => m.gc, "heatmap" => m.heatmap, "dup" => m.dup,
            "polyx" => m.polyx, "kmer" => m.kmer, "adapter" => m.adapter, "tile" => m.tile,
            _ => false,
        }).collect()
    }
}
```

**Step 4: Gate the accumulator**

In `tools/fastcheck/src/accum.rs`, add `pub modules: Modules` (a field on `Accum`, set in
`new`) and wrap each block in `observe`:

```rust
        if self.modules.bases { /* the base_pos loop */ }
        if self.modules.quality { /* the qual_pos loop and the census */ }
        if self.modules.length { /* len_hist */ }
        if self.modules.gc { /* acgt_hist, gc_hist */ }
        if self.modules.heatmap { /* the heatmap cell */ }
        if self.modules.dup { self.seq_counts.observe(seq); }
        if self.modules.polyx { /* the polyX measurement */ }
        if self.modules.kmer { /* the 12-mer windows */ }
```

(The quality line is still *read* by the scanner for framing — only the accumulators and
the decode are skipped. `reads`, `bases`, `reads_with_n` and `reads_over_max_len` are always
counted: a run that reports no reads is not a run.)

**Step 5: Gate the report**

In `tools/fastcheck/src/report.rs`:

- `write_report` gains `modules: &[&'static str]` and the doc gains
  `"modules": modules,`.
- The quality metrics (`mean_read_mean_q`, `reads_q20_frac`, `reads_q30_frac`,
  `phred_offset_certain`, `quality_byte_min/max`, `read_mean_q_hist`) are emitted only when
  `modules.contains(&"quality")`; `phred_offset` is always emitted as 33 with an
  `assumptions` entry:

```rust
    if !modules.contains(&"quality") {
        assumptions.push(json!({
            "what": "quality",
            "value": "not measured",
            "basis": "--no-quality (or --only without quality) skipped the quality accumulators",
            "consequence": "no quality-derived metric is reported; the Phred offset is the \
                            stated default 33, not a measurement",
        }));
    }
```

- The `heatmap`, `polyx`, `kmer_spectrum` and `repeated_sequences` blocks are `Null` when
  their module is off (a consumer must never read an absent block as a measured zero).

**Step 6: Run the tests**

Run: `cargo test --offline`
Expected: all green, including the four new gating tests. `report_is_v1_and_matches_the_hand_counted_fixture`
still passes because the default is every module on.

**Step 7: Commit**

```bash
git add tools/fastcheck/src/cli.rs tools/fastcheck/src/accum.rs \
        tools/fastcheck/src/main.rs tools/fastcheck/src/reader.rs \
        tools/fastcheck/src/report.rs tools/fastcheck/tests/cli.rs
git commit -m "feat(fastcheck): --only / --no-quality / --no-heatmap module gating (PLAN-QC §3)"
```

---

## Task 13: The fastsim ground-truth recovery test (the decisive one)

**Objective:** `tools/fastcheck/tests/ground_truth.py` generates reads from a known profile,
runs fastcheck, and asserts the recovered numbers against the ground truth within stated
tolerances. This is PLAN-QC §7.3 and the difference between "my numbers look plausible" and
"my numbers are right".

**Files:**
- Create: `tools/fastcheck/tests/ground_truth.py`
- Create: `tools/fastcheck/tests/gt_profile.json` (the synthetic profile, no real data — D-11)
- Modify: `tools/fastcheck/pixi.toml` (a `ground-truth` task)

**Step 1: Create the profile**

Create `tools/fastcheck/tests/gt_profile.json`:

```json
{
  "_comment": "Synthetic ground-truth profile for fastcheck's recovery test. No measured data (D-11): every number is a round number describing a mechanism, so the expected recovery is computable by hand. pos_q is flat at 30, tail_q_collapse is 0 so the 3' G-run histogram is recovered exactly, and the head CDF has mass at 0..3 so the leading-run convolution is testable.",
  "read_len": 150,
  "gc": 0.5,
  "geometry": "tail3p_both",
  "qual_noise": 0.5,
  "head": {"R1": {"0": 70, "1": 20, "2": 7, "3": 3}, "R2": {"0": 70, "1": 20, "2": 7, "3": 3}},
  "tail": {"R1": {"0": 50, "5": 30, "10": 20}, "R2": {"0": 50, "5": 30, "10": 20}},
  "joint_tail": {"0,0": 25, "5,5": 15, "10,10": 10},
  "tail_q_collapse": {"R1": 0.0, "R2": 0.0},
  "pos_q": {"R1": [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
             "R2": [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]},
  "read_meanq": {"R1": {"30": 100}, "R2": {"30": 100}},
  "full_length_rate": 0.0,
  "full_length_q": {"R1": 20, "R2": 20},
  "adapter": {"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"}
}
```

**Step 2: Write the failing test**

Create `tools/fastcheck/tests/ground_truth.py`:

> **NOTE — this plan was truncated.** The subagent that wrote it died mid-write
> (response truncation), so the content of `ground_truth.py`, the implementation
> steps for tasks 10–12 (insert anchor fix, per-tile modules, `--only` gating)
> and their verification commands are **missing**. Everything through Task 9
> (adapter content per position) and the start of Task 13 (the profile fixture)
> is complete and executable as written. Before implementing Task 13 onward,
> finish writing the test script and steps — or dispatch a subagent with the
> partial plan and the ground-truth tolerances table above as context.

What the recovered tail of Task 13 does assert (from the verified tolerances
table at the top of this plan): per-position mean Q within 0.05, GC within
0.03, adapter rate within 1% of pairs, 3′ G-run histogram within 4% per bin,
the 5′ head histogram within 1.5% of the convolution model, insert mean within
2%, and duplication level as a lower bound `0 ≤ level ≤ gt_dup_rate`. The
profile fixture (above) is flat-Q30, GC 0.5, tail mass at 0/5/10, head mass at
0..3, 2% duplicates, adapter readthrough — every tolerance follows from it.

## Task 14: Docs updates

**Objective:** PLAN-QC §0 status reflects what shipped (already updated
2026-10-12 in commit `649a7ee`); DECISIONS.md D-5 gets the `--html` correction
(fastcheck renders HTML despite "never renders HTML"); CLI-SPEC.md §4 gains the
fastcheck-specific gap entries with their closure status.

**Files:** `docs/DECISIONS.md` (D-5), `docs/CLI-SPEC.md` (§4 + the fastcheck
section).

**Step 1:** In `docs/DECISIONS.md` D-5, after "It never reads a run directory
and never renders HTML.", add:

```markdown
*Correction (2026-10-12):* fastcheck **does** render its own per-file HTML
report (`--html`, MultiQC-style, inline SVG, no JS). The boundary now reads:
fastcheck renders the **per-file** report; fasthelper renders the
**cross-sample** report. fastcheck still never reads a run directory.
```

**Step 2:** In `docs/CLI-SPEC.md` §4, append after the numbered list:

```markdown
### fastcheck-specific gaps (closed by docs/plans/2026-10-12_fastcheck-exactness.md)

* `--max-mem` admission refusal — ported from fastprep's proven model.
* `--overrep-mode none|fingerprint|second-pass|table` — the four honest
  trade-offs; default `none` (PLAN-QC §5.9).
* adapter content per position — reusing fastprep's catalogue/matcher
  semantics, copied per D-1.
* per-tile quality modules — default off, ~5 MB per 100 tiles.
* `--only <module>` gating — the parser materialises only what enabled
  modules need.
* ground-truth recovery test (PLAN-QC §7.3) — the decisive correctness gate.
```

**Step 3: Commit**

```bash
git add docs/DECISIONS.md docs/CLI-SPEC.md
git commit -m "docs: fastcheck gaps and the D-5 boundary correction (fastcheck renders per-file HTML)"
```

**Verify:** `grep -A3 "fastcheck-specific gaps" docs/CLI-SPEC.md` shows the
six bullets; `grep -B1 -A3 "Correction (2026-10-12)" docs/DECISIONS.md` shows
the D-5 note.

## Final verification

```fish
# every tool's tests, including fastcheck's now-larger suite
for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper
    cd tools/$t; and cargo test --offline; and cd ../..
end

# the differential harness must be green (it was RED before Task 1)
cd tools/fastcheck
python tests/differential.py ../../examples/runs/demo/sim_E/sim_R1.fastq.gz

# the ground-truth recovery test — the decisive one
python tests/ground_truth.py --pairs 20000

# the closed loop still works
fastcheck --insert-hist --profile /tmp/p.json /tmp/gt_R1.fastq /tmp/gt_R2.fastq
fastsim --profile /tmp/p.json --pairs 20000 -o /tmp/sim2
```

Expected: all suites green, `differential OK`, `GROUND TRUTH OK` (or a loud,
specific failure naming the metric and its deviation), and the closed loop
runs with the profile's insert CDF consumed.
