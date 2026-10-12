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

/// cgroup v2. The limit is frequently set on a parent slice rather than on the leaf, so the
/// whole ancestor chain is walked and the minimum taken.
fn cgroup_v2(root: &Path, phys: u64) -> Option<Limit> {
    // `/proc/self/cgroup` under v2 has a single line of the form `0::/path`.
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

/// cgroup v1. The controller path comes from the `:memory:` line.
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

    // The sentinel means unlimited; so, in practice, does anything at or above physical
    // memory.
    (bytes != V1_UNLIMITED && bytes < phys)
        .then(|| Limit { bytes, source: LimitSource::CgroupV1(rel) })
}

/// The scheduler's grant. Worth reading even under cgroup enforcement: some sites set it
/// without cgroups.
fn slurm(env: &[(String, String)]) -> Option<Limit> {
    let get = |k: &str| env.iter().find(|(n, _)| n == k).map(|(_, v)| v.as_str());

    if let Some(mb) = get("SLURM_MEM_PER_NODE").and_then(|v| v.parse::<u64>().ok()) {
        return Some(Limit {
            bytes: mb * 1024 * 1024,
            source: LimitSource::Slurm("SLURM_MEM_PER_NODE"),
        });
    }
    let per_cpu = get("SLURM_MEM_PER_CPU")?.parse::<u64>().ok()?;
    let cpus = get("SLURM_CPUS_ON_NODE")?.parse::<u64>().ok()?;
    Some(Limit {
        bytes: per_cpu * cpus * 1024 * 1024,
        source: LimitSource::Slurm("SLURM_MEM_PER_CPU x SLURM_CPUS_ON_NODE"),
    })
}

/// `RLIMIT_AS` bounds *virtual address space*, not resident set, so it is respected only
/// when set below physical memory and is deliberately the weakest signal.
#[cfg(target_os = "linux")]
fn rlimit_as(phys: u64) -> Option<Limit> {
    let mut rl = libc::rlimit { rlim_cur: 0, rlim_max: 0 };
    // SAFETY: `rl` is a valid, fully initialised `rlimit`; getrlimit only writes through the
    // pointer.
    if unsafe { libc::getrlimit(libc::RLIMIT_AS, &mut rl) } != 0 {
        return None;
    }
    let cur = rl.rlim_cur as u64;
    (rl.rlim_cur != libc::RLIM_INFINITY && cur < phys)
        .then_some(Limit { bytes: cur, source: LimitSource::RlimitAs })
}

#[cfg(not(target_os = "linux"))]
fn rlimit_as(_phys: u64) -> Option<Limit> { None }

/// `MemAvailable`, not `MemTotal`: it accounts for reclaimable page cache and is the honest
/// answer to "how much can this process actually take".
fn mem_available(root: &Path) -> Option<u64> {
    meminfo_field(root, "MemAvailable:")
}

fn mem_total(root: &Path) -> Option<u64> {
    meminfo_field(root, "MemTotal:")
}

fn meminfo_field(root: &Path, key: &str) -> Option<u64> {
    read(root.join("proc/meminfo"))?
        .lines()
        .find(|l| l.starts_with(key))
        .and_then(|l| l.split_whitespace().nth(1)?.parse::<u64>().ok())
        .map(|kb| kb * 1024)
}

fn read(p: impl AsRef<Path>) -> Option<String> {
    std::fs::read_to_string(p).ok()
}

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

pub fn project(capacity: usize, fingerprint_bytes: usize) -> Projection {
    if capacity == 0 {
        return Projection { buckets: 0, bytes: 0 };
    }
    let need = (capacity as u64 * 8).div_ceil(7);
    let buckets = need.next_power_of_two();
    Projection { buckets, bytes: buckets * (fingerprint_bytes as u64 + 1) }
}

/// Largest capacity that still fits in `2^k` buckets, i.e. `7 * 2^k / 8`.
///
/// Bucket counts double at these points, **not** at powers of two: 5e8 sits 6% above the
/// 2^29 boundary and so projects 9216 MB where 4608 MB would do.
pub const fn capacity_ceiling(pow2_buckets: u32) -> usize { (7usize << pow2_buckets) / 8 }

/// Replacement for a round-number preallocation clamp: the largest capacity that fits in
/// 2^29 buckets, i.e. 4608 MB at 64-bit fingerprints.
pub const MAX_PREALLOC_ENTRIES: usize = capacity_ceiling(29); // 469_762_048

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
pub fn reserve_for(limit: u64) -> u64 { (512 * 1024 * 1024).max(limit / 10) }

/// Admission rule. `limit: None` means undiscoverable, which degrades to
/// [`Admission::Proceed`]: refusing on an unknown bound would make the tool unusable
/// wherever discovery fails.
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
