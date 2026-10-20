use crate::hasher::HashType;
use anyhow::{Context, Result, bail};
use std::fs;
use std::path::{Path, PathBuf};

/// Maximum number of entries to pre-allocate in the hash set.
/// Prevents multi-GB allocations from inaccurate file-size estimates.
const MAX_PREALLOC_ENTRIES: usize = 500_000_000; // 500M

pub fn birthday_problem_square_approximation(x: usize, n: &HashType) -> f64 {
    (x as f64).powi(2) / 2.0_f64.powi(n.to_num() as i32 + 1)
}

/// Returns the appropriate `HashType` given the dataset size and collision probability threshold.
pub fn get_hash_method(size: usize, threshold: f64) -> HashType {
    if (2f64 * 2.0f64.powi(64) * threshold).sqrt() < size as f64 {
        HashType::XXH3_128
    } else {
        HashType::XXH3_64
    }
}

/// Computes the write buffer size from the expected read length.
///
/// A FASTQ record is roughly:
///   - 1 ID line:      ~40 bytes (conservative estimate)
///   - 1 sequence:     read_length bytes
///   - 1 '+' line:     2 bytes
///   - 1 quality line: read_length bytes
///
/// Total ≈ `2 * read_length + 42` bytes.
///
/// We target holding ~512 records per flush as a heuristic — large enough
/// to amortize syscall and gzip overhead, small enough to avoid memory pressure
/// when many instances run in parallel on an HPC node.
/// The result is clamped between 64 KB (floor, never regress below OS page granularity)
/// and 64 MB (ceiling, avoids OOM on ultra-long reads at high parallelism).
pub fn compute_write_buffer_size(read_length: usize) -> usize {
    let bytes_per_record = 2 * read_length + 42;
    let target = bytes_per_record * 512;
    target.clamp(64 * 1024, 64 * 1024 * 1024)
}

/// Computes the bytes-per-read divisor for capacity pre-allocation.
///
/// Same record size formula as above, with a compression factor of ~4x for gzip.
/// Clamped to a minimum of 1 to prevent division by zero.
pub fn compute_bytes_per_read(read_length: usize, is_gz: bool) -> usize {
    let plain = (2 * read_length + 42).max(1);
    if is_gz {
        (plain / 4).max(1) // ~3.5–4x compression on DNA is typical
    } else {
        plain
    }
}

/// Estimates the number of sequences in a file to pre-allocate memory.
/// The result is capped at `MAX_PREALLOC_ENTRIES` to avoid over-allocation.
pub fn estimate_sequence_capacity<P: AsRef<Path>>(path: P, read_length: usize) -> Result<usize> {
    let path_ref = path.as_ref();
    if !path_ref.exists() {
        return Ok(0);
    }

    let file_size_bytes = fs::metadata(path_ref)?.len() as usize;

    let is_gz = path_ref.to_string_lossy().to_lowercase().ends_with(".gz");

    let divisor = compute_bytes_per_read(read_length, is_gz);

    Ok((file_size_bytes / divisor).clamp(1, MAX_PREALLOC_ENTRIES))
}

/// Resolves a CLI path to a filesystem identity, for "is this the same file?".
///
/// Comparing the strings is not comparing the files: `-o ./coll2.fq -p
/// coll2.fq` are two spellings of one path, so a string check let two writers
/// interleave into it, `--json j.fq -o j.fq` overwrote a FASTQ with JSON, and
/// `-1 f.fq -o f.fq --force` truncated its own input. Canonicalizing resolves
/// `.`/`..`, symlinks and absolute-vs-relative spellings.
///
/// A path that does not exist yet cannot be canonicalized; resolve its parent
/// directory instead and keep the file name, which covers the usual "output
/// does not exist yet" case. Anything we cannot resolve at all is compared as
/// written — a best effort, never a silent skip.
pub fn path_identity(path: &str) -> PathBuf {
    let raw = Path::new(path);
    if let Ok(resolved) = fs::canonicalize(raw) {
        return resolved;
    }
    let parent = match raw.parent() {
        // `Path::new("out.fq").parent()` is `Some("")`, which is not a directory.
        Some(p) if p.as_os_str().is_empty() => Path::new("."),
        Some(p) => p,
        None => return raw.to_path_buf(),
    };
    match (fs::canonicalize(parent), raw.file_name()) {
        (Ok(dir), Some(name)) => dir.join(name),
        _ => raw.to_path_buf(),
    }
}

/// Performs a fast first-pass count of sequences in a FASTX file.
/// Uses needletail's streaming parser — no storage, just counting.
pub fn count_sequences(path: &str) -> Result<usize> {
    use needletail::parse_fastx_file;
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("Failed to open input for pre-count: {path}"))?;
    let mut count = 0;
    while let Some(r) = reader.next() {
        r.with_context(|| format!("Invalid sequence at record #{} during pre-count", count + 1))?;
        count += 1;
    }
    Ok(count)
}

/// Validates that the output format is compatible with the input data.
/// FASTA input cannot be written as FASTQ (no quality scores).
pub fn validate_format_compatibility(
    has_quality: bool,
    output_format: crate::writer::OutputFormat,
) -> Result<()> {
    if !has_quality && output_format == crate::writer::OutputFormat::Fastq {
        bail!(
            "FASTA → FASTQ conversion not supported: the input file does not contain \
             quality scores. Use a .fasta/.fa/.fna extension for output."
        );
    }
    Ok(())
}

/// Extracts the base read ID, stripping both the description (anything after
/// the first whitespace) and the `/1` or `/2` mate suffix used by pre-CASAVA-1.8
/// Illumina headers and some SRA `fastq-dump --split-files` outputs.
#[inline]
pub fn base_read_id(id: &[u8]) -> &[u8] {
    // Strip description: "ID desc" → "ID"
    let head = id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id);
    // Strip /1 or /2 mate suffix if present
    if head.len() >= 2 {
        let tail = &head[head.len() - 2..];
        if tail == b"/1" || tail == b"/2" {
            return &head[..head.len() - 2];
        }
    }
    head
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hasher::HashType;

    // --- birthday_problem_square_approximation ---

    #[test]
    fn test_birthday_64bit_at_2_pow_32() {
        let n = (u32::MAX as usize) + 1;
        let prob = birthday_problem_square_approximation(n, &HashType::XXH3_64);
        assert!(
            (prob - 0.5).abs() < 1e-9,
            "Expected probability 0.5, got {:.12}",
            prob
        );
    }

    #[test]
    fn test_birthday_128bit_at_2_pow_32() {
        let n = (u32::MAX as usize) + 1;
        let prob = birthday_problem_square_approximation(n, &HashType::XXH3_128);
        let expected = (n as f64).powi(2) / 2.0_f64.powi(129);
        assert!(
            (prob - expected).abs() < 1e-30,
            "Expected probability {:.3e}, got {:.3e}",
            expected,
            prob
        );
    }

    #[test]
    fn test_birthday_zero_sequences() {
        let prob = birthday_problem_square_approximation(0, &HashType::XXH3_64);
        assert_eq!(prob, 0.0);
    }

    #[test]
    fn test_birthday_one_sequence() {
        let prob = birthday_problem_square_approximation(1, &HashType::XXH3_64);
        assert!(prob > 0.0 && prob < 1e-15);
    }

    // --- get_hash_method ---

    #[test]
    fn test_get_hash_method_large_dataset_selects_128() {
        let large_size = 200_000_000_usize;
        let method = get_hash_method(large_size, 0.001);
        assert!(
            matches!(method, HashType::XXH3_128),
            "A dataset of {} sequences should use XXH3_128",
            large_size
        );
    }

    #[test]
    fn test_get_hash_method_small_dataset_selects_64() {
        let small_size = 1_000_usize;
        let method = get_hash_method(small_size, 0.001);
        assert!(
            matches!(method, HashType::XXH3_64),
            "A dataset of {} sequences should use XXH3_64",
            small_size
        );
    }

    #[test]
    fn test_get_hash_method_threshold_zero_forces_128() {
        let method = get_hash_method(1, 0.0);
        assert!(
            matches!(method, HashType::XXH3_128),
            "A threshold of 0 should always select XXH3_128"
        );
    }

    #[test]
    fn test_get_hash_method_threshold_one_forces_64() {
        let method = get_hash_method(1_000_000, 1.0);
        assert!(
            matches!(method, HashType::XXH3_64),
            "A threshold of 1.0 should select XXH3_64 for a normal-sized dataset"
        );
    }

    #[test]
    fn test_get_hash_method_tipping_point() {
        let boundary = 4_294_967_296_usize;
        let just_below = get_hash_method(boundary - 1, 0.5);
        let just_above = get_hash_method(boundary + 1, 0.5);
        assert!(
            matches!(just_below, HashType::XXH3_64),
            "Just below the boundary should be XXH3_64"
        );
        assert!(
            matches!(just_above, HashType::XXH3_128),
            "Just above the boundary should be XXH3_128"
        );
    }

    // --- compute_write_buffer_size ---

    #[test]
    fn test_buffer_size_illumina_150bp() {
        let buf = compute_write_buffer_size(150);
        assert_eq!(buf, 175_104);
        assert!(buf >= 64 * 1024);
        assert!(buf <= 64 * 1024 * 1024);
    }

    #[test]
    fn test_buffer_size_floor_clamping() {
        let buf = compute_write_buffer_size(1);
        assert_eq!(buf, 64 * 1024);
    }

    #[test]
    fn test_buffer_size_ceiling_clamping() {
        let buf = compute_write_buffer_size(100_000);
        assert_eq!(buf, 64 * 1024 * 1024);
    }

    // --- compute_bytes_per_read ---

    #[test]
    fn test_bytes_per_read_plain_150bp() {
        assert_eq!(compute_bytes_per_read(150, false), 342);
    }

    #[test]
    fn test_bytes_per_read_gz_150bp() {
        assert_eq!(compute_bytes_per_read(150, true), 85);
    }

    #[test]
    fn test_bytes_per_read_zero_length_no_panic() {
        assert!(compute_bytes_per_read(0, false) >= 1);
        assert!(compute_bytes_per_read(0, true) >= 1);
    }

    // --- path_identity ---

    #[test]
    fn test_path_identity_ignores_spelling() {
        let dir = std::env::temp_dir().join(format!("fastdedup_identity_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let file = dir.join("out.fq");
        fs::write(&file, b"@r\nA\n+\nI\n").unwrap();

        let plain = file.to_string_lossy().to_string();
        let dotted = dir.join(".").join("out.fq").to_string_lossy().to_string();
        assert_eq!(path_identity(&plain), path_identity(&dotted));

        // Un fichier qui n'existe pas encore : l'identité passe par le parent.
        let missing = dir.join("nouveau.fq").to_string_lossy().to_string();
        let missing_dotted = dir
            .join(".")
            .join("nouveau.fq")
            .to_string_lossy()
            .to_string();
        assert_eq!(path_identity(&missing), path_identity(&missing_dotted));
        assert_ne!(path_identity(&plain), path_identity(&missing));

        let _ = fs::remove_dir_all(&dir);
    }

    // --- base_read_id ---

    #[test]
    fn test_base_read_id_strips_mate_suffix() {
        assert_eq!(
            base_read_id(b"HWUSI-EAS100R:6:73:941:1973#0/1"),
            b"HWUSI-EAS100R:6:73:941:1973#0"
        );
        assert_eq!(
            base_read_id(b"HWUSI-EAS100R:6:73:941:1973#0/2"),
            b"HWUSI-EAS100R:6:73:941:1973#0"
        );
        // CASAVA 1.8+ format (no /N suffix, metadata after space)
        assert_eq!(
            base_read_id(b"HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT"),
            b"HWI-ST1276:71:C1162ACXX:1:1101:1208:2458"
        );
        // SRA split-files variant
        assert_eq!(
            base_read_id(b"SRR001666.1/1 HWI:1:1:82:1089 length=36"),
            b"SRR001666.1"
        );
        // No suffix, no description
        assert_eq!(base_read_id(b"read42"), b"read42");
    }
}
