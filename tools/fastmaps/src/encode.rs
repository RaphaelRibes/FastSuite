use std::fs::File;
use std::io::{BufRead, BufReader};

use flate2::read::MultiGzDecoder;

/// Encode a read sequence as `(x = T+C count, y = T+G count, n = ACGT bases)`.
///
/// Takes bytes, not `&str`: the FASTQ loop reads with `read_until` into reused
/// buffers, so a read never has to be UTF-8-validated (or allocated) just to be
/// counted. Every non-ACGT byte — `N`, a stray `\r`, an arbitrary binary byte —
/// is skipped identically.
///
/// Non-ACGT characters (e.g. N) are skipped from the coordinates, so a read
/// with no ACGT base at all encodes to `(0, 0)` — the same pixel as an all-A
/// read. `n` is what lets the caller tell those apart: an `n == 0` read has no
/// measurable composition and must not be counted as polyA-rich.
pub fn encode_read(seq: &[u8]) -> (u32, u32, u32) {
    let mut x = 0u32;
    let mut y = 0u32;
    let mut n = 0u32;
    for &b in seq {
        let (bx, by) = match b {
            b'A' => (0, 0),
            b'T' => (1, 1),
            b'G' => (0, 1),
            b'C' => (1, 0),
            _ => continue,
        };
        x += bx;
        y += by;
        n += 1;
    }
    (x, y, n)
}

/// `(sum of Phred scores, number of quality bytes)` — the integer pair whose
/// ratio is the mean quality.
///
/// Integer accumulation is what makes the per-pixel merge order-independent
/// (FM-12): summing `f64` partial means in whatever order worker threads
/// happen to finish is not associative, summing `u64` is exactly so.
pub fn quality_sums(qual: &[u8]) -> (u64, u64) {
    if qual.is_empty() {
        return (0, 0);
    }
    let sum: u64 = qual.iter().map(|&b| b.saturating_sub(33) as u64).sum();
    (sum, qual.len() as u64)
}

/// Open a FASTQ file for reading, transparently decompressing `.gz` files.
pub fn get_reader(path: &str) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Escape a string for use as a JSON string value (including the surrounding
/// double-quotes). Used to avoid pulling in serde.
pub fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_a() {
        assert_eq!(encode_read(b"AAAA"), (0, 0, 4));
    }

    #[test]
    fn all_t() {
        assert_eq!(encode_read(b"TTTT"), (4, 4, 4));
    }

    #[test]
    fn mixed() {
        assert_eq!(encode_read(b"ATGC"), (2, 2, 4));
    }

    #[test]
    fn skips_n() {
        assert_eq!(encode_read(b"ANNN"), (0, 0, 1));
    }

    #[test]
    fn all_n_has_no_acgt() {
        assert_eq!(encode_read(b"NNNN"), (0, 0, 0));
    }

    /// FM-10: a byte that is not valid UTF-8 is skipped like any other
    /// non-ACGT byte instead of aborting the file. `lines()` would have failed
    /// on this input and silently ended the run at that read.
    #[test]
    fn non_utf8_bytes_are_skipped_not_fatal() {
        assert_eq!(encode_read(&[0xff, 0xfe, b'A', b'C']), (1, 0, 2));
        assert_eq!(encode_read(&[0x80, 0x81]), (0, 0, 0));
    }

    fn mq(qual: &[u8]) -> f64 {
        let (sum, n) = quality_sums(qual);
        if n == 0 { 0.0 } else { sum as f64 / n as f64 }
    }

    #[test]
    fn quality() {
        assert_eq!(quality_sums(b"!"), (0, 1)); // ASCII 33 → Q0
        assert_eq!(quality_sums(b"I"), (40, 1)); // ASCII 73 → Q40
        assert!((mq(b"!I") - 20.0).abs() < 1e-9);
        assert_eq!(quality_sums(b""), (0, 0));
    }

    /// The integer path must reproduce the old `f64` sum exactly for anything a
    /// FASTQ file can hold: scores are ≤ 61 and a read is ≪ 2^53 bytes.
    #[test]
    fn quality_sums_match_the_float_mean() {
        let q = b"!5I?J+;";
        let by_float: f64 = q.iter().map(|&b| b.saturating_sub(33) as f64).sum::<f64>() / q.len() as f64;
        assert_eq!(mq(q), by_float);
        assert_eq!(quality_sums(q).0, q.iter().map(|&b| b.saturating_sub(33) as u64).sum::<u64>());
    }
}
