use std::fs::File;
use std::io::{BufRead, BufReader};

use flate2::read::MultiGzDecoder;

/// Encode a read sequence as `(x = T+C count, y = T+G count)`.
/// Non-ACGT characters (e.g. N) are ignored.
pub fn encode_read(seq: &str) -> (u32, u32) {
    let mut x = 0u32;
    let mut y = 0u32;
    for b in seq.bytes() {
        let (bx, by) = match b {
            b'A' => (0, 0),
            b'T' => (1, 1),
            b'G' => (0, 1),
            b'C' => (1, 0),
            _ => continue,
        };
        x += bx;
        y += by;
    }
    (x, y)
}

/// Mean Phred quality score for a quality string (ASCII-33 encoded).
pub fn mean_quality(qual: &str) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: f64 = qual.bytes().map(|b| b.saturating_sub(33) as f64).sum();
    sum / qual.len() as f64
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
        assert_eq!(encode_read("AAAA"), (0, 0));
    }

    #[test]
    fn all_t() {
        assert_eq!(encode_read("TTTT"), (4, 4));
    }

    #[test]
    fn mixed() {
        assert_eq!(encode_read("ATGC"), (2, 2));
    }

    #[test]
    fn skips_n() {
        assert_eq!(encode_read("ANNN"), (0, 0));
    }

    #[test]
    fn quality() {
        assert_eq!(mean_quality("!"), 0.0);      // ASCII 33 → Q0
        assert_eq!(mean_quality("I"), 40.0);     // ASCII 73 → Q40
        assert!((mean_quality("!I") - 20.0).abs() < 1e-9);
    }
}
