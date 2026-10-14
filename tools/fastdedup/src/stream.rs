//! D-4 standard streams: `-` means stdin (input) or stdout (output).
//!
//! A stream has no filename, so nothing may be inferred from one:
//!   * compression on INPUT is sniffed from magic bytes (needletail's
//!     auto-detection: gzip/bzip2/xz/zstd, plain as fallback),
//!   * format and compression on OUTPUT are explicit flags
//!     (`--out-format`, `--out-compression`), defaulting to plain FASTQ.

use anyhow::Result;

/// The path that means "the stream" (docs/DECISIONS.md, D-4).
pub const STREAM: &str = "-";

/// True when `path` designates stdin/stdout rather than a file.
pub fn is_stream(path: &str) -> bool {
    path == STREAM
}

/// Open a FASTX input: stdin for `-`, the file otherwise.
///
/// `parse_fastx_reader` sniffs the first two bytes and selects the decoder,
/// so gzip on stdin behaves exactly like gzip in a mislabelled file.
pub fn open_fastx(path: &str) -> Result<Box<dyn needletail::FastxReader>> {
    if is_stream(path) {
        needletail::parse_fastx_reader(std::io::stdin())
            .map_err(|e| anyhow::anyhow!("reading FASTX from stdin: {e}"))
    } else {
        needletail::parse_fastx_file(path)
            .map_err(|e| anyhow::anyhow!("Failed to read input file {path}: {e}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dash_is_the_stream() {
        assert!(is_stream("-"));
        assert!(!is_stream("--"));
        assert!(!is_stream("./-"));
        assert!(!is_stream("in.fastq"));
    }

    #[test]
    fn a_stream_is_just_a_reader() {
        // A stream is a `Read`; the same sniffing path must apply to it.
        let body = b"@r1\nACGT\n+\nIIII\n";
        let mut rdr = needletail::parse_fastx_reader(&body[..]).unwrap();
        let rec = rdr.next().expect("one record").unwrap();
        assert_eq!(rec.id(), b"r1");
        assert!(rdr.next().is_none());
    }
}
