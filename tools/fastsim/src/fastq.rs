//! FASTQ writer: plain or gzip, with NovaSeq-style instrument headers.
//!
//! Header format mimics the real data:
//! `@LH00916:<flowcell>:<lane>:<tile>:<x>:<y> <mate>:N:0:<index>`
//! where `<x>:<y>` is a unique coordinate per pair (derived from the pair
//! index), so R1 and R2 of a pair share an ID and a trimmed file can be
//! re-paired by coordinate. Phred+33 encoding.

use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::quality::MateRead;

pub const PHIRED_OFFSET: u8 = 33;
const INSTRUMENT: &str = "LH00916";
const FLOWCELL: &str = "SIMFLOW1";
const INDEX: &str = "ATCGATCG"; // fake dual-index barcode

pub struct FastqWriter {
    out: BufWriter<Box<dyn Write>>,
    run: u32,
    flowcell: String,
    lane: u8,
    tile: u32,
    index: String,
}

impl FastqWriter {
    /// Open a gzip FASTQ at `path` (creates/truncates).
    pub fn open_gz(path: &Path) -> Result<Self> {
        let f =
            std::fs::File::create(path).with_context(|| format!("create {}", path.display()))?;
        let gz = GzEncoder::new(Box::new(f) as Box<dyn Write>, Compression::new(6));
        Ok(FastqWriter::new(BufWriter::new(
            Box::new(gz) as Box<dyn Write>
        )))
    }

    /// Open a plain (uncompressed) FASTQ at `path`.
    pub fn open_plain(path: &Path) -> Result<Self> {
        let f =
            std::fs::File::create(path).with_context(|| format!("create {}", path.display()))?;
        Ok(FastqWriter::new(BufWriter::new(
            Box::new(f) as Box<dyn Write>
        )))
    }

    fn new(out: BufWriter<Box<dyn Write>>) -> Self {
        FastqWriter {
            out,
            run: 1,
            flowcell: FLOWCELL.to_string(),
            lane: 1,
            tile: 1101,
            index: INDEX.to_string(),
        }
    }

    /// Header for pair `pair_idx`, mate 1 or 2.
    pub fn header(&self, pair_idx: u64, mate: u8) -> String {
        // coordinates: x = pair_idx % 10_000, y = pair_idx / 10_000
        let x = pair_idx % 10_000;
        let y = pair_idx / 10_000;
        format!(
            "@{}:{}:{}:{}:{}:{}:{} {}:N:0:{}",
            INSTRUMENT, self.run, self.flowcell, self.lane, self.tile, x, y, mate, self.index
        )
    }

    /// Write one record.
    pub fn write_record(&mut self, name: &str, mr: &MateRead) -> Result<()> {
        self.out.write_all(name.as_bytes())?;
        self.out.write_all(b"\n")?;
        self.out.write_all(&mr.seq)?;
        self.out.write_all(b"\n+\n")?;
        // Phred+33: ASCII = phred + 33
        for &q in &mr.qual {
            self.out.write_all(&[q + PHIRED_OFFSET])?;
        }
        self.out.write_all(b"\n")?;
        Ok(())
    }

    /// Write a mate pair.
    pub fn write_pair(&mut self, pair_idx: u64, r1: &MateRead, r2: &MateRead) -> Result<()> {
        self.write_record(&self.header(pair_idx, 1), r1)?;
        self.write_record(&self.header(pair_idx, 2), r2)
    }

    /// Flush and finish (finishes gzip stream).
    pub fn finish(mut self) -> Result<()> {
        self.out.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn mr(seq: &[u8], qual: &[u8]) -> MateRead {
        MateRead {
            seq: seq.to_vec(),
            qual: qual.to_vec(),
        }
    }

    #[test]
    fn header_format_matches_novaseq_style() {
        let path = std::env::temp_dir().join("fastsim_hdr.fastq");
        let w = FastqWriter::open_plain(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let h = w.header(123456, 1);
        assert!(h.starts_with("@LH00916:1:SIMFLOW1:1:1101:"));
        assert!(h.ends_with(" 1:N:0:ATCGATCG"));
        // unique coordinate per pair
        assert_ne!(w.header(0, 1), w.header(1, 1));
        // mates share the coordinate, differ only in the mate field
        assert_eq!(
            w.header(42, 1).split(' ').next(),
            w.header(42, 2).split(' ').next()
        );
    }

    #[test]
    fn write_pair_roundtrip_plain() {
        let path = std::env::temp_dir().join("fastsim_r1.fastq");
        let mut w = FastqWriter::open_plain(&path).unwrap();
        let r1 = mr(b"ACGT", &[30, 30, 30, 30]);
        let r2 = mr(b"TGCA", &[20, 20, 20, 20]);
        w.write_pair(7, &r1, &r2).unwrap();
        w.finish().unwrap();

        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 8, "2 records x 4 lines");
        assert!(lines[0].starts_with("@LH00916:"));
        assert_eq!(lines[1], "ACGT");
        assert_eq!(lines[3], "????"); // Q30+33 = 63 = '?'
        assert_eq!(lines[5], "TGCA"); // R2 seq
        assert_eq!(lines[7], "5555"); // Q20+33 = 53 = '5'
    }

    #[test]
    fn gzip_output_decompresses() {
        let path = std::env::temp_dir().join("fastsim_gz.fastq.gz");
        let mut w = FastqWriter::open_gz(&path).unwrap();
        let r1 = mr(b"AAAA", &[20, 21, 22, 23]);
        let r2 = mr(b"TTTT", &[24, 25, 26, 27]);
        w.write_pair(1, &r1, &r2).unwrap();
        w.finish().unwrap();

        let f = std::fs::File::open(&path).unwrap();
        let mut gz = flate2::read::GzDecoder::new(f);
        let mut s = String::new();
        gz.read_to_string(&mut s).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert!(s.contains("AAAA"));
        assert!(s.contains("TTTT"));
        assert_eq!(s.lines().count(), 8);
    }

    #[test]
    fn deterministic_bytes_for_same_input() {
        let p1 = std::env::temp_dir().join("det_a.fastq");
        let p2 = std::env::temp_dir().join("det_b.fastq");
        let r1 = mr(b"ACGTACGT", &[30; 8]);
        let r2 = mr(b"TGCATGCA", &[25; 8]);
        for p in [&p1, &p2] {
            let mut w = FastqWriter::open_plain(p).unwrap();
            w.write_pair(3, &r1, &r2).unwrap();
            w.finish().unwrap();
        }
        let a = std::fs::read(&p1).unwrap();
        let b = std::fs::read(&p2).unwrap();
        std::fs::remove_file(&p1).unwrap();
        std::fs::remove_file(&p2).unwrap();
        assert_eq!(a, b);
    }
}
