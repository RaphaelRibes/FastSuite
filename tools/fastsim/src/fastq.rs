//! FASTQ writer: plain or gzip, with NovaSeq-style instrument headers.
//!
//! Header format mimics the real data:
//! `@LH00916:<flowcell>:<lane>:<tile>:<x>:<y> <mate>:N:0:<index>`
//! where `<x>:<y>` is a unique coordinate per pair (derived from the pair
//! index), so R1 and R2 of a pair share an ID and a trimmed file can be
//! re-paired by coordinate. Phred+33 encoding.
//!
//! Three properties this writer is built around, because its bytes are the
//! corpus generator every downstream bench and paper is validated against:
//!
//! 1. **The payload is fixed.** `--gz-level` changes the DEFLATE encoding
//!    only: decompressed bytes are identical at every level, so the flag can
//!    never invalidate an existing campaign.
//! 2. **The record is emitted in six `write_all` calls, not 300+.** The
//!    quality line is encoded once into a reused buffer (FS-6) instead of one
//!    `write_all` per base.
//! 3. **The sink is an enum, not `Box<dyn Write>`** (FS-7), so `put` inlines
//!    into the copy loop and there is no vtable hop per record. The
//!    `BufWriter` capacity and layering are unchanged, so the DEFLATE encoder
//!    still sees the same 8 KiB chunks and produces the same `.gz` bytes.
//!
//! ## The gzip level trade (measured, 200 000 pairs, `taskset -c 0`)
//!
//! Deflate at level 6 is ~75 % of fastsim's wall time: grep the wall cells in
//! the tool's own benchmark and the compressor is the run. FASTQ quality
//! strings are near-incompressible noise, so the hash-chain search that level 6
//! pays for mostly fails. Two measurements, because they disagree and the
//! difference matters:
//!
//! | backend | level | input rate | bytes written (both mates) |
//! |---|---|---|---|
//! | flate2/miniz_oxide (what fastsim uses) | 6 | ~11 MB/s | 38.75 MB (27.4 % of plain) |
//! | flate2/miniz_oxide | 1 | ~200 MB/s | 52.63 MB (+35.8 %) |
//! | system `gzip` CLI (review proxy) | 6 | 9.2 MB/s | 19.54 MB/mate (27.7 %) |
//! | system `gzip` CLI (review proxy) | 1 | 109 MB/s | 23.04 MB/mate (+17.9 %) |
//!
//! The proxy *understates* what level 1 costs in fastsim's own backend (miniz's
//! fast mode is worse at packing than C zlib's), and understates its speed as
//! well — end to end, level 1 takes 200 000 pairs from ~14.8 s to ~2.3 s on the
//! reference laptop, i.e. ~6x, of which the generator's ~1.8 s is the floor.
//! Quote +35.8 %, not +17.9 %, when sizing a campaign.
//!
//! Level 6 stays the default so the bytes on disk do not change for anyone who
//! does not ask; `--gz-level 1` (and 4, the balanced choice) is the opt-in.
//! Nothing downstream parses `.gz` size — fastcheck and fastprep re-inflate,
//! and inflate is ~4× cheaper than deflate.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::fmtint::push_u64;
use crate::quality::MateRead;

pub const PHIRED_OFFSET: u8 = 33;
/// gzip level used when `--gz-level` is not given: today's 0.1.0 bytes.
pub const DEFAULT_GZ_LEVEL: u32 = 6;
/// Smallest accepted level (deflate's fastest).
pub const MIN_GZ_LEVEL: u32 = 1;
/// Largest accepted level (deflate's smallest output).
pub const MAX_GZ_LEVEL: u32 = 9;

const INSTRUMENT: &str = "LH00916";
const FLOWCELL: &str = "SIMFLOW1";
const INDEX: &str = "ATCGATCG"; // fake dual-index barcode

/// Phred score → Phred+33 ASCII, for every value a `u8` can hold.
///
/// A 256-entry table (not 42) so a profile with an out-of-range
/// `full_length_q` still emits `q + 33` instead of panicking — the same bytes
/// the per-byte `q + PHIRED_OFFSET` loop produced.
const QUAL_ASCII: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = (i as u8).wrapping_add(PHIRED_OFFSET);
        i += 1;
    }
    t
};

/// The output end of a [`FastqWriter`]: plain or gzip, dispatch resolved by an
/// enum the compiler can inline instead of a `dyn Write` vtable.
pub enum Sink {
    Plain(BufWriter<File>),
    Gz(BufWriter<GzEncoder<File>>),
}

impl Sink {
    #[inline]
    fn put(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        match self {
            Sink::Plain(w) => w.write_all(bytes),
            Sink::Gz(w) => w.write_all(bytes),
        }
    }

    /// Flush the buffer into the encoder (which finishes its stream on drop,
    /// exactly as the boxed writer did).
    fn finish(mut self) -> std::io::Result<()> {
        match &mut self {
            Sink::Plain(w) => w.flush(),
            Sink::Gz(w) => w.flush(),
        }
    }
}

pub struct FastqWriter {
    out: Sink,
    /// Constant part of every header, built once (`@INSTRUMENT:run:flowcell:lane:tile:`).
    prefix: Vec<u8>,
    /// Reused header buffer, one per record (FS-10).
    hdr: Vec<u8>,
    /// Reused quality-line buffer (FS-6).
    qbuf: Vec<u8>,
    run: u32,
    flowcell: String,
    lane: u8,
    tile: u32,
}

impl FastqWriter {
    /// Open a gzip FASTQ at `path` (creates/truncates).
    ///
    /// `level` is the DEFLATE level, 1..=9; it changes only the encoding, never
    /// the payload. See the module docs for the measured rate/size trade and
    /// use [`DEFAULT_GZ_LEVEL`] for 0.1.0's bytes.
    pub fn open_gz(path: &Path, level: u32) -> Result<Self> {
        let f = File::create(path).with_context(|| format!("create {}", path.display()))?;
        let gz = GzEncoder::new(f, Compression::new(level));
        Ok(FastqWriter::new(Sink::Gz(BufWriter::new(gz))))
    }

    /// Open a gzip FASTQ at the default level ([`DEFAULT_GZ_LEVEL`]).
    pub fn open_gz_default(path: &Path) -> Result<Self> {
        FastqWriter::open_gz(path, DEFAULT_GZ_LEVEL)
    }

    /// Open a plain (uncompressed) FASTQ at `path`.
    pub fn open_plain(path: &Path) -> Result<Self> {
        let f = File::create(path).with_context(|| format!("create {}", path.display()))?;
        Ok(FastqWriter::new(Sink::Plain(BufWriter::new(f))))
    }

    fn new(out: Sink) -> Self {
        let mut w = FastqWriter {
            out,
            prefix: Vec::with_capacity(48),
            hdr: Vec::with_capacity(64),
            qbuf: Vec::with_capacity(320),
            run: 1,
            flowcell: FLOWCELL.to_string(),
            lane: 1,
            tile: 1101,
        };
        // The header prefix is constant for the whole run; only x/y/mate vary.
        w.prefix.push(b'@');
        w.prefix.extend_from_slice(INSTRUMENT.as_bytes());
        w.prefix.push(b':');
        push_u64(&mut w.prefix, w.run as u64);
        w.prefix.push(b':');
        w.prefix.extend_from_slice(w.flowcell.as_bytes());
        w.prefix.push(b':');
        push_u64(&mut w.prefix, w.lane as u64);
        w.prefix.push(b':');
        push_u64(&mut w.prefix, w.tile as u64);
        w.prefix.push(b':');
        w
    }

    /// Append the header for pair `pair_idx`, mate 1 or 2, into `buf`.
    #[inline]
    fn header_into(prefix: &[u8], buf: &mut Vec<u8>, pair_idx: u64, mate: u8) {
        // coordinates: x = pair_idx % 10_000, y = pair_idx / 10_000
        buf.clear();
        buf.extend_from_slice(prefix);
        push_u64(buf, pair_idx % 10_000);
        buf.push(b':');
        push_u64(buf, pair_idx / 10_000);
        buf.push(b' ');
        push_u64(buf, mate as u64);
        buf.extend_from_slice(b":N:0:");
        buf.extend_from_slice(INDEX.as_bytes());
    }

    /// Header for pair `pair_idx`, mate 1 or 2.
    pub fn header(&self, pair_idx: u64, mate: u8) -> String {
        let mut buf = Vec::with_capacity(64);
        Self::header_into(&self.prefix, &mut buf, pair_idx, mate);
        // ASCII by construction: instrument/flowcell/index are literals or
        // decimal digits.
        String::from_utf8(buf).expect("header is ASCII")
    }

    /// Encode the Phred line into the reused buffer: one pass, no per-base
    /// `write_all`.
    #[inline]
    fn encode_qual(&mut self, qual: &[u8]) {
        self.qbuf.clear();
        if self.qbuf.capacity() < qual.len() {
            self.qbuf.reserve(qual.len());
        }
        self.qbuf.extend(qual.iter().map(|&q| QUAL_ASCII[q as usize]));
    }

    /// Write one record under an already formatted name.
    pub fn write_record(&mut self, name: &str, mr: &MateRead) -> Result<()> {
        self.out.put(name.as_bytes())?;
        self.out.put(b"\n")?;
        self.out.put(&mr.seq)?;
        self.out.put(b"\n+\n")?;
        self.encode_qual(&mr.qual);
        self.out.put(&self.qbuf)?;
        self.out.put(b"\n")?;
        Ok(())
    }

    /// Write one mate of pair `pair_idx`, formatting the header into a reused
    /// buffer (no `String` per record).
    pub fn write_mate(&mut self, pair_idx: u64, mate: u8, mr: &MateRead) -> Result<()> {
        Self::header_into(&self.prefix, &mut self.hdr, pair_idx, mate);
        // Split the borrow: `hdr`/`qbuf` are read while `out` is written.
        self.out.put(&self.hdr)?;
        self.out.put(b"\n")?;
        self.out.put(&mr.seq)?;
        self.out.put(b"\n+\n")?;
        self.encode_qual(&mr.qual);
        self.out.put(&self.qbuf)?;
        self.out.put(b"\n")?;
        Ok(())
    }

    /// Write a mate pair.
    pub fn write_pair(&mut self, pair_idx: u64, r1: &MateRead, r2: &MateRead) -> Result<()> {
        self.write_mate(pair_idx, 1, r1)?;
        self.write_mate(pair_idx, 2, r2)
    }

    /// Flush and finish (finishes gzip stream).
    pub fn finish(self) -> Result<()> {
        self.out.finish()?;
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

    fn seq_of(n: usize) -> Vec<u8> {
        (0..n).map(|i| b"ACGT"[i % 4]).collect()
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

    /// FS-10: the header is appended from integers, so it must still equal the
    /// documented format string at every boundary of the x/y split and for
    /// pair indices that overflow a 32-bit coordinate.
    #[test]
    fn header_append_matches_the_format_macro() {
        let path = std::env::temp_dir().join("fastsim_hdr2.fastq");
        let w = FastqWriter::open_plain(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        for pair_idx in [
            0u64,
            1,
            9,
            9_999,
            10_000,
            12_345,
            123_456,
            4_294_967_296,
            18_446_744_073_709_551_615,
        ] {
            for mate in [1u8, 2] {
                let x = pair_idx % 10_000;
                let y = pair_idx / 10_000;
                let want = format!(
                    "@{}:{}:{}:{}:{}:{}:{} {}:N:0:{}",
                    "LH00916", 1, "SIMFLOW1", 1, 1101, x, y, mate, "ATCGATCG"
                );
                assert_eq!(w.header(pair_idx, mate), want, "pair {pair_idx} mate {mate}");
            }
        }
    }

    /// FS-6/FS-10 guard: the record layout is a contract — header, sequence,
    /// `+`, the Phred+33 quality line, newline — including the extremes of the
    /// Phred range.
    #[test]
    fn record_layout_is_exact() {
        let path = std::env::temp_dir().join("fastsim_layout.fastq");
        let mut w = FastqWriter::open_plain(&path).unwrap();
        let qual: Vec<u8> = (0..42u8).collect();
        let seq = seq_of(42);
        w.write_pair(9, &mr(&seq, &qual), &mr(&seq, &qual)).unwrap();
        w.finish().unwrap();
        let got = std::fs::read(&path).unwrap();
        std::fs::remove_file(&path).unwrap();

        let qline: String = (0..42u8).map(|q| (q + PHIRED_OFFSET) as char).collect();
        let sline = String::from_utf8(seq).unwrap();
        let want = format!(
            "@LH00916:1:SIMFLOW1:1:1101:9:0 1:N:0:ATCGATCG\n{sline}\n+\n{qline}\n\
             @LH00916:1:SIMFLOW1:1:1101:9:0 2:N:0:ATCGATCG\n{sline}\n+\n{qline}\n"
        );
        assert_eq!(String::from_utf8(got).unwrap(), want);
    }

    /// FS-1: the level is parameterised, and it may not change the payload.
    #[test]
    fn gz_level_is_parameterised_and_the_payload_is_fixed() {
        assert_eq!(DEFAULT_GZ_LEVEL, 6, "the default must stay 0.1.0's level");

        let mut payloads = Vec::new();
        let mut sizes = Vec::new();
        for level in [MIN_GZ_LEVEL, DEFAULT_GZ_LEVEL, MAX_GZ_LEVEL] {
            let path = std::env::temp_dir().join(format!("fastsim_lvl{level}.fastq.gz"));
            let mut w = FastqWriter::open_gz(&path, level).unwrap();
            // Entropy like real quality strings, so the level difference is
            // visible rather than hidden by a tiny record.
            for i in 0..400u64 {
                let qual: Vec<u8> = (0..150).map(|j| ((i as usize * 7 + j * 13) % 42) as u8).collect();
                let seq: Vec<u8> = (0..150).map(|j| b"ACGT"[(i as usize + j) % 4]).collect();
                w.write_pair(i, &mr(&seq, &qual), &mr(&seq, &qual)).unwrap();
            }
            w.finish().unwrap();

            let f = File::open(&path).unwrap();
            let mut plain = Vec::new();
            flate2::read::GzDecoder::new(f)
                .read_to_end(&mut plain)
                .unwrap();
            payloads.push(plain);
            sizes.push(std::fs::metadata(&path).unwrap().len());
            std::fs::remove_file(&path).unwrap();
        }
        assert_eq!(payloads[0], payloads[1], "level 1 must inflate to level 6's bytes");
        assert_eq!(payloads[1], payloads[2], "level 9 must inflate to level 6's bytes");
        assert!(
            sizes[0] > sizes[1],
            "level 1 trades size for speed: {sizes:?}"
        );
        assert!(
            sizes[2] <= sizes[1] + 64,
            "level 9 must not be materially larger than level 6: {sizes:?}"
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
        let mut w = FastqWriter::open_gz_default(&path).unwrap();
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

    /// `write_mate` and `write_record` are two entry points to the same bytes:
    /// the duplicate path uses the first, the format test the second.
    #[test]
    fn write_mate_and_write_record_agree() {
        let p1 = std::env::temp_dir().join("fastsim_mate_a.fastq");
        let p2 = std::env::temp_dir().join("fastsim_mate_b.fastq");
        let r = mr(&seq_of(150), &[30; 150]);
        let mut a = FastqWriter::open_plain(&p1).unwrap();
        a.write_mate(5, 2, &r).unwrap();
        a.finish().unwrap();
        let mut b = FastqWriter::open_plain(&p2).unwrap();
        let name = b.header(5, 2);
        b.write_record(&name, &r).unwrap();
        b.finish().unwrap();
        let (x, y) = (std::fs::read(&p1).unwrap(), std::fs::read(&p2).unwrap());
        std::fs::remove_file(&p1).unwrap();
        std::fs::remove_file(&p2).unwrap();
        assert_eq!(x, y);
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
