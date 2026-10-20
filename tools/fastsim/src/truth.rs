//! Ground-truth writer: `prefix.gt.tsv`, one row per simulated pair.
//!
//! The GT defines the exact correct answer for a trimmer, so a scorer can
//! compute trim accuracy absolutely (not relative to another tool):
//!
//! ```text
//! pair_id  insert  r1_h  r1_t  r2_h  r2_t  r1_ad  r2_ad  full_flag  is_dup  r1_n  r2_n
//! ```
//!
//! - `r1_n` / `r2_n`: ambiguous (`N`) bases the sequencing overlay injected
//!   into that mate (0 when the profile carries no N model). An injected `N`
//!   is an instrument no-call recorded INSIDE the keep window — it never
//!   touches `h`/`t` — so the GT states it explicitly rather than letting a
//!   reader mistake it for fragment biology.
//!
//! - `pair_id`: unique integer (matches the FASTQ coordinate x:y).
//! - `insert`: true fragment length I.
//! - `r1_h` / `r1_t`: 5′/3′ G-run lengths of R1; same for R2.
//! - `full_flag`: 1 = entirely artefact (full-length spike or biology
//!   window < 1 bp) — the correct answer is to DROP the pair; 0 = trim.
//! - `is_dup`: 1 = PCR duplicate of a previously emitted pair: the whole
//!   stored read (template, artefact, substitution errors and quality) is
//!   replayed byte-for-byte.
//!
//! The correct trimmed window per mate is `[h, min(I, L−t))`:
//! - `[0, h)`     — 5′ G-run (artefact)
//! - `[h, min(I, L−t))` — biology (keep)
//! - `[min(I, L−t), L−t)` — adapter readthrough (artefact, only when I < L−t)
//! - `[L−t, L)`   — 3′ G-run (artefact)
//!
//! So `h`, `t`, `I` and `L` (from the profile) fully determine truth; no
//! per-base annotation needed.

use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::artefact::Artefact;
use crate::fmtint::push_u64;

pub struct GtWriter<W: Write> {
    out: W,
    read_len: u32,
    single: bool,
    /// Reused row buffer: one `write_all` per pair, no `fmt::Display` dispatch
    /// for the twelve integer fields (FS-10).
    row: Vec<u8>,
}

impl GtWriter<BufWriter<std::fs::File>> {
    pub fn open(path: &Path, read_len: u32) -> Result<Self> {
        Self::open_with(path, read_len, false)
    }

    pub fn open_single(path: &Path, read_len: u32) -> Result<Self> {
        Self::open_with(path, read_len, true)
    }

    fn open_with(path: &Path, read_len: u32, single: bool) -> Result<Self> {
        let f =
            std::fs::File::create(path).with_context(|| format!("create {}", path.display()))?;
        let mut w = GtWriter {
            out: BufWriter::new(f),
            read_len,
            single,
            row: Vec::with_capacity(96),
        };
        w.out.write_all(if single {
            b"pair_id\tinsert\tr1_h\tr1_t\tr1_adapter\tfull_flag\tis_dup\tr1_n\n"
        } else {
            b"pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter\tfull_flag\tis_dup\tr1_n\tr2_n\n"
        })?;
        Ok(w)
    }
}

impl<W: Write> GtWriter<W> {
    /// `n1`/`n2` are the numbers of ambiguous (`N`) bases the sequencing
    /// overlay injected into each mate. They are recorded so the GT names the
    /// injected artefact instead of leaving a scorer to infer it from the
    /// FASTQ: an injected N is an instrument no-call, not fragment biology, and
    /// it sits inside the keep window (the overlay never touches `h`/`t`), so
    /// it is visible here and nowhere else in the row.
    pub fn write_pair(
        &mut self,
        pair_idx: u64,
        insert_len: u32,
        a: &Artefact,
        is_dup: bool,
        n1: u32,
        n2: u32,
    ) -> Result<()> {
        Self::row_into(
            &mut self.row,
            self.read_len,
            self.single,
            pair_idx,
            insert_len,
            a,
            is_dup,
            n1,
            n2,
        );
        self.out.write_all(&self.row)?;
        Ok(())
    }

    /// Append one GT row (with its trailing newline) into `buf`.
    ///
    /// The twelve fields are integers, so they are appended directly —
    /// `writeln!` routed each one through `fmt::Display`, which is ~14
    /// formatter dispatches per pair for 6.1 MB per 200 000-pair unit of
    /// output. The bytes are identical to the `writeln!` form: same order,
    /// same tabs, same clamping to `read_len`.
    #[allow(clippy::too_many_arguments)]
    pub fn row_into(
        buf: &mut Vec<u8>,
        read_len: u32,
        single: bool,
        pair_idx: u64,
        insert_len: u32,
        a: &Artefact,
        is_dup: bool,
        n1: u32,
        n2: u32,
    ) {
        // Clamp artefact lengths to read_len for display (full pairs may
        // have t = L).
        let l = read_len;
        let full = u8::from(a.full);
        let dup = u8::from(is_dup);
        buf.clear();
        push_u64(buf, pair_idx);
        buf.push(b'\t');
        push_u64(buf, insert_len as u64);
        buf.push(b'\t');
        push_u64(buf, a.h1.min(l) as u64);
        buf.push(b'\t');
        push_u64(buf, a.t1.min(l) as u64);
        buf.push(b'\t');
        if single {
            push_u64(buf, a.a1.min(l) as u64);
        } else {
            push_u64(buf, a.h2.min(l) as u64);
            buf.push(b'\t');
            push_u64(buf, a.t2.min(l) as u64);
            buf.push(b'\t');
            push_u64(buf, a.a1.min(l) as u64);
            buf.push(b'\t');
            push_u64(buf, a.a2.min(l) as u64);
        }
        buf.push(b'\t');
        push_u64(buf, full as u64);
        buf.push(b'\t');
        push_u64(buf, dup as u64);
        buf.push(b'\t');
        push_u64(buf, n1 as u64);
        if !single {
            buf.push(b'\t');
            push_u64(buf, n2 as u64);
        }
        buf.push(b'\n');
    }

    pub fn finish(mut self) -> Result<()> {
        self.out.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn artefact(h1: u32, t1: u32, h2: u32, t2: u32, full: bool) -> Artefact {
        Artefact {
            h1,
            t1,
            h2,
            t2,
            a1: 0,
            a2: 0,
            full,
        }
    }

    #[test]
    fn gt_row_matches_injected_values() {
        let path = std::env::temp_dir().join("fastsim_gt.tsv");
        let mut w = GtWriter::open(&path, 150).unwrap();
        let a = artefact(2, 5, 0, 3, false);
        // an injected N overlay on R1 only: the r1_n column must carry it
        w.write_pair(7, 300, &a, false, 3, 0).unwrap();
        w.write_pair(8, 120, &artefact(0, 0, 0, 0, true), true, 0, 0)
            .unwrap();
        w.finish().unwrap();

        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(
            lines[0],
            "pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter\tfull_flag\tis_dup\tr1_n\tr2_n"
        );
        assert_eq!(lines[1], "7\t300\t2\t5\t0\t3\t0\t0\t0\t0\t3\t0");
        assert_eq!(lines[2], "8\t120\t0\t0\t0\t0\t0\t0\t1\t1\t0\t0");
    }

    #[test]
    fn full_pair_lengths_clamped_to_read_len() {
        let path = std::env::temp_dir().join("fastsim_gt2.tsv");
        let mut w = GtWriter::open(&path, 150).unwrap();
        // full spike: t = L
        w.write_pair(0, 300, &artefact(0, 150, 0, 150, true), false, 0, 0)
            .unwrap();
        w.finish().unwrap();
        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(
            s.lines().nth(1).unwrap(),
            "0\t300\t0\t150\t0\t150\t0\t0\t1\t0\t0\t0"
        );
    }

    /// FS-10: the appended row must be byte-identical to what `writeln!` with
    /// twelve `Display` arguments produced, including the `min(read_len)`
    /// clamping and both column layouts.
    #[test]
    fn row_append_matches_the_format_macro() {
        let l = 150u32;
        let cases: Vec<(u64, u32, Artefact, bool, u32, u32)> = vec![
            (0, 300, artefact(2, 5, 0, 3, false), false, 3, 0),
            (7, 120, artefact(0, 0, 0, 0, true), true, 0, 0),
            (9_999, 20, artefact(150, 150, 149, 3, false), false, 0, 12),
            (123_456, 288, artefact(200, 151, 0, 64, true), true, 41, 41),
            (u64::MAX, 1_000_000, artefact(0, 150, 0, 150, false), false, 41, 42),
        ];
        for (pair_idx, insert_len, a, is_dup, n1, n2) in cases {
            let full = u8::from(a.full);
            let dup = u8::from(is_dup);
            let want = format!(
                "{pair_idx}\t{insert_len}\t{}\t{}\t{}\t{}\t{}\t{}\t{full}\t{dup}\t{n1}\t{n2}\n",
                a.h1.min(l), a.t1.min(l), a.h2.min(l), a.t2.min(l), a.a1.min(l), a.a2.min(l)
            );
            let mut buf = Vec::new();
            GtWriter::<Vec<u8>>::row_into(
                &mut buf, l, false, pair_idx, insert_len, &a, is_dup, n1, n2,
            );
            assert_eq!(String::from_utf8(buf).unwrap(), want, "paired row {pair_idx}");

            let want_single = format!(
                "{pair_idx}\t{insert_len}\t{}\t{}\t{}\t{full}\t{dup}\t{n1}\n",
                a.h1.min(l), a.t1.min(l), a.a1.min(l)
            );
            let mut buf = Vec::new();
            GtWriter::<Vec<u8>>::row_into(&mut buf, l, true, pair_idx, insert_len, &a, is_dup, n1, n2);
            assert_eq!(String::from_utf8(buf).unwrap(), want_single, "single row {pair_idx}");
        }
    }

    #[test]
    fn row_count_matches_fastq_pairs() {
        let path = std::env::temp_dir().join("fastsim_gt3.tsv");
        let mut w = GtWriter::open(&path, 150).unwrap();
        let n = 1000u64;
        for i in 0..n {
            w.write_pair(i, 300, &artefact(0, 0, 0, 0, i % 7 == 0), false, 0, 0)
                .unwrap();
        }
        w.finish().unwrap();
        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(s.lines().count(), n as usize + 1, "header + n rows");
    }
}
