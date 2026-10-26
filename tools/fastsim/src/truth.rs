//! Ground-truth writer: `prefix.gt.tsv`, one row per simulated pair.
//!
//! The GT defines the exact correct answer for a trimmer, so a scorer can
//! compute trim accuracy absolutely (not relative to another tool):
//!
//! ```text
//! pair_id  insert  r1_h  r1_t  r2_h  r2_t  r1_ad  r2_ad  full_flag  is_dup  r1_n  r2_n  r1_full  r2_full
//! ```
//!
//! - `r1_n` / `r2_n`: ambiguous (`N`) bases in that mate's emitted record,
//!   EXCLUDING the no-adapter length filler. Two sources produce them: the
//!   sequencing overlay's no-call model (an instrument no-call recorded INSIDE
//!   the keep window — it never touches `h`/`t`) and a template `N` inherited
//!   from a `--reference` contig. Both are in the read the scorer sees, so both
//!   are recorded here. The `N` filler of the no-adapter mode is NOT a no-call
//!   and is not counted here (it is reported separately in the run's
//!   `adapter_pad_bases` metric).
//!
//! - `pair_id`: unique integer (matches the FASTQ coordinate x:y).
//! - `insert`: true fragment length I.
//! - `r1_h` / `r1_t`: 5′/3′ G-run lengths of R1 **as drawn** — the fragment
//!   truth. Sequencing errors are applied after the artefact is decided, so a
//!   substitution inside an injected run can make the emitted run differ from
//!   `h`; scoring that difference is the trimmer's error, not a GT mismatch.
//!   Same for R2.
//! - `r1_full` / `r2_full`: 1 = that mate is entirely artefact
//!   (`h + t >= L`, from the drawn extents). A trimmer must drop it.
//! - `full_flag`: 1 = the PAIR is entirely artefact (`r1_full && r2_full`,
//!   identical to `r1_full` in `--single` mode) — the correct answer is to
//!   DROP the pair; 0 = trim. Per-mate drops read `r1_full`/`r2_full`.
//! - `is_dup`: 1 = PCR duplicate of a previously emitted pair: the whole
//!   stored read (template, artefact, substitution errors and quality) is
//!   replayed byte-for-byte.
//!
//! The correct trimmed window per mate is `[h, h + min(I, L−h−t))`:
//! - `[0, h)`     — 5′ G-run (artefact)
//! - `[h, h + min(I, L−h−t))` — biology (keep)
//! - `[h + min(I, L−h−t), L−t)` — adapter readthrough (artefact, only when
//!   the insert is shorter than the window `L−h−t`)
//! - `[L−t, L)`   — 3′ G-run (artefact)
//!
//! So `h`, `t`, `I` and `L` (from the profile) fully determine truth; no
//! per-base annotation needed.

use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::artefact::Artefact;
use crate::fmtint::push_u64;
use crate::profile::Mate;

/// Is one mate's emitted record entirely artefact (`h + t >= L`)?
///
/// The single definition of "full" for a mate: the GT's `r1_full`/`r2_full`
/// columns and `full_flag` are computed from it, and the run's
/// `full_artefact_pairs` counters call it too, so a published count can never
/// disagree with the rows a scorer reads.
pub fn mate_is_full(a: &Artefact, read_len: u32, mate: Mate) -> bool {
    let (h, t) = match mate {
        Mate::R1 => (a.h1, a.t1),
        Mate::R2 => (a.h2, a.t2),
    };
    h.min(read_len) + t.min(read_len) >= read_len
}

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
            b"pair_id\tinsert\tr1_h\tr1_t\tr1_adapter\tfull_flag\tis_dup\tr1_n\tr1_full\n"
        } else {
            b"pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter\tfull_flag\tis_dup\tr1_n\tr2_n\tr1_full\tr2_full\n"
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
        // "Entirely artefact" is measured per mate from the EMITTED record the
        // GT names: `a.h1`/`a.t1` are the emitted runs after substitutions, so
        // a row can never claim a run the FASTQ does not contain.
        let f1: u8 = u8::from(mate_is_full(a, l, Mate::R1));
        let f2: u8 = u8::from(mate_is_full(a, l, Mate::R2));
        // The pair must be dropped only when BOTH mates are useless; a pair
        // with one full mate still has valid biology in the other.
        let full: u8 = if single { f1 } else { f1 & f2 };
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
        if single {
            buf.push(b'\t');
            push_u64(buf, f1 as u64);
        } else {
            buf.push(b'\t');
            push_u64(buf, n2 as u64);
            buf.push(b'\t');
            push_u64(buf, f1 as u64);
            buf.push(b'\t');
            push_u64(buf, f2 as u64);
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
        w.write_pair(8, 120, &artefact(150, 0, 0, 150, true), true, 0, 0)
            .unwrap();
        w.finish().unwrap();

        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(
            lines[0],
            "pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter\tfull_flag\tis_dup\tr1_n\tr2_n\tr1_full\tr2_full"
        );
        assert_eq!(lines[1], "7\t300\t2\t5\t0\t3\t0\t0\t0\t0\t3\t0\t0\t0");
        assert_eq!(lines[2], "8\t120\t150\t0\t0\t150\t0\t0\t1\t1\t0\t0\t1\t1");
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
            "0\t300\t0\t150\t0\t150\t0\t0\t1\t0\t0\t0\t1\t1"
        );
    }

    /// FS-10: the appended row must be byte-identical to the `writeln!` form,
    /// including the `min(read_len)` clamping and both column layouts. The
    /// expected strings are literal so the per-mate `full` flags are pinned
    /// independently of the writer's arithmetic.
    #[test]
    fn row_append_matches_the_format_macro() {
        let l = 150u32;
        let cases: Vec<(u64, u32, Artefact, bool, u32, u32, &str, &str)> = vec![
            (
                0,
                300,
                artefact(2, 5, 0, 3, false),
                false,
                3,
                0,
                "0\t300\t2\t5\t0\t3\t0\t0\t0\t0\t3\t0\t0\t0\n",
                "0\t300\t2\t5\t0\t0\t0\t3\t0\n",
            ),
            (
                7,
                120,
                artefact(150, 0, 0, 150, true),
                true,
                0,
                0,
                "7\t120\t150\t0\t0\t150\t0\t0\t1\t1\t0\t0\t1\t1\n",
                "7\t120\t150\t0\t0\t1\t1\t0\t1\n",
            ),
            (
                9_999,
                20,
                artefact(150, 150, 149, 3, false),
                false,
                0,
                12,
                "9999\t20\t150\t150\t149\t3\t0\t0\t1\t0\t0\t12\t1\t1\n",
                "9999\t20\t150\t150\t0\t1\t0\t0\t1\n",
            ),
            (
                123_456,
                288,
                artefact(200, 151, 0, 64, true),
                true,
                41,
                41,
                "123456\t288\t150\t150\t0\t64\t0\t0\t0\t1\t41\t41\t1\t0\n",
                "123456\t288\t150\t150\t0\t1\t1\t41\t1\n",
            ),
            (
                u64::MAX,
                1_000_000,
                artefact(0, 150, 0, 150, false),
                false,
                41,
                42,
                "18446744073709551615\t1000000\t0\t150\t0\t150\t0\t0\t1\t0\t41\t42\t1\t1\n",
                "18446744073709551615\t1000000\t0\t150\t0\t1\t0\t41\t1\n",
            ),
        ];
        for (pair_idx, insert_len, a, is_dup, n1, n2, want, want_single) in cases {
            let mut buf = Vec::new();
            GtWriter::<Vec<u8>>::row_into(
                &mut buf, l, false, pair_idx, insert_len, &a, is_dup, n1, n2,
            );
            assert_eq!(String::from_utf8(buf).unwrap(), want, "paired row {pair_idx}");

            let mut buf = Vec::new();
            GtWriter::<Vec<u8>>::row_into(&mut buf, l, true, pair_idx, insert_len, &a, is_dup, n1, n2);
            assert_eq!(
                String::from_utf8(buf).unwrap(),
                want_single,
                "single row {pair_idx}"
            );
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
