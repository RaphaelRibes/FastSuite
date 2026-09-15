//! Ground-truth writer: `prefix.gt.tsv`, one row per simulated pair.
//!
//! The GT defines the exact correct answer for a trimmer, so a scorer can
//! compute trim accuracy absolutely (not relative to another tool):
//!
//! ```text
//! pair_id  insert  r1_h  r1_t  r2_h  r2_t  full_flag  is_dup
//! ```
//!
//! - `pair_id`: unique integer (matches the FASTQ coordinate x:y).
//! - `insert`: true fragment length I.
//! - `r1_h` / `r1_t`: 5′/3′ G-run lengths of R1; same for R2.
//! - `full_flag`: 1 = entirely artefact (full-length spike or biology
//!   window < 1 bp) — the correct answer is to DROP the pair; 0 = trim.
//! - `is_dup`: 1 = PCR duplicate of a previously emitted fragment (same
//!   template, independently re-sequenced: fresh errors, quality and
//!   polyG artefacts).
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

pub struct GtWriter<W: Write> {
    out: W,
    read_len: u32,
}

impl GtWriter<BufWriter<std::fs::File>> {
    pub fn open(path: &Path, read_len: u32) -> Result<Self> {
        let f =
            std::fs::File::create(path).with_context(|| format!("create {}", path.display()))?;
        let mut w = GtWriter {
            out: BufWriter::new(f),
            read_len,
        };
        w.out.write_all(
            b"pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tfull_flag\tis_dup\n",
        )?;
        Ok(w)
    }
}

impl<W: Write> GtWriter<W> {
    pub fn write_pair(
        &mut self,
        pair_idx: u64,
        insert_len: u32,
        a: &Artefact,
        is_dup: bool,
    ) -> Result<()> {
        // Clamp artefact lengths to read_len for display (full pairs may
        // have t = L).
        let l = self.read_len;
        let full = u8::from(a.full);
        let dup = u8::from(is_dup);
        writeln!(
            self.out,
            "{pair_idx}\t{insert_len}\t{}\t{}\t{}\t{}\t{full}\t{dup}",
            a.h1.min(l),
            a.t1.min(l),
            a.h2.min(l),
            a.t2.min(l)
        )?;
        Ok(())
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
            full,
        }
    }

    #[test]
    fn gt_row_matches_injected_values() {
        let path = std::env::temp_dir().join("fastsim_gt.tsv");
        let mut w = GtWriter::open(&path, 150).unwrap();
        let a = artefact(2, 5, 0, 3, false);
        w.write_pair(7, 300, &a, false).unwrap();
        w.write_pair(8, 120, &artefact(0, 0, 0, 0, true), true)
            .unwrap();
        w.finish().unwrap();

        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(
            lines[0],
            "pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tfull_flag\tis_dup"
        );
        assert_eq!(lines[1], "7\t300\t2\t5\t0\t3\t0\t0");
        assert_eq!(lines[2], "8\t120\t0\t0\t0\t0\t1\t1");
    }

    #[test]
    fn full_pair_lengths_clamped_to_read_len() {
        let path = std::env::temp_dir().join("fastsim_gt2.tsv");
        let mut w = GtWriter::open(&path, 150).unwrap();
        // full spike: t = L
        w.write_pair(0, 300, &artefact(0, 150, 0, 150, true), false)
            .unwrap();
        w.finish().unwrap();
        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(s.lines().nth(1).unwrap(), "0\t300\t0\t150\t0\t150\t1\t0");
    }

    #[test]
    fn row_count_matches_fastq_pairs() {
        let path = std::env::temp_dir().join("fastsim_gt3.tsv");
        let mut w = GtWriter::open(&path, 150).unwrap();
        let n = 1000u64;
        for i in 0..n {
            w.write_pair(i, 300, &artefact(0, 0, 0, 0, i % 7 == 0), false)
                .unwrap();
        }
        w.finish().unwrap();
        let s = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(s.lines().count(), n as usize + 1, "header + n rows");
    }
}
