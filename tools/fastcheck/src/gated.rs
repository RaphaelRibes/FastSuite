//! Module gating for the per-read measurement (`--only`, `--no-quality`,
//! `--no-heatmap`).
//!
//! Gating means the work is *not done*, not that its answer is hidden: with
//! `--only dup` nothing walks the quality bytes, nothing grows a position row and
//! nothing touches the 256x256 grid. `Accum::observe` remains the ungated,
//! all-modules path and is still what a run with no gating flag takes, so the
//! default report is unchanged byte for byte; [`observe_gated`] is its gated twin,
//! used only when a gating flag was given.
//!
//! Why this lives outside `accum.rs`: that file is the dense per-read accumulator
//! and the gated twin is far easier to hold to it from beside it than from inside
//! it. Rust allows an inherent `impl` for a type from another module of the
//! same crate, and every `Accum` field is `pub`, so the gated path can be written
//! beside the type without bloating it. The handful of helpers copied from
//! `accum.rs` (`grow`, `record_n_run`, the adapter-prefix search) are private
//! there; `gated_all_modules_is_exactly_observe` below holds this copy to the
//! original counter for counter — including the `tile` module — so a divergence
//! cannot pass unnoticed.
//!
//! [`observe_gated`]: Accum::observe_gated

use crate::accum::{contains_ci, qbin, Accum, QBYTES, QBINS};
use crate::cli::Modules;
use crate::heatmap;
use crate::polyx;

/// `[u64; QBINS]` has no `Default`, so rows are grown with an explicit literal
/// (the same one `accum.rs` uses).
const QUAL_ROW: [u64; QBINS] = [0; QBINS];

/// The ambiguous-base code (A=0 C=1 G=2 T=3, everything else 4) — the same
/// mapping `Accum::observe` applies.
fn base_code(raw: u8) -> usize {
    match raw.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // N and every other IUPAC code: counted, never folded into ACGT
    }
}

fn grow<T: Clone + Default>(v: &mut Vec<T>, n: usize) {
    if v.len() < n {
        v.resize(n, T::default());
    }
}

/// Record one ambiguous run `[start, end)` (copy of `Accum::record_n_run`, which
/// is private to `accum`).
fn record_n_run(a: &mut Accum, start: usize, end: usize, read_len: usize) {
    let len = end - start;
    grow(&mut a.n_run_hist, len + 1);
    a.n_run_hist[len] += 1;
    if start == 0 {
        a.n_run_leading += 1;
    }
    if end >= read_len {
        a.n_run_trailing += 1;
    }
}

// Case-insensitive substring search: the gated twin calls `accum::contains_ci` rather
// than carrying a copy, so the two paths cannot drift (the copy that lived here had to
// be kept byte-identical by hand).

/// What a gated pass actually measured, in units of work rather than of output.
///
/// The claim that gating saves work is made against these counts, never against a
/// wall-clock time: every field is incremented where the operation is performed,
/// so a counter that reads zero is a loop that did not run.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Work {
    /// Reads handed to the gated accumulator.
    pub reads: u64,
    /// `base_pos` cells written (the `bases` module).
    pub base_cells: u64,
    /// `qual_pos` cells written (the `quality` module).
    pub qual_cells: u64,
    /// Quality bytes decoded into the census and the per-read mean (the `quality` module).
    pub quality_bytes: u64,
    /// Reads whose polyX runs were measured (the `polyx` module).
    pub polyx_reads: u64,
    /// Heatmap cells incremented (the `heatmap` module).
    pub heatmap_cells: u64,
    /// Reads fingerprinted (the `dup` module).
    pub fingerprints: u64,
    /// Distinct 12-mer windows inserted (the `kmer` module).
    pub kmer_windows: u64,
}

impl Work {
    pub fn add(&mut self, o: &Work) {
        self.reads += o.reads;
        self.base_cells += o.base_cells;
        self.qual_cells += o.qual_cells;
        self.quality_bytes += o.quality_bytes;
        self.polyx_reads += o.polyx_reads;
        self.heatmap_cells += o.heatmap_cells;
        self.fingerprints += o.fingerprints;
        self.kmer_windows += o.kmer_windows;
    }

    /// One line, for the run's diagnostic output.
    pub fn summary(&self) -> String {
        format!(
            "{} reads measured: {} base cells, {} quality cells, {} quality bytes, \
             {} fingerprints, {} polyX reads, {} heatmap cells, {} 12-mer windows",
            self.reads,
            self.base_cells,
            self.qual_cells,
            self.quality_bytes,
            self.fingerprints,
            self.polyx_reads,
            self.heatmap_cells,
            self.kmer_windows,
        )
    }
}

impl Accum {
    /// One record, measuring only what `m` enables — `observe`, with every block
    /// behind its module.
    ///
    /// `reads`, `bases` and `reads_over_max_len` are always counted, and so is
    /// `reads_with_n` (a run that reports no reads is not a run); everything else
    /// is conditional. The quality *line* is still delivered here by the scanner
    /// for framing, but with `quality` off its bytes are never touched. `header`
    /// carries the tile coordinates the `tile` module needs.
    pub fn observe_gated(&mut self, header: &[u8], seq: &[u8], qual: &[u8], m: &Modules) -> Work {
        let mut w = Work { reads: 1, ..Work::default() };
        let len = seq.len();
        self.reads += 1;
        self.bases += len as u64;

        if m.length {
            grow(&mut self.len_hist, len + 1);
            self.len_hist[len] += 1;
        }

        let per_pos = len.min(self.max_pos);
        if len > self.max_pos {
            self.reads_over_max_len += 1;
        }
        if m.bases {
            grow(&mut self.base_pos, per_pos);
        }
        if m.quality && self.qual_pos.len() < per_pos {
            self.qual_pos.resize(per_pos, QUAL_ROW);
        }

        // Composition and the ambiguous-run model share one pass over the bases. With
        // neither module on, the only thing the sequence is still needed for is the
        // "does this read contain an N" flag, which is one cheap scan.
        if m.bases || m.gc {
            let (mut acgt, mut gc) = (0u64, 0u64);
            let mut has_n = false;
            let mut run_start = 0usize;
            let mut in_run = false;
            for (i, &raw) in seq.iter().enumerate() {
                let code = base_code(raw);
                if code == 4 {
                    if !in_run {
                        in_run = true;
                        run_start = i;
                    }
                    has_n = true;
                } else {
                    if in_run {
                        if m.bases {
                            record_n_run(self, run_start, i, len);
                        }
                        in_run = false;
                    }
                    acgt += 1;
                    if code == 1 || code == 2 {
                        gc += 1;
                    }
                }
                if i < per_pos && m.bases {
                    self.base_pos[i][code] += 1;
                    w.base_cells += 1;
                }
            }
            if in_run && m.bases {
                record_n_run(self, run_start, len, len);
            }
            if has_n {
                self.reads_with_n += 1;
            }
            if m.gc {
                grow(&mut self.acgt_hist, acgt as usize + 1);
                self.acgt_hist[acgt as usize] += 1;
                grow(&mut self.gc_hist, gc as usize + 1);
                self.gc_hist[gc as usize] += 1;
            }
        } else if seq.iter().any(|&raw| base_code(raw) == 4) {
            self.reads_with_n += 1;
        }

        if m.quality {
            let mut qsum = 0u64;
            for (i, &q) in qual.iter().enumerate() {
                self.qual_byte[q as usize] += 1;
                qsum += q as u64;
                match qbin(q) {
                    Some(b) => {
                        if i < per_pos {
                            self.qual_pos[i][b] += 1;
                            w.qual_cells += 1;
                        }
                    }
                    // Not a Phred byte: in the census, never in a per-position row.
                    None => self.qual_out_of_range += 1,
                }
            }
            w.quality_bytes = qual.len() as u64;
            let raw_mean = if qual.is_empty() {
                0
            } else {
                (qsum / qual.len() as u64) as usize
            };
            self.meanq_raw_hist[raw_mean.min(QBYTES - 1)] += 1;
        }

        // The --profile extras are not a module: they are the price of the profile
        // the user asked for, and main.rs refuses --profile when a module it needs
        // is off, so they stay ungated here.
        if self.track_profile {
            let mut hit = false;
            for (i, (_, pat)) in crate::accum::ADAPTER_PREFIXES.iter().enumerate() {
                if contains_ci(seq, pat.as_bytes()) {
                    self.adapter_hits[i] += 1;
                    hit = true;
                }
            }
            if !hit {
                let start = seq
                    .iter()
                    .position(|&b| b.to_ascii_uppercase() != b'G')
                    .unwrap_or(seq.len());
                for &b in &seq[start..] {
                    match b.to_ascii_uppercase() {
                        b'G' | b'C' => {
                            self.frag_gc += 1;
                            self.frag_bases += 1;
                        }
                        b'A' | b'T' => {
                            self.frag_bases += 1;
                        }
                        _ => {}
                    }
                }
            }
        }

        if m.polyx {
            let (lb, ll, tb, tl) = polyx::measure_polyx(seq);
            w.polyx_reads = 1;
            if lb != usize::MAX {
                grow(&mut self.polyx5[lb], ll + 1);
                self.polyx5[lb][ll] += 1;
            }
            if tb != usize::MAX {
                grow(&mut self.polyx3[tb], tl + 1);
                self.polyx3[tb][tl] += 1;
            }
            // Quality INSIDE each G-run, per mate: the profile's run_q_mean reads
            // these, and they are the polyX measurement's own output.
            const G: usize = 2; // polyx::base_code(b'G')
            if lb == G {
                let n = ll.min(qual.len());
                self.run_q_head.0 += qual[..n].iter().map(|&b| b as u64).sum::<u64>();
                self.run_q_head.1 += n as u64;
            }
            if tb == G {
                let start = qual.len().saturating_sub(tl);
                self.run_q_tail.0 += qual[start..].iter().map(|&b| b as u64).sum::<u64>();
                self.run_q_tail.1 += (qual.len() - start) as u64;
            }
        }

        if m.heatmap {
            let (hx, hy) = heatmap::encode(seq);
            self.heatmap[heatmap::index(hx, hy)] += 1;
            w.heatmap_cells = 1;
        }

        if m.dup {
            self.seq_counts.observe(seq);
            w.fingerprints = 1;
            if self.track_overrep_table {
                // See accum::observe: the scan refuses a sequence that table mode cannot
                // pack (reader::scan_file_gated), so this only ever sees ACGT and observe
                // always packs it. It is not a silent drop.
                if !self.overrep_table.observe(seq) {
                    debug_assert!(
                        false,
                        "reader::scan_file_gated must refuse a sequence table mode cannot pack"
                    );
                }
            }
        }

        if m.adapter {
            if let Some(ap) = self.adapter_pos.as_mut() {
                ap.observe(seq, &self.adapter_set);
            }
        }

        // Per-tile quality: gated on its own module, so it is skipped when --tiles is
        // off even under an --only that selected other modules.
        if m.tile {
            if let Some(t) = self.tile.as_mut() {
                t.observe(header, qual);
            }
        }

        if m.kmer && seq.len() >= 12 {
            // The SAME rolling register the fused path uses (accum::KmerRoll). The
            // pre-fusion loop that used to live here (`seq.windows(12)` +
            // `pack_kmer`) decoded each base 11 more times, which made `--only kmer`
            // 2.6x slower than running every module at once.
            let mut roll = crate::accum::KmerRoll::default();
            for &raw in seq {
                if let Some(k) = roll.push(crate::accum::BASE_CODE[raw as usize]) {
                    self.kmer_counts.add(k);
                    self.kmer_windows += 1;
                    w.kmer_windows += 1;
                }
            }
        }

        w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Modules;

    /// The gated path with every module on must be the ungated path, counter for
    /// counter — otherwise "no gating flag = today's behaviour" would not hold the
    /// moment a flag is added on the command line.
    #[test]
    fn gated_all_modules_is_exactly_observe() {
        let reads: [(&[u8], &[u8], &[u8]); 4] = [
            (b"@i:1:f:1:1101:1:1 1:N:0:x", b"ACGTACGTAC", b"IIIIIIIIII"),
            (b"@i:1:f:1:1101:2:2 1:N:0:x", b"AANNGGCCTT", b"!!!!55555."),
            (b"@i:1:f:1:2202:3:3 1:N:0:x", b"GGGGNNNNNN", b"##########"),
            (b"@no_coords", b"ACGTACGTAC", b"JJJJJJJJJJ"),
        ];
        let mut plain = Accum::new(64);
        let mut gated = Accum::new(64);
        gated.track_profile = true;
        gated.adapter_set = crate::adapter::AdapterSet::illumina_default();
        gated.adapter_pos = Some(crate::adapter_pos::AdapterPos::new(&gated.adapter_set));
        plain.adapter_set = gated.adapter_set.clone();
        plain.adapter_pos = Some(crate::adapter_pos::AdapterPos::new(&plain.adapter_set));
        plain.track_profile = true;
        plain.tile = Some(crate::tile::TileQual::default());
        gated.tile = Some(crate::tile::TileQual::default());
        let all = Modules {
            adapter: true,
            tile: true,
            ..Modules::defaults(true)
        };
        for (header, seq, qual) in reads {
            plain.observe(header, seq, qual);
            let w = gated.observe_gated(header, seq, qual, &all);
            assert_eq!(w.reads, 1);
        }
        assert_eq!(plain.reads, gated.reads);
        assert_eq!(plain.bases, gated.bases);
        assert_eq!(plain.reads_with_n, gated.reads_with_n);
        assert_eq!(plain.reads_over_max_len, gated.reads_over_max_len);
        assert_eq!(plain.n_run_hist, gated.n_run_hist);
        assert_eq!(plain.n_run_leading, gated.n_run_leading);
        assert_eq!(plain.n_run_trailing, gated.n_run_trailing);
        assert_eq!(plain.len_hist, gated.len_hist);
        assert_eq!(plain.acgt_hist, gated.acgt_hist);
        assert_eq!(plain.gc_hist, gated.gc_hist);
        assert_eq!(plain.meanq_raw_hist, gated.meanq_raw_hist);
        assert_eq!(plain.qual_byte, gated.qual_byte);
        assert_eq!(plain.base_pos, gated.base_pos);
        assert_eq!(plain.qual_pos, gated.qual_pos);
        assert_eq!(plain.polyx3, gated.polyx3);
        assert_eq!(plain.polyx5, gated.polyx5);
        assert_eq!(plain.run_q_head, gated.run_q_head);
        assert_eq!(plain.run_q_tail, gated.run_q_tail);
        assert_eq!(plain.heatmap, gated.heatmap);
        assert_eq!(plain.adapter_hits, gated.adapter_hits);
        assert_eq!(plain.frag_gc, gated.frag_gc);
        assert_eq!(plain.frag_bases, gated.frag_bases);
        assert_eq!(plain.seq_counts.distinct(), gated.seq_counts.distinct());
        assert_eq!(plain.kmer_counts.report(1), gated.kmer_counts.report(1));
        assert_eq!(
            plain.adapter_pos.as_ref().unwrap().to_json(0.2, 14),
            gated.adapter_pos.as_ref().unwrap().to_json(0.2, 14)
        );
        assert_eq!(
            plain.tile.as_ref().unwrap().to_json(33),
            gated.tile.as_ref().unwrap().to_json(33),
            "the tile module must be the same gated and ungated"
        );
    }

    /// The gated path measures what its modules say and nothing else: the counters
    /// are incremented where the work happens, so a zero here is a loop that never
    /// ran. This is the measurement the "gating saves work" test rests on.
    #[test]
    fn gating_only_measures_the_selected_modules() {
        let (seq, qual) = (&b"ACGTACGTAC"[..], &b"IIIIIIIIII"[..]);
        let only = |set: fn(&mut Modules)| {
            let m = {
                let mut m = Modules::none();
                set(&mut m);
                m
            };
            let mut a = Accum::new(64);
            let w = a.observe_gated(b"@i:1:f:1:1101:1:1 1:N:0:x", seq, qual, &m);
            (a, w)
        };

        let (a, w) = only(|m| m.dup = true);
        assert_eq!(w.fingerprints, 1, "dup fingerprints every read");
        assert_eq!((w.base_cells, w.qual_cells), (0, 0), "no position row was grown");
        assert_eq!(w.quality_bytes, 0, "no quality byte was decoded");
        assert_eq!(w.kmer_windows, 0);
        assert_eq!(a.base_pos.len(), 0, "bases off: no row allocated");
        assert_eq!(a.qual_pos.len(), 0, "quality off: no row allocated");
        assert_eq!(a.heatmap.iter().sum::<u64>(), 0, "heatmap off: no cell touched");
        assert_eq!(a.reads, 1, "reads and bases are always counted");
        assert_eq!(a.bases, 10);

        let (a, w) = only(|m| m.bases = true);
        assert_eq!(w.base_cells, 10, "one cell per base");
        assert_eq!(w.fingerprints, 0);
        assert_eq!(a.base_pos.len(), 10);
        assert_eq!(a.qual_pos.len(), 0);

        let (_, w) = only(|m| m.quality = true);
        assert_eq!(w.qual_cells, 10);
        assert_eq!(w.quality_bytes, 10);
        assert_eq!(w.base_cells, 0);

        let (_, w) = only(|m| m.kmer = true);
        assert_eq!(w.kmer_windows, 0, "10 bp is shorter than a 12-mer");
        let (_, w) = only(|m| {
            m.kmer = true;
            m.bases = true;
        });
        assert_eq!(w.kmer_windows, 0);
        let long = &b"ACGTACGTACGTACGT"[..];
        let mut a = Accum::new(64);
        let w = a.observe_gated(b"@r", long, &long[..], &Modules { kmer: true, ..Modules::none() });
        assert_eq!(w.kmer_windows, 5, "16 bp carries 5 distinct windows");
    }
}
