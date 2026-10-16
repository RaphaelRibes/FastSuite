//! Exact, one-pass accumulators.
//!
//! Every counter is an integer count. Nothing here estimates, samples or
//! thresholds. Memory depends on read length and the modules enabled, never on
//! the number of reads.
//!
//! Layout: position-major (`pos_row[pos]`), grown lazily to the longest read
//! actually seen. A plane-major layout (`plane[code * max_pos + pos]`) was
//! tried, on the theory that a read walks positions sequentially and would
//! then touch few cache lines. Measured on 2 M reads x 150 bp it was **25 %
//! slower** (pure scan 5.14 s vs 4.10 s at `--max-read-len 1000`): the
//! quality plane spans `256 * max_read_len * 8 B` = 2 MB, so consecutive bases
//! land in different 8 KB rows and each read walks 2 MB instead of the
//! position-major working set of `read_len * 2 KB` = 300 KB. The measurement
//! wins; the theory is recorded here so nobody re-tries it by accident.
//!
//! Sizing at R = 150, per file:
//!   base_pos  R x 5 x 8 B   = 6 KB
//!   qual_pos  R x 256 x 8 B = 307 KB
//! `qual_pos` is indexed by the raw quality byte rather than by Phred, so the
//! accumulator never needs to know the encoding: the offset is decided once,
//! at report time, from the exact byte census. A run therefore cannot be
//! silently wrong about +33/+64.

pub const NT: usize = 5; // A, C, G, T, other/N
pub const QBINS: usize = 256;

fn grow<T: Clone + Default>(v: &mut Vec<T>, n: usize) {
    if v.len() < n {
        v.resize(n, T::default());
    }
}

// [u64; 256] has no `Default`, so rows are grown with an explicit literal.
const QUAL_ROW: [u64; QBINS] = [0; QBINS];


use crate::heatmap;
use crate::polyx;
use crate::repeated;
use crate::repeated::RepeatedCounter;
#[derive(Clone)]
pub struct Accum {
    /// Position rows are grown to the longest read seen; reads longer than
    /// `--max-read-len` are counted and flagged, never misattributed.
    pub max_pos: usize,
    pub reads: u64,
    pub bases: u64,
    pub reads_with_n: u64,
    pub reads_over_max_len: u64,
    /// Ambiguous-base (N/IUPAC) run-length histogram: `n_run_hist[len]` = number
    /// of maximal ambiguous runs of that length. `n_run_leading` /
    /// `n_run_trailing` count the runs touching the read's first / last base.
    /// Together with `reads_with_n` these are the three measured terms the
    /// simulator's N overlay is fitted from; the per-position `base_pos[..][4]`
    /// plane cannot recover a run length, because a run is not a per-position
    /// event.
    pub n_run_hist: Vec<u64>,
    pub n_run_leading: u64,
    pub n_run_trailing: u64,
    /// read length -> count
    pub len_hist: Vec<u64>,
    /// ACGT count per read -> count (N and other codes excluded; the stated
    /// denominator for %GC)
    pub acgt_hist: Vec<u64>,
    /// G+C count per read -> count
    pub gc_hist: Vec<u64>,
    /// floor(raw byte mean) of each read -> count. For an integer offset `k`,
    /// `floor(x) - k == floor(x - k)`, so the Phred offset is applied at report
    /// time without touching the data again.
    pub meanq_raw_hist: [u64; QBINS],
    /// Exact census of quality bytes; this is what decides +33 vs +64.
    pub qual_byte: [u64; QBINS],
    /// `base_pos[pos][code]`, position-major.
    pub base_pos: Vec<[u64; NT]>,
    /// `qual_pos[pos][byte]`, position-major.
    pub qual_pos: Vec<[u64; QBINS]>,
    /// Adapter-prefix hits per catalogue entry (TruSeq, Nextera), and the G+C of
    /// the adapter-free reads with the 5' G-run stripped. Both feed the fitter's
    /// readthrough rate and FRAGMENT gc; they are only tracked when a profile is
    /// requested, so the common path pays nothing.
    pub adapter_hits: [u64; 2],
    pub frag_gc: u64,
    pub frag_bases: u64,
    /// Whether the per-read profile scan (adapter prefix + fragment GC) runs.
    pub track_profile: bool,
    /// polyX run-length histograms: `polyx3[base][len]` (trailing 3′ run) and
    /// `polyx5[base][len]` (leading 5′ run), base 0=A 1=C 2=G 3=T.
    pub polyx3: [Vec<u64>; 4],
    pub polyx5: [Vec<u64>; 4],
    /// `(sum of raw quality bytes, base count)` for the bases lying INSIDE a
    /// leading (5′) / trailing (3′) G-run. The profile's flat in-run quality is
    /// the mean of these bytes minus the Phred offset — the same "store bytes,
    /// apply the offset at report time" rule the other quality counters here
    /// follow. Without them the fitted per-position curve cannot be split into
    /// its out-of-run and in-run parts, which is exactly what the simulator's
    /// quality model consumes (`pos_q` + `run_q_mean`).
    pub run_q_head: (u64, u64),
    pub run_q_tail: (u64, u64),
    /// The fastmaps-style heatmap grid (SIZE*SIZE, row-major y*SIZE+x).
    pub heatmap: Vec<u64>,
    /// Exact duplication counter over 64-bit fingerprints of the whole read. Quality is
    /// ignored (stated in the report). This is the baseline every --overrep-mode builds on
    /// (PLAN-QC §5.9): it gives the exact multiplicity histogram with no sequence text.
    pub seq_counts: crate::fingerprint::FingerprintTable,
    /// Bounded 12-mer counter (packed u64 keys).
    pub kmer_counts: RepeatedCounter<u64>,
    /// 2-bit packed sequence table, filled only in `--overrep-mode table` (the memory is
    /// real: ~46 B per unique 150 bp read, admission-checked in main.rs).
    pub overrep_table: crate::overrep::SequenceTable,
    /// Whether the packed sequence table is filled (set in main.rs for `table` mode).
    pub track_overrep_table: bool,
    /// Adapter content per position, only when --adapter-content is on.
    pub adapter_pos: Option<crate::adapter_pos::AdapterPos>,
    /// The adapter set for --adapter-content (only consulted when adapter_pos is Some).
    pub adapter_set: crate::adapter::AdapterSet,
}

/// Adapter prefixes the fitter's readthrough scan looks for.
pub const ADAPTER_PREFIXES: [(&str, &str); 2] =
    [("TruSeq", "AGATCGGAAGAGC"), ("Nextera", "CTGTCTCTTATA")];

/// Case-insensitive substring search (the adapter prefix in a read).
fn contains_ci(hay: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() || hay.len() < needle.len() {
        return false;
    }
    hay.windows(needle.len()).any(|w| {
        w.iter().zip(needle).all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
    })
}

impl Accum {
    pub fn new(max_pos: usize) -> Self {
        Self {
            max_pos,
            reads: 0,
            bases: 0,
            reads_with_n: 0,
            n_run_hist: Vec::new(),
            n_run_leading: 0,
            n_run_trailing: 0,
            reads_over_max_len: 0,
            len_hist: Vec::new(),
            acgt_hist: Vec::new(),
            gc_hist: Vec::new(),
            meanq_raw_hist: [0; QBINS],
            qual_byte: [0; QBINS],
            base_pos: Vec::new(),
            qual_pos: Vec::new(),
            adapter_hits: [0; 2],
            frag_gc: 0,
            frag_bases: 0,
            track_profile: false,
            polyx3: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            polyx5: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            run_q_head: (0, 0),
            run_q_tail: (0, 0),
            heatmap: vec![0; heatmap::SIZE * heatmap::SIZE],
            seq_counts: crate::fingerprint::FingerprintTable::with_capacity(1024),
            kmer_counts: RepeatedCounter::new(repeated::DEFAULT_CAP),
            overrep_table: crate::overrep::SequenceTable::default(),
            track_overrep_table: false,
            adapter_pos: None,
            adapter_set: crate::adapter::AdapterSet::illumina_default(),
        }
    }

    /// Accumulator bytes for the rows currently allocated — printed before the
    /// numbers so a run's analysis memory is known, not discovered.
    pub fn resident_bytes(&self) -> usize {
        self.base_pos.len() * std::mem::size_of::<[u64; NT]>()
            + self.qual_pos.len() * std::mem::size_of::<[u64; QBINS]>()
            + (self.len_hist.len() + self.acgt_hist.len() + self.gc_hist.len()) * 8
            + self.n_run_hist.len() * 8
            + 2 * QBINS * 8
    }

    /// Record one ambiguous run `[start, end)` of a read of length `read_len`.
    fn record_n_run(&mut self, start: usize, end: usize, read_len: usize) {
        let len = end - start;
        grow(&mut self.n_run_hist, len + 1);
        self.n_run_hist[len] += 1;
        if start == 0 {
            self.n_run_leading += 1;
        }
        if end >= read_len {
            self.n_run_trailing += 1;
        }
    }

    /// One record. `seq` and `qual` must have equal length (the scanner
    /// enforces it: a mismatch is an error, not an assumption).
    pub fn observe(&mut self, seq: &[u8], qual: &[u8]) {
        let len = seq.len();
        self.reads += 1;
        self.bases += len as u64;
        grow(&mut self.len_hist, len + 1);
        self.len_hist[len] += 1;

        let per_pos = len.min(self.max_pos);
        if len > self.max_pos {
            self.reads_over_max_len += 1;
        }
        grow(&mut self.base_pos, per_pos);
        if self.qual_pos.len() < per_pos {
            self.qual_pos.resize(per_pos, QUAL_ROW);
        }

        let (mut acgt, mut gc, mut has_n) = (0u64, 0u64, false);
        let mut run_start: usize = 0; // first base of the ambiguous run in progress
        let mut in_run = false;
        for (i, &raw) in seq.iter().enumerate() {
            let code = match raw.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4, // N and every other IUPAC code: counted, never folded into ACGT
            };
            if code == 4 {
                if !in_run {
                    in_run = true;
                    run_start = i;
                }
                has_n = true;
            } else {
                if in_run {
                    self.record_n_run(run_start, i, len);
                    in_run = false;
                }
                acgt += 1;
                if code == 1 || code == 2 {
                    gc += 1;
                }
            }
            if i < per_pos {
                self.base_pos[i][code] += 1;
            }
        }
        if in_run {
            self.record_n_run(run_start, len, len);
        }
        if has_n {
            self.reads_with_n += 1;
        }
        grow(&mut self.acgt_hist, acgt as usize + 1);
        self.acgt_hist[acgt as usize] += 1;
        grow(&mut self.gc_hist, gc as usize + 1);
        self.gc_hist[gc as usize] += 1;

        let mut qsum = 0u64;
        for (i, &q) in qual.iter().enumerate() {
            self.qual_byte[q as usize] += 1;
            qsum += q as u64;
            if i < per_pos {
                self.qual_pos[i][q as usize] += 1;
            }
        }
        let raw_mean = if qual.is_empty() {
            0
        } else {
            (qsum / qual.len() as u64) as usize
        };
        self.meanq_raw_hist[raw_mean.min(QBINS - 1)] += 1;

        // Profile extras (only when a profile was requested): the adapter
        // readthrough rate and the FRAGMENT gc. Readthrough is the adapter
        // prefix appearing anywhere in the read; the fragment gc is taken from
        // reads WITHOUT one, with the 5' G-run stripped — those are pure
        // fragment, while adapter-carrying reads are exactly the short-insert
        // ones and would weight the estimate by short-fragment composition.
        if self.track_profile {
            let mut hit = false;
            for (i, (_, pat)) in ADAPTER_PREFIXES.iter().enumerate() {
                if contains_ci(seq, pat.as_bytes()) {
                    self.adapter_hits[i] += 1;
                    hit = true;
                }
            }
            if !hit {
                let start = seq.iter().position(|&b| b.to_ascii_uppercase() != b'G').unwrap_or(seq.len());
                for &b in &seq[start..] {
                    match b.to_ascii_uppercase() {
                        b'G' | b'C' => { self.frag_gc += 1; self.frag_bases += 1; }
                        b'A' | b'T' => { self.frag_bases += 1; }
                        _ => {}
                    }
                }
            }
        }

        // polyX run lengths (always: cheap, and the report always wants them).
        let (lb, ll, tb, tl) = polyx::measure_polyx(seq);
        if lb != usize::MAX {
            grow(&mut self.polyx5[lb], ll + 1);
            self.polyx5[lb][ll] += 1;
        }
        if tb != usize::MAX {
            grow(&mut self.polyx3[tb], tl + 1);
            self.polyx3[tb][tl] += 1;
        }
        // Quality INSIDE each G-run, per mate. The fitter subtracts this flat
        // in-run mean out of the per-position curve to recover the out-of-run
        // curve the simulator wants, so it has to be measured, not inferred.
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

        // Heatmap cell.
        let (hx, hy) = heatmap::encode(seq);
        self.heatmap[heatmap::index(hx, hy)] += 1;

        // Duplication: every read is fingerprinted, exactly, with no sample and no cap.
        // The table is pre-sized from the projected record count in main.rs, so this does
        // not rehash. Over-representation naming is a separate, opt-in concern (§5.9).
        self.seq_counts.observe(seq);
        if self.track_overrep_table {
            // A non-ACGT base cannot be packed; refusing is the caller's job (main.rs
            // checks the first read), and here the sequence is simply skipped rather than
            // silently mis-counted — the report's unique_sequences count shows the gap.
            let _ = self.overrep_table.observe(seq);
        }
        if let Some(ap) = self.adapter_pos.as_mut() {
            ap.observe(seq, &self.adapter_set);
        }
        if seq.len() >= 12 {
            for w in seq.windows(12) {
                if let Some(pk) = repeated::pack_kmer(w) {
                    self.kmer_counts.insert(pk);
                }
            }
        }
    }

    /// Exact sum of two accumulators.
    pub fn merge(&mut self, o: &Accum) {
        self.max_pos = self.max_pos.max(o.max_pos);
        self.reads += o.reads;
        self.bases += o.bases;
        self.reads_with_n += o.reads_with_n;
        self.n_run_leading += o.n_run_leading;
        self.n_run_trailing += o.n_run_trailing;
        add_hist(&mut self.n_run_hist, &o.n_run_hist);
        self.reads_over_max_len += o.reads_over_max_len;
        self.adapter_hits[0] += o.adapter_hits[0];
        self.adapter_hits[1] += o.adapter_hits[1];
        self.frag_gc += o.frag_gc;
        self.frag_bases += o.frag_bases;
        for b in 0..4 {
            add_hist(&mut self.polyx3[b], &o.polyx3[b]);
            add_hist(&mut self.polyx5[b], &o.polyx5[b]);
        }
        self.run_q_head.0 += o.run_q_head.0;
        self.run_q_head.1 += o.run_q_head.1;
        self.run_q_tail.0 += o.run_q_tail.0;
        self.run_q_tail.1 += o.run_q_tail.1;
        for (i, &c) in o.heatmap.iter().enumerate() {
            self.heatmap[i] += c;
        }
        // The fingerprint table merges exactly, like every other accumulator here.
        self.seq_counts.merge(&o.seq_counts);
        match (self.adapter_pos.as_mut(), o.adapter_pos.as_ref()) {
            (Some(a), Some(b)) => a.merge(b),
            // The run total starts without the optional counter: the first file that
            // tracked it is ADOPTED whole (a clone of nothing-plus-itself), so its
            // counts are not added twice.
            (None, Some(b)) => self.adapter_pos = Some(b.clone()),
            _ => {}
        }
        if self.track_overrep_table || o.track_overrep_table {
            self.track_overrep_table = true;
            self.overrep_table.merge(&o.overrep_table);
        }
        add_hist(&mut self.len_hist, &o.len_hist);
        add_hist(&mut self.acgt_hist, &o.acgt_hist);
        add_hist(&mut self.gc_hist, &o.gc_hist);
        for i in 0..QBINS {
            self.meanq_raw_hist[i] += o.meanq_raw_hist[i];
            self.qual_byte[i] += o.qual_byte[i];
        }
        grow(&mut self.base_pos, o.base_pos.len());
        for (i, row) in o.base_pos.iter().enumerate() {
            for (d, s) in self.base_pos[i].iter_mut().zip(row.iter()) {
                *d += s;
            }
        }
        if self.qual_pos.len() < o.qual_pos.len() {
            self.qual_pos.resize(o.qual_pos.len(), QUAL_ROW);
        }
        for (i, row) in o.qual_pos.iter().enumerate() {
            for (d, s) in self.qual_pos[i].iter_mut().zip(row.iter()) {
                *d += s;
            }
        }
    }

    pub fn base_at(&self, pos: usize, code: usize) -> u64 {
        self.base_pos.get(pos).map_or(0, |r| r[code])
    }

    /// All 256 per-position quality bins for `pos`.
    pub fn qual_at(&self, pos: usize) -> [u64; QBINS] {
        self.qual_pos.get(pos).copied().unwrap_or(QUAL_ROW)
    }

    /// Positions carrying at least one quality observation.
    pub fn covered_positions(&self) -> usize {
        self.qual_pos
            .iter()
            .filter(|row| row.iter().any(|&c| c > 0))
            .count()
    }
}

fn add_hist(dst: &mut Vec<u64>, src: &[u64]) {
    grow(dst, src.len());
    for (d, s) in dst.iter_mut().zip(src.iter()) {
        *d += s;
    }
}

/// Nearest-rank quantile read off an exact histogram: the smallest value whose
/// cumulative count reaches `ceil(q * n)`. No interpolation, no binning.
pub fn quantile(hist: &[u64], q: f64) -> f64 {
    let n: u64 = hist.iter().sum();
    if n == 0 {
        return 0.0;
    }
    let rank = (q * n as f64).ceil().max(1.0) as u64;
    let mut seen = 0u64;
    for (v, &c) in hist.iter().enumerate() {
        seen += c;
        if seen >= rank {
            return v as f64;
        }
    }
    (hist.len().saturating_sub(1)) as f64
}

#[derive(Debug, Clone, PartialEq)]
pub struct OffsetVerdict {
    /// Phred offset to apply (33 or 64).
    pub offset: u8,
    /// True only when the census makes the other encoding impossible.
    pub certain: bool,
    pub reason: String,
    /// Quality bytes outside the printable Phred range were seen.
    pub malformed: bool,
}

/// Decide the Phred offset from the exact byte census.
///
/// Phred+33 uses bytes 33..=126, Phred+64 uses 64..=126. A byte below 64 is
/// arithmetically impossible under +64, so its presence *proves* +33. Its
/// absence proves nothing: a +33 file whose lowest quality is Q31 also sits
/// entirely above 64. In that case we report the ambiguity instead of
/// asserting a guess, and `--phred-offset` settles it.
pub fn offset_verdict(census: &[u64; QBINS], requested: Option<u8>) -> OffsetVerdict {
    let min = (0..QBINS).find(|&b| census[b] > 0).map(|b| b as u8);
    let malformed = min.map_or(false, |m| m < 33);

    match requested {
        Some(o) => {
            let certain = match (o, min) {
                (33, Some(m)) => m < 64,
                (64, Some(m)) => m >= 64,
                _ => false,
            };
            let reason = if certain {
                format!("--phred-offset {o}, consistent with the byte census")
            } else if o == 64 {
                "forced +64, but the census contains bytes below 64 (impossible under +64)"
                    .to_string()
            } else {
                "forced +33; the census cannot exclude +64 either".to_string()
            };
            OffsetVerdict {
                offset: o,
                certain,
                reason,
                malformed,
            }
        }
        None => match min {
            None => OffsetVerdict {
                offset: 33,
                certain: false,
                reason: "no quality bytes observed".to_string(),
                malformed,
            },
            Some(m) if m < 64 => OffsetVerdict {
                offset: 33,
                certain: true,
                reason: format!("census minimum {m} is below 64, which Phred+64 cannot represent"),
                malformed,
            },
            Some(m) => OffsetVerdict {
                offset: 33,
                certain: false,
                reason: format!(
                    "census minimum {m} >= 64: +33 and +64 are both arithmetically possible; \
                     pass --phred-offset to decide"
                ),
                malformed,
            },
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-computed fixture. Three reads, all 4 bp, so every counter below can
    /// be checked by eye.
    ///
    /// read 1: ACGT / IIII  -> acgt 4, gc 2, mean byte 73
    /// read 2: AANN / !!!!  -> acgt 2, gc 0, has_n, mean byte 33
    /// read 3: GGCC / 55..  -> acgt 4, gc 4, mean byte 53
    fn fixture() -> Accum {
        let mut a = Accum::new(8);
        a.observe(b"ACGT", b"IIII");
        a.observe(b"AANN", b"!!!!");
        a.observe(b"GGCC", b"55..");
        a
    }

    #[test]
    fn totals_and_lengths() {
        let a = fixture();
        assert_eq!(a.reads, 3);
        assert_eq!(a.bases, 12);
        assert_eq!(a.len_hist[4], 3);
        assert_eq!(a.reads_with_n, 1);
        assert_eq!(a.reads_over_max_len, 0);
    }

    #[test]
    fn composition_is_exact_and_n_is_not_acgt() {
        let a = fixture();
        // position 1 (index 0): read1 A, read2 A, read3 G -> A=2, G=1
        assert_eq!(a.base_at(0, 0), 2, "A at pos 1");
        assert_eq!(a.base_at(0, 2), 1, "G at pos 1");
        // position 3 (index 2): read1 G, read2 N, read3 C
        assert_eq!(a.base_at(2, 1), 1, "C at pos 3");
        assert_eq!(a.base_at(2, 2), 1, "G at pos 3");
        assert_eq!(a.base_at(2, 4), 1, "N at pos 3");
        // ACGT denominators per read exclude the Ns
        assert_eq!(a.acgt_hist[4], 2);
        assert_eq!(a.acgt_hist[2], 1);
        assert_eq!(a.gc_hist[4], 1);
        assert_eq!(a.gc_hist[2], 1);
        assert_eq!(a.gc_hist[0], 1);
    }

    #[test]
    fn quality_byte_census_and_offset_verdict() {
        let a = fixture();
        // '!' = 33 appears 4x (read 2), '5' = 53 twice, '.' = 46 twice, 'I' = 73 4x
        assert_eq!(a.qual_byte[33], 4);
        assert_eq!(a.qual_byte[53], 2);
        assert_eq!(a.qual_byte[46], 2);
        assert_eq!(a.qual_byte[73], 4);
        let v = offset_verdict(&a.qual_byte, None);
        assert_eq!(v.offset, 33);
        assert!(v.certain, "byte 33 proves +33");
        assert!(!v.malformed);
    }

    #[test]
    fn offset_stays_unasserted_when_ambiguous() {
        let mut census = [0u64; QBINS];
        census[70] = 5; // 'F' = Q37 under +33, Q6 under +64
        let v = offset_verdict(&census, None);
        assert_eq!(v.offset, 33, "default, but not claimed as certain");
        assert!(!v.certain);
        assert!(v.reason.contains("--phred-offset"));
        let v64 = offset_verdict(&census, Some(64));
        assert_eq!(v64.offset, 64);
        assert!(v64.certain);
        // Forcing +64 on a census that cannot be +64 is reported, not hidden.
        let mut low = [0u64; QBINS];
        low[33] = 1;
        let forced = offset_verdict(&low, Some(64));
        assert!(!forced.certain);
        assert!(forced.reason.contains("impossible"));
    }

    #[test]
    fn malformed_quality_bytes_are_flagged() {
        let mut census = [0u64; QBINS];
        census[10] = 1;
        let v = offset_verdict(&census, None);
        assert!(v.malformed);
    }

    #[test]
    fn quantiles_are_nearest_rank() {
        // lengths 1,2,3,4 -> median at rank ceil(0.5*4)=2 -> 2; q90 -> 4
        let mut h = vec![0u64; 5];
        h[1] = 1;
        h[2] = 1;
        h[3] = 1;
        h[4] = 1;
        assert_eq!(quantile(&h, 0.5), 2.0);
        assert_eq!(quantile(&h, 0.9), 4.0);
        assert_eq!(quantile(&h, 0.0), 1.0);
        assert_eq!(quantile(&[0u64; 4], 0.5), 0.0);
    }

    #[test]
    fn per_position_quality_is_binned_by_byte() {
        let a = fixture();
        let row = a.qual_at(0);
        assert_eq!(row[73], 1); // 'I'
        assert_eq!(row[33], 1); // '!'
        assert_eq!(row[53], 1); // '5'
        assert_eq!(a.covered_positions(), 4);
    }

    #[test]
    fn positions_beyond_max_read_len_are_excluded_not_misattributed() {
        let mut a = Accum::new(4);
        a.observe(b"ACGTAACC", b"IIIIIIII"); // 8 bp read, max_pos 4
        assert_eq!(a.reads, 1);
        assert_eq!(a.bases, 8);
        assert_eq!(a.reads_over_max_len, 1);
        assert_eq!(a.base_at(0, 0), 1, "'A' at position 1");
        assert_eq!(a.base_at(1, 1), 1, "'C' at position 2");
        assert_eq!(a.base_at(3, 3), 1, "'T' at position 4, the last covered one");
        assert_eq!(a.covered_positions(), 4, "nothing beyond max_pos");
        assert_eq!(a.base_at(4, 0), 0, "position 5 must not be counted");
        assert_eq!(a.base_pos.len(), 4, "no rows beyond max_pos are allocated");
    }

    #[test]
    fn merge_is_exact() {
        let mut x = Accum::new(8);
        x.observe(b"ACGT", b"IIII");
        let mut y = Accum::new(8);
        y.observe(b"AANN", b"!!!!");
        y.observe(b"GGCC", b"55..");
        let mut m = x.clone();
        m.merge(&y);
        let f = fixture();
        assert_eq!(m.reads, f.reads);
        assert_eq!(m.bases, f.bases);
        assert_eq!(m.gc_hist, f.gc_hist);
        assert_eq!(m.qual_byte, f.qual_byte);
        for p in 0..8 {
            for code in 0..NT {
                assert_eq!(m.base_at(p, code), f.base_at(p, code), "base {p}/{code}");
            }
            assert_eq!(m.qual_at(p), f.qual_at(p), "qual at {p}");
        }
    }
}
