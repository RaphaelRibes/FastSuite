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
//!   qual_pos  R x 128 x 8 B = 154 KB
//! `qual_pos` is indexed by the raw quality byte (`byte - 33`) rather than by
//! Phred, so the accumulator never needs to know the encoding: the offset is
//! decided once, at report time, from the exact byte census. A run therefore
//! cannot be silently wrong about +33/+64.
//!
//! The plane is 128 cells wide because the only bytes a Phred quality can be are
//! 33..=126 — `byte - 33` is the cell. The census below is still 256 bytes wide:
//! it is what proves +33 vs +64 and what reports the byte range, so a byte
//! outside 33..=126 must be *counted* there even though it is excluded from the
//! per-position rows (FC-8). That exclusion is the only place where a malformed
//! file and a valid one can report different per-position quality.

pub const NT: usize = 5; // A, C, G, T, other/N
/// Width of the exact quality-byte census: one cell per possible byte.
pub const QBYTES: usize = 256;
/// Lowest and highest byte a Phred quality can be (Phred+33, Q0..Q93).
pub const QMIN: u8 = 33;
pub const QMAX: u8 = 126;
/// Cells in a per-position quality row: indexed by `byte - QMIN`, so only 94 of them
/// can ever be touched (33..=126). The row is padded to 128 cells so it stays a power
/// of two — 1 024 B, 16 cache lines — which is half the 2 048 B the 256-cell row cost.
pub const QBINS: usize = 128;

/// The per-position quality cell for a raw byte, or `None` when the byte is not a
/// Phred quality at all (it goes to the census, not to the rows).
#[inline]
pub fn qbin(b: u8) -> Option<usize> {
    if b >= QMIN && b <= QMAX {
        Some((b - QMIN) as usize)
    } else {
        None
    }
}

/// 256-entry base LUT: A=0 C=1 G=2 T=3 (either case), N and every other IUPAC code 4.
/// One table lookup replaces the `to_ascii_uppercase()` + five-arm match the loop used
/// to run per base, per pass.
pub const BASE_CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    let mut i = 0usize;
    while i < 256 {
        t[i] = match i as u8 {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        };
        i += 1;
    }
    t
};

/// The 12-mer rolling register's mask: 12 bases x 2 bits.
const KMER_MASK: u32 = (1 << (2 * repeated::KMER_K)) - 1;

/// The 12-mer rolling 2-bit register, shared by the fused path (`Accum::observe`)
/// and the gated path (`Accum::observe_gated`).
///
/// Both paths push every base code through this one definition of the window rule,
/// so `windows(12)` + `pack_kmer` (the pre-fusion loop the gated twin kept, which
/// decoded each base ~11 times) can never disagree with the fused register about
/// what a counted window is: 12 consecutive A/C/G/T codes, and any other code
/// (N, IUPAC) resets the run.
#[derive(Default)]
pub struct KmerRoll {
    win: u32,
    valid: u32,
}

impl KmerRoll {
    /// Push one base code. Returns the packed k-mer when the last 12 codes are all
    /// A/C/G/T, i.e. exactly when `pack_kmer` would have accepted that window.
    #[inline]
    pub fn push(&mut self, code: u8) -> Option<u32> {
        self.win = ((self.win << 2) | code as u32) & KMER_MASK;
        if code == 4 {
            self.valid = 0;
        } else {
            self.valid = (self.valid + 1).min(repeated::KMER_K as u32);
        }
        if self.valid == repeated::KMER_K as u32 {
            Some(self.win)
        } else {
            None
        }
    }
}

fn grow<T: Clone + Default>(v: &mut Vec<T>, n: usize) {
    if v.len() < n {
        v.resize(n, T::default());
    }
}

// [u64; QBINS] has no `Default`, so rows are grown with an explicit literal.
const QUAL_ROW: [u64; QBINS] = [0; QBINS];


use crate::heatmap;
use crate::polyx;
use crate::repeated;
use crate::repeated::KmerTable;

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
    pub meanq_raw_hist: [u64; QBYTES],
    /// Exact census of quality bytes; this is what decides +33 vs +64.
    pub qual_byte: [u64; QBYTES],
    /// Quality bytes outside `33..=126`: counted here, in the census, and NOT folded
    /// into the per-position rows (which are indexed by `byte - 33`). Zero on a valid
    /// FASTQ, and the report says so only when it is not.
    pub qual_out_of_range: u64,
    /// `base_pos[pos][code]`, position-major.
    pub base_pos: Vec<[u64; NT]>,
    /// `qual_pos[pos][byte - 33]`, position-major; 128 cells, 33..=126.
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
    /// Exact 12-mer counter: a dense 4^12 `u16` table plus an exact width-overflow map
    /// (`crate::repeated::KmerTable`). No cap, no eviction, no approximation.
    pub kmer_counts: KmerTable,
    /// Windows actually counted into `kmer_counts`: one per 12 consecutive A/C/G/T
    /// bases. This is the denominator a k-mer share must divide by — a read of `L`
    /// bases carries `L - 11` windows, not one, so dividing a count by the read count
    /// is unbounded (a 12-mer repeating inside every read reaches `L - 11` per read).
    pub kmer_windows: u64,
    /// Bases of each code over EVERY base read, including those past `--max-read-len`
    /// (which the per-position `base_pos` rows do not reach). This is what makes the
    /// reported composition sum to `bases` rather than to the covered window only.
    pub base_totals: [u64; NT],
    /// 2-bit packed sequence table, filled only in `--overrep-mode table` (the memory is
    /// real: ~46 B per unique 150 bp read, admission-checked in main.rs).
    pub overrep_table: crate::overrep::SequenceTable,
    /// Whether the packed sequence table is filled (set in main.rs for `table` mode).
    pub track_overrep_table: bool,
    /// Adapter content per position, only when --adapter-content is on.
    pub adapter_pos: Option<crate::adapter_pos::AdapterPos>,
    /// The adapter set for --adapter-content (only consulted when adapter_pos is Some).
    pub adapter_set: crate::adapter::AdapterSet,
    /// Per-tile quality, only when --tiles is on. `None` means the flag was off, the work
    /// was not done, and no tile block is emitted.
    pub tile: Option<crate::tile::TileQual>,
}

/// Adapter prefixes the fitter's readthrough scan looks for.
pub const ADAPTER_PREFIXES: [(&str, &str); 2] =
    [("TruSeq", "AGATCGGAAGAGC"), ("Nextera", "CTGTCTCTTATA")];

/// Case-insensitive substring search (the adapter prefix in a read).
///
/// Word-at-a-time prefilter, no dependency: eight candidate start offsets are tested
/// per `u64` load by folding the haystack word and the needle's first byte with
/// `& 0xDF` and looking for a zero byte (`x - 0x0101..` & `!x` & `0x8080..`). The
/// 13-byte pattern is then compared — still folded, word at a time — only at the
/// candidates, and confirmed byte-exactly with `to_ascii_uppercase`.
///
/// `& 0xDF` is a *superset* of case-insensitive equality: `'@'` and ``'`'`` fold
/// together without being the same letter case, so a folded match is a candidate, never
/// an answer. Folded *inequality* is conclusive, which is what makes the reject cheap:
/// the 2 x 138 x 13 `to_ascii_uppercase` calls the window scan performed per read
/// become ~19 `u64` loads plus one exact verification at each real candidate.
pub(crate) fn contains_ci(hay: &[u8], needle: &[u8]) -> bool {
    const FOLD: u64 = 0xDFDF_DFDF_DFDF_DFDF;
    const ONES: u64 = 0x0101_0101_0101_0101;
    const HIGH: u64 = 0x8080_8080_8080_8080;
    let n = needle.len();
    if n == 0 || hay.len() < n {
        return false;
    }
    #[inline]
    fn word(s: &[u8]) -> u64 {
        u64::from_le_bytes([s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]])
    }
    // Reject fast on folded words, then confirm exactly. Only ever called at a
    // candidate, so its cost is paid once per genuine first-byte match.
    #[inline]
    fn verify(h: &[u8], needle: &[u8]) -> bool {
        let n = needle.len();
        let mut i = 0;
        while i + 8 <= n {
            if (word(&h[i..]) & FOLD) != (word(&needle[i..]) & FOLD) {
                return false;
            }
            i += 8;
        }
        while i < n {
            if (h[i] & 0xDF) != (needle[i] & 0xDF) {
                return false;
            }
            i += 1;
        }
        h.iter().zip(needle).all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
    }

    let first = (needle[0] & 0xDF) as u64 * ONES;
    let last = hay.len() - n;
    let mut i = 0usize;
    while i <= last {
        if last - i + 1 >= 8 {
            // Eight candidate starts in one load; the lanes are all <= `last` because
            // eight candidates remain, so a match here is a legal window.
            let x = (word(&hay[i..]) & FOLD) ^ first;
            let z = x.wrapping_sub(ONES) & !x & HIGH;
            if z == 0 {
                i += 8;
                continue;
            }
            let mut bits = z;
            while bits != 0 {
                let lane = (bits.trailing_zeros() / 8) as usize;
                if verify(&hay[i + lane..], needle) {
                    return true;
                }
                bits &= bits - 1;
            }
            i += 8;
        } else {
            if (hay[i] & 0xDF) == (needle[0] & 0xDF) && verify(&hay[i..], needle) {
                return true;
            }
            i += 1;
        }
    }
    false
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
            meanq_raw_hist: [0; QBYTES],
            qual_byte: [0; QBYTES],
            qual_out_of_range: 0,
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
            kmer_counts: KmerTable::new(),
            kmer_windows: 0,
            base_totals: [0; NT],
            overrep_table: crate::overrep::SequenceTable::default(),
            track_overrep_table: false,
            adapter_pos: None,
            adapter_set: crate::adapter::AdapterSet::illumina_default(),
            tile: None,
        }
    }

    /// Accumulator bytes for the rows currently allocated — printed before the
    /// numbers so a run's analysis memory is known, not discovered.
    pub fn resident_bytes(&self) -> usize {
        self.base_pos.len() * std::mem::size_of::<[u64; NT]>()
            + self.qual_pos.len() * std::mem::size_of::<[u64; QBINS]>()
            + (self.len_hist.len() + self.acgt_hist.len() + self.gc_hist.len()) * 8
            + self.n_run_hist.len() * 8
            + 2 * QBYTES * 8
            + self.kmer_counts.resident_bytes()
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

    /// One record, in one pass over the sequence.
    ///
    /// Every per-base consumer reads the same decoded base: the composition cell and
    /// the N-run state machine, the heatmap's `(x = #(T+C), y = #(T+G))`, the rolling
    /// 12-mer register, and the fragment start the `--profile` fit needs. Before this
    /// the sequence was walked four times (composition, heatmap, xxh3, 12-mer windows)
    /// and the 12-mer loop decoded each base 11.1x; now it is decoded once through the
    /// 256-entry [`BASE_CODE`] LUT.
    ///
    /// The N-run state machine and the polyX measurement are unchanged: a run is
    /// `record_n_run`ed where it ends, exactly as before, and `polyx::measure_polyx`
    /// is still its own read (it breaks at the first transition, so it costs 2–3 byte
    /// compares per read, not per base).
    ///
    /// `header` carries the tile coordinates for `--tiles`; `seq` and `qual` must have
    /// equal length (the scanner enforces it: a mismatch is an error, not an
    /// assumption).
    pub fn observe(&mut self, header: &[u8], seq: &[u8], qual: &[u8]) {
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
        // Heatmap counters: x = #(T+C), y = #(T+G), capped at the grid edge at the end
        // (`heatmap::encode`'s own rule, so both paths agree on a saturating read).
        let (mut hx, mut hy) = (0u64, 0u64);
        // Rolling 12-mer register (shared with the gated path). Any other code (N,
        // IUPAC) resets the run, which is exactly `pack_kmer`'s rule.
        let mut roll = KmerRoll::default();
        // First non-G base: the fragment start the profile fitter measures from.
        let mut first_non_g: Option<usize> = None;
        for (i, &raw) in seq.iter().enumerate() {
            let code = BASE_CODE[raw as usize];
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
                match code {
                    1 => {
                        gc += 1;
                        hx += 1;
                    }
                    2 => {
                        gc += 1;
                        hy += 1;
                    }
                    3 => {
                        hx += 1;
                        hy += 1;
                    }
                    _ => {}
                }
            }
            self.base_totals[code as usize] += 1;
            if first_non_g.is_none() && code != 2 {
                first_non_g = Some(i);
            }
            if i < per_pos {
                self.base_pos[i][code as usize] += 1;
            }
            if let Some(k) = roll.push(code) {
                self.kmer_counts.add(k);
                self.kmer_windows += 1;
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
            match qbin(q) {
                Some(b) => {
                    if i < per_pos {
                        self.qual_pos[i][b] += 1;
                    }
                }
                // Not a Phred byte: in the census above, never in a per-position row.
                None => self.qual_out_of_range += 1,
            }
        }
        let raw_mean = if qual.is_empty() {
            0
        } else {
            (qsum / qual.len() as u64) as usize
        };
        self.meanq_raw_hist[raw_mean.min(QBYTES - 1)] += 1;

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
                // The fragment starts at the first non-G base. Every base before it is a
                // G, so it contributes exactly `lead` to both ACGT and GC: the fragment's
                // counts are the read's own counts minus `lead`, with no second pass.
                let lead = first_non_g.unwrap_or(len) as u64;
                debug_assert!(lead <= gc && lead <= acgt);
                self.frag_gc += gc - lead;
                self.frag_bases += acgt - lead;
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

        // Heatmap cell, from the counters the base loop already kept.
        self.heatmap[heatmap::index(
            (hx as usize).min(heatmap::SIZE - 1),
            (hy as usize).min(heatmap::SIZE - 1),
        )] += 1;

        // Duplication: every read is fingerprinted, exactly, with no sample and no cap.
        // The table starts at 1024 entries and grows (main.rs admission-checks the size it
        // reaches at the projected record count). Over-representation naming is a separate,
        // opt-in concern (§5.9).
        self.seq_counts.observe(seq);
        if self.track_overrep_table {
            // --overrep-mode table packs sequences 2 bits per base, and a read with an N or
            // an IUPAC code cannot be packed. The scan refuses such a read outright
            // (reader::scan_file_gated, naming the record and the base's position), so an
            // unpackable sequence never reaches here — it is not silently dropped from the
            // unique-sequence table, which is exactly what the old `let _ =` did.
            if !self.overrep_table.observe(seq) {
                debug_assert!(
                    false,
                    "reader::scan_file_gated must refuse a sequence table mode cannot pack"
                );
            }
        }
        if let Some(ap) = self.adapter_pos.as_mut() {
            ap.observe(seq, &self.adapter_set);
        }
        // Per-tile quality: only when --tiles is on (the accumulator is Some).
        if let Some(t) = self.tile.as_mut() {
            t.observe(header, qual);
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
        // The 12-mer table merges exactly, cell by cell, with the width overflow carried
        // through `hot`: the merged spectrum is the spectrum of the concatenated input.
        // Without this the merged accumulator carries an empty spectrum and every report
        // says kmer_top_count 0 however over-represented the k-mers are.
        self.kmer_counts.merge(&o.kmer_counts);
        self.kmer_windows += o.kmer_windows;
        match (self.adapter_pos.as_mut(), o.adapter_pos.as_ref()) {
            (Some(a), Some(b)) => a.merge(b),
            // The run total starts without the optional counter: the first file that
            // tracked it is ADOPTED whole (a clone of nothing-plus-itself), so its
            // counts are not added twice.
            (None, Some(b)) => self.adapter_pos = Some(b.clone()),
            _ => {}
        }
        // The per-tile accumulator merges the same way: adopted whole by the first file
        // that tracked it, summed exactly thereafter.
        match (self.tile.as_mut(), o.tile.as_ref()) {
            (Some(a), Some(b)) => a.merge(b),
            (None, Some(b)) => self.tile = Some(b.clone()),
            _ => {}
        }
        if self.track_overrep_table || o.track_overrep_table {
            self.track_overrep_table = true;
            self.overrep_table.merge(&o.overrep_table);
        }
        add_hist(&mut self.len_hist, &o.len_hist);
        add_hist(&mut self.acgt_hist, &o.acgt_hist);
        add_hist(&mut self.gc_hist, &o.gc_hist);
        for i in 0..QBYTES {
            self.meanq_raw_hist[i] += o.meanq_raw_hist[i];
            self.qual_byte[i] += o.qual_byte[i];
        }
        self.qual_out_of_range += o.qual_out_of_range;
        for (d, s) in self.base_totals.iter_mut().zip(o.base_totals.iter()) {
            *d += s;
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

    /// The 128 per-position quality cells for `pos`, indexed by `byte - 33`.
    pub fn qual_at(&self, pos: usize) -> [u64; QBINS] {
        self.qual_pos.get(pos).copied().unwrap_or(QUAL_ROW)
    }

    /// Bases the per-position planes hold: the first `--max-read-len` windows of every read.
    ///
    /// This is the denominator for a fraction of the per-position G sum. `bases` counts
    /// EVERY base in the file, so dividing by it made the head-dropout baseline and the
    /// profile's `g_mean` scale with the cap (120 bp reads at `--max-read-len 20`: a G
    /// fraction six times too small).
    pub fn per_pos_bases(&self) -> u64 {
        self.base_pos.iter().map(|r| r.iter().sum::<u64>()).sum()
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
pub fn offset_verdict(census: &[u64; QBYTES], requested: Option<u8>) -> OffsetVerdict {
    let min = (0..QBYTES).find(|&b| census[b] > 0).map(|b| b as u8);
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
    use serde_json::json;

    /// Hand-computed fixture. Three reads, all 4 bp, so every counter below can
    /// be checked by eye.
    ///
    /// read 1: ACGT / IIII  -> acgt 4, gc 2, mean byte 73
    /// read 2: AANN / !!!!  -> acgt 2, gc 0, has_n, mean byte 33
    /// read 3: GGCC / 55..  -> acgt 4, gc 4, mean byte 53
    fn fixture() -> Accum {
        let mut a = Accum::new(8);
        a.observe(b"@r", b"ACGT", b"IIII");
        a.observe(b"@r", b"AANN", b"!!!!");
        a.observe(b"@r", b"GGCC", b"55..");
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
        let mut census = [0u64; QBYTES];
        census[70] = 5; // 'F' = Q37 under +33, Q6 under +64
        let v = offset_verdict(&census, None);
        assert_eq!(v.offset, 33, "default, but not claimed as certain");
        assert!(!v.certain);
        assert!(v.reason.contains("--phred-offset"));
        let v64 = offset_verdict(&census, Some(64));
        assert_eq!(v64.offset, 64);
        assert!(v64.certain);
        // Forcing +64 on a census that cannot be +64 is reported, not hidden.
        let mut low = [0u64; QBYTES];
        low[33] = 1;
        let forced = offset_verdict(&low, Some(64));
        assert!(!forced.certain);
        assert!(forced.reason.contains("impossible"));
    }

    #[test]
    fn malformed_quality_bytes_are_flagged() {
        let mut census = [0u64; QBYTES];
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
        assert_eq!(row[73 - 33], 1); // 'I' = 73 -> cell 40
        assert_eq!(row[0], 1); // '!' = 33 -> cell 0
        assert_eq!(row[53 - 33], 1); // '5' = 53 -> cell 20
        assert_eq!(a.covered_positions(), 4);
    }

    /// FC-8: the per-position plane is 128 cells, indexed by `byte - 33`, because the
    /// only bytes a valid FASTQ can carry are 33..=126. A byte outside that range is
    /// still in the exact byte census (which decides the offset and reports the range),
    /// but it is not folded into a per-position row where it would silently shift the
    /// position's mean and quantiles.
    #[test]
    fn per_position_quality_plane_covers_only_valid_phred_bytes() {
        let mut a = Accum::new(8);
        // One read per interesting byte, so every position carries the same quality and
        // the row can be read cell by cell.
        a.observe(b"@r", b"ACGT", &[33u8, 33, 33, 33]);
        a.observe(b"@r", b"ACGT", &[73u8, 73, 73, 73]);
        a.observe(b"@r", b"ACGT", &[126u8, 126, 126, 126]);
        a.observe(b"@r", b"ACGT", &[200u8, 200, 200, 200]);
        let row = a.qual_at(0);
        assert_eq!(row.len(), 128, "the plane is 128 cells, not 256");
        assert_eq!(row[0], 1, "'!' = 33 -> cell 0");
        assert_eq!(row[73 - 33], 1, "'I' = 73 -> cell 40");
        assert_eq!(row[126 - 33], 1, "the highest valid byte, 126, is the last used cell");
        assert_eq!(row.iter().sum::<u64>(), 3, "byte 200 is not a per-position quality");
        assert_eq!(a.qual_byte[200], 4, "but it is still in the exact byte census");
        // The census, and therefore the offset verdict and the reported byte range, are
        // untouched: a malformed file is still named as one.
        let v = offset_verdict(&a.qual_byte, None);
        assert!(v.certain, "byte 33 proves +33, exactly as before");
        assert!(!v.malformed, "the malformed flag is about bytes BELOW 33");
    }

    /// FC-8, the other half: a byte above 126 is out of range for a Phred byte too, and
    /// it is counted where the report can name it rather than dropped.
    #[test]
    fn quality_bytes_above_the_phred_range_are_counted_not_dropped() {
        let mut a = Accum::new(8);
        a.observe(b"@r", b"ACGT", &[33u8, 73, 126, 200]);
        assert_eq!(a.qual_out_of_range, 1, "byte 200 is outside 33..=126");
        let mut ok = Accum::new(8);
        ok.observe(b"@r", b"ACGT", &[33u8, 73, 126, 40]);
        assert_eq!(ok.qual_out_of_range, 0, "a valid file pays nothing and reports nothing");
        // Both merge exactly, like every other counter here.
        let mut m = Accum::new(8);
        m.merge(&a);
        m.merge(&ok);
        assert_eq!(m.qual_out_of_range, 1);
    }

    #[test]
    fn positions_beyond_max_read_len_are_excluded_not_misattributed() {
        let mut a = Accum::new(4);
        a.observe(b"@r", b"ACGTAACC", b"IIIIIIII"); // 8 bp read, max_pos 4
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
        x.observe(b"@r", b"ACGT", b"IIII");
        let mut y = Accum::new(8);
        y.observe(b"@r", b"AANN", b"!!!!");
        y.observe(b"@r", b"GGCC", b"55..");
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

    /// C1: the composition totals cover EVERY base, including the bases past
    /// `--max-read-len` that no per-position row can hold. Ten reads of 8 bp at
    /// max_pos 4 used to report 40 bases of composition against `bases` 80.
    #[test]
    fn base_totals_cover_every_base_including_past_max_read_len() {
        let mut a = Accum::new(4);
        for _ in 0..10 {
            a.observe(b"@r", b"ACGTACGT", b"IIIIIIII");
        }
        assert_eq!(a.bases, 80);
        assert_eq!(a.base_totals.iter().sum::<u64>(), 80, "every base is accounted for");
        assert_eq!(a.base_totals[0], 20, "4 As per read x 10 reads");
        assert_eq!(a.base_totals[4], 0, "no ambiguous base in this fixture");
        // ... and it merges exactly, so a multi-file report sums the same way.
        let mut b = Accum::new(4);
        for _ in 0..10 {
            b.observe(b"@r", b"NNNNACGT", b"IIIIIIII");
        }
        a.merge(&b);
        assert_eq!(a.base_totals.iter().sum::<u64>(), 160);
        assert_eq!(a.base_totals[4], 40);
    }

    /// Regression: `Accum::merge` summed every counter except the 12-mer one, so the merged
    /// accumulator — the only one any report renders — always had an empty spectrum:
    /// `kmer_top_count` 0 and `kmer_spectrum.top` [] even when every read carried the same
    /// 12-mer. 2000 reads of 40 bp, split over two accumulators and merged.
    #[test]
    fn kmer_spectrum_survives_merge() {
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bp, 29 12-mer windows
        let qual = [b'I'; 40];
        let mut a = Accum::new(40);
        let mut b = Accum::new(40);
        for _ in 0..1000 {
            a.observe(b"@r", read, &qual);
            b.observe(b"@r", read, &qual);
        }
        let mut total = Accum::new(40);
        total.merge(&a);
        total.merge(&b);
        assert_eq!(total.reads, 2000);

        // Through the report path, not the raw counter: this is what the JSON emits.
        // 29 windows per read over a period-4 sequence give four distinct 12-mers, the
        // phase-0 one (ACGTACGTACGT) at 8 windows per read: 8 * 2000 = 16000.
        let (m, block) = crate::rich::kmer_metrics(&total);
        assert_eq!(m["kmer_top_count"], json!(16_000));
        // C2: the share is of the 29 windows per read (58 000 windows), not of the
        // 2 000 reads: 16000/58000 = 0.2759, where the read denominator said 8.0.
        assert_eq!(m["kmer_top_frac"], json!(0.2759));
        assert_eq!(block["top"].as_array().unwrap().len(), 4);
        assert_eq!(block["top"][0]["kmer"], json!("ACGTACGTACGT"));
        assert_eq!(block["top"][0]["count"], json!(16_000));
    }

    /// FC-1: the spectrum is exact now, and it says so. A 16 bp read carries five
    /// 12-mer windows; every one of them is counted, `distinct` comes off the non-zero
    /// cells, and the block no longer claims an approximation.
    #[test]
    fn kmer_spectrum_is_exact_and_reports_distinct() {
        let mut a = Accum::new(32);
        for _ in 0..20 {
            // 16 bp, five windows, all five different. 20 reads puts the top count (20)
            // above the report's 0.1%-of-reads / 10-count threshold.
            a.observe(b"@r", b"ACGTACGTACGTACGA", &[b'I'; 16]);
        }
        let (m, block) = crate::rich::kmer_metrics(&a);
        assert_eq!(block["approximate"], json!(false), "the dense table is exact");
        assert!(
            block.get("approximation").is_none(),
            "there is no bound left to state: {block}"
        );
        assert_eq!(block["distinct"], json!(5), "5 windows, all distinct");
        assert_eq!(m["kmer_top_count"], json!(20));
        assert_eq!(block["top"].as_array().unwrap().len(), 5);
        let named: Vec<&str> = block["top"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| e["kmer"].as_str().unwrap())
            .collect();
        assert!(named.contains(&"ACGTACGTACGT"), "{named:?}");
    }

    // ---------------------------------------------------------------------------
    // FC-2: the fused single pass, held to the per-function path it replaced.
    // ---------------------------------------------------------------------------

    /// Deterministic xorshift64 — no dev-dependency, and the same reads every run.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    /// The pre-FC-2 path, kept verbatim as the reference the fused loop is held to:
    /// an inline composition/N-run pass, then `heatmap::encode`, then `pack_kmer` per
    /// window, then `contains_ci` per adapter prefix plus the fragment-GC loop, then
    /// `polyx::measure_polyx`. Everything here is a separate function call over the
    /// whole read, which is exactly what `observe` no longer does.
    fn reference(reads: &[(Vec<u8>, Vec<u8>)], max_pos: usize) -> (Accum, Vec<(u64, u32)>) {
        let mut a = Accum::new(max_pos);
        a.track_profile = true;
        let mut kmers: std::collections::BTreeMap<u64, u32> = std::collections::BTreeMap::new();
        for (seq, qual) in reads {
            let len = seq.len();
            a.reads += 1;
            a.bases += len as u64;
            grow(&mut a.len_hist, len + 1);
            a.len_hist[len] += 1;
            let per_pos = len.min(max_pos);
            if len > max_pos {
                a.reads_over_max_len += 1;
            }
            grow(&mut a.base_pos, per_pos);
            if a.qual_pos.len() < per_pos {
                a.qual_pos.resize(per_pos, QUAL_ROW);
            }

            let (mut acgt, mut gc, mut has_n) = (0u64, 0u64, false);
            let mut run_start = 0usize;
            let mut in_run = false;
            for (i, &raw) in seq.iter().enumerate() {
                let code = match raw.to_ascii_uppercase() {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 4,
                };
                if code == 4 {
                    if !in_run {
                        in_run = true;
                        run_start = i;
                    }
                    has_n = true;
                } else {
                    if in_run {
                        a.record_n_run(run_start, i, len);
                        in_run = false;
                    }
                    acgt += 1;
                    if code == 1 || code == 2 {
                        gc += 1;
                    }
                }
                if i < per_pos {
                    a.base_pos[i][code] += 1;
                }
            }
            if in_run {
                a.record_n_run(run_start, len, len);
            }
            if has_n {
                a.reads_with_n += 1;
            }
            grow(&mut a.acgt_hist, acgt as usize + 1);
            a.acgt_hist[acgt as usize] += 1;
            grow(&mut a.gc_hist, gc as usize + 1);
            a.gc_hist[gc as usize] += 1;

            let mut qsum = 0u64;
            for (i, &q) in qual.iter().enumerate() {
                a.qual_byte[q as usize] += 1;
                qsum += q as u64;
                if (33..=126).contains(&q) {
                    if i < per_pos {
                        a.qual_pos[i][(q - 33) as usize] += 1;
                    }
                } else {
                    a.qual_out_of_range += 1;
                }
            }
            let raw_mean = if qual.is_empty() { 0 } else { (qsum / qual.len() as u64) as usize };
            a.meanq_raw_hist[raw_mean.min(QBYTES - 1)] += 1;

            let mut hit = false;
            for (i, (_, pat)) in ADAPTER_PREFIXES.iter().enumerate() {
                if contains_ci(seq, pat.as_bytes()) {
                    a.adapter_hits[i] += 1;
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
                            a.frag_gc += 1;
                            a.frag_bases += 1;
                        }
                        b'A' | b'T' => a.frag_bases += 1,
                        _ => {}
                    }
                }
            }

            let (lb, ll, tb, tl) = polyx::measure_polyx(seq);
            if lb != usize::MAX {
                grow(&mut a.polyx5[lb], ll + 1);
                a.polyx5[lb][ll] += 1;
            }
            if tb != usize::MAX {
                grow(&mut a.polyx3[tb], tl + 1);
                a.polyx3[tb][tl] += 1;
            }
            const G: usize = 2;
            if lb == G {
                let n = ll.min(qual.len());
                a.run_q_head.0 += qual[..n].iter().map(|&b| b as u64).sum::<u64>();
                a.run_q_head.1 += n as u64;
            }
            if tb == G {
                let start = qual.len().saturating_sub(tl);
                a.run_q_tail.0 += qual[start..].iter().map(|&b| b as u64).sum::<u64>();
                a.run_q_tail.1 += (qual.len() - start) as u64;
            }

            let (hx, hy) = heatmap::encode(seq);
            a.heatmap[heatmap::index(hx, hy)] += 1;

            if seq.len() >= 12 {
                for w in seq.windows(12) {
                    if let Some(pk) = repeated::pack_kmer(w) {
                        *kmers.entry(pk).or_insert(0) += 1;
                    }
                }
            }
        }
        let mut kv: Vec<(u64, u32)> = kmers.into_iter().collect();
        kv.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        (a, kv)
    }

    /// Reads with everything the loop has to get right: both cases, N and other IUPAC
    /// codes, empty reads, reads shorter and longer than 12 bp, and a max_pos that
    /// truncates the per-position planes.
    fn random_reads(seed: u64, n: usize, max_len: usize) -> Vec<(Vec<u8>, Vec<u8>)> {
        const ALPHABET: &[u8] = b"ACGTACGTacgtNNnRYry.-";
        let mut rng = Rng(seed | 1);
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            let len = rng.below(max_len + 1);
            let mut seq = Vec::with_capacity(len);
            for _ in 0..len {
                seq.push(ALPHABET[rng.below(ALPHABET.len())]);
            }
            // Every 7th read carries a real adapter prefix, in a random case, so both
            // the hit path and the fragment path are exercised.
            if rng.below(7) == 0 && seq.len() >= 13 {
                let pat: &[u8] = ADAPTER_PREFIXES[rng.below(2)].1.as_bytes();
                let at = rng.below(seq.len() - 12);
                for (j, &b) in pat.iter().enumerate() {
                    seq[at + j] = if j % 2 == 0 { b } else { b.to_ascii_lowercase() };
                }
            }
            let qual: Vec<u8> = (0..len).map(|_| (33 + rng.below(94)) as u8).collect();
            out.push((seq, qual));
        }
        out
    }

    /// The fused loop must be the per-function path, counter for counter. If a base is
    /// decoded once instead of four times, the four answers must not move.
    #[test]
    fn the_fused_pass_is_the_per_function_path() {
        let reads = random_reads(0x5EED_1234, 400, 40);
        let mut fused = Accum::new(24); // shorter than the longest read: truncation too
        fused.track_profile = true;
        for (seq, qual) in &reads {
            fused.observe(b"@r", seq, qual);
        }
        let (want, want_kmers) = reference(&reads, 24);

        assert_eq!(fused.reads, want.reads);
        assert_eq!(fused.bases, want.bases);
        assert_eq!(fused.reads_with_n, want.reads_with_n);
        assert_eq!(fused.reads_over_max_len, want.reads_over_max_len);
        assert_eq!(fused.n_run_hist, want.n_run_hist, "N-run lengths");
        assert_eq!(fused.n_run_leading, want.n_run_leading);
        assert_eq!(fused.n_run_trailing, want.n_run_trailing);
        assert_eq!(fused.len_hist, want.len_hist);
        assert_eq!(fused.acgt_hist, want.acgt_hist);
        assert_eq!(fused.gc_hist, want.gc_hist);
        assert_eq!(fused.meanq_raw_hist, want.meanq_raw_hist);
        assert_eq!(fused.qual_byte, want.qual_byte);
        assert_eq!(fused.qual_out_of_range, want.qual_out_of_range);
        assert_eq!(fused.base_pos, want.base_pos, "composition plane");
        assert_eq!(fused.qual_pos, want.qual_pos, "quality plane");
        assert_eq!(fused.polyx3, want.polyx3);
        assert_eq!(fused.polyx5, want.polyx5);
        assert_eq!(fused.run_q_head, want.run_q_head);
        assert_eq!(fused.run_q_tail, want.run_q_tail);
        assert_eq!(fused.heatmap, want.heatmap, "heatmap grid");
        assert_eq!(fused.adapter_hits, want.adapter_hits, "adapter readthrough rate");
        assert_eq!(fused.frag_gc, want.frag_gc, "fragment GC");
        assert_eq!(fused.frag_bases, want.frag_bases);
        assert_eq!(fused.kmer_counts.report(1), want_kmers, "12-mer spectrum");
    }

    /// FC-6: the word-at-a-time prefilter must answer exactly what the double-fold
    /// window scan answered, including the bytes that survive `& 0xDF` together
    /// ('@' and '`' fold to the same value but are not the same letter case).
    #[test]
    fn contains_ci_agrees_with_the_window_reference() {
        fn reference(hay: &[u8], needle: &[u8]) -> bool {
            if needle.is_empty() || hay.len() < needle.len() {
                return false;
            }
            hay.windows(needle.len())
                .any(|w| w.iter().zip(needle).all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase()))
        }
        let mut cases: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        // the two real prefixes, present and absent, in every case mix
        for pat in ["AGATCGGAAGAGC", "CTGTCTCTTATA"] {
            cases.push((pat.as_bytes().to_vec(), pat.as_bytes().to_vec()));
            cases.push((pat.to_ascii_lowercase().into_bytes(), pat.as_bytes().to_vec()));
            cases.push((b"ACGTACGTACGTACGT".to_vec(), pat.as_bytes().to_vec()));
            let mut h = vec![b'g'; 40];
            h.extend_from_slice(pat.as_bytes());
            cases.push((h, pat.as_bytes().to_vec()));
        }
        // fold collisions: '@'/'`' and '['/'{' differ by 0x20 but are not one letter
        cases.push((b"@@@AGATCGGAAGAGC".to_vec(), b"AGATCGGAAGAGC".to_vec()));
        cases.push((b"@@@`GATCGGAAGAGC".to_vec(), b"AGATCGGAAGAGC".to_vec()));
        cases.push((b"[{][{][{][{][{][{][{".to_vec(), b"[{][{][{][{][{".to_vec()));
        // needles shorter than a word and exactly a word
        cases.push((b"aaAAaa".to_vec(), b"AAAA".to_vec()));
        cases.push((b"xACGTACGTx".to_vec(), b"acgtacgt".to_vec()));
        // random haystacks over a nasty alphabet, and a random needle from each
        const HAY_ALPHABET: &[u8] = b"ACGTacgt@`[{N. nAGATCG";
        const NEEDLE_ALPHABET: &[u8] = b"ACGTacgt@`[";
        let mut rng = Rng(0xC0FF_EE11);
        for _ in 0..300 {
            let hl = rng.below(60);
            let hay: Vec<u8> = (0..hl).map(|_| HAY_ALPHABET[rng.below(HAY_ALPHABET.len())]).collect();
            let nl = 1 + rng.below(13);
            let needle: Vec<u8> =
                (0..nl).map(|_| NEEDLE_ALPHABET[rng.below(NEEDLE_ALPHABET.len())]).collect();
            cases.push((hay, needle));
        }
        for (hay, needle) in &cases {
            assert_eq!(
                contains_ci(hay, needle),
                reference(hay, needle),
                "hay {hay:?} needle {needle:?}"
            );
        }
        // and the semantics the fitter relies on: an empty needle never matches
        assert!(!contains_ci(b"ACGT", b""));
        assert!(contains_ci(b"xxagatcggaagagcxx", b"AGATCGGAAGAGC"));
    }
}
