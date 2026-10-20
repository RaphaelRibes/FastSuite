use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::sync::{mpsc, Arc, Mutex};

use crate::cli::DataMode;
use crate::encode::{encode_read, get_reader, quality_sums};
use crate::file_kind::FileKind;

// ── Types shared with the spiller ──────────────────────────────────────────

/// A single read kept for the `.preview` file.
#[derive(Clone)]
pub struct PreviewRead {
    pub x: u16,
    pub y: u16,
    pub file_id: u16,
    pub read_num: u32,
    pub header: String,
    pub seq: String,
    pub qual: String,
}

#[derive(Clone)]
pub struct SpillChunk {
    pub indices: Vec<(u16, u16, u16, u32)>, // (x, y, file_id, read_num)
    pub previews: Vec<PreviewRead>,
}

/// Per-pixel accumulation: read count plus the *integer* quality sums.
///
/// Storing `(sum_q, sum_len)` instead of an accumulated `f64` mean is what
/// makes the merge order-independent (FM-12): worker threads finish in whatever
/// order the OS schedules them, and `f64` addition is not associative, so
/// summing floating-point partial means in `Done` arrival order made
/// `heatmap_3d.html` differ run to run. `u64` addition is exactly associative,
/// and `sum_q as f64 / sum_len as f64` is a single rounded division — the same
/// value in every order and every process.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PixelAcc {
    pub count: u64,
    pub sum_q: u64,
    pub sum_len: u64,
}

impl PixelAcc {
    /// Mean Phred quality over the reads that carried a quality string.
    pub fn mean_quality(&self) -> f64 {
        if self.sum_len == 0 {
            0.0
        } else {
            self.sum_q as f64 / self.sum_len as f64
        }
    }
}

#[derive(Clone)]
pub enum SpillerMsg {
    Chunk(SpillChunk),
    Done {
        counts: HashMap<(u16, u16), PixelAcc>,
        counts_r1: HashMap<(u16, u16), u64>,
        counts_r2: HashMap<(u16, u16), u64>,
        per_file_counts: HashMap<u16, u64>,
        no_call: u64,
        no_call_r1: u64,
        no_call_r2: u64,
        min_seq_len: Option<u32>,
        max_seq_len: u32,
        n_reads: u64,
    },
}

// ── Worker configuration ────────────────────────────────────────────────────

/// Lightweight config copied into every worker thread.
#[derive(Clone, Copy)]
pub struct WorkerCfg {
    pub preview_n: usize,
    /// Byte threshold after which the worker flushes to the spiller.
    pub threshold_bytes: usize,
    /// Byte threshold for the *preview* buffer (see `preview_flush_bytes`).
    pub preview_flush_bytes: usize,
    /// Whether to keep the per-read index stream at all (FM-1).
    pub data_mode: DataMode,
}

/// Upper bound on the preview flush threshold.
///
/// Previews carry read payloads (header+seq+qual), so they are flushed on their
/// own, much sooner than the index buffer: the spiller streams each one to its
/// row's temp file, which keeps preview memory flat in library size instead of
/// `preview_n x pixels` payloads resident for the whole run.
const PREVIEW_FLUSH_MAX: usize = 256 * 1024;

/// Floor for the preview flush threshold: below this a flush would be a message
/// per few reads, which costs more than the payload it saves.
const PREVIEW_FLUSH_MIN: usize = 32 * 1024;

/// Preview flush threshold for a worker whose index share is `threshold_bytes`.
///
/// FM-2: this used to be a hard-coded 4 MiB *per worker* that `--buffer-mb`
/// could not see, so the preview buffer was a memory floor the budget could not
/// move (`--buffer-mb 1` still held up to 4 MiB per worker). Deriving it from
/// the worker's share puts it inside the budget; the 256 KiB cap keeps it well
/// under the old value while still batching a few thousand reads per message.
/// Output-neutral because previews are capped per `(pixel, kind)` and re-sorted
/// at encode time, so flush boundaries cannot change a byte of `heatmap.preview`.
pub fn preview_flush_bytes(threshold_bytes: usize) -> usize {
    threshold_bytes.clamp(PREVIEW_FLUSH_MIN, PREVIEW_FLUSH_MAX)
}

/// A corner extraction region. Each instance writes matching reads to its
/// own output file.
#[derive(Clone)]
pub struct CornerExtractor {
    pub x1: u32,
    pub y1: u32,
    pub x2: u32,
    pub y2: u32,
    pub writer: Arc<Mutex<BufWriter<File>>>,
}

// ── Worker-local state ──────────────────────────────────────────────────────

/// State accumulated by one worker thread.
///
/// `indices` and `previews` are cleared on every spill; everything else
/// persists to the final `Done` message.
pub struct WorkerLocal {
    pub indices: Vec<(u16, u16, u16, u32)>,
    pub previews: Vec<PreviewRead>,
    /// Per-`(pixel, kind, file)` preview cap — persists across spills.
    /// Keyed by file as well as pixel and kind so the quota is the *first*
    /// `preview_n` of every file, not of whichever file a worker happened to
    /// meet first: the encode-time rule picks the first `preview_n` by
    /// `(file_id, read_num)`, and a quota that merged files gave that rule a
    /// different candidate set at `-j 1` than at `-j 3` (review 2026-10-21,
    /// item 7). RAM stays bounded by the preview flush threshold.
    pub preview_count: HashMap<((u16, u16), FileKind, u16), u8>,
    /// Per-`((x, y), file_id)` count of index records kept, for
    /// `DataMode::Capped` (FM-1). Persists across spills so the cap is a cap on
    /// the whole run, not per chunk. Empty in `Exact`/`Off` mode.
    pub data_count: HashMap<((u16, u16), u16), u8>,
    /// Per-pixel (count, quality sums) — for PNG/HTML.
    pub counts: HashMap<(u16, u16), PixelAcc>,
    /// Per-pixel R1-only count.
    pub counts_r1: HashMap<(u16, u16), u64>,
    /// Per-pixel R2-only count.
    pub counts_r2: HashMap<(u16, u16), u64>,
    /// Reads with zero ACGT bases (all-N / no-call). They land on (0,0) like an
    /// all-A read but have no measurable composition, so the corner metrics
    /// keep them out of their denominators and report them separately.
    pub no_call: u64,
    pub no_call_r1: u64,
    pub no_call_r2: u64,
    /// Per-file total read count.
    pub per_file_counts: HashMap<u16, u64>,
    /// Shortest read length seen (`None` until a read arrives).
    ///
    /// Published next to `max_seq_len` because `size` is derived from the
    /// *longest* read: a read shorter than that cannot reach a corner pixel.
    pub min_seq_len: Option<u32>,
    pub max_seq_len: u32,
    pub total_reads: u64,
    pub bytes_estimate: usize,
    /// Bytes currently held in `previews` (flushed on its own threshold).
    pub preview_bytes: usize,
}

impl WorkerLocal {
    pub fn new() -> Self {
        Self {
            indices: Vec::new(),
            previews: Vec::new(),
            preview_count: HashMap::new(),
            data_count: HashMap::new(),
            counts: HashMap::new(),
            counts_r1: HashMap::new(),
            counts_r2: HashMap::new(),
            no_call: 0,
            no_call_r1: 0,
            no_call_r2: 0,
            per_file_counts: HashMap::new(),
            min_seq_len: None,
            max_seq_len: 0,
            total_reads: 0,
            bytes_estimate: 0,
            preview_bytes: 0,
        }
    }

    /// Record one read in the per-read index stream, honouring `--data-mode`.
    ///
    /// * `Exact` — always, as before.
    /// * `Capped(n)` — only while this `(pixel, file)` has fewer than `n`.
    /// * `Off` — never: no push, so nothing can spill and the render pass has
    ///   nothing to read.
    ///
    /// Returns whether the record was kept. Pixel *counts* are updated by the
    /// caller regardless, so this only ever affects `heatmap.data`.
    #[inline]
    pub fn push_index(
        &mut self,
        xy: (u16, u16),
        file_id: u16,
        read_num: u32,
        mode: DataMode,
    ) -> bool {
        let keep = match mode {
            DataMode::Exact => true,
            DataMode::Off => false,
            DataMode::Capped(n) => {
                let c = self.data_count.entry((xy, file_id)).or_insert(0);
                if *c < n {
                    *c += 1;
                    true
                } else {
                    false
                }
            }
        };
        if keep {
            self.indices.push((xy.0, xy.1, file_id, read_num));
            // 16 B charged per record against a 12 B in-memory record; the
            // budget is intentionally conservative here (FM-5, not taken).
            self.bytes_estimate += 16;
        }
        keep
    }
}

// ── Helpers ─────────────────────────────────────────────────────────────────

/// Move accumulated indices/previews into a `Chunk` and send to the spiller.
pub fn flush_local(state: &mut WorkerLocal, tx: &mpsc::Sender<SpillerMsg>) {
    if state.indices.is_empty() && state.previews.is_empty() {
        return;
    }
    let chunk = SpillChunk {
        indices: std::mem::take(&mut state.indices),
        previews: std::mem::take(&mut state.previews),
    };
    state.bytes_estimate = 0;
    state.preview_bytes = 0;
    tx.send(SpillerMsg::Chunk(chunk))
        .expect("spiller died while workers still running");
}

/// Flush pending previews alone (indices stay buffered). Keeps the read
/// payloads from accumulating: the spiller writes them straight to the row's
/// preview temp file, so they never become `preview_n x pixels` resident.
pub fn flush_previews(state: &mut WorkerLocal, tx: &mpsc::Sender<SpillerMsg>) {
    if state.previews.is_empty() {
        return;
    }
    let chunk = SpillChunk {
        indices: Vec::new(),
        previews: std::mem::take(&mut state.previews),
    };
    state.preview_bytes = 0;
    tx.send(SpillerMsg::Chunk(chunk))
        .expect("spiller died while workers still running");
}

/// Read one FASTQ line into `buf` (cleared first), returning the byte count
/// `read_until` reported (0 = EOF).
///
/// FM-10: `BufRead::lines()` allocated a `String` per line and validated UTF-8
/// on every base — four allocations and two full UTF-8 scans per read, for data
/// that is only ever counted byte by byte. This reuses one `Vec<u8>` per line
/// slot instead.
///
/// It must reproduce `lines()`'s line ending handling exactly: strip `\n`, then
/// strip a `\r` left by CRLF input. A `\r` left in `seq` would add a byte to
/// the read length and move every downstream number (`size`, the (x, y) plane,
/// `data_row_offsets`).
fn read_line(reader: &mut dyn BufRead, buf: &mut Vec<u8>) -> std::io::Result<usize> {
    buf.clear();
    let n = reader.read_until(b'\n', buf)?;
    if buf.last() == Some(&b'\n') {
        buf.pop();
        if buf.last() == Some(&b'\r') {
            buf.pop();
        }
    }
    Ok(n)
}

// ── FASTQ processing ─────────────────────────────────────────────────────────

/// Process one FASTQ file, updating `state` and writing corner extractions.
pub fn process_fastq(
    path: &str,
    file_id: u16,
    state: &mut WorkerLocal,
    tx: &mpsc::Sender<SpillerMsg>,
    cfg: &WorkerCfg,
    corners: &[CornerExtractor],
    file_kinds: &[FileKind],
) -> std::io::Result<()> {
    let mut reader = get_reader(path)?;
    // Four reused line buffers for the whole file (FM-10).
    let mut header: Vec<u8> = Vec::with_capacity(256);
    let mut seq: Vec<u8> = Vec::with_capacity(256);
    let mut plus: Vec<u8> = Vec::with_capacity(16);
    let mut qual: Vec<u8> = Vec::with_capacity(256);
    let mut read_num: u32 = 0;
    let kind = file_kinds
        .get(file_id as usize)
        .copied()
        .unwrap_or(FileKind::Unknown);

    loop {
        // A truncated record ends the file exactly as `lines()` running dry did.
        // A genuine I/O error, in contrast, used to end the loop silently
        // (`lines()` yielded `Some(Err(_))`, the pattern did not match) and
        // leave a plausible-looking partial heatmap; `?` now aborts loudly.
        if read_line(&mut *reader, &mut header)? == 0 {
            break;
        }
        if read_line(&mut *reader, &mut seq)? == 0 {
            break;
        }
        if read_line(&mut *reader, &mut plus)? == 0 {
            break;
        }
        if read_line(&mut *reader, &mut qual)? == 0 {
            break;
        }

        // A FASTQ record's qual string is exactly as long as its seq string.
        // The preview wire format stores one length for both, so a record that
        // breaks the invariant used to abort in the *decoder* (slice panic,
        // exit 101) — or, with no previews kept, not at all. Validate here,
        // where the file and the read number are still known (review
        // 2026-10-21, item 3).
        if seq.len() != qual.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "read {}: {} seq bases but {} quality bytes — a FASTQ record's \
                     seq and qual must be the same length",
                    read_num as u64 + 1,
                    seq.len(),
                    qual.len()
                ),
            ));
        }

        let (cx, cy, n_acgt) = encode_read(&seq);
        let (sum_q, sum_len) = quality_sums(&qual);

        let slen = seq.len() as u32;
        if slen > u16::MAX as u32 {
            eprintln!(
                "Error: read length {} in {} exceeds u16::MAX, aborting",
                slen, path
            );
            std::process::exit(1);
        }
        if slen > state.max_seq_len {
            state.max_seq_len = slen;
        }
        state.min_seq_len = Some(state.min_seq_len.map_or(slen, |m| m.min(slen)));

        let xy = (cx as u16, cy as u16);

        // Per-pixel counts (combined — for PNG / HTML).
        let e = state.counts.entry(xy).or_default();
        e.count += 1;
        e.sum_q += sum_q;
        e.sum_len += sum_len;
        // Per-pixel split counts.
        match kind {
            FileKind::R1 => {
                *state.counts_r1.entry(xy).or_insert(0) += 1;
            }
            FileKind::R2 => {
                *state.counts_r2.entry(xy).or_insert(0) += 1;
            }
            FileKind::Unknown => {}
        }
        // A read with no ACGT base has no measurable composition: it keeps its
        // (0,0) pixel (the plane is unchanged) but is counted here so the
        // corner metrics can leave it out of the polyA denominators.
        if n_acgt == 0 {
            state.no_call += 1;
            match kind {
                FileKind::R1 => state.no_call_r1 += 1,
                FileKind::R2 => state.no_call_r2 += 1,
                FileKind::Unknown => {}
            }
        }
        state.total_reads += 1;
        *state.per_file_counts.entry(file_id).or_insert(0) += 1;

        // Index entry for heatmap.data (unless --data-mode says otherwise).
        state.push_index(xy, file_id, read_num, cfg.data_mode);

        // Preview accumulation (cap per (pixel, kind, file), so R1/R2 each get
        // their own preview_n quota and so does each file — the encode-time
        // rule then sees every file's prefix, whatever -j is).
        if cfg.preview_n > 0 {
            let pc = state.preview_count.entry((xy, kind, file_id)).or_insert(0);
            if (*pc as usize) < cfg.preview_n {
                *pc += 1;
                let cost = header.len() + seq.len() + qual.len() + 64;
                // Only a kept preview pays for its Strings (FM-10).
                state.previews.push(PreviewRead {
                    x: xy.0,
                    y: xy.1,
                    file_id,
                    read_num,
                    header: String::from_utf8_lossy(&header).into_owned(),
                    seq: String::from_utf8_lossy(&seq).into_owned(),
                    qual: String::from_utf8_lossy(&qual).into_owned(),
                });
                state.preview_bytes += cost;
                if state.preview_bytes >= cfg.preview_flush_bytes {
                    flush_previews(state, tx);
                }
            }
        }

        // Corner extractions — check every registered box.
        for c in corners {
            if cx >= c.x1 && cx <= c.x2 && cy >= c.y1 && cy <= c.y2 {
                let mut w = c.writer.lock().unwrap();
                w.write_all(&header)?;
                w.write_all(b"\n")?;
                w.write_all(&seq)?;
                w.write_all(b"\n")?;
                w.write_all(&plus)?;
                w.write_all(b"\n")?;
                w.write_all(&qual)?;
                w.write_all(b"\n")?;
            }
        }

        read_num = read_num
            .checked_add(1)
            .expect("read_num overflowed u32");

        if state.bytes_estimate >= cfg.threshold_bytes {
            flush_local(state, tx);
        }
    }
    Ok(())
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn st() -> WorkerLocal {
        WorkerLocal::new()
    }

    // ── FM-1: the index stream is optional ──────────────────────────────────

    #[test]
    fn index_push_exact_keeps_every_read() {
        let mut s = st();
        for rn in 0..5 {
            assert!(s.push_index((3, 4), 0, rn, DataMode::Exact));
        }
        assert_eq!(s.indices.len(), 5);
        assert_eq!(s.bytes_estimate, 5 * 16);
        assert!(s.data_count.is_empty(), "exact mode must not track a cap");
    }

    #[test]
    fn index_push_capped_limits_per_pixel_and_file() {
        let mut s = st();
        // 3 reads at one pixel from file 0: only the first 2 are kept.
        assert!(s.push_index((1, 1), 0, 0, DataMode::Capped(2)));
        assert!(s.push_index((1, 1), 0, 1, DataMode::Capped(2)));
        assert!(!s.push_index((1, 1), 0, 2, DataMode::Capped(2)));
        // Same pixel, different file: its own quota.
        assert!(s.push_index((1, 1), 1, 0, DataMode::Capped(2)));
        // Different pixel: its own quota.
        assert!(s.push_index((2, 1), 0, 0, DataMode::Capped(2)));
        assert_eq!(s.indices.len(), 4);
        assert_eq!(
            s.indices,
            vec![(1, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0), (2, 1, 0, 0)]
        );
        // The budget only ever sees what is really stored.
        assert_eq!(s.bytes_estimate, 4 * 16);
    }

    #[test]
    fn index_push_off_stores_nothing() {
        let mut s = st();
        for rn in 0..10 {
            assert!(!s.push_index((0, 0), 0, rn, DataMode::Off));
        }
        assert!(s.indices.is_empty());
        assert_eq!(s.bytes_estimate, 0, "off: nothing buffered ⇒ nothing spills");
        assert!(s.data_count.is_empty());
    }

    /// The cap must survive a spill: `data_count` is not part of the chunk.
    #[test]
    fn capped_counter_survives_a_flush() {
        let mut s = st();
        let (tx, _rx) = mpsc::channel();
        for rn in 0..3 {
            s.push_index((5, 5), 0, rn, DataMode::Capped(2));
        }
        flush_local(&mut s, &tx);
        assert!(s.indices.is_empty());
        assert!(!s.push_index((5, 5), 0, 3, DataMode::Capped(2)));
        assert!(s.data_count.contains_key(&((5, 5), 0)));
    }

    // ── FM-2: the preview flush lives inside the budget ─────────────────────

    #[test]
    fn preview_flush_is_budget_aware() {
        // The old constant was 4 MiB *per worker*, invisible to --buffer-mb.
        assert_eq!(preview_flush_bytes(4 * 1024 * 1024), 256 * 1024);
        assert_eq!(preview_flush_bytes(64 * 1024 * 1024), 256 * 1024);
        // A worker share below the cap is honoured: never exceed the budget.
        assert_eq!(preview_flush_bytes(512 * 1024), 256 * 1024);
        assert_eq!(preview_flush_bytes(64 * 1024), 64 * 1024);
        // Degenerate tiny shares keep a usable batch size.
        assert_eq!(preview_flush_bytes(0), PREVIEW_FLUSH_MIN);
        assert_eq!(preview_flush_bytes(1), PREVIEW_FLUSH_MIN);
        assert!(preview_flush_bytes(8 * 1024 * 1024) <= 8 * 1024 * 1024);
    }

    // ── FM-10: byte-level line reading ──────────────────────────────────────

    #[test]
    fn read_line_strips_lf_and_crlf() {
        let data = b"@h\r\nACGT\r\n+\r\nIIII\r\n@z\nAC\n+\nII\n";
        let mut r = Cursor::new(&data[..]);
        let mut buf = Vec::new();
        let mut got = Vec::new();
        while read_line(&mut r, &mut buf).unwrap() > 0 {
            got.push(String::from_utf8(buf.clone()).unwrap());
        }
        assert_eq!(got, vec!["@h", "ACGT", "+", "IIII", "@z", "AC", "+", "II"]);
    }

    #[test]
    fn read_line_reports_eof_and_keeps_an_unterminated_last_line() {
        let mut r = Cursor::new(&b"abc"[..]);
        let mut buf = Vec::new();
        assert_eq!(read_line(&mut r, &mut buf).unwrap(), 3);
        assert_eq!(buf, b"abc");
        assert_eq!(read_line(&mut r, &mut buf).unwrap(), 0);
        assert!(buf.is_empty(), "buffer is cleared on every call");

        // A blank line is a line, not EOF.
        let mut r = Cursor::new(&b"\n\n"[..]);
        assert_eq!(read_line(&mut r, &mut buf).unwrap(), 1);
        assert!(buf.is_empty());
        assert_eq!(read_line(&mut r, &mut buf).unwrap(), 1);
        assert_eq!(read_line(&mut r, &mut buf).unwrap(), 0);
    }

    // ── FM-12: integer quality accumulation is order-independent ────────────

    #[test]
    fn pixel_acc_mean_is_order_independent() {
        let parts: Vec<(u64, u64)> = vec![(40, 1), (20, 1), (61, 1), (0, 1)];
        let mut fwd = PixelAcc::default();
        for &(q, n) in &parts {
            fwd.sum_q += q;
            fwd.sum_len += n;
            fwd.count += 1;
        }
        let mut rev = PixelAcc::default();
        for &(q, n) in parts.iter().rev() {
            rev.sum_q += q;
            rev.sum_len += n;
            rev.count += 1;
        }
        assert_eq!(fwd, rev);
        assert_eq!(fwd.mean_quality(), rev.mean_quality());
        assert_eq!(fwd.mean_quality(), 121.0 / 4.0);

        // Two workers' partials merged in either order agree bit for bit.
        let a = PixelAcc { count: 3, sum_q: 100, sum_len: 3 };
        let b = PixelAcc { count: 2, sum_q: 90, sum_len: 2 };
        let mut ab = a;
        ab.count += b.count;
        ab.sum_q += b.sum_q;
        ab.sum_len += b.sum_len;
        let mut ba = b;
        ba.count += a.count;
        ba.sum_q += a.sum_q;
        ba.sum_len += a.sum_len;
        assert_eq!(ab, ba);
        assert_eq!(ab.mean_quality().to_bits(), ba.mean_quality().to_bits());

        assert_eq!(PixelAcc::default().mean_quality(), 0.0);
    }
}
