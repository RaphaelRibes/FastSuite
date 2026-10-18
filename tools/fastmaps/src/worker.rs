use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::sync::{mpsc, Arc, Mutex};

use crate::encode::{encode_read, get_reader, mean_quality};
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

pub struct SpillChunk {
    pub indices: Vec<(u16, u16, u16, u32)>, // (x, y, file_id, read_num)
    pub previews: Vec<PreviewRead>,
}

pub enum SpillerMsg {
    Chunk(SpillChunk),
    Done {
        counts: HashMap<(u16, u16), (u64, f64)>,
        counts_r1: HashMap<(u16, u16), u64>,
        counts_r2: HashMap<(u16, u16), u64>,
        per_file_counts: HashMap<u16, u64>,
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
}

/// Byte threshold for the *preview* buffer. Previews carry read payloads
/// (header+seq+qual), so they are flushed on their own, much sooner than the
/// index buffer: the spiller streams each one to its row's temp file, which
/// keeps preview memory flat in library size instead of `preview_n x pixels`
/// read payloads resident for the whole run.
const PREVIEW_FLUSH_BYTES: usize = 4 * 1024 * 1024;

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
    /// Per-pixel-per-kind preview cap — persists across spills.
    /// Keyed by ((x, y), kind) so each kind gets its own preview_n quota.
    pub preview_count: HashMap<((u16, u16), FileKind), u8>,
    /// Per-pixel (count, cumulative_mean_quality) — for PNG/HTML.
    pub counts: HashMap<(u16, u16), (u64, f64)>,
    /// Per-pixel R1-only count.
    pub counts_r1: HashMap<(u16, u16), u64>,
    /// Per-pixel R2-only count.
    pub counts_r2: HashMap<(u16, u16), u64>,
    /// Per-file total read count.
    pub per_file_counts: HashMap<u16, u64>,
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
            counts: HashMap::new(),
            counts_r1: HashMap::new(),
            counts_r2: HashMap::new(),
            per_file_counts: HashMap::new(),
            max_seq_len: 0,
            total_reads: 0,
            bytes_estimate: 0,
            preview_bytes: 0,
        }
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
    let reader = get_reader(path)?;
    let mut lines = reader.lines();
    let mut read_num: u32 = 0;
    let kind = file_kinds
        .get(file_id as usize)
        .copied()
        .unwrap_or(FileKind::Unknown);

    while let (Some(Ok(header)), Some(Ok(seq)), Some(Ok(plus)), Some(Ok(qual))) =
        (lines.next(), lines.next(), lines.next(), lines.next())
    {
        let (cx, cy) = encode_read(&seq);
        let mq = mean_quality(&qual);

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

        let xy = (cx as u16, cy as u16);

        // Per-pixel counts (combined — for PNG / HTML).
        let e = state.counts.entry(xy).or_insert((0, 0.0));
        e.0 += 1;
        e.1 += mq;
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
        state.total_reads += 1;
        *state.per_file_counts.entry(file_id).or_insert(0) += 1;

        // Index entry for heatmap.data.
        state.indices.push((xy.0, xy.1, file_id, read_num));
        state.bytes_estimate += 16;

        // Preview accumulation (cap per (pixel, kind), so R1/R2 each get
        // their own preview_n quota — otherwise alphabetical file_id order
        // could starve R2 entirely).
        if cfg.preview_n > 0 {
            let pc = state.preview_count.entry((xy, kind)).or_insert(0);
            if (*pc as usize) < cfg.preview_n {
                *pc += 1;
                let cost = header.len() + seq.len() + qual.len() + 64;
                state.previews.push(PreviewRead {
                    x: xy.0,
                    y: xy.1,
                    file_id,
                    read_num,
                    header: header.clone(),
                    seq: seq.clone(),
                    qual: qual.clone(),
                });
                state.preview_bytes += cost;
                if state.preview_bytes >= PREVIEW_FLUSH_BYTES {
                    flush_previews(state, tx);
                }
            }
        }

        // Corner extractions — check every registered box.
        for c in corners {
            if cx >= c.x1 && cx <= c.x2 && cy >= c.y1 && cy <= c.y2 {
                let mut w = c.writer.lock().unwrap();
                writeln!(w, "{}\n{}\n{}\n{}", header, seq, plus, qual)?;
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
