use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use crate::format::write_preview_entry;
use crate::worker::{PixelAcc, SpillerMsg};

// ── Types ────────────────────────────────────────────────────────────────────

/// All data collected by the spiller thread, handed back to main after
/// all workers have finished.
pub struct SpillerState {
    /// Combined per-pixel (count, quality sums).
    pub counts: HashMap<(u16, u16), PixelAcc>,
    /// Per-pixel R1 counts.
    pub counts_r1: HashMap<(u16, u16), u64>,
    /// Per-pixel R2 counts.
    pub counts_r2: HashMap<(u16, u16), u64>,
    /// Per-file total read count.
    pub per_file_counts: HashMap<u16, u64>,
    /// Reads with zero ACGT bases (no-call): unplaceable on the composition
    /// plane, so the corner metrics exclude them instead of calling them polyA.
    pub no_call: u64,
    pub no_call_r1: u64,
    pub no_call_r2: u64,
    /// Shortest / longest read length over the whole run. `min_seq_len` is
    /// `None` only when no read was seen at all.
    pub min_seq_len: Option<u32>,
    pub max_seq_len: u32,
    pub total_reads: u64,
    pub temp_dir: PathBuf,
}

/// One spilled index record: `x:u16 LE | file_id:u16 LE | read_num:u32 LE`.
/// `y` is implicit in the row's file name.
pub type IndexRecord = [u8; 8];

/// Size of one spilled index record on disk.
pub const INDEX_RECORD_BYTES: usize = 8;

/// Per-row `BufWriter` capacity.
///
/// FM-3: this was `BufWriter::new` (8 KiB) on 152 row files *plus* 152 preview
/// row files — 2.5 MB of buffers, most of it idle. 4 KiB still holds 512
/// records, which is far more than the per-flush row run a worker produces.
pub const ROW_WRITER_BUF: usize = 4 * 1024;

// ── Helpers ──────────────────────────────────────────────────────────────────

pub fn row_temp_path(temp_dir: &Path, y: u16) -> PathBuf {
    temp_dir.join(format!("heatmap_row_{y:05}.tmp"))
}

/// Temp file holding one row's preview reads (header+seq+qual payloads).
/// Streamed to disk by the spiller so preview memory never grows with the
/// library; format.rs reads a row, then deletes the file.
pub fn preview_row_temp_path(temp_dir: &Path, y: u16) -> PathBuf {
    temp_dir.join(format!("preview_row_{y:05}.tmp"))
}

fn create_row_writer(temp_dir: &Path, y: u16, tag: &str) -> BufWriter<File> {
    let path = if tag.is_empty() {
        row_temp_path(temp_dir, y)
    } else {
        preview_row_temp_path(temp_dir, y)
    };
    let f = File::create(&path).unwrap_or_else(|e| {
        eprintln!("Cannot create temp file {}: {e}", path.display());
        std::process::exit(1);
    });
    BufWriter::with_capacity(ROW_WRITER_BUF, f)
}

// ── Spiller thread ───────────────────────────────────────────────────────────

/// Receive chunks from workers, write index entries directly to per-y temp
/// files (8 bytes each: x:u16 LE | fid:u16 LE | rn:u32 LE), stream preview
/// reads to per-y preview temp files, and aggregate counts in memory.
///
/// Peak RAM ≈ the largest single row's index/preview buffer per open temp
/// writer + the mpsc queue; neither the index stream nor the preview payloads
/// are retained for the whole run.
///
/// FM-3: the writer tables are `HashMap`s keyed by row, not a
/// `Vec<Option<_>>` sized to `u16::MAX` (2.1 MB of `None`s allocated before a
/// single row is seen). A short-read run has `max_seq_len + 1` ≈ 152 rows, so
/// the map holds ~152 entries whatever the read length.
pub fn spiller_thread(rx: mpsc::Receiver<SpillerMsg>, temp_dir: PathBuf) -> SpillerState {
    // Lazily-opened writers, one per y value actually seen.
    let mut row_writers: HashMap<u16, BufWriter<File>> = HashMap::new();
    let mut preview_writers: HashMap<u16, BufWriter<File>> = HashMap::new();

    let mut counts: HashMap<(u16, u16), PixelAcc> = HashMap::new();
    let mut counts_r1: HashMap<(u16, u16), u64> = HashMap::new();
    let mut counts_r2: HashMap<(u16, u16), u64> = HashMap::new();
    let mut per_file_counts: HashMap<u16, u64> = HashMap::new();
    let mut no_call = 0u64;
    let mut no_call_r1 = 0u64;
    let mut no_call_r2 = 0u64;
    let mut min_seq_len: Option<u32> = None;
    let mut max_seq_len = 0u32;
    let mut total_reads = 0u64;

    while let Ok(msg) = rx.recv() {
        match msg {
            SpillerMsg::Chunk(chunk) => {
                for (x, y, fid, rn) in chunk.indices {
                    let w = row_writers
                        .entry(y)
                        .or_insert_with(|| create_row_writer(&temp_dir, y, ""));
                    // FM-9: one packed 8-byte store instead of three
                    // `write_all` calls per record (30 M → 10 M per 10 M reads).
                    let mut rec: IndexRecord = [0u8; INDEX_RECORD_BYTES];
                    rec[0..2].copy_from_slice(&x.to_le_bytes());
                    rec[2..4].copy_from_slice(&fid.to_le_bytes());
                    rec[4..8].copy_from_slice(&rn.to_le_bytes());
                    w.write_all(&rec).unwrap();
                }
                for p in chunk.previews {
                    let y = p.y;
                    let w = preview_writers
                        .entry(y)
                        .or_insert_with(|| create_row_writer(&temp_dir, y, "preview"));
                    write_preview_entry(w, &p).unwrap_or_else(|e| {
                        eprintln!("Error writing preview temp row: {e}");
                        std::process::exit(1);
                    });
                }
            }
            SpillerMsg::Done {
                counts: c,
                counts_r1: cr1,
                counts_r2: cr2,
                per_file_counts: pfc,
                no_call: nc,
                no_call_r1: nc1,
                no_call_r2: nc2,
                min_seq_len: mn,
                max_seq_len: m,
                n_reads: n,
            } => {
                if m > max_seq_len {
                    max_seq_len = m;
                }
                min_seq_len = match (min_seq_len, mn) {
                    (None, x) => x,
                    (a, None) => a,
                    (Some(a), Some(b)) => Some(a.min(b)),
                };
                total_reads += n;
                no_call += nc;
                no_call_r1 += nc1;
                no_call_r2 += nc2;
                // Integer accumulation: exact in any `Done` arrival order
                // (FM-12). The old `f64` partial means were not.
                for (k, p) in c {
                    let e = counts.entry(k).or_default();
                    e.count += p.count;
                    e.sum_q += p.sum_q;
                    e.sum_len += p.sum_len;
                }
                for (k, cnt) in cr1 {
                    *counts_r1.entry(k).or_insert(0) += cnt;
                }
                for (k, cnt) in cr2 {
                    *counts_r2.entry(k).or_insert(0) += cnt;
                }
                for (fid, cnt) in pfc {
                    *per_file_counts.entry(fid).or_insert(0) += cnt;
                }
            }
        }
    }

    // Flush all open writers before write_format reads them.
    for w in row_writers.values_mut() {
        w.flush().unwrap_or_else(|e| {
            eprintln!("Error flushing spill file: {e}");
            std::process::exit(1);
        });
    }
    for w in preview_writers.values_mut() {
        w.flush().unwrap_or_else(|e| {
            eprintln!("Error flushing preview spill file: {e}");
            std::process::exit(1);
        });
    }
    drop(row_writers); // close all handles
    drop(preview_writers);

    SpillerState {
        counts,
        counts_r1,
        counts_r2,
        per_file_counts,
        no_call,
        no_call_r1,
        no_call_r2,
        min_seq_len,
        max_seq_len,
        total_reads,
        temp_dir,
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::worker::{PixelAcc, SpillChunk};

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("maps_spill_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn done(
        counts: &[((u16, u16), PixelAcc)],
        n_reads: u64,
        max_seq_len: u32,
    ) -> SpillerMsg {
        SpillerMsg::Done {
            counts: counts.iter().copied().collect(),
            counts_r1: HashMap::new(),
            counts_r2: HashMap::new(),
            per_file_counts: HashMap::new(),
            no_call: 0,
            no_call_r1: 0,
            no_call_r2: 0,
            min_seq_len: Some(max_seq_len),
            max_seq_len,
            n_reads,
        }
    }

    /// FM-3 + FM-9: one temp file per *seen* row (not per possible row), each
    /// holding packed 8-byte records in arrival order.
    #[test]
    fn one_file_per_seen_row_with_packed_records() {
        let dir = tmpdir("rows");
        let (tx, rx) = mpsc::channel();
        tx.send(SpillerMsg::Chunk(SpillChunk {
            indices: vec![(7, 3, 0, 11), (9, 3, 1, 2)],
            previews: vec![],
        }))
        .unwrap();
        // A second chunk for an already-open row appends; a new row opens.
        tx.send(SpillerMsg::Chunk(SpillChunk {
            indices: vec![(2, 3, 0, 12), (1, 7, 0, 0)],
            previews: vec![],
        }))
        .unwrap();
        tx.send(done(&[], 4, 9)).unwrap();
        drop(tx);

        let state = spiller_thread(rx, dir.clone());
        assert_eq!(state.total_reads, 4);
        assert_eq!(state.max_seq_len, 9);

        // Exactly two row files exist: rows 3 and 7.
        let mut names: Vec<String> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
            .collect();
        names.sort();
        assert_eq!(names.len(), 2, "writer map must only open seen rows: {names:?}");
        assert_eq!(names, vec!["heatmap_row_00003.tmp", "heatmap_row_00007.tmp"]);

        let row3 = std::fs::read(row_temp_path(&dir, 3)).unwrap();
        assert_eq!(row3.len(), 3 * INDEX_RECORD_BYTES);
        // Byte layout: x LE | fid LE | rn LE, in arrival order across chunks.
        let rec = |i: usize| -> (u16, u16, u32) {
            let r = &row3[i * 8..i * 8 + 8];
            (
                u16::from_le_bytes(r[0..2].try_into().unwrap()),
                u16::from_le_bytes(r[2..4].try_into().unwrap()),
                u32::from_le_bytes(r[4..8].try_into().unwrap()),
            )
        };
        assert_eq!(rec(0), (7, 0, 11));
        assert_eq!(rec(1), (9, 1, 2));
        assert_eq!(rec(2), (2, 0, 12));

        let row7 = std::fs::read(row_temp_path(&dir, 7)).unwrap();
        assert_eq!(row7, vec![1, 0, 0, 0, 0, 0, 0, 0]);

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// FM-12: two workers' quality partials merge to the same value in either
    /// `Done` order — that is the whole point of the integer accumulator.
    #[test]
    fn done_merge_is_order_independent() {
        let a = ((4u16, 5u16), PixelAcc { count: 3, sum_q: 100, sum_len: 3 });
        let b = ((4u16, 5u16), PixelAcc { count: 2, sum_q: 71, sum_len: 2 });

        let merge = |first: &SpillerMsg, second: &SpillerMsg| {
            let dir = tmpdir("order");
            let (tx, rx) = mpsc::channel();
            tx.send(first.clone()).unwrap();
            tx.send(second.clone()).unwrap();
            drop(tx);
            let s = spiller_thread(rx, dir.clone());
            let _ = std::fs::remove_dir_all(&dir);
            s.counts[&(4, 5)]
        };

        // Two workers, one Done message each for the same pixel: the merge
        // order is the message arrival order, which races.
        let ab = merge(&done(&[a], 3, 6), &done(&[b], 2, 0));
        let ba = merge(&done(&[b], 2, 0), &done(&[a], 3, 6));
        assert_eq!(ab, ba);
        assert_eq!(ab.mean_quality().to_bits(), ba.mean_quality().to_bits());
        assert_eq!(ab.count, 5);
        assert_eq!(ab.sum_q, 171);
        assert_eq!(ab.sum_len, 5);
        assert_eq!(ab.mean_quality(), 171.0 / 5.0);
    }
}
