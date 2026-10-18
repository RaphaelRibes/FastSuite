use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use crate::format::write_preview_entry;
use crate::worker::SpillerMsg;

// ── Types ────────────────────────────────────────────────────────────────────

/// All data collected by the spiller thread, handed back to main after
/// all workers have finished.
pub struct SpillerState {
    /// Combined per-pixel (count, cumulative_mean_quality).
    pub counts: HashMap<(u16, u16), (u64, f64)>,
    /// Per-pixel R1 counts.
    pub counts_r1: HashMap<(u16, u16), u64>,
    /// Per-pixel R2 counts.
    pub counts_r2: HashMap<(u16, u16), u64>,
    /// Per-file total read count.
    pub per_file_counts: HashMap<u16, u64>,
    pub max_seq_len: u32,
    pub total_reads: u64,
    pub temp_dir: PathBuf,
}

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

// ── Spiller thread ───────────────────────────────────────────────────────────

/// Receive chunks from workers, write index entries directly to per-y temp
/// files (8 bytes each: x:u16 LE | fid:u16 LE | rn:u32 LE), stream preview
/// reads to per-y preview temp files, and aggregate counts in memory.
///
/// Peak RAM ≈ the largest single row's index/preview buffer per open temp
/// writer + the mpsc queue; neither the index stream nor the preview payloads
/// are retained for the whole run.
pub fn spiller_thread(rx: mpsc::Receiver<SpillerMsg>, temp_dir: PathBuf) -> SpillerState {
    // Lazily-opened writers, one per y value (at most max_seq_len+1 open).
    let mut row_writers: Vec<Option<BufWriter<File>>> =
        (0..=u16::MAX as usize).map(|_| None).collect();

    let mut preview_writers: HashMap<u16, BufWriter<File>> = HashMap::new();
    let mut counts: HashMap<(u16, u16), (u64, f64)> = HashMap::new();
    let mut counts_r1: HashMap<(u16, u16), u64> = HashMap::new();
    let mut counts_r2: HashMap<(u16, u16), u64> = HashMap::new();
    let mut per_file_counts: HashMap<u16, u64> = HashMap::new();
    let mut max_seq_len = 0u32;
    let mut total_reads = 0u64;

    while let Ok(msg) = rx.recv() {
        match msg {
            SpillerMsg::Chunk(chunk) => {
                for (x, y, fid, rn) in chunk.indices {
                    let slot = &mut row_writers[y as usize];
                    if slot.is_none() {
                        let path = row_temp_path(&temp_dir, y);
                        let f = File::create(&path).unwrap_or_else(|e| {
                            eprintln!("Cannot create temp file {}: {e}", path.display());
                            std::process::exit(1);
                        });
                        *slot = Some(BufWriter::new(f));
                    }
                    let w = slot.as_mut().unwrap();
                    w.write_all(&x.to_le_bytes()).unwrap();
                    w.write_all(&fid.to_le_bytes()).unwrap();
                    w.write_all(&rn.to_le_bytes()).unwrap();
                }
                for p in chunk.previews {
                    let y = p.y;
                    let w = preview_writers.entry(y).or_insert_with(|| {
                        let path = preview_row_temp_path(&temp_dir, y);
                        let f = File::create(&path).unwrap_or_else(|e| {
                            eprintln!("Cannot create temp file {}: {e}", path.display());
                            std::process::exit(1);
                        });
                        BufWriter::new(f)
                    });
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
                max_seq_len: m,
                n_reads: n,
            } => {
                if m > max_seq_len {
                    max_seq_len = m;
                }
                total_reads += n;
                for (k, (cnt, q)) in c {
                    let e = counts.entry(k).or_insert((0, 0.0));
                    e.0 += cnt;
                    e.1 += q;
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
    for slot in row_writers.iter_mut().flatten() {
        slot.flush().unwrap_or_else(|e| {
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
        max_seq_len,
        total_reads,
        temp_dir,
    }
}
