use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::cli::DataMode;
use crate::encode::json_escape;
use crate::file_kind::FileKind;
use crate::spiller::{preview_row_temp_path, row_temp_path, SpillerState};
use crate::varint::varint_write;
use crate::worker::PreviewRead;

/// Everything `write_format` needs besides the spiller's state. A struct
/// rather than ten positional arguments: the data-mode knob (FM-1) was the
/// eleventh.
pub struct FormatOpts<'a> {
    pub files: &'a [String],
    pub file_kinds: &'a [FileKind],
    pub preview_n: usize,
    pub data_mode: DataMode,
    pub output_dir: &'a Path,
    /// `--json PATH`: where the `fastsuite.report.v1` document goes. The grid
    /// JSON is *not* redirected — it always lands at
    /// `<output_dir>/heatmap.json` (review 2026-10-21, item 1).
    pub report_path: Option<&'a Path>,
    pub elapsed_s: f64,
    pub quiet: bool,
}

/// Name of the grid JSON inside the output directory. One definition: the
/// writer, the report's `grid` block and the clobber guard all use it.
pub const GRID_JSON_NAME: &str = "heatmap.json";

/// Bytes charged per plane cell when deciding whether a run is allowed.
///
/// The end-of-run plane is three dense `u64` grids of `size²` cells (24 B/cell)
/// plus their JSON rendering and the `Vec<String>` join temporaries. Measured:
/// 61 B/cell at a 4 001×4 001 grid (976 MB). 64 B/cell rounds that up, so the
/// cap below bounds peak RSS rather than only the arrays.
pub const PLANE_BYTES_PER_CELL: usize = 64;

/// Ceiling on the output plane, published in `--help`.
///
/// 256 MiB at 64 B/cell ⇒ `MAX_PLANE_SIDE² = 2048²` cells ⇒ reads of at most
/// `MAX_PLANE_SIDE - 1 = 2047` bp. `--buffer-mb` sizes the read-index buffer
/// and cannot bound the plane, so without a ceiling a 65 kb read asks for
/// ~100 GB and the process is OOM-killed with no message.
pub const MAX_PLANE_BYTES: usize = 256 * 1024 * 1024;

/// Largest accepted `size` (= longest read + 1).
pub const MAX_PLANE_SIDE: usize = 2048;

// ── Row encoding ─────────────────────────────────────────────────────────────

/// Encode one row of index entries into the `.data` wire format.
///
/// Format per row:
/// ```text
/// row    := varint(n_pixels) pixel*
/// pixel  := varint(x) varint(n_groups) group*
/// group  := varint(file_id) varint(n_reads) varint(idx_0) varint(delta_1)…
/// ```
/// Read indices within a group are delta-encoded (ascending order).
pub fn encode_data_row(indices: Vec<(u16, u16, u32)>) -> Vec<u8> {
    // Group by x.
    let mut by_x: HashMap<u16, Vec<(u16, u32)>> = HashMap::new();
    for (x, fid, rn) in indices {
        by_x.entry(x).or_default().push((fid, rn));
    }
    let mut xs: Vec<u16> = by_x.keys().copied().collect();
    xs.sort_unstable();

    let mut buf = Vec::new();
    varint_write(&mut buf, xs.len() as u64);

    for x in xs {
        let mut entries = by_x.remove(&x).unwrap();
        entries.sort_unstable(); // (file_id, read_num)

        // Group consecutive runs of the same file_id.
        let mut groups: Vec<(u16, Vec<u32>)> = Vec::new();
        for (fid, rn) in entries {
            match groups.last_mut() {
                Some(g) if g.0 == fid => g.1.push(rn),
                _ => groups.push((fid, vec![rn])),
            }
        }

        varint_write(&mut buf, x as u64);
        varint_write(&mut buf, groups.len() as u64);
        for (fid, rns) in groups {
            varint_write(&mut buf, fid as u64);
            varint_write(&mut buf, rns.len() as u64);
            let mut prev = 0u32;
            for (i, &rn) in rns.iter().enumerate() {
                varint_write(&mut buf, if i == 0 { rn as u64 } else { (rn - prev) as u64 });
                prev = rn;
            }
        }
    }
    buf
}

/// Encode one row of preview reads into the `.preview` wire format.
///
/// Format per row:
/// ```text
/// row          := varint(n_pixels) pixel_preview*
/// pixel_preview := varint(x) varint(n_stored) read*
/// read          := varint(file_id) varint(read_num) varint(hdr_len) hdr
///                  varint(seq_len) seq qual
/// ```
/// `qual_len == seq_len` (not stored separately).
/// Reads per (pixel, kind) are kept as the first `preview_n` by
/// `(file_id, read_num)`, so R1 and R2 each get their own quota. The worker
/// keeps the first `preview_n` of *each file*, so this rule sees the same
/// candidates at any `-j`.
pub fn encode_preview_row(
    previews: Vec<PreviewRead>,
    preview_n: usize,
    file_kinds: &[FileKind],
) -> Vec<u8> {
    // Group by x.
    let mut by_x: HashMap<u16, Vec<PreviewRead>> = HashMap::new();
    for p in previews {
        by_x.entry(p.x).or_default().push(p);
    }

    // Per pixel: split by kind, take preview_n of each, merge.
    for v in by_x.values_mut() {
        let mut by_kind: HashMap<FileKind, Vec<PreviewRead>> = HashMap::new();
        for p in v.drain(..) {
            let kind = file_kinds
                .get(p.file_id as usize)
                .copied()
                .unwrap_or(FileKind::Unknown);
            by_kind.entry(kind).or_default().push(p);
        }
        let mut out = Vec::new();
        for (_kind, mut group) in by_kind {
            group.sort_unstable_by_key(|p| (p.file_id, p.read_num));
            group.truncate(preview_n);
            out.extend(group);
        }
        // Final order: by (file_id, read_num) so the wire layout is stable.
        out.sort_unstable_by_key(|p| (p.file_id, p.read_num));
        *v = out;
    }

    let mut xs: Vec<u16> = by_x.iter().filter(|(_, v)| !v.is_empty()).map(|(x, _)| *x).collect();
    xs.sort_unstable();

    let mut buf = Vec::new();
    varint_write(&mut buf, xs.len() as u64);
    for x in xs {
        let reads = by_x.get(&x).unwrap();
        varint_write(&mut buf, x as u64);
        varint_write(&mut buf, reads.len() as u64);
        for r in reads {
            varint_write(&mut buf, r.file_id as u64);
            varint_write(&mut buf, r.read_num as u64);
            varint_write(&mut buf, r.header.len() as u64);
            buf.extend_from_slice(r.header.as_bytes());
            varint_write(&mut buf, r.seq.len() as u64);
            buf.extend_from_slice(r.seq.as_bytes());
            buf.extend_from_slice(r.qual.as_bytes());
        }
    }
    buf
}

// ── Preview temp-file codec ──────────────────────────────────────────────────
//
// The spiller streams preview reads to one temp file per grid row so their
// payloads (header+seq+qual) never stay resident. Entry layout, little-endian:
//   x:u16 | file_id:u16 | read_num:u32 | hlen:u32 | header | slen:u32 | seq | qual
// `y` is implicit in the row's file name; `qual` has `slen` bytes.

/// Append one preview read to a row's preview temp file.
pub fn write_preview_entry(w: &mut impl Write, p: &PreviewRead) -> std::io::Result<()> {
    w.write_all(&p.x.to_le_bytes())?;
    w.write_all(&p.file_id.to_le_bytes())?;
    w.write_all(&p.read_num.to_le_bytes())?;
    w.write_all(&(p.header.len() as u32).to_le_bytes())?;
    w.write_all(p.header.as_bytes())?;
    w.write_all(&(p.seq.len() as u32).to_le_bytes())?;
    w.write_all(p.seq.as_bytes())?;
    w.write_all(p.qual.as_bytes())
}

/// Decode a whole preview temp row (written by `write_preview_entry`).
///
/// Every field is bounds-checked: a truncated or corrupt row is an error, not
/// a slice panic. The qual string carries no length of its own — the format
/// relies on `seq.len() == qual.len()`, which `process_fastq` enforces for
/// every record it reads (review 2026-10-21, item 3).
pub fn decode_preview_row(raw: &[u8], y: u16) -> std::io::Result<Vec<PreviewRead>> {
    fn corrupt(why: String) -> std::io::Error {
        std::io::Error::new(std::io::ErrorKind::InvalidData, why)
    }
    /// Take the next `n` bytes, or fail without moving `pos`.
    fn take<'a>(raw: &'a [u8], pos: &mut usize, n: usize) -> std::io::Result<&'a [u8]> {
        let end = pos
            .checked_add(n)
            .ok_or_else(|| corrupt(format!("field length {n} overflows at offset {pos}")))?;
        if end > raw.len() {
            return Err(corrupt(format!(
                "preview row truncated: {n} B needed at offset {pos}, row is {} B",
                raw.len()
            )));
        }
        let out = &raw[*pos..end];
        *pos = end;
        Ok(out)
    }
    fn u16_at(raw: &[u8], pos: &mut usize) -> std::io::Result<u16> {
        Ok(u16::from_le_bytes(take(raw, pos, 2)?.try_into().unwrap()))
    }
    fn u32_at(raw: &[u8], pos: &mut usize) -> std::io::Result<u32> {
        Ok(u32::from_le_bytes(take(raw, pos, 4)?.try_into().unwrap()))
    }
    let mut out = Vec::new();
    let mut pos = 0usize;
    while pos < raw.len() {
        let x = u16_at(raw, &mut pos)?;
        let file_id = u16_at(raw, &mut pos)?;
        let read_num = u32_at(raw, &mut pos)?;
        let hlen = u32_at(raw, &mut pos)? as usize;
        let header = String::from_utf8_lossy(take(raw, &mut pos, hlen)?).into_owned();
        let slen = u32_at(raw, &mut pos)? as usize;
        let seq = String::from_utf8_lossy(take(raw, &mut pos, slen)?).into_owned();
        let qual = String::from_utf8_lossy(take(raw, &mut pos, slen)?).into_owned();
        out.push(PreviewRead {
            x,
            y,
            file_id,
            read_num,
            header,
            seq,
            qual,
        });
    }
    Ok(out)
}

// ── Output ───────────────────────────────────────────────────────────────────

/// Write `heatmap.data`, `heatmap.preview`, and `heatmap.json` to `output_dir`.
///
/// Index entries are read one row at a time from the per-y temp files written
/// by the spiller, then deleted, so peak RAM ≈ the largest single row's data.
/// Preview payloads are read back the same way (one row file at a time).
///
/// With `--data-mode off` there is no render pass at all: no row file is read,
/// no decoder runs, and each row is written as the empty pixel list so the file
/// and its offsets stay structurally valid. `capped:N` renders whatever the
/// worker kept (at most N reads per (pixel, file)).
pub fn write_format(state: SpillerState, opts: &FormatOpts<'_>) -> std::io::Result<()> {
    let FormatOpts {
        files,
        file_kinds,
        preview_n,
        data_mode,
        output_dir,
        report_path,
        elapsed_s,
        quiet,
    } = *opts;

    let size = (state.max_seq_len + 1) as usize;
    // The plane is O(max_read_len²): three dense u64 grids of `size²` cells
    // plus their JSON rendering. `--buffer-mb` sizes the read-index buffer and
    // cannot bound it, so the plane gets its own published ceiling and a
    // refusal that names both (review 2026-10-21, item 4). 4 000 bp used to
    // peak at 976 MB and 8 000 bp at 3 788 MB, identically at --buffer-mb 1 and
    // 8; the u16 guard admitted ~65 kb reads (~100 GB).
    let cells = size.saturating_mul(size);
    let plane_bytes = cells.saturating_mul(PLANE_BYTES_PER_CELL);
    if plane_bytes > MAX_PLANE_BYTES {
        eprintln!(
            "Error: max read length {} bp needs a {size}x{size} plane ({cells} cells x \
             {PLANE_BYTES_PER_CELL} B/cell = {} MiB), over the {} MiB plane cap. The plane \
             grows as O(max read length^2); --buffer-mb bounds the read-index buffer only and \
             cannot bound it. The cap allows reads up to {} bp: trim reads to <= {} bp upstream \
             and re-run.",
            state.max_seq_len,
            plane_bytes / (1024 * 1024),
            MAX_PLANE_BYTES / (1024 * 1024),
            MAX_PLANE_SIDE - 1,
            MAX_PLANE_SIDE - 1,
        );
        std::process::exit(1);
    }
    if size > u16::MAX as usize {
        eprintln!("Error: max read length {} exceeds u16::MAX", state.max_seq_len);
        std::process::exit(1);
    }

    // The grid is sized for the *longest* read, so a shorter read cannot reach
    // a corner: on a mixed-length run the corner fractions silently describe
    // only the longest reads. The corner definition is unchanged here; instead
    // the run says so out loud, and both outputs publish the lengths (and the
    // window) the fractions were computed from, so a biased run is visible
    // after the fact even if nobody read stderr.
    let min_seq_len = state.min_seq_len.unwrap_or(state.max_seq_len);
    let corner_window = crate::report::corner_window(size as u32);
    let spread = state.max_seq_len - min_seq_len;
    if spread > corner_window {
        eprintln!(
            "fastmaps: warning: read lengths span {min_seq_len}..{max} bp (spread {spread} > \
             corner window {corner_window}); corner fractions are length-biased: a read \
             shorter than {max} bp cannot reach a corner",
            max = state.max_seq_len,
        );
    }

    let mut data_offsets: Vec<u64>    = Vec::with_capacity(size + 1);
    let mut preview_offsets: Vec<u64> = Vec::with_capacity(size + 1);
    let mut data_written:    u64 = 0;
    let mut preview_written: u64 = 0;
    data_offsets.push(0);
    preview_offsets.push(0);

    // The grid JSON keeps its own home: `--json` names the *report*, not the
    // grid (review 2026-10-21, item 1). One name for the writer, the report's
    // `grid` block and the clobber guard.
    let grid_path = output_dir.join(GRID_JSON_NAME);
    let grid_name = GRID_JSON_NAME;

    let mut data_file    = BufWriter::new(File::create(output_dir.join("heatmap.data"))?);
    let mut preview_file = BufWriter::new(File::create(output_dir.join("heatmap.preview"))?);

    // `varint(0)`: a row with no pixels. Off mode emits this for every row.
    const EMPTY_ROW: [u8; 1] = [0];

    for y in 0..size as u16 {
        let drow: Vec<u8> = if data_mode == DataMode::Off {
            EMPTY_ROW.to_vec()
        } else {
            // Read temp file for this row, then delete it immediately.
            let temp = row_temp_path(&state.temp_dir, y);
            let indices: Vec<(u16, u16, u32)> = if temp.exists() {
                let raw = std::fs::read(&temp)?;
                if raw.len() % 8 != 0 {
                    eprintln!("Corrupt temp file {} (len {})", temp.display(), raw.len());
                    std::process::exit(1);
                }
                let mut v = Vec::with_capacity(raw.len() / 8);
                for chunk in raw.chunks_exact(8) {
                    let x   = u16::from_le_bytes(chunk[0..2].try_into().unwrap());
                    let fid = u16::from_le_bytes(chunk[2..4].try_into().unwrap());
                    let rn  = u32::from_le_bytes(chunk[4..8].try_into().unwrap());
                    v.push((x, fid, rn));
                }
                std::fs::remove_file(&temp)?;
                v
            } else {
                vec![]
            };
            encode_data_row(indices)
        };

        let previews = {
            // Preview payloads were streamed to a per-row temp file by the
            // spiller; read (and delete) this row's file instead of holding
            // them in RAM for the whole run.
            let path = preview_row_temp_path(&state.temp_dir, y);
            if path.exists() {
                let raw = std::fs::read(&path)?;
                std::fs::remove_file(&path)?;
                decode_preview_row(&raw, y)?
            } else {
                Vec::new()
            }
        };

        let prow = encode_preview_row(previews, preview_n, file_kinds);

        data_file.write_all(&drow)?;
        preview_file.write_all(&prow)?;
        data_written    += drow.len() as u64;
        preview_written += prow.len() as u64;

        data_offsets.push(data_written);
        preview_offsets.push(preview_written);
    }

    data_file.flush()?;
    preview_file.flush()?;

    // Flat row-major counts grid: index = y*size + x.
    let mut counts_arr = vec![0u64; size * size];
    for (&(x, y), acc) in &state.counts {
        counts_arr[y as usize * size + x as usize] = acc.count;
    }
    let mut counts_r1_arr = vec![0u64; size * size];
    for (&(x, y), &c) in &state.counts_r1 {
        counts_r1_arr[y as usize * size + x as usize] = c;
    }
    let mut counts_r2_arr = vec![0u64; size * size];
    for (&(x, y), &c) in &state.counts_r2 {
        counts_r2_arr[y as usize * size + x as usize] = c;
    }

    // Corner metrics (formerly derived by fasthelper from this very file).
    // No-call reads (zero ACGT bases) sit on (0,0) in the grid like an all-A
    // read; each window passes its own count so the fractions can leave them
    // out instead of reporting them as polyA.
    let mut metrics = crate::report::corner_metrics(
        size as u32,
        &[
            ("", &counts_arr, state.no_call),
            ("r1", &counts_r1_arr, state.no_call_r1),
            ("r2", &counts_r2_arr, state.no_call_r2),
        ],
    );
    // The lengths the corner window was derived from: a consumer can see a
    // length-biased run from the report alone (`max - min > corner_window`).
    metrics.insert("min_seq_len".into(), min_seq_len as f64);
    metrics.insert("max_seq_len".into(), state.max_seq_len as f64);
    metrics.insert("corner_window".into(), corner_window as f64);
    // FM-1: a non-exact data mode is a labelled approximation of
    // `heatmap.data` alone. Both the grid JSON and the v1 wrapper carry the
    // label, so no consumer can mistake a capped/absent stream for the exact
    // one. `data_complete` is a boolean in the grid JSON and 0.0 in the metric
    // map (which is `f64`-valued by schema).
    let data_labels = if data_mode.is_exact() {
        String::new()
    } else {
        metrics.insert("data_complete".into(), 0.0);
        metrics.insert(
            "data_reads_per_pixel".into(),
            data_mode.reads_per_pixel() as f64,
        );
        format!(
            "  \"data_reads_per_pixel\": {},\n  \"data_complete\": false,\n",
            data_mode.reads_per_pixel()
        )
    };
    // A dense Vec indexed by file id. Out-of-range ids are a gap, not an
    // error, so they serialise as 0 exactly as the old HashMap's `.get()` did.
    let per_file_json = (0..files.len())
        .map(|i| state.per_file_counts.get(i).copied().unwrap_or(0).to_string())
        .collect::<Vec<_>>()
        .join(",");

    // File kinds, ordered by file index.
    let file_kinds_json = file_kinds
        .iter()
        .map(|k| format!("\"{}\"", k.as_str()))
        .collect::<Vec<_>>()
        .join(",");

    let counts_join =
        |arr: &[u64]| arr.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(",");

    let json = format!(
        concat!(
            "{{\n",
            "  \"version\": 2,\n",
            "  \"size\": {size},\n",
            "  \"min_seq_len\": {min_seq_len},\n",
            "  \"max_seq_len\": {max_len},\n",
            "  \"corner_window\": {corner_window},\n",
            "  \"files\": [{files}],\n",
            "  \"file_kinds\": [{kinds}],\n",
            "  \"total_reads\": {total},\n",
            "  \"per_file_counts\": [{pfc}],\n",
            "  \"preview_reads_per_pixel\": {pn},\n",
            "{labels}",
            "  \"counts\": [{counts}],\n",
            "  \"counts_r1\": [{counts_r1}],\n",
            "  \"counts_r2\": [{counts_r2}],\n",
            "  \"data_row_offsets\": [{doff}],\n",
            "  \"preview_row_offsets\": [{poff}]\n",
            "}}\n"
        ),
        size      = size,
        min_seq_len = min_seq_len,
        max_len   = state.max_seq_len,
        corner_window = corner_window,
        files     = files.iter().map(|f| json_escape(f)).collect::<Vec<_>>().join(","),
        kinds     = file_kinds_json,
        total     = state.total_reads,
        pfc       = per_file_json,
        pn        = preview_n,
        labels    = data_labels,
        counts    = counts_join(&counts_arr),
        counts_r1 = counts_join(&counts_r1_arr),
        counts_r2 = counts_join(&counts_r2_arr),
        doff      = data_offsets.iter().map(|o| o.to_string()).collect::<Vec<_>>().join(","),
        poff      = preview_offsets.iter().map(|o| o.to_string()).collect::<Vec<_>>().join(","),
    );

    let mut jf = BufWriter::new(File::create(grid_path)?);
    jf.write_all(json.as_bytes())?;
    jf.flush()?;

    // `--json PATH` writes the v1 report at PATH. The grid JSON is never
    // redirected: it is `<output-dir>/heatmap.json`, so `-o` always has it
    // (review 2026-10-21, item 1). The report's `grid` block names that file.
    let mut report_suffix = String::new();
    if let Some(report_path) = report_path {
        crate::report::write_report_v1(
            report_path,
            &std::env::args().collect::<Vec<_>>().join(" "),
            elapsed_s,
            files,
            &metrics,
            grid_name,
        )?;
        report_suffix = format!(
            " + {} (fastsuite.report.v1)",
            report_path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| report_path.display().to_string())
        );
    }

    if !quiet {
        println!(
            "Saved heatmap.data ({data_written} B), \
             heatmap.preview ({preview_written} B), \
             heatmap.json (size={size}, {} reads){report_suffix}",
            state.total_reads,
        );
    }
    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::varint::varint_read;
    #[test]
    fn data_row_empty() {
        let buf = encode_data_row(vec![]);
        assert_eq!(buf, vec![0]); // varint(0)
    }

    #[test]
    fn data_row_single() {
        let buf = encode_data_row(vec![(3, 2, 42)]);
        let mut pos = 0;
        assert_eq!(varint_read(&buf, &mut pos), 1);  // n_pixels
        assert_eq!(varint_read(&buf, &mut pos), 3);  // x
        assert_eq!(varint_read(&buf, &mut pos), 1);  // n_groups
        assert_eq!(varint_read(&buf, &mut pos), 2);  // file_id
        assert_eq!(varint_read(&buf, &mut pos), 1);  // n_reads
        assert_eq!(varint_read(&buf, &mut pos), 42); // read_num
        assert_eq!(pos, buf.len());
    }

    #[test]
    fn data_row_deltas() {
        let buf = encode_data_row(vec![(5, 1, 100), (5, 1, 110), (5, 1, 115)]);
        let mut pos = 0;
        assert_eq!(varint_read(&buf, &mut pos), 1);   // n_pixels
        assert_eq!(varint_read(&buf, &mut pos), 5);   // x
        assert_eq!(varint_read(&buf, &mut pos), 1);   // n_groups
        assert_eq!(varint_read(&buf, &mut pos), 1);   // file_id
        assert_eq!(varint_read(&buf, &mut pos), 3);   // n_reads
        assert_eq!(varint_read(&buf, &mut pos), 100); // abs
        assert_eq!(varint_read(&buf, &mut pos), 10);  // delta
        assert_eq!(varint_read(&buf, &mut pos), 5);   // delta
    }

    fn make_preview(fid: u16, rn: u32) -> PreviewRead {
        PreviewRead {
            x: 0, y: 0, file_id: fid, read_num: rn,
            header: format!("@r{rn}"),
            seq: "ACGT".into(),
            qual: "IIII".into(),
        }
    }

    #[test]
    fn preview_row_truncates() {
        // All file_ids are 0/1/2, all unknown kind → single bucket truncates to 3.
        let kinds = vec![FileKind::Unknown; 4];
        let buf = encode_preview_row(
            vec![
                make_preview(1, 50), make_preview(0, 10), make_preview(0, 5),
                make_preview(2, 1),  make_preview(0, 20), make_preview(1, 30),
            ],
            3,
            &kinds,
        );
        let mut pos = 0;
        assert_eq!(varint_read(&buf, &mut pos), 1); // n_pixels
        assert_eq!(varint_read(&buf, &mut pos), 0); // x
        assert_eq!(varint_read(&buf, &mut pos), 3); // n_stored — truncated
    }

    #[test]
    fn preview_row_truncates_per_kind() {
        // 4 R1 reads + 4 R2 reads at one pixel, preview_n=2 → 2 R1 + 2 R2 = 4.
        let kinds = vec![
            FileKind::R1, FileKind::R1, FileKind::R1, FileKind::R1,
            FileKind::R2, FileKind::R2, FileKind::R2, FileKind::R2,
        ];
        let mut reads = Vec::new();
        for fid in 0..8u16 {
            reads.push(make_preview(fid, 10));
        }
        let buf = encode_preview_row(reads, 2, &kinds);
        let mut pos = 0;
        assert_eq!(varint_read(&buf, &mut pos), 1); // n_pixels
        assert_eq!(varint_read(&buf, &mut pos), 0); // x
        assert_eq!(varint_read(&buf, &mut pos), 4); // 2 per kind × 2 kinds
    }

    #[test]
    fn preview_temp_entry_roundtrip() {
        let mut p = make_preview(2, 42);
        p.x = 7;
        p.y = 9;
        p.header = "@NS500:1:FC".into();
        p.seq = "ACGTACGTAC".into();
        p.qual = "IIIIIIIIII".into();
        let mut buf = Vec::new();
        write_preview_entry(&mut buf, &p).unwrap();
        let back = decode_preview_row(&buf, 9).unwrap();
        assert_eq!(back.len(), 1);
        let q = &back[0];
        assert_eq!(
            (q.x, q.y, q.file_id, q.read_num),
            (7, 9, 2, 42)
        );
        assert_eq!(q.header, p.header);
        assert_eq!(q.seq, p.seq);
        assert_eq!(q.qual, p.qual);
    }

    #[test]
    fn preview_temp_row_roundtrip_multiple() {
        let mut reads = Vec::new();
        for i in 0..5u32 {
            let mut p = make_preview(i as u16, i);
            p.x = (i * 3) as u16;
            p.seq = "ACGT".repeat(1 + i as usize);
            p.qual = "I".repeat(p.seq.len());
            reads.push(p);
        }
        let mut buf = Vec::new();
        for p in &reads {
            write_preview_entry(&mut buf, p).unwrap();
        }
        assert!(decode_preview_row(&[], 3).unwrap().is_empty());
        let back = decode_preview_row(&buf, 3).unwrap();
        assert_eq!(back.len(), reads.len());
        for (a, b) in back.iter().zip(&reads) {
            assert_eq!((a.x, a.file_id, a.read_num, &a.seq, &a.qual),
                       (b.x, b.file_id, b.read_num, &b.seq, &b.qual));
        }
    }

    /// A truncated row is an error, never a slice panic — every cut point of a
    /// valid entry, which is the shape of the short-qual record that used to
    /// abort the process with exit 101 (review 2026-10-21, item 3).
    #[test]
    fn preview_row_truncated_is_an_error_not_a_panic() {
        let mut p = make_preview(0, 1);
        p.seq = "ACGTACGT".into();
        p.qual = "IIIIIIII".into();
        let mut buf = Vec::new();
        write_preview_entry(&mut buf, &p).unwrap();
        for cut in 1..buf.len() {
            assert!(
                decode_preview_row(&buf[..cut], 0).is_err(),
                "a row cut at {cut} of {} B must be reported, not decoded",
                buf.len()
            );
        }
        assert_eq!(decode_preview_row(&buf, 0).unwrap().len(), 1);
    }

    /// The plane constants are one equation: the accepted side is exactly the
    /// largest square the byte cap allows.
    #[test]
    fn plane_constants_agree() {
        assert_eq!(
            MAX_PLANE_SIDE * MAX_PLANE_SIDE * PLANE_BYTES_PER_CELL,
            MAX_PLANE_BYTES
        );
        assert_eq!(MAX_PLANE_SIDE, 2048);
    }
}
