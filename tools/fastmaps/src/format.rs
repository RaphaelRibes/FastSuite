use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::encode::json_escape;
use crate::file_kind::FileKind;
use crate::spiller::{row_temp_path, SpillerState};
use crate::varint::varint_write;
use crate::worker::PreviewRead;

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
/// `(file_id, read_num)`, so R1 and R2 each get their own quota.
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

// ── Output ───────────────────────────────────────────────────────────────────

/// Write `heatmap.data`, `heatmap.preview`, and `heatmap.json` to `output_dir`.
///
/// Index entries are read one row at a time from the per-y temp files written
/// by the spiller, then deleted, so peak RAM ≈ the largest single row's data.
pub fn write_format(
    mut state: SpillerState,
    files: &[String],
    file_kinds: &[FileKind],
    preview_n: usize,
    output_dir: &Path,
) -> std::io::Result<()> {
    let size = (state.max_seq_len + 1) as usize;
    if size > u16::MAX as usize {
        eprintln!("Error: max read length {} exceeds u16::MAX", state.max_seq_len);
        std::process::exit(1);
    }

    let mut data_offsets: Vec<u64>    = Vec::with_capacity(size + 1);
    let mut preview_offsets: Vec<u64> = Vec::with_capacity(size + 1);
    let mut data_written:    u64 = 0;
    let mut preview_written: u64 = 0;
    data_offsets.push(0);
    preview_offsets.push(0);

    let mut data_file    = BufWriter::new(File::create(output_dir.join("heatmap.data"))?);
    let mut preview_file = BufWriter::new(File::create(output_dir.join("heatmap.preview"))?);

    for y in 0..size as u16 {
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

        let previews = state.preview_rows.remove(&y).unwrap_or_default();

        let drow = encode_data_row(indices);
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
    for (&(x, y), &(c, _)) in &state.counts {
        counts_arr[y as usize * size + x as usize] = c;
    }
    let mut counts_r1_arr = vec![0u64; size * size];
    for (&(x, y), &c) in &state.counts_r1 {
        counts_r1_arr[y as usize * size + x as usize] = c;
    }
    let mut counts_r2_arr = vec![0u64; size * size];
    for (&(x, y), &c) in &state.counts_r2 {
        counts_r2_arr[y as usize * size + x as usize] = c;
    }

    // Per-file counts, ordered by file index.
    let per_file_json = (0..files.len())
        .map(|i| {
            state
                .per_file_counts
                .get(&(i as u16))
                .copied()
                .unwrap_or(0)
                .to_string()
        })
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
            "  \"files\": [{files}],\n",
            "  \"file_kinds\": [{kinds}],\n",
            "  \"total_reads\": {total},\n",
            "  \"per_file_counts\": [{pfc}],\n",
            "  \"preview_reads_per_pixel\": {pn},\n",
            "  \"counts\": [{counts}],\n",
            "  \"counts_r1\": [{counts_r1}],\n",
            "  \"counts_r2\": [{counts_r2}],\n",
            "  \"data_row_offsets\": [{doff}],\n",
            "  \"preview_row_offsets\": [{poff}]\n",
            "}}\n"
        ),
        size      = size,
        files     = files.iter().map(|f| json_escape(f)).collect::<Vec<_>>().join(","),
        kinds     = file_kinds_json,
        total     = state.total_reads,
        pfc       = per_file_json,
        pn        = preview_n,
        counts    = counts_join(&counts_arr),
        counts_r1 = counts_join(&counts_r1_arr),
        counts_r2 = counts_join(&counts_r2_arr),
        doff      = data_offsets.iter().map(|o| o.to_string()).collect::<Vec<_>>().join(","),
        poff      = preview_offsets.iter().map(|o| o.to_string()).collect::<Vec<_>>().join(","),
    );

    let mut jf = BufWriter::new(File::create(output_dir.join("heatmap.json"))?);
    jf.write_all(json.as_bytes())?;
    jf.flush()?;

    println!(
        "Saved heatmap.data ({data_written} B), \
         heatmap.preview ({preview_written} B), \
         heatmap.json (size={size}, {} reads)",
        state.total_reads
    );
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
}
