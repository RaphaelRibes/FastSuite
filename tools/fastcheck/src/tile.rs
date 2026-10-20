//! Per-tile quality (PLAN-QC module 7).
//!
//! FastQC's per-tile modules answer "is one region of the flowcell failing", which a
//! whole-file quality curve cannot see: a bad tile is invisible once its reads are averaged
//! with the good ones. The accumulator is `tile x position x quality cell`, ~2.5 MB per 100
//! tiles at R=150, so it is off by default (`--tiles`) and the projection is printed when
//! it is on.
//!
//! Only 7-field Illumina headers carry a tile. Reads with any other header shape are
//! counted in `reads_without_coordinates` rather than dropped: a run whose headers were
//! rewritten must not silently lose its reads.

use std::collections::BTreeMap;

use crate::accum::{qbin, quantile, QBINS};

/// Extract the tile number from a FASTQ header, if it has the 7-field Illumina form.
///
/// `@instrument:run:flowcell:lane:tile:x:y <rest>` — the tile is field 5 (index 4) of the
/// part before the first space. Anything else (a plain name, a rewritten header, a
/// different vendor's layout) returns `None`.
///
/// The seven field starts are walked with an index and held in a fixed array: the
/// `split(':').collect::<Vec<_>>()` this replaces allocated a `Vec` per read — 20 M
/// allocations on the 20 M-record cell — to read one field.
pub fn tile_of(header: &[u8]) -> Option<u32> {
    let end = header.iter().position(|&b| b == b' ').unwrap_or(header.len());
    let f = &header[..end];
    let mut starts = [0usize; 7];
    let mut nfields = 0usize;
    let mut field_start = 0usize;
    for i in 0..=f.len() {
        if i == f.len() || f[i] == b':' {
            if nfields == 7 {
                return None; // more than seven fields
            }
            starts[nfields] = field_start;
            nfields += 1;
            field_start = i + 1;
        }
    }
    if nfields != 7 {
        return None;
    }
    // Field 4 runs from its start to the colon before field 5.
    let tile = &f[starts[4]..starts[5] - 1];
    std::str::from_utf8(tile).ok()?.parse::<u32>().ok()
}

/// `tile -> position -> quality-byte counts` (128 cells per position, `byte - 33`).
#[derive(Clone, Default)]
pub struct TileQual {
    pub tiles: BTreeMap<u32, Vec<[u64; QBINS]>>,
    pub reads_with_coordinates: u64,
    pub reads_without_coordinates: u64,
}

impl TileQual {
    #[inline]
    pub fn observe(&mut self, header: &[u8], qual: &[u8]) {
        match tile_of(header) {
            Some(t) => {
                self.reads_with_coordinates += 1;
                let rows = self.tiles.entry(t).or_default();
                // One resize for the whole read, not one length test and possible resize
                // per base: after this `rows[i]` is in bounds for every base of the read.
                if rows.len() < qual.len() {
                    rows.resize(qual.len(), [0; QBINS]);
                }
                for (i, &q) in qual.iter().enumerate() {
                    if let Some(b) = qbin(q) {
                        rows[i][b] += 1;
                    }
                }
            }
            // No coordinates: counted, never dropped.
            None => self.reads_without_coordinates += 1,
        }
    }

    /// Exact sum of two accumulators.
    pub fn merge(&mut self, o: &TileQual) {
        for (&t, rows) in &o.tiles {
            let dst = self.tiles.entry(t).or_default();
            if dst.len() < rows.len() {
                dst.resize(rows.len(), [0; QBINS]);
            }
            for (d, s) in dst.iter_mut().zip(rows.iter()) {
                for i in 0..QBINS {
                    d[i] += s[i];
                }
            }
        }
        self.reads_with_coordinates += o.reads_with_coordinates;
        self.reads_without_coordinates += o.reads_without_coordinates;
    }

    /// Resident bytes, for the projection printed before the numbers.
    pub fn resident_bytes(&self) -> usize {
        self.tiles
            .values()
            .map(|rows| rows.len() * std::mem::size_of::<[u64; QBINS]>())
            .sum()
    }

    /// The report block: one entry per tile with its mean quality and per-position means.
    pub fn to_json(&self, offset: u8) -> serde_json::Value {
        let tiles: Vec<serde_json::Value> = self
            .tiles
            .iter()
            .map(|(&t, rows)| {
                let mut sum = 0.0f64;
                let mut n = 0u64;
                let mut per_pos = Vec::with_capacity(rows.len());
                for row in rows {
                    let depth: u64 = row.iter().sum();
                    if depth == 0 {
                        per_pos.push(serde_json::json!(0.0));
                        continue;
                    }
                    let (mean, _q10, _q50, _q90) = pos_mean(row, offset);
                    per_pos.push(serde_json::json!((mean * 100.0).round() / 100.0));
                    sum += mean * depth as f64;
                    n += depth;
                }
                serde_json::json!({
                    "tile": t.to_string(),
                    "reads": rows.first().map(|r| r.iter().sum::<u64>()).unwrap_or(0),
                    "mean_q": if n > 0 { ((sum / n as f64) * 100.0).round() / 100.0 } else { 0.0 },
                    "per_position_mean_q": per_pos,
                })
            })
            .collect();
        serde_json::json!({
            "tiles": tiles,
            "tile_count": self.tiles.len(),
            "reads_with_coordinates": self.reads_with_coordinates,
            "reads_without_coordinates": self.reads_without_coordinates,
            "note": "reads whose header has no 7-field Illumina tile are counted in \
                     reads_without_coordinates, never dropped",
        })
    }
}

/// Mean quality of one position's exact histogram, in Phred space. The row is indexed
/// by `byte - 33`, so the byte is recovered before the offset is applied.
fn pos_mean(row: &[u64; QBINS], offset: u8) -> (f64, f64, f64, f64) {
    let depth: u64 = row.iter().sum();
    if depth == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut q = vec![0u64; QBINS];
    let mut sum = 0.0;
    for (b, &c) in row.iter().enumerate() {
        let byte = b + crate::accum::QMIN as usize;
        let phred = ((byte as i64 - offset as i64).max(0) as usize).min(QBINS - 1);
        q[phred] += c;
        sum += phred as f64 * c as f64;
    }
    (sum / depth as f64, quantile(&q, 0.10), quantile(&q, 0.50), quantile(&q, 0.90))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tile_is_field_five_of_a_seven_field_header() {
        assert_eq!(tile_of(b"@LH00916:1:SIMFLOW1:1:1101:1:1 1:N:0:ATCGATCG"), Some(1101));
        assert_eq!(tile_of(b"@LH00916:1:SIMFLOW1:1:2202:9:9 2:N:0:ATCGATCG"), Some(2202));
    }

    #[test]
    fn a_header_without_coordinates_has_no_tile() {
        assert_eq!(tile_of(b"@simple_read"), None);
        assert_eq!(tile_of(b"@a:b:c 1:N:0:x"), None, "only 3 fields");
        assert_eq!(tile_of(b"@LH00916:1:SIMFLOW1:1:notanumber:1:1 1:N:0:x"), None);
    }

    /// The index walk must accept and reject exactly what `split(':')` did: an empty
    /// leading field still counts, an eighth field still disqualifies, and an empty tile
    /// field is not a number.
    #[test]
    fn the_index_walk_agrees_with_split_on_odd_headers() {
        assert_eq!(tile_of(b"@:1:2:3:77:5:6 1:N:0:x"), Some(77), "empty first field");
        assert_eq!(tile_of(b"@a:b:c:d:1234:f:g:x 1:N:0:x"), None, "eight fields");
        assert_eq!(tile_of(b"@a:b:c:d::f:g 1:N:0:x"), None, "empty tile field");
        assert_eq!(tile_of(b":1:2:3:9:5:6"), Some(9), "no '@', no space");
        assert_eq!(tile_of(b""), None, "the empty header splits into one field");
        // Same answer as the Vec-based reference for every prefix of a real header.
        let h: &[u8] = b"@LH00916:1:SIMFLOW1:1:1101:1:1 1:N:0:ATCGATCG";
        for cut in 0..h.len() {
            // The reference: cut at the first space, split on every colon, take field 5.
            let f = &h[..cut];
            let end = f.iter().position(|&b| b == b' ').unwrap_or(f.len());
            let fields: Vec<&[u8]> = f[..end].split(|&b| b == b':').collect();
            let want = if fields.len() == 7 {
                std::str::from_utf8(fields[4]).ok().and_then(|s| s.parse::<u32>().ok())
            } else {
                None
            };
            assert_eq!(tile_of(&h[..cut]), want, "prefix of length {cut}");
        }
    }

    #[test]
    fn quality_is_separated_by_tile() {
        let mut t = TileQual::default();
        t.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"IIII");
        t.observe(b"@i:1:f:1:1101:2:2 1:N:0:x", b"IIII");
        t.observe(b"@i:1:f:1:1102:3:3 1:N:0:x", b"####");
        t.observe(b"@no_coords", b"IIII");
        assert_eq!(t.reads_with_coordinates, 3);
        assert_eq!(t.reads_without_coordinates, 1);
        assert_eq!(t.tiles.len(), 2);
        let good = &t.tiles[&1101];
        // 'I' = 73; offset 33 -> Q40, two reads x four bases.
        let (mean, ..) = pos_mean(&good[0], 33);
        assert!((mean - 40.0).abs() < 1e-9);
        let bad = &t.tiles[&1102];
        let (mean, ..) = pos_mean(&bad[0], 33);
        assert!((mean - 2.0).abs() < 1e-9, "'#' = 35 -> Q2, got {mean}");
    }

    /// A read longer than any seen before grows the rows once, and a later read reuses
    /// them: the resize is not inside the per-base loop any more, so this is the state it
    /// has to produce.
    #[test]
    fn rows_grow_once_to_the_longest_read() {
        let mut t = TileQual::default();
        t.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"II");
        assert_eq!(t.tiles[&1101].len(), 2);
        t.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"IIIIII");
        assert_eq!(t.tiles[&1101].len(), 6);
        t.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"IIII");
        assert_eq!(t.tiles[&1101].len(), 6, "a shorter read does not shrink the rows");
        assert_eq!(t.tiles[&1101][0][73 - 33], 3);
        assert_eq!(t.tiles[&1101][5][73 - 33], 1);
    }

    #[test]
    fn merge_is_exact() {
        let mut a = TileQual::default();
        a.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"IIII");
        let mut b = TileQual::default();
        b.observe(b"@i:1:f:1:1101:2:2 1:N:0:x", b"IIII");
        a.merge(&b);
        assert_eq!(a.tiles[&1101][0][73 - 33], 2);
        assert_eq!(a.reads_with_coordinates, 2);
    }
}
