//! Per-tile quality (PLAN-QC module 7).
//!
//! FastQC's per-tile modules answer "is one region of the flowcell failing", which a
//! whole-file quality curve cannot see: a bad tile is invisible once its reads are averaged
//! with the good ones. The accumulator is `tile x position x quality bin`, ~5 MB per 100
//! tiles at R=150, so it is off by default (`--tiles`) and the projection is printed when
//! it is on.
//!
//! Only 7-field Illumina headers carry a tile. Reads with any other header shape are
//! counted in `reads_without_coordinates` rather than dropped: a run whose headers were
//! rewritten must not silently lose its reads.

use std::collections::BTreeMap;

use crate::accum::{quantile, QBINS};

/// Extract the tile number from a FASTQ header, if it has the 7-field Illumina form.
///
/// `@instrument:run:flowcell:lane:tile:x:y <rest>` — the tile is field 5 (index 4) of the
/// part before the first space. Anything else (a plain name, a rewritten header, a
/// different vendor's layout) returns `None`.
pub fn tile_of(header: &[u8]) -> Option<u32> {
    let end = header.iter().position(|&b| b == b' ').unwrap_or(header.len());
    let fields: Vec<&[u8]> = header[..end].split(|&b| b == b':').collect();
    if fields.len() != 7 {
        return None;
    }
    std::str::from_utf8(fields[4]).ok()?.parse::<u32>().ok()
}

/// `tile -> position -> quality-byte counts`.
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
                for (i, &q) in qual.iter().enumerate() {
                    if rows.len() <= i {
                        rows.resize(i + 1, [0; QBINS]);
                    }
                    rows[i][q as usize] += 1;
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

/// Mean quality of one position's exact histogram, in Phred space.
fn pos_mean(row: &[u64; QBINS], offset: u8) -> (f64, f64, f64, f64) {
    let depth: u64 = row.iter().sum();
    if depth == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut q = vec![0u64; QBINS];
    let mut sum = 0.0;
    for (b, &c) in row.iter().enumerate() {
        let phred = ((b as i64 - offset as i64).max(0) as usize).min(QBINS - 1);
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

    #[test]
    fn merge_is_exact() {
        let mut a = TileQual::default();
        a.observe(b"@i:1:f:1:1101:1:1 1:N:0:x", b"IIII");
        let mut b = TileQual::default();
        b.observe(b"@i:1:f:1:1101:2:2 1:N:0:x", b"IIII");
        a.merge(&b);
        assert_eq!(a.tiles[&1101][0][73], 2);
        assert_eq!(a.reads_with_coordinates, 2);
    }
}
