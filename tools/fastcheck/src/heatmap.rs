//! The fastmaps-style heatmap: each read encoded as `(x = T+C count, y = T+G
//! count)` into a `SIZE x SIZE` grid, row-major `index = y*SIZE + x`. This is
//! the same encoding fastmaps uses (`encode.rs`), so the two tools' heatmaps
//! are directly comparable.

/// Grid edge length. 64 matches fastmaps' default and keeps the JSON small
/// (~32 KB per mate).
pub const SIZE: usize = 64;

/// Encode a read as `(x, y)` where `x = #(T+C)` and `y = #(T+G)`.
/// Non-ACGT characters (e.g. N) are ignored, matching fastmaps.
pub fn encode(seq: &[u8]) -> (usize, usize) {
    let (mut x, mut y) = (0usize, 0usize);
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'T' => {
                x += 1;
                y += 1;
            }
            b'C' => x += 1,
            b'G' => y += 1,
            _ => {}
        }
    }
    (x.min(SIZE - 1), y.min(SIZE - 1))
}

/// The flat row-major index for a grid cell.
pub fn index(x: usize, y: usize) -> usize {
    y * SIZE + x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_fastmaps_encoding() {
        // fastmaps encode_read: A=(0,0) T=(1,1) G=(0,1) C=(1,0).
        assert_eq!(encode(b"AAAA"), (0, 0));
        assert_eq!(encode(b"TTTT"), (4, 4));
        assert_eq!(encode(b"ATGC"), (2, 2));
        assert_eq!(encode(b"ANNN"), (0, 0), "N is ignored");
    }

    #[test]
    fn grid_index_is_row_major() {
        assert_eq!(index(4, 4), 4 * SIZE + 4);
        assert_eq!(index(0, 1), SIZE);
        assert_eq!(index(SIZE - 1, SIZE - 1), SIZE * SIZE - 1);
    }

    #[test]
    fn counts_are_capped_at_the_grid_edge() {
        // A read longer than SIZE with all T saturates at the corner.
        let long = vec![b'T'; SIZE * 2];
        assert_eq!(encode(&long), (SIZE - 1, SIZE - 1));
    }
}
