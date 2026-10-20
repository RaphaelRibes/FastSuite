//! Decimal integer appending without `core::fmt`.
//!
//! The FASTQ header and the GT row are 6–12 integer fields per pair, and
//! `format!`/`writeln!` route each one through `fmt::Display` (trait object
//! dispatch + a formatter adapter) plus an allocation for the header. The hot
//! loop appends into reused buffers instead, so the only work is the digits.

/// Append `v` in decimal to `buf` (no allocation, no `fmt`).
#[inline]
pub fn push_u64(buf: &mut Vec<u8>, v: u64) {
    // u64::MAX has 20 digits.
    let mut tmp = [0u8; 20];
    let mut i = tmp.len();
    let mut n = v;
    loop {
        i -= 1;
        // `n % 10` / `n / 10` on a u64 — the compiler turns both into one
        // multiply + shift pair.
        tmp[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_the_format_macro() {
        for v in [
            0u64,
            1,
            9,
            10,
            99,
            100,
            9_999,
            10_000,
            12_345,
            123_456,
            4_294_967_295,
            18_446_744_073_709_551_615,
        ] {
            let mut buf = Vec::new();
            push_u64(&mut buf, v);
            assert_eq!(String::from_utf8(buf).unwrap(), format!("{v}"));
        }
    }

    #[test]
    fn appends_after_existing_bytes() {
        let mut buf = b"x\t".to_vec();
        push_u64(&mut buf, 300);
        push_u64(&mut buf, 0);
        assert_eq!(buf, b"x\t3000");
    }
}
