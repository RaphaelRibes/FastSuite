/// Write `val` as a LEB128 unsigned varint into `buf`.
pub fn varint_write(buf: &mut Vec<u8>, mut val: u64) {
    while val >= 0x80 {
        buf.push((val as u8 & 0x7F) | 0x80);
        val >>= 7;
    }
    buf.push(val as u8);
}

/// Read one LEB128 unsigned varint from `buf` starting at `*pos`.
/// Only needed by tests; not part of the runtime decode path (the notebook
/// decoder is in Python).
#[cfg(test)]
pub fn varint_read(buf: &[u8], pos: &mut usize) -> u64 {
    let mut result = 0u64;
    let mut shift = 0u32;
    loop {
        let b = buf[*pos];
        *pos += 1;
        result |= ((b & 0x7F) as u64) << shift;
        if b & 0x80 == 0 {
            return result;
        }
        shift += 7;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let cases = [0u64, 1, 127, 128, 255, 256, 16383, 16384, 0xFFFF_FFFF, u64::MAX];
        for &v in &cases {
            let mut buf = Vec::new();
            varint_write(&mut buf, v);
            let mut pos = 0;
            assert_eq!(varint_read(&buf, &mut pos), v, "failed for {v}");
            assert_eq!(pos, buf.len());
        }
    }

    #[test]
    fn sequence() {
        let mut buf = Vec::new();
        for &v in &[0u64, 127, 128, 12345] {
            varint_write(&mut buf, v);
        }
        let mut pos = 0;
        assert_eq!(varint_read(&buf, &mut pos), 0);
        assert_eq!(varint_read(&buf, &mut pos), 127);
        assert_eq!(varint_read(&buf, &mut pos), 128);
        assert_eq!(varint_read(&buf, &mut pos), 12345);
        assert_eq!(pos, buf.len());
    }
}