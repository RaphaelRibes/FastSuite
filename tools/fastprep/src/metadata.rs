//! Binary per-read trim metadata.
//!
//! Every read is recorded in the metadata file, not just dropped or
//! trimmed reads, enabling full downstream inspection.
//!
//! # File format
//!
//! ```text
//! HEADER  4 bytes  b"PXMT"
//!         1 byte   version (0x03)
//!
//! RECORD (repeated, variable length, little-endian):
//!   u16  id_len          read-ID byte count (capped at u16::MAX)
//!   [u8] id              raw ID bytes (no null terminator)
//!   u8   mate            0 = R1 / single-end,  1 = R2
//!   u8   status          0 = kept
//!                        1 = drop:hard
//!                        2 = drop:poly5
//!                        3 = drop:poly3
//!                        4 = drop:poly_int
//!                        5 = drop:meanq
//!                        6 = drop:adapter
//!                        7 = drop:duplicate
//!   u32  t5_adapter      adapter bases trimmed from 5′ end
//!   u32  t3_adapter      adapter bases trimmed from 3′ end
//!   u32  t5_poly         polyX bases trimmed from 5′ end
//!   u32  t3_poly         polyX bases trimmed from 3′ end
//!   u32  read_len        original read length (before any trimming)
//!   u8   n_runs          number of interior-run entries that follow
//!   [n_runs ×
//!     u32  start         0-based start position in the original read
//!     u32  len           run length in bases
//!   ]
//! ```text
//!
//! Drop records can carry interior runs (e.g. `PolyInt`). `t5_adapter`/`t3_adapter` carry the
//! adapter contribution even for dropped reads.

use std::io::{self, Write};

use crate::algorithm::{DropReason, TrimResult};

pub const MAGIC: &[u8; 4] = b"PXMT";
pub const VERSION: u8 = 3;

/// Which mate this record belongs to.
#[derive(Clone, Copy)]
pub enum Mate {
    /// Single-end read, or R1 in a paired run.
    R1,
    /// R2 in a paired run.
    R2,
}

/// Write the 5-byte file header. Call once, before any records.
pub fn write_file_header(w: &mut dyn Write) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&[VERSION])
}

/// Write one metadata record for a touched read.
pub fn write_record(
    w: &mut dyn Write,
    id: &[u8],
    mate: Mate,
    res: &TrimResult,
    read_len: usize,
) -> io::Result<()> {
    let id_len = id.len().min(u16::MAX as usize) as u16;
    w.write_all(&id_len.to_le_bytes())?;
    w.write_all(&id[..id_len as usize])?;
    w.write_all(&[match mate { Mate::R1 => 0u8, Mate::R2 => 1 }])?;

    match res {
        TrimResult::Drop { reason, t5_adapter, t3_adapter, t5_poly, t3_poly, runs } => {
            let status: u8 = match reason {
                DropReason::Hard    => 1,
                DropReason::Poly5   => 2,
                DropReason::Poly3   => 3,
                DropReason::PolyInt => 4,
                DropReason::MeanQ   => 5,
                DropReason::Adapter => 6,
                DropReason::Duplicate => 7,
            };
            w.write_all(&[status])?;
            w.write_all(&(*t5_adapter as u32).to_le_bytes())?;
            w.write_all(&(*t3_adapter as u32).to_le_bytes())?;
            w.write_all(&(*t5_poly as u32).to_le_bytes())?;
            w.write_all(&(*t3_poly as u32).to_le_bytes())?;
            w.write_all(&(read_len as u32).to_le_bytes())?;
            let n = runs.len().min(u8::MAX as usize) as u8;
            w.write_all(&[n])?;
            for &(start, len) in &runs[..n as usize] {
                w.write_all(&(start as u32).to_le_bytes())?;
                w.write_all(&(len as u32).to_le_bytes())?;
            }
        }
        TrimResult::Emit { t5_adapter, t3_adapter, t5_poly, t3_poly, runs } => {
            w.write_all(&[0u8])?; // status = kept
            w.write_all(&(*t5_adapter as u32).to_le_bytes())?;
            w.write_all(&(*t3_adapter as u32).to_le_bytes())?;
            w.write_all(&(*t5_poly as u32).to_le_bytes())?;
            w.write_all(&(*t3_poly as u32).to_le_bytes())?;
            w.write_all(&(read_len as u32).to_le_bytes())?;
            let n = runs.len().min(u8::MAX as usize) as u8;
            w.write_all(&[n])?;
            for &(start, len) in &runs[..n as usize] {
                w.write_all(&(start as u32).to_le_bytes())?;
                w.write_all(&(len as u32).to_le_bytes())?;
            }
        }
    }
    Ok(())
}

/// Returns `true` for reads that should be recorded in the metadata file.
#[inline]
pub fn is_touched(res: &TrimResult) -> bool {
    match res {
        TrimResult::Drop { .. } => true,
        TrimResult::Emit { t5_adapter, t3_adapter, t5_poly, t3_poly, runs } => {
            *t5_adapter > 0 || *t3_adapter > 0 || *t5_poly > 0 || *t3_poly > 0 || !runs.is_empty()
        }
    }
}
