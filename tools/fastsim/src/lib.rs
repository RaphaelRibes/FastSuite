//! fastsim: NovaSeq-like paired-end read simulator with ground truth.
//!
//! Library modules (used by the `fastsim` binary and by `tests/`):
//! - [`profile`]: typed profile JSON (artefact CDFs, quality curves, adapters)
//! - [`randutil`]: CDF + joint-distribution sampling

pub mod artefact;
pub mod biology;
pub mod fastq;
pub mod fmtint;
pub mod profile;
pub mod quality;
pub mod randutil;
pub mod stats;
pub mod truth;
