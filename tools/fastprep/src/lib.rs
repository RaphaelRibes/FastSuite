//! PolyX trimming library.
//!
//! See the accompanying specification for semantics. The core entry
//! points are [`algorithm::polyx_trim`] (polyX only) and [`algorithm::trim`]
//! (adapter + polyX). Other modules build the CLI binary on top of these.

pub mod adapter;
pub mod algorithm;
pub mod cli;
pub mod correct;
pub mod io;
pub mod metadata;
pub mod pe_overlap;
pub mod pipeline;
pub mod stats;
pub mod dedup;
pub mod memory;

pub use algorithm::{polyx_trim, trim, DropReason, MatchRule, MinLen, Params, TrimOutput, TrimResult};
pub use pe_overlap::{detect_pe_overlap, detect_pe_overlap_trim};
