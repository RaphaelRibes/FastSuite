//! fasthelper library surface.
//!
//! The binary (`src/main.rs`) is a thin dispatcher over these modules; they are
//! public so the report adapters can be regression-tested directly
//! (`tests/adapters.rs`).

pub mod cli;
pub mod config;
pub mod report;
pub mod runner;
