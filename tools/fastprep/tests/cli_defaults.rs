//! The shipped CLI defaults, asserted against the clap definition itself.
//!
//! `Params::default()` is the LIBRARY's default and deliberately differs from
//! what the binary runs with in six knobs:
//!
//! | knob                | `Params::default()` | shipped CLI   |
//! |---------------------|---------------------|---------------|
//! | `min_len`           | `0.6L`              | `0.67L`       |
//! | `min_run`           | 6                   | 12            |
//! | `head_dropout`      | `None`              | `Some(4)`     |
//! | `scan_interior`     | `true`              | `false`       |
//! | `boundary_anchored` | `false`             | `true`        |
//! | `mean_q`            | `None`              | `Some(25.0)`  |
//!
//! The library default stays the library's: it is the undecorated algorithm
//! (no quality gate, no dropout prior, no interior scan) that the 300+ unit
//! tests exercise, and it is a stable base for callers that configure their own
//! policy. The CLI default is the measured shipped configuration
//! (`docs/DECISIONS.md` D-16, `docs/TRIM-ACCURACY.md`). Flipping one to match
//! the other would silently repoint every existing test at the shipped policy,
//! so this test is the drift guard instead: it parses the CLI with NO flags and
//! pins the documented set.

use clap::Parser;
use fastprep::algorithm::MinLen;
use fastprep::cli::Cli;

#[test]
fn no_flags_parses_the_documented_shipped_defaults() {
    // `-1`/`-o` are Options checked in `finalize`, so a bare command line
    // parses: this is literally the CLI with nothing typed.
    let c = Cli::try_parse_from(["fastprep"]).expect("a bare command line parses");

    assert_eq!(c.min_len, MinLen::Frac(0.67), "--min-len");
    assert_eq!(c.mercy, 3.0, "--mercy");
    assert_eq!(c.min_run, 12, "--min-run");
    assert_eq!(c.head_dropout, Some(4), "--head-dropout");
    assert_eq!(c.mean_q, Some(25.0), "--mean-q");
    assert!(!c.no_mean_q, "--mean-q is not disabled by default");
    assert_eq!(c.adapter_error_rate, 0.2, "--adapter-error-rate");
    assert_eq!(c.min_adapter_overlap, 14, "--min-adapter-overlap");
    assert_eq!(c.phred_sat, 20, "--phred-sat");
    assert!(c.auto_phred && !c.no_auto_phred, "auto Phred census ON");
    assert!(c.dedup && !c.no_dedup, "--dedup ON");
    // Boundary anchoring is ON through the ABSENCE of the flag: `finalize`'s
    // `(None, false)` arm resolves it to true.
    assert_eq!(c.boundary_anchored, None, "boundary-anchored ON by default");
    assert!(!c.no_boundary_anchored, "--no-boundary-anchored absent");
    // The interior scan is OFF: `finalize` computes
    // `interior || !no_interior`, and the default `no_interior` is true.
    assert!(!c.interior && c.no_interior, "interior scan OFF");
    assert_eq!(c.overtrim_margin, 0, "--overtrim-margin");
    assert!(!c.no_5 && !c.no_3, "both end scans ON");
}
