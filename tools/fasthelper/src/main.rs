//! fasthelper — FastSuite orchestration and report aggregation.
//!
//! `fasthelper run -c fastsuite.toml`   run the chain, then aggregate
//! `fasthelper concatenate <dir> -o out`  concatenate existing reports
//!                                       (`report` is kept as an alias)
//! `fasthelper doctor`                  check the install
//! `fasthelper schema`                  print the report contract

use anyhow::Result;
use clap::Parser;

use fasthelper::{cli, report, runner};

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    match cli.command {
        cli::Cmd::Run(a) => runner::run(&a),
        cli::Cmd::Concatenate(a) => {
            // Two flavours of the same job: reports -> one table, or
            // per-position TSVs -> one positions x samples matrix.
            if a.positions {
                fasthelper::positions::concatenate_positions(
                    &a.inputs,
                    &a.outdir,
                    &a.metric,
                    a.all_metrics,
                    a.stdout,
                )
            } else {
                report::cmd_concatenate(&a)
            }
        }
        cli::Cmd::Compare(a) => fasthelper::compare::cmd_compare(&a.before, &a.after, &a.outdir),
        cli::Cmd::Doctor(a) => runner::doctor(&a),
        cli::Cmd::Check(a) => fasthelper::check::cmd_check(&a),
        cli::Cmd::Schema => {
            print!("{}", report::SCHEMA_DOC);
            Ok(())
        }
    }
}
