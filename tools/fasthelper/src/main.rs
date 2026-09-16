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
        cli::Cmd::Concatenate(a) => report::cmd_concatenate(&a),
        cli::Cmd::Doctor(a) => runner::doctor(&a),
        cli::Cmd::Schema => {
            print!("{}", report::SCHEMA_DOC);
            Ok(())
        }
    }
}
