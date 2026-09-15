//! fasthelper — FastSuite orchestration and report aggregation.
//!
//! `fasthelper run -c fastsuite.toml`   run the chain, then aggregate
//! `fasthelper report <dir> -o out`     aggregate existing reports only
//! `fasthelper doctor`                  check the install
//! `fasthelper schema`                  print the report contract

use anyhow::Result;
use clap::Parser;

use fasthelper::{cli, report, runner};

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    match cli.command {
        cli::Cmd::Run(a) => runner::run(&a),
        cli::Cmd::Report(a) => report::cmd_report(&a),
        cli::Cmd::Doctor(a) => runner::doctor(&a),
        cli::Cmd::Schema => {
            print!("{}", report::SCHEMA_DOC);
            Ok(())
        }
    }
}
