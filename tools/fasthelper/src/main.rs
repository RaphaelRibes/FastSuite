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
                // Neither flag has a meaning in this flavour, and both used to
                // be dropped without a word: --fail-on exited 0 whatever the
                // metric said, and --format produced none of the named files.
                if !a.fail_on.is_empty() {
                    anyhow::bail!(
                        "--fail-on is not implemented for --positions (it gates report metrics, \
                         and positions mode writes positions_matrix.tsv) — drop the flag or run \
                         without --positions"
                    );
                }
                if a.format != cli::DEFAULT_FORMAT {
                    anyhow::bail!(
                        "--format is not used with --positions (positions mode always writes \
                         positions_matrix.tsv) — drop --format or run without --positions"
                    );
                }
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
        cli::Cmd::Install(a) => fasthelper::install::cmd_install(&a),
        cli::Cmd::Guard(a) => fasthelper::adapter_guard::cmd_guard(&a),
        cli::Cmd::Schema => {
            print!("{}", report::SCHEMA_DOC);
            Ok(())
        }
    }
}
