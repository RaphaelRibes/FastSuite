//! fastsim binary: simulate NovaSeq-like paired-end reads with ground
//! truth, driven by a measured profile.
//!
//! Writes:
//! - `PREFIX_R1.fastq[.gz]`
//! - `PREFIX_R2.fastq[.gz]`
//! - `PREFIX.gt.tsv` (ground truth: insert, h/t per mate, full flag)

mod cli;

use std::path::Path;

use anyhow::{bail, Context, Result};
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use fastsim::artefact::{build_read, decide, PooledPair};
use fastsim::biology::{insert_from_reference, random_dna_bias, read_fasta, History};
use fastsim::fastq::FastqWriter;
use fastsim::profile::{AdapterMode, Mate, Profile};
use fastsim::quality::{qualify, MateRead};
use fastsim::truth::GtWriter;

fn main() -> Result<()> {
    let cfg = cli::parse()?;
    run(&cfg)
}

/// Ambiguous (`N`) bases in one emitted mate: what the sequencing overlay
/// (and, for the unreachable pad path, the adapter filler) put there. Counted
/// from the emitted bytes so the GT records exactly what a scorer will see.
fn n_bases(mr: &MateRead) -> u32 {
    mr.seq.iter().filter(|&&b| b == b'N').count() as u32
}

fn run(cfg: &cli::Config) -> Result<()> {
    // --force: a generator always writes, so the clobber check is the whole
    // safety net. Probe the artefacts the run will create.
    if !cfg.force {
        for suffix in [".gt.tsv", "_R1.fastq.gz", "_R2.fastq.gz", "_R1.fastq", "_R2.fastq"] {
            let probe = format!("{}{}", cfg.out, suffix);
            if std::path::Path::new(&probe).exists() {
                bail!("{probe} exists; pass --force to overwrite");
            }
        }
    }
    let t0 = std::time::Instant::now();
    // Load profile and apply CLI overrides.
    let mut profile = Profile::from_json(Path::new(&cfg.profile))?;
    if let Some(seq) = &cfg.adapter_seq {
        if seq.is_empty() {
            profile.adapter.mode = AdapterMode::None;
        } else {
            profile.adapter.seq = seq.as_bytes().to_vec();
            profile.adapter.mode = AdapterMode::Readthrough;
        }
    }
    // The profile is the authority on read length: an explicit --read-len that
    // disagrees is a user error worth refusing, but omitting it must just work,
    // so a profile goes straight in.
    if let Some(l) = cfg.read_len {
        if l != profile.read_len {
            bail!("--read-len {} does not match profile read_len {}", l, profile.read_len);
        }
    }
    // Duplicate-class rate: the profile's measured dup_rate is the authority
    // when present — the closed loop compares the sim's QC against the real
    // library's QC, and a sim pinned to --dup-rate's 0.02 default sits at
    // 0.015-0.018 whatever the real library's rate is (measured dup_ratio up to
    // 20x — the other GATE-2 cluster). --dup-rate is ONLY the fallback for a
    // profile with no measured rate: clap cannot tell an explicit flag from its
    // default, so once the profile carries a rate the flag is never consulted
    // (see `profile_dup_rate_beats_the_cli_flag` in tests/cli_test.rs).
    if profile.dup_rate.is_none() {
        profile.dup_rate = Some(cfg.dup_rate);
    }
    // Insert floor: the profile's read length is the hard ceiling (a floor
    // above it is unrepresentable), so the default is min(150, read_len) —
    // the documented `fastcheck --profile p.json` -> `fastsim --profile p.json`
    // loop must work unmodified for a 100 bp or 75 bp library. An explicit
    // --insert-min still wins, and is still validated against the read length.
    let insert_min = cfg.insert_min.unwrap_or(profile.read_len.min(150));
    if insert_min > profile.read_len {
        bail!(
            "--insert-min ({}) should be <= the profile read_len ({})",
            insert_min,
            profile.read_len
        );
    }

    // Reference contigs (optional).
    let contigs = match &cfg.reference {
        Some(p) => Some(read_fasta(Path::new(p))?),
        None => None,
    };

    // Open outputs.
    let (r1_path, r2_path) = if cfg.gzip {
        (
            format!("{}_R1.fastq.gz", cfg.out),
            format!("{}_R2.fastq.gz", cfg.out),
        )
    } else {
        (
            format!("{}_R1.fastq", cfg.out),
            format!("{}_R2.fastq", cfg.out),
        )
    };
    let gt_path = format!("{}.gt.tsv", cfg.out);

    let mut w1 = if cfg.gzip {
        FastqWriter::open_gz(Path::new(&r1_path))?
    } else {
        FastqWriter::open_plain(Path::new(&r1_path))?
    };
    let mut w2 = if cfg.single {
        None
    } else if cfg.gzip {
        Some(FastqWriter::open_gz(Path::new(&r2_path))?)
    } else {
        Some(FastqWriter::open_plain(Path::new(&r2_path))?)
    };
    let mut gt = if cfg.single {
        GtWriter::open_single(Path::new(&gt_path), profile.read_len)?
    } else {
        GtWriter::open(Path::new(&gt_path), profile.read_len)?
    };

    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed);

    // Duplicate pool: recently emitted pairs. A duplicate replays a stored
    // pair byte-for-byte — same template, artefact, errors and quality.
    let mut hist = History::new(cfg.dup_pool as usize);

    let mut n_full = 0u64;
    let mut n_adapter = 0u64;
    let mut n_dup = 0u64;
    let mut n_drop = 0u64;
    let mut n_ambig = 0u64; // ambiguous (N) bases emitted across both mates
    let mut n_ambig_reads = 0u64; // mates carrying at least one N
    for pair_idx in 0..cfg.pairs {
        // 1. Insert. Two orthogonal forced classes, then the ordinary draw:
        //   - PCR duplicate: re-sequence a stored fragment (only fresh
        //     inserts enter the pool, so a copy always has an earlier
        //     original that was itself sequenced once).
        //   - must-drop: fragment too short to survive trimming. With
        //     insert < drop_min_len the GT keep window [h, min(I, L-t)) is
        //     below min_length on BOTH mates, so a correct trimmer drops
        //     the pair. The read is still read_len bases (adapter
        //     readthrough fills the rest) — we are simulating sequencer
        //     output, not the fragment.
        let dup_rate = profile.dup_rate.unwrap_or(cfg.dup_rate);
        let is_dup = dup_rate > 0.0 && rng.gen::<f64>() < dup_rate;
        let is_drop = !is_dup && cfg.drop_rate > 0.0 && rng.gen::<f64>() < cfg.drop_rate;

        // A duplicate REPLAYS a stored pair verbatim: same bases, same artefact,
        // same substitution errors. This is what a real optical duplicate is
        // (same tile, same cycle), and it is what makes the copy detectable by
        // an exact sequence fingerprint — redrawing the artefact and the errors
        // for the copy left the measured duplicate rate at ~0.35 of the
        // requested rate, independent of pool size.
        if is_dup {
            if let Some(src) = hist.sample(&mut rng) {
                let a = src.artefact;
                w1.write_record(&w1.header(pair_idx, 1), &src.r1)?;
                if let (Some(w2), Some(m2)) = (w2.as_mut(), src.r2.as_ref()) {
                    w2.write_record(&w2.header(pair_idx, 2), m2)?;
                }
                let (d1, d2) = (n_bases(&src.r1), src.r2.as_ref().map_or(0, n_bases));
                n_ambig += (d1 + d2) as u64;
                n_ambig_reads += (d1 > 0) as u64 + (d2 > 0) as u64;
                gt.write_pair(pair_idx, src.insert_len, &a, true, d1, d2)?;
                n_dup += 1;
                continue;
            }
            // Pool not primed yet (first pairs): fall through and emit fresh.
        }

        let (insert, insert_len): (Vec<u8>, u32) = if is_drop {
            // Forced short fragment: I ~ U[1, drop_min_len). Short enough
            // that even with zero artefacts the kept window is < min_len.
            let l = rng.gen_range(1..=cfg.drop_min_len.saturating_sub(1).max(1));
            let i = match &contigs {
                Some(cs) => insert_from_reference(cs, l, &mut rng)?.1,
                None => random_dna_bias(l as usize, cfg.gc, Some(&profile.head_bias), Some(&profile.tail_bias), &mut rng),
            };
            (i, l)
        } else {
            let l = profile.sample_insert(
                cfg.insert_mean,
                cfg.insert_sd,
                insert_min,
                &mut rng,
            );
            let i = match &contigs {
                Some(cs) => insert_from_reference(cs, l, &mut rng)?.1,
                None => random_dna_bias(l as usize, cfg.gc, Some(&profile.head_bias), Some(&profile.tail_bias), &mut rng),
            };
            (i, l)
        };
        if is_dup {
            n_dup += 1;
        }
        if is_drop {
            n_drop += 1;
        }

        // 2. Artefact decision (correlated tails, heads). Drawn fresh for every
        // non-duplicate pair; a duplicate never reaches here — it replays the
        // stored artefact verbatim in the is_dup branch above.
        let mut a = decide(&profile, &mut rng);
        if a.full {
            n_full += 1;
        }

        // 3. Build reads.
        let r1 = build_read(&profile, &insert, insert_len, Mate::R1, &a, &mut rng);
        let r2 = if cfg.single {
            None
        } else {
            Some(build_read(&profile, &insert, insert_len, Mate::R2, &a, &mut rng))
        };
        // Record where the adapter actually landed, so the GT labels it explicitly.
        a.a1 = r1.adapter_len;
        a.a2 = r2.as_ref().map_or(0, |r| r.adapter_len);
        // Sequencer output is always exactly read_len bases — every branch
        // of build_read pads to L (adapter readthrough or Ns). This is the
        // contract the bench relies on, so assert it per pair.
        debug_assert_eq!(r1.seq.len() as u32, profile.read_len, "R1 length");
        if let Some(r2) = &r2 {
            debug_assert_eq!(r2.seq.len() as u32, profile.read_len, "R2 length");
        }
        if r1.adapter_len > 0 {
            n_adapter += 1;
        }

        // 4. Quality + errors.
        let m1 = qualify(&profile, Mate::R1, &r1, &mut rng);
        let m2 = r2.as_ref().map(|r| qualify(&profile, Mate::R2, r, &mut rng));

        // 5. Write.
        w1.write_record(&w1.header(pair_idx, 1), &m1)?;
        if let (Some(w2), Some(m2)) = (w2.as_mut(), m2.as_ref()) {
            w2.write_record(&w2.header(pair_idx, 2), m2)?;
        }
        let (n1, n2) = (n_bases(&m1), m2.as_ref().map_or(0, n_bases));
        n_ambig += (n1 + n2) as u64;
        n_ambig_reads += (n1 > 0) as u64 + (n2 > 0) as u64;
        gt.write_pair(pair_idx, insert_len, &a, is_dup, n1, n2)?;
        // Pool the finished pair so a later duplicate replays these exact bytes.
        // The forced-drop class is excluded: it is a contaminant model, not a
        // library member, and must not seed the duplicate pool.
        if !is_drop {
            hist.push(PooledPair {
                insert: insert.clone(),
                insert_len,
                artefact: a,
                r1: m1,
                r2: m2,
            });
        }
    }

    w1.finish().context("finish R1")?;
    if let Some(w2) = w2 {
        w2.finish().context("finish R2")?;
    }
    gt.finish().context("finish GT")?;

    if !cfg.quiet {
        eprintln!(
            "fastsim: {} pairs, read_len {}, full-artefact {} ({:.4}%), adapter-readthrough {} ({:.4}%), duplicates {} ({:.4}%), must-drop {} ({:.4}%), ambiguous {} reads ({:.4}%) / {} N bases",
            cfg.pairs,
            profile.read_len,
            n_full,
            n_full as f64 / cfg.pairs as f64 * 100.0,
            n_adapter,
            n_adapter as f64 / cfg.pairs as f64 * 100.0,
            n_dup,
            n_dup as f64 / cfg.pairs as f64 * 100.0,
            n_drop,
            n_drop as f64 / cfg.pairs as f64 * 100.0,
            n_ambig_reads,
            n_ambig_reads as f64 / (cfg.pairs as f64 * 2.0) * 100.0,
            n_ambig,
        );
    }
    // Citation request (MIT stays permissive; this is the community-norm ask).
    eprintln!(
        "If you use fastsim in published research, please cite: \
         Ribes R. FastSuite — fastsim. https://github.com/RaphaelRibes/FastSuite"
    );

    // FastSuite report (schema v1) so the synthetic control is aggregatable
    // alongside the real samples.
    if let Some(path) = &cfg.json {
        let report = serde_json::json!({
            "schema": "fastsuite.report.v1",
            "tool": "fastsim",
            "tool_version": env!("CARGO_PKG_VERSION"),
            "command": std::env::args().collect::<Vec<_>>().join(" "),
            "mode": if cfg.single { "single" } else { "paired" },
            "dry_run": false,
            "inputs": { "profile": cfg.profile },
            "outputs": {
                "r1": format!("{}_R1.fastq{}", cfg.out, if cfg.gzip { ".gz" } else { "" }),
                "r2": format!("{}_R2.fastq{}", cfg.out, if cfg.gzip { ".gz" } else { "" }),
                "gt": format!("{}.gt.tsv", cfg.out),
            },
            "elapsed_s": t0.elapsed().as_secs_f64(),
            "metrics": {
                "reads_in": 2 * cfg.pairs,
                "pairs_simulated": cfg.pairs,
                "read_len": profile.read_len,
                "seed": cfg.seed,
                "full_artefact_pairs": n_full,
                "adapter_readthrough_pairs": n_adapter,
                "duplicate_pairs": n_dup,
                "must_drop_pairs": n_drop,
                "full_artefact_rate": n_full as f64 / cfg.pairs as f64,
                "adapter_readthrough_rate": n_adapter as f64 / cfg.pairs as f64,
                "duplicate_rate_observed": n_dup as f64 / cfg.pairs as f64,
                "must_drop_rate": n_drop as f64 / cfg.pairs as f64,
                "ambiguous_bases": n_ambig,
                "reads_with_ambiguous": n_ambig_reads,
                "reads_with_ambiguous_rate": n_ambig_reads as f64 / (cfg.pairs as f64 * 2.0)
            }
        });
        std::fs::write(path, serde_json::to_string_pretty(&report)? + "\n")
            .with_context(|| format!("writing JSON report to {path}"))?;
    }
    Ok(())
}
