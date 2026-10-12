//! CLI-level tests: what a user actually gets.
//!
//! These run the built binary, so they cover argument wiring, the JSON
//! contract and the promise that nothing is written to disk except what was
//! asked for (the old fastmaps spilled 8 bytes per read to temp files).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fastcheck() -> &'static str {
    env!("CARGO_BIN_EXE_fastcheck")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fastcheck_cli_{}_{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// 6 reads: two identical, one with Ns, one with a bad-quality tail.
const FASTQ: &[u8] = b"@r1\nACGTACGTAC\n+\nIIIIIIIIII\n\
@r2\nACGTACGTAC\n+\nIIIIIIIIII\n\
@r3\nNNNNACGTAC\n+\n##########\n\
@r4\nGGGGCCCCGG\n+\n5555555555\n\
@r5\nACGTACGTAC\n+\n!!!!!!!!!!\n\
@r6\nTTTTAAAACG\n+\nJJJJJJJJJJ\n";

fn write_fixture(dir: &Path) -> PathBuf {
    let p = dir.join("reads.fastq");
    std::fs::write(&p, FASTQ).unwrap();
    p
}

fn run(dir: &Path, args: &[&str]) -> (String, String, i32) {
    let out = Command::new(fastcheck())
        .args(args)
        .current_dir(dir)
        .output()
        .expect("spawn fastcheck");
    (
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

fn read_json(p: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

#[test]
fn report_is_v1_and_matches_the_hand_counted_fixture() {
    let d = tmpdir("report");
    let fq = write_fixture(&d);
    let json = d.join("r.json");
    let (stdout, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--json", json.to_str().unwrap()],
    );
    assert_eq!(code, 0, "stderr: {stderr}");

    let v = read_json(&json);
    assert_eq!(v["schema"], "fastsuite.report.v1");
    assert_eq!(v["tool"], "fastcheck");
    let m = &v["metrics"];
    assert_eq!(m["reads"], 6);
    assert_eq!(m["bases"], 60);
    assert_eq!(m["mean_len"], 10.0);
    assert_eq!(m["min_len"], 10.0);
    assert_eq!(m["max_len"], 10.0);
    assert_eq!(m["median_len"], 10.0);
    assert_eq!(m["reads_with_non_acgt"], 1);
    // 60 bases total, 4 of them N (r3) -> 56 ACGT bases, 30 of which are G/C.
    assert_eq!(m["bases_acgt"], 56);
    assert_eq!(m["gc_fraction"], 30.0 / 56.0);
    assert_eq!(m["phred_offset"], 33);
    assert_eq!(m["phred_offset_certain"], 1, "'#' = 35 < 64 proves +33");
    assert_eq!(m["quality_byte_min"], 33, "'!' = 33");
    assert_eq!(m["quality_byte_max"], 74, "'J' = 74");
    assert_eq!(m["partial"], 0);
    assert_eq!(v["partial"]["is_partial"], false);

    // Per-file block carries the same numbers for the single file.
    assert_eq!(v["per_file"][0]["reads"], 6);
    assert_eq!(v["per_file"][0]["decompression"], "plain");

    // The summary line is one line, starts with the tool name, and reports
    // the same read count.
    let line = stdout.lines().next().unwrap();
    assert!(line.starts_with("fastcheck: "), "{line}");
    assert!(line.contains("6 reads"), "{line}");
}

#[test]
fn identical_bytes_only_in_the_input_give_identical_reports() {
    let d = tmpdir("determinism");
    let fq = write_fixture(&d);
    let a = d.join("a.json");
    let b = d.join("b.json");
    run(&d, &[fq.to_str().unwrap(), "--json", a.to_str().unwrap(), "-j", "1", "-q"]);
    run(&d, &[fq.to_str().unwrap(), "--json", b.to_str().unwrap(), "-j", "8", "-q"]);

    let mut va = read_json(&a);
    let mut vb = read_json(&b);
    // Only timing and the recorded command line may differ.
    for v in [&mut va, &mut vb] {
        v["elapsed_s"] = serde_json::json!(0.0);
        v["command"] = serde_json::json!("");
        v["metrics"]["elapsed_s"] = serde_json::json!(0.0);
    }
    assert_eq!(va, vb, "-j must not change any reported number");
}

#[test]
fn multi_file_reports_are_identical_whatever_the_thread_split() {
    let d = tmpdir("multifile");
    let a = write_fixture(&d);
    let b = d.join("second.fastq");
    std::fs::write(&b, &FASTQ[..FASTQ.len() / 2]).unwrap();

    let j1 = d.join("j1.json");
    let j8 = d.join("j8.json");
    let args1 = [a.to_str().unwrap(), b.to_str().unwrap(), "--json", j1.to_str().unwrap(), "-j", "1", "-q"];
    let args8 = [a.to_str().unwrap(), b.to_str().unwrap(), "--json", j8.to_str().unwrap(), "-j", "8", "-q"];
    run(&d, &args1);
    run(&d, &args8);

    let mut v1 = read_json(&j1);
    let mut v8 = read_json(&j8);
    for v in [&mut v1, &mut v8] {
        v["elapsed_s"] = serde_json::json!(0.0);
        v["command"] = serde_json::json!("");
        v["metrics"]["elapsed_s"] = serde_json::json!(0.0);
    }
    assert_eq!(v1, v8);
    assert_eq!(v1["metrics"]["files"], 2);
}

#[test]
fn read_budget_marks_the_report_partial() {
    let d = tmpdir("budget");
    let fq = write_fixture(&d);
    let json = d.join("r.json");
    let (_, _, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "--reads", "2", "-q"],
    );
    assert_eq!(code, 0);
    let v = read_json(&json);
    assert_eq!(v["metrics"]["reads"], 2);
    assert_eq!(v["metrics"]["partial"], 1);
    assert_eq!(v["partial"]["is_partial"], true);
    assert_eq!(v["partial"]["reason"], "reads budget");
}

#[test]
fn positions_tsv_has_one_row_per_position_and_exact_depths() {
    let d = tmpdir("positions");
    let fq = write_fixture(&d);
    let tsv = d.join("pos.tsv");
    run(&d, &[fq.to_str().unwrap(), "--positions-tsv", tsv.to_str().unwrap(), "-q"]);
    let text = std::fs::read_to_string(&tsv).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 11, "header + 10 positions");
    assert_eq!(lines[0], "pos\ta\tc\tg\tt\tn\tdepth\tmean_q\tq10\tq50\tq90\tg_frac\tpolyx3_frac\tpolyx5_frac");
    assert!(lines[1].starts_with("1\t"), "1-based positions");
    // every position is covered by all 6 reads
    for l in &lines[1..] {
        assert_eq!(l.split('\t').nth(6).unwrap(), "6", "depth in {l}");
    }
}

#[test]
fn nothing_is_written_to_disk_unless_asked() {
    let d = tmpdir("notemp");
    let fq = write_fixture(&d);
    let (_, _, code) = run(&d, &[fq.to_str().unwrap(), "-q"]);
    assert_eq!(code, 0);
    let mut left: Vec<String> = std::fs::read_dir(&d)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
        .collect();
    left.sort();
    assert_eq!(left, vec!["reads.fastq".to_string()], "no temp files, no spill");
}

#[test]
fn mixed_read_lengths_report_the_true_minimum() {
    // Reads of 4, 8 and 12 bp. `min_len` must be the first non-empty bin (4),
    // not a quantile: q=1.0 would hand back the maximum, which is a bug that
    // a uniform-length fixture cannot see.
    let d = tmpdir("mixed");
    let p = d.join("mixed.fastq");
    std::fs::write(
        &p,
        b"@a\nACGT\n+\nIIII\n@b\nACGTACGT\n+\nIIIIIIII\n@c\nACGTACGTACGT\n+\nIIIIIIIIIIII\n",
    )
    .unwrap();
    let json = d.join("m.json");
    let (_, stderr, code) = run(&d, &[p.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
    assert_eq!(code, 0, "{stderr}");
    let m = read_json(&json)["metrics"].clone();
    assert_eq!(m["min_len"], 4.0);
    assert_eq!(m["max_len"], 12.0);
    assert_eq!(m["median_len"], 8.0, "nearest rank, 3 reads -> rank 2");
    assert_eq!(m["mean_len"], 8.0);
    assert_eq!(m["reads"], 3);
    assert_eq!(m["bases"], 24);
}

#[test]
fn stdin_is_read_exactly_like_a_file() {
    // D-4: `-` must behave like the file it replaces, in every decompression
    // path, and produce the same numbers.
    let d = tmpdir("stdin");
    let fq = write_fixture(&d);
    let gz = d.join("reads.fastq.gz");
    {
        let mut enc = flate2::write::GzEncoder::new(
            std::fs::File::create(&gz).unwrap(),
            flate2::Compression::default(),
        );
        use std::io::Write;
        enc.write_all(FASTQ).unwrap();
        enc.finish().unwrap();
    }

    // plain file vs plain stdin
    let a = d.join("file.json");
    run(&d, &[fq.to_str().unwrap(), "--json", a.to_str().unwrap(), "-q"]);
    let from_pipe = Command::new(fastcheck())
        .args(["-", "--json", "-", "-q"])
        .current_dir(&d)
        .stdin(std::fs::File::open(&fq).unwrap())
        .output()
        .unwrap();
    assert!(from_pipe.status.success());
    let piped = String::from_utf8_lossy(&from_pipe.stdout);
    let file_text = std::fs::read_to_string(&a).unwrap();
    assert!(
        piped.contains("\"reads\": 6"),
        "stdin/metrics: {piped}"
    );
    // The metric block is what must not change; the file name and the
    // decompression-mode label legitimately differ for a stream.
    let metrics_of = |s: &str| -> serde_json::Value {
        let mut v: serde_json::Value = serde_json::from_str(s).unwrap();
        v["metrics"]["elapsed_s"] = serde_json::json!(0.0);
        v["metrics"].clone()
    };
    assert_eq!(
        metrics_of(&piped),
        metrics_of(&file_text),
        "stdin must not change any number"
    );
    let v: serde_json::Value = serde_json::from_str(&piped).unwrap();
    assert_eq!(v["per_file"][0]["name"], "-");
    assert_eq!(
        v["per_file"][0]["decompression"], "stdin-plain",
        "the report must say the data came from a stream"
    );

    // gzip stdin: magic-byte detection, not the filename. No `-q` here, so the
    // summary line is what lands on stdout.
    let gz_pipe = Command::new(fastcheck())
        .args(["-"])
        .current_dir(&d)
        .stdin(std::fs::File::open(&gz).unwrap())
        .output()
        .unwrap();
    assert!(gz_pipe.status.success());
    let gz_line = String::from_utf8_lossy(&gz_pipe.stdout);
    assert!(gz_line.contains("6 reads"), "gzip stdin: {gz_line:?}");
}

#[test]
fn assumptions_are_reported_as_data() {
    // D-6: an unprovable Phred offset becomes an explicit assumption entry,
    // not a silent default.
    let d = tmpdir("assumptions");
    let p = d.join("ambiguous.fastq");
    // Every quality byte >= 64, so +33 cannot be proven: 'h' = 104 (Q40 under
    // +33, Q71 under +64 — implausible, but arithmetically possible).
    std::fs::write(
        &p,
        b"@r1\nACGTACGT\n+\nhhhhhhhh\n@r2\nACGTACGT\n+\nhhhhhhhh\n",
    )
    .unwrap();
    let json = d.join("a.json");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    assert_eq!(v["metrics"]["phred_offset_certain"], 0);
    let a = &v["assumptions"][0];
    assert_eq!(a["what"], "phred_offset");
    assert!(
        a["consequence"].as_str().unwrap().contains("Phred+64"),
        "the consequence must be stated: {a}"
    );

    // A file that proves +33 has no assumptions at all.
    let d2 = tmpdir("assumptions_none");
    let fq = write_fixture(&d2);
    let json2 = d2.join("b.json");
    run(&d2, &[fq.to_str().unwrap(), "--json", json2.to_str().unwrap(), "-q"]);
    assert_eq!(read_json(&json2)["assumptions"].as_array().unwrap().len(), 0);
}

#[test]
fn fasta_input_fails_loudly() {
    let d = tmpdir("fasta");
    let p = d.join("ref.fastq");
    std::fs::write(&p, b">contig1\nACGTACGT\n").unwrap();
    let (_, stderr, code) = run(&d, &[p.to_str().unwrap(), "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("FASTA"), "{stderr}");
}


#[test]
fn rich_report_has_polyx_repeats_and_heatmap() {
    let d = tmpdir("rich");
    let fq = write_fixture(&d);
    let html = d.join("r.html");
    let json = d.join("r.json");
    let (_, _, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--html", html.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0);
    let h = std::fs::read_to_string(&html).unwrap();
    for sec in [
        "polyX run lengths",
        "Repeated sequences",
        "Heatmap",
        "tabs-polyx",
        "class=\"heatmap\"",
    ] {
        assert!(h.contains(sec), "HTML missing section: {sec}");
    }
    let j = read_json(&json);
    let m = &j["metrics"];
    for k in [
        "polyx3_max", "polyx5_max", "head_g_frac", "head_dropout",
        "repeated_seq_count", "distinct_seqs", "duplication_level",
        "kmer_top_count", "heatmap_peak",
    ] {
        assert!(m.get(k).is_some(), "missing metric {k}");
    }
    for b in ["polyx", "repeated_sequences", "kmer_spectrum", "heatmap"] {
        assert!(j.get(b).is_some(), "missing block {b}");
    }
    // The 10-copy floor means a 6-read fixture cannot trigger a repeat; the
    // over-representation logic is unit-tested in repeated.rs. But the polyX
    // measurement is exact: r4 = GGGGCCCCGG has a leading G run of 4.
    assert!(m["repeated_seq_count"].is_u64());
    assert!(m["polyx5_max"].as_u64().unwrap() >= 4);
}


#[test]
fn max_mem_flag_is_accepted_and_the_projection_is_printed_before_the_scan() {
    let d = tmpdir("maxmem");
    let fq = write_fixture(&d);
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--max-mem", "4096", "-q"],
    );
    assert_eq!(code, 0, "stderr: {stderr}");
    // The projection line must appear, and must name the fingerprint table.
    let proj = stderr
        .lines()
        .find(|l| l.contains("fingerprint table"))
        .unwrap_or_else(|| panic!("no projection line in stderr:\n{stderr}"));
    assert!(proj.contains("MB"), "{proj}");
    // The existing accumulator-state line is untouched.
    assert!(stderr.contains("accumulator state"), "{stderr}");
}

#[test]
fn processed_files_produce_a_before_after_report() {
    let d = tmpdir("cmp");
    let fq = write_fixture(&d);
    let html = d.join("c.html");
    let json = d.join("c.json");
    let (_, _, code) = run(
        &d,
        &[
            fq.to_str().unwrap(),
            "--processed-1", fq.to_str().unwrap(),
            "--html", html.to_str().unwrap(),
            "--json", json.to_str().unwrap(),
            "-q",
        ],
    );
    assert_eq!(code, 0);
    let j = read_json(&json);
    assert!(j.get("processed").is_some() && j["processed"].is_object(), "processed block missing");
    assert_eq!(j["processed"]["reads"].as_u64().unwrap(), 6);
    let h = std::fs::read_to_string(&html).unwrap();
    assert!(h.contains("Before / After"), "missing Before/After section");
    assert!(h.contains("Per-base sequence content"), "missing per-base section");
    assert!(h.contains("ppolyx"), "missing processed polyX tabs");
}
