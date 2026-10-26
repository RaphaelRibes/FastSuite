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
    run_env(dir, args, &[])
}

/// `run` with extra environment variables: `$FASTSUITE_THREADS` has to be set for the
/// process under test, not for this one, or the test file would leak into every case.
fn run_env(dir: &Path, args: &[&str], envs: &[(&str, &str)]) -> (String, String, i32) {
    let mut cmd = Command::new(fastcheck());
    cmd.args(args).current_dir(dir);
    cmd.env_remove("FASTSUITE_THREADS");
    cmd.envs(envs.iter().copied());
    let out = cmd.output().expect("spawn fastcheck");
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
    // Only the run parameters may differ (timing, the recorded command line, and the
    // resolved thread count the report now states — see `blank_run_variance`).
    blank_run_variance(&mut va);
    blank_run_variance(&mut vb);
    assert_eq!(va, vb, "-j must not change any MEASURED number");
    // ... and the parameter it does change is stated, not silently different.
    assert_eq!(read_json(&b)["metrics"]["threads"], 8);
}

#[test]
fn metrics_record_the_resolved_thread_count() {
    let d = tmpdir("threads_metric");
    let fq = write_fixture(&d);
    for j in ["1", "3"] {
        let json = d.join(format!("t{j}.json"));
        let (_, stderr, code) = run(
            &d,
            &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-j", j, "-q"],
        );
        assert_eq!(code, 0, "{stderr}");
        assert_eq!(
            read_json(&json)["metrics"]["threads"],
            j.parse::<u64>().unwrap(),
            "the report must state the -j it ran with"
        );
    }
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
    blank_run_variance(&mut v1);
    blank_run_variance(&mut v8);
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
        // `busy_cores` (CPU-s / wall-s) describes the run, like `elapsed_s`: it is not a
        // function of the bytes, so two runs cannot agree on it.
        v["metrics"]["busy_cores"] = serde_json::json!(0.0);
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
    // The exact fingerprint table reports multiplicities, never names (PLAN-QC §5.9).
    assert_eq!(j["repeated_sequences"]["names_available"], false);
    assert!(j["repeated_sequences"]["top"].is_array());
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
fn a_projection_above_max_mem_refuses_before_reading() {
    let d = tmpdir("refuse");
    let fq = write_fixture(&d);
    // 1 MB ceiling against a ~6 MB projection on the 6-read fixture: it must refuse.
    let (stdout, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--max-mem", "1", "-q"],
    );
    assert_ne!(code, 0, "must refuse: {stdout}{stderr}");
    assert!(stdout.is_empty(), "no report may be written on refusal: {stdout}");
    assert!(stderr.contains("fingerprint table"), "{stderr}");
    assert!(stderr.contains("limited to"), "{stderr}");
    assert!(stderr.contains("--overrep-mode none"), "{stderr}"); // a way out
    // The refusal must happen before any numbers are printed.
    assert!(!stderr.contains("accumulator state"), "{stderr}");
}

#[test]
fn the_kmer_spectrum_is_exact_and_says_so() {
    let d = tmpdir("kmer_exact");
    // 32 bp reads: 21 windows each, so the block has something to count.
    let fq = write_overrep(&d);
    let json = d.join("k.json");
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
    assert_eq!(code, 0, "stderr: {stderr}");
    let v = read_json(&json);
    let k = &v["kmer_spectrum"];
    assert_eq!(k["approximate"], false, "the dense 4^12 table is exact");
    assert!(
        k.get("approximation").is_none(),
        "there is no bound left to state, so nothing may claim one: {k}"
    );
    // `distinct` is counted from the non-zero cells, and the packed fixture is dominated
    // by one 32 bp sequence appearing 4 times: its 21 windows are its own.
    let distinct = k["distinct"].as_u64().unwrap();
    assert!(distinct >= 21, "distinct = {distinct}, expected at least one read's windows");
    // The block is a function of the data, not of the thread split (FC-1's bonus: the
    // old counter evicted in read order, so -j could change the top list).
    let json8 = d.join("k8.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json8.to_str().unwrap(), "-j", "8", "-q"]);
    let v8 = read_json(&json8);
    assert_eq!(v8["kmer_spectrum"], *k, "-j must not move a single k-mer count");
}

/// A fixture where one sequence appears 4 times out of 12 reads.
///
/// NOTE: the plan's original fixture had seq/qual length mismatches (which fastcheck
/// correctly refuses as wrapped FASTQ) and a repeated sequence containing N (which the
/// 2-bit packer of --overrep-mode table cannot represent). Quality lengths match their
/// sequences here and the repeated sequence is pure ACGT, so every mode can name or
/// count it.
const OVERREP_FASTQ: &[u8] = b"@a\n\
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@b\n\
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@c\n\
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@d\n\
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n\
+\n\
#################################\n\
@e\n\
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGT\n\
+\n\
555555555555555555555555555555555\n\
@f\n\
ACGTACGTACGTACGTACGTACGTACGTACGT\n\
+\n\
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n\
@g\n\
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@h\n\
CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT\n\
+\n\
55555555555555555555555555555555\n\
@i\n\
AACCAACCAACCAACCAACCAACCAACCAACC\n\
+\n\
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n\
@j\n\
TGCATGCATGCATGCATGCATGCATGCATGCA\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@k\n\
GGATCCTTAGGCCAATTGGCCAATTGGCCAA\n\
+\n\
5555555555555555555555555555555\n\
@l\n\
CATCATCATCATCATCATCATCATCATCATCA\n\
+\n\
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n\
";

fn write_overrep(dir: &Path) -> std::path::PathBuf {
    let p = dir.join("over.fastq");
    std::fs::write(&p, OVERREP_FASTQ).unwrap();
    p
}

/// A sequence repeated three times whose LAST base is an N — exactly the case the 2-bit
/// packer of `--overrep-mode table` cannot represent. The old code dropped every such read
/// from the unique-sequence table, silently under-counting the library; the mode must refuse.
const OVERREP_N_FASTQ: &[u8] = b"@a\n\
AGATCGGAAGAGCACACGTCTGAACTCCAGTCN\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@b\n\
AGATCGGAAGAGCACACGTCTGAACTCCAGTCN\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@c\n\
AGATCGGAAGAGCACACGTCTGAACTCCAGTCN\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@d\n\
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGT\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
";

#[test]
fn overrep_mode_none_is_the_default_and_reports_no_names() {
    let d = tmpdir("overrep_none");
    let fq = write_overrep(&d);
    let json = d.join("r.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    assert_eq!(v["overrep"]["mode"], "none");
    assert_eq!(v["overrep"]["names_available"], false);
    // Duplication levels are still exact in this mode.
    assert!(v["metrics"]["distinct_seqs"].is_u64());
    assert_eq!(v["overrep"]["min_count"], 10, "the default 0.1%-of-reads floor");
    assert!(v["overrep"]["sequences"].as_array().unwrap().is_empty());
}

#[test]
fn overrep_mode_fingerprint_reports_exact_multiplicity_without_names() {
    let d = tmpdir("overrep_fp");
    let fq = write_overrep(&d);
    let json = d.join("r.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(),
              "--overrep-mode", "fingerprint", "-q"]);
    let v = read_json(&json);
    assert_eq!(v["overrep"]["mode"], "fingerprint");
    assert_eq!(v["overrep"]["names_available"], false);
    // 4 copies of the first sequence, threshold max(10, 0.1% of 12) = 10 -> none reported.
    assert!(v["overrep"]["sequences"].as_array().unwrap().is_empty());
    // The threshold is a flag, so it can be lowered to see the 4x sequence.
    let json2 = d.join("r2.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json2.to_str().unwrap(),
              "--overrep-mode", "fingerprint", "--overrep-min-count", "3", "-q"]);
    let v2 = read_json(&json2);
    let seqs = v2["overrep"]["sequences"].as_array().unwrap();
    assert_eq!(seqs.len(), 1, "only the 4x sequence clears a threshold of 3");
    assert_eq!(seqs[0]["count"], 4);
    assert!(seqs[0].get("seq").is_none(), "fingerprint mode cannot name the sequence");
    assert_eq!(v2["overrep"]["names_available"], false);
}

#[test]
fn overrep_mode_table_names_the_sequences_and_second_pass_records_two_passes() {
    let d = tmpdir("overrep_tbl");
    let fq = write_overrep(&d);
    let json = d.join("t.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(),
              "--overrep-mode", "table", "--overrep-min-count", "3", "-q"]);
    let v = read_json(&json);
    assert_eq!(v["overrep"]["mode"], "table");
    assert_eq!(v["overrep"]["names_available"], true);
    assert_eq!(v["overrep"]["unique_sequences_retained"], true);
    let seqs = v["overrep"]["sequences"].as_array().unwrap();
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0]["seq"], "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    assert_eq!(seqs[0]["count"], 4);

    let json2 = d.join("s.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json2.to_str().unwrap(),
              "--overrep-mode", "second-pass", "--overrep-min-count", "3", "-q"]);
    let v2 = read_json(&json2);
    assert_eq!(v2["overrep"]["mode"], "second-pass");
    assert_eq!(v2["overrep"]["names_available"], true);
    assert_eq!(v2["overrep"]["passes"], 2, "the second read must be declared");
    let seqs2 = v2["overrep"]["sequences"].as_array().unwrap();
    assert_eq!(seqs2.len(), 1);
    assert_eq!(seqs2[0]["seq"], "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
}

#[test]
fn overrep_second_pass_is_thread_count_invariant() {
    let d = tmpdir("overrep_threads");
    // BGZF: the one input whose second pass actually builds an inflate pool, so the
    // resolved -j reaches `overrep::second_pass` -> `reader::open` (FC-12).
    let plain = write_overrep(&d);
    let mut src = std::fs::File::open(&plain).unwrap();
    let bgzf = d.join("over.fastq.gz");
    let mut w = gzp::par::compress::ParCompressBuilder::<gzp::deflate::Bgzf>::new()
        .num_threads(2)
        .unwrap()
        .from_writer(std::fs::File::create(&bgzf).unwrap());
    std::io::copy(&mut src, &mut w).unwrap();
    drop(w);

    let mut reports = Vec::new();
    // -j 2 and -j 8 both inflate BGZF in parallel, so the reports must be identical;
    // -j 1 would legitimately report "bgzf-serial" (a different decompression MODE, not
    // a different number) and is covered by the plain-file identity tests.
    for j in ["2", "8"] {
        let json = d.join(format!("s{j}.json"));
        let (_, stderr, code) = run(
            &d,
            &[bgzf.to_str().unwrap(), "--json", json.to_str().unwrap(),
              "--overrep-mode", "second-pass", "--overrep-min-count", "3", "-j", j, "-q"],
        );
        assert_eq!(code, 0, "-j {j}: {stderr}");
        let mut v = read_json(&json);
        assert_eq!(v["overrep"]["passes"], 2, "the second read must still be declared");
        blank_run_variance(&mut v);
        reports.push(v);
    }
    if reports[0] != reports[1] {
        let a = reports[0].as_object().unwrap();
        let b = reports[1].as_object().unwrap();
        for (k, v) in a {
            if b.get(k) != Some(v) {
                panic!("-j changed {k}:\n 1: {v}\n 8: {}", b.get(k).unwrap());
            }
        }
        panic!("-j changed the report in a way the top-level keys do not name");
    }
    // The pass still finds what it found with one thread.
    let seqs = reports[0]["overrep"]["sequences"].as_array().unwrap();
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0]["count"], 4);
}

#[test]
fn overrep_mode_table_refuses_above_max_mem_like_the_fingerprint_table() {
    let d = tmpdir("overrep_refuse");
    let fq = write_overrep(&d);
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--overrep-mode", "table", "--max-mem", "1", "-q"],
    );
    assert_ne!(code, 0, "must refuse: {stderr}");
    assert!(
        stderr.contains("unique sequences") || stderr.contains("fingerprint table"),
        "the refusal must name the table: {stderr}"
    );
    assert!(stderr.contains("--overrep-mode fingerprint"), "must offer the cheaper mode: {stderr}");
}

/// Under `--overrep-mode table` every sequence is packed 2 bits per base, so an N cannot be
/// represented. The old code dropped every such read from the unique-sequence table, silently
/// under-counting exactly the libraries the mode exists to name — this must be a refusal that
/// names the record and the offending base's position, not a silent gap.
#[test]
fn overrep_mode_table_refuses_a_read_it_cannot_pack_instead_of_dropping_it() {
    let d = tmpdir("overrep_table_n");
    let fq = d.join("n.fastq");
    std::fs::write(&fq, OVERREP_N_FASTQ).unwrap();

    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--overrep-mode", "table", "--overrep-min-count", "3", "-q"],
    );
    assert_ne!(code, 0, "table mode must refuse a sequence it cannot pack: {stderr}");
    assert!(stderr.contains("cannot pack"), "the refusal must say why: {stderr}");
    assert!(stderr.contains("@a"), "and name the offending read: {stderr}");
    assert!(stderr.contains("position 32"), "and the offending base's position: {stderr}");

    // The modes that do not pack read the same file without complaint.
    let json = d.join("fp.json");
    let (_, stderr2, code2) = run(
        &d,
        &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(),
          "--overrep-mode", "fingerprint", "--overrep-min-count", "3", "-q"],
    );
    assert_eq!(code2, 0, "fingerprint mode must accept N-bearing reads: {stderr2}");
    assert_eq!(read_json(&json)["overrep"]["mode"], "fingerprint");
}

#[test]
fn overrep_mode_is_rejected_when_unknown() {
    let d = tmpdir("overrep_bad");
    let fq = write_fixture(&d);
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--overrep-mode", "sampled", "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("none") && stderr.contains("second-pass"), "{stderr}");
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

/// 6 reads of 150 bp: two carry the full 19-base Nextera adapter starting at 1-based
/// position 121, one carries a truncated 16-base Nextera tail at 135, three are clean.
const ADAPTER_FASTQ: &[u8] = b"@a\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTGTCTCTTATACACATCTAAAAAAAAAAA\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@b\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTGTCTCTTATACACATCTAAAAAAAAAAA\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@c\n\
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTGTCTCTTATACACA\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@d\n\
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n\
+\n\
555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555\n\
@e\n\
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n\
+\n\
555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555\n\
@f\n\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n\
+\n\
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n\
";

#[test]
fn adapter_content_per_position_counts_every_match_start() {
    let d = tmpdir("adapter_pos");
    let p = d.join("ad.fastq");
    std::fs::write(&p, ADAPTER_FASTQ).unwrap();
    let json = d.join("a.json");
    let tsv = d.join("pos.tsv");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--adapter-content", "--adapter-preset", "nextera",
          "--json", json.to_str().unwrap(), "--positions-tsv", tsv.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let blk = &v["adapter_content"];
    assert_eq!(blk["adapters"].as_array().unwrap().len(), 1, "nextera preset = one adapter");
    assert_eq!(blk["adapters"][0]["name"], "Nextera");
    // Per-position counts: 2 full adapters at position 121 (0-based 120), one truncated
    // 16-base match at position 135 (0-based 134).
    let per_pos = blk["adapters"][0]["per_position"].as_object().unwrap();
    assert_eq!(per_pos["121"], 2, "two reads carry the full adapter at 121");
    assert_eq!(per_pos["135"], 1, "one read carries a 16-base truncated adapter at 135");
    // Totals: 3 of 6 reads carry an adapter somewhere.
    assert_eq!(blk["adapters"][0]["reads_with_match"], 3);
    assert_eq!(blk["error_rate"], 0.2);
    assert_eq!(blk["min_overlap"], 14);

    // The TSV gains one column per adapter.
    let text = std::fs::read_to_string(&tsv).unwrap();
    let header = text.lines().next().unwrap();
    assert!(header.contains("adapter_Nextera"), "header: {header}");
    let rows: Vec<&str> = text.lines().skip(1).collect();
    let idx = header.split('\t').position(|c| c == "adapter_Nextera").unwrap();
    assert_eq!(rows[120].split('\t').nth(idx).unwrap(), "2");
    assert_eq!(rows[134].split('\t').nth(idx).unwrap(), "1");
    assert_eq!(rows[0].split('\t').nth(idx).unwrap(), "0", "clean at position 1");
}

#[test]
fn adapter_content_is_off_by_default_and_the_columns_are_absent() {
    let d = tmpdir("adapter_off");
    let fq = write_fixture(&d);
    let json = d.join("a.json");
    let tsv = d.join("pos.tsv");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(),
              "--positions-tsv", tsv.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    assert!(v.get("adapter_content").is_none(), "must be opt-in");
    let header = std::fs::read_to_string(&tsv).unwrap().lines().next().unwrap().to_string();
    assert!(!header.contains("adapter_"), "no adapter columns by default: {header}");
}

#[test]
fn adapter_content_honours_the_error_rate_and_overlap_flags() {
    let d = tmpdir("adapter_flags");
    let p = d.join("ad.fastq");
    std::fs::write(&p, ADAPTER_FASTQ).unwrap();
    let json = d.join("a.json");
    // An 18-base floor rejects the 16-base truncated match but keeps the full 19-base
    // ones (the plan's 30-base floor would reject both: Nextera is 19 bases long).
    run(&d, &[p.to_str().unwrap(), "--adapter-content", "--adapter-preset", "nextera",
              "--min-adapter-overlap", "18", "--json", json.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    let per_pos = v["adapter_content"]["adapters"][0]["per_position"].as_object().unwrap();
    assert_eq!(per_pos.get("135"), None, "16 < 18: below the floor");
    assert_eq!(per_pos["121"], 2);
    assert_eq!(v["adapter_content"]["min_overlap"], 18);
}

// ---------------------------------------------------------------------------
// Task 12: --only / --no-quality / --no-heatmap module gating.
//
// The gating contract has three parts, and each is tested below: (1) with no
// gating flag the report is exactly what it was before --only existed; (2) a
// selected module's fields are present and the others' are absent, never zero;
// (3) the work those modules own is not performed at all.
// ---------------------------------------------------------------------------

/// Every module that a bare run measures. `adapter` is excluded: it is opt-in
/// behind --adapter-content, and stays so under --only.
const EVERY_MODULE: &str = "bases,quality,length,gc,heatmap,dup,polyx,kmer";

/// The blocks and metrics a run with every module measures. Listing them here (rather
/// than anywhere near the report writer) keeps the "default is unchanged" test honest:
/// if the default report ever stops carrying one of them, this fails.
const ALL_BLOCKS: &[&str] = &[
    "base_composition",
    "base_position",
    "read_length_hist",
    "read_mean_q_hist",
    "polyx",
    "repeated_sequences",
    "overrep",
    "kmer_spectrum",
    "heatmap",
];

const ALL_METRICS: &[&str] = &[
    "min_len", "max_len", "median_len", "bases_acgt", "bases_non_acgt", "gc_fraction",
    "mean_read_mean_q", "reads_q20_frac", "reads_q30_frac", "phred_offset_certain",
    "quality_byte_min", "quality_byte_max", "head_g_frac", "head_dropout", "polyx3_max",
    "polyx5_max", "distinct_seqs", "duplication_level", "kmer_top_count", "heatmap_peak",
];

/// Blank the fields that are run parameters rather than measurements: the elapsed time,
/// the recorded command line, and the resolved thread count.
///
/// REVERSAL, named deliberately: this helper used not to blank `metrics.threads` because
/// the report did not carry one, and the two identity tests below asserted that `-j`
/// changes no reported number at all. The report now STATES the -j it ran with (a run
/// parameter a consumer needs, and the one thing that legitimately differs between two
/// runs of the same bytes), so the assertion is narrowed to what it was really about:
/// every MEASURED number is thread-count invariant. `metrics_record_the_resolved_thread_count`
/// pins the field itself.
fn blank_run_variance(v: &mut serde_json::Value) {
    v["elapsed_s"] = serde_json::json!(0.0);
    v["command"] = serde_json::json!("");
    v["metrics"]["elapsed_s"] = serde_json::json!(0.0);
    v["metrics"]["threads"] = serde_json::json!(0);
    // `busy_cores` is a measurement of the RUN (CPU-s / wall-s), not of the input: two runs
    // of the same bytes on a different machine, or at a different -j, cannot agree on it.
    v["metrics"]["busy_cores"] = serde_json::json!(0.0);
}

#[test]
fn without_a_gating_flag_the_report_is_todays_report() {
    let d = tmpdir("only_absent");
    let fq = write_fixture(&d);
    let json = d.join("r.json");
    run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
    let v = read_json(&json);
    assert!(v.get("modules").is_none(), "no --only flag, no modules array");
    for b in ALL_BLOCKS {
        assert!(v.get(b).is_some(), "default report lost block {b}");
    }
    for k in ALL_METRICS {
        assert!(v["metrics"].get(k).is_some(), "default report lost metric {k}");
    }
    assert!(v["metrics"]["phred_offset_certain"].is_u64());
    assert!(v["assumptions"].as_array().unwrap().is_empty());

    // And the whole-module path measures exactly the default path: listing every module
    // changes nothing but the declaration that it was listed.
    let listed = d.join("l.json");
    run(&d, &[fq.to_str().unwrap(), "--only", EVERY_MODULE, "--json", listed.to_str().unwrap(), "-q"]);
    let mut vl = read_json(&listed);
    assert_eq!(
        vl["modules"].as_array().unwrap().len(),
        8,
        "the array must name the modules that ran"
    );
    vl.as_object_mut().unwrap().remove("modules");
    let mut vd = v.clone();
    blank_run_variance(&mut vl);
    blank_run_variance(&mut vd);
    assert_eq!(vl, vd, "--only with every module must measure today's report exactly");
}

#[test]
fn only_dup_runs_the_duplication_module_and_nothing_else() {
    let d = tmpdir("only_dup");
    let fq = write_fixture(&d);
    let json = d.join("d.json");
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--only", "dup", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let mods = v["modules"].as_array().unwrap();
    assert!(mods.iter().any(|m| m == "dup"), "dup must be enabled: {mods:?}");
    assert!(!mods.iter().any(|m| m == "quality"), "{mods:?}");
    assert!(!mods.iter().any(|m| m == "heatmap"), "{mods:?}");
    // Duplication is still exact, and the gated-out modules are absent, not zero.
    assert!(v["metrics"]["distinct_seqs"].is_u64());
    assert!(v.get("heatmap").is_none(), "the heatmap must not be built");
    assert!(v.get("base_position").is_none(), "no per-position matrices were materialised");
    assert!(v.get("polyx").is_none());
    assert!(v.get("kmer_spectrum").is_none());
    assert!(v.get("read_length_hist").is_none());
    // The metrics that belong to a gated-out module are absent rather than zero, which
    // is what distinguishes "not measured" from "measured and empty".
    for k in ["mean_read_mean_q", "gc_fraction", "min_len", "heatmap_peak", "polyx3_max", "kmer_top_count"] {
        assert!(v["metrics"].get(k).is_none(), "gated-out metric {k} must be absent");
    }
    // Reads and bases are always counted: a run that reports no reads is not a run.
    assert_eq!(v["metrics"]["reads"], 6);
    assert_eq!(v["metrics"]["bases"], 60);
    assert_eq!(v["metrics"]["phred_offset"], 33, "a stated default, not a measurement");
    assert!(v["assumptions"].as_array().unwrap().iter().any(|a| a["what"] == "quality"));
}

#[test]
fn only_bases_keeps_the_composition_and_drops_the_quality() {
    let d = tmpdir("only_bases");
    let fq = write_fixture(&d);
    let json = d.join("b.json");
    let tsv = d.join("b.tsv");
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--only", "bases", "--json", json.to_str().unwrap(),
          "--positions-tsv", tsv.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    assert_eq!(v["modules"].as_array().unwrap(), &vec![serde_json::json!("bases")]);
    // The composition is measured, with one row per position even though the quality
    // rows (which normally define the coverage) were never grown.
    assert_eq!(v["base_position"].as_object().unwrap().len(), 10, "10 positions");
    assert_eq!(v["base_position"]["1"][0], 3, "3 reads start with A");
    // The composition block is summed from the every-base totals, which the gated path
    // fills exactly as the ungated one does: a gated run must not report zeros here.
    let comp = &v["base_composition"];
    let comp_sum: u64 = ["a", "c", "g", "t", "other"]
        .iter()
        .map(|k| comp[*k].as_u64().unwrap())
        .sum();
    assert_eq!(comp_sum, v["metrics"]["bases"].as_u64().unwrap(), "{comp}");
    assert_eq!(comp["a"], 15, "hand count from the fixture: {comp}");
    assert!(v["metrics"].get("head_g_frac").is_some());
    assert!(v.get("read_mean_q_hist").is_none());
    assert!(v["metrics"].get("quality_byte_min").is_none());
    assert!(v["metrics"].get("phred_offset_certain").is_none());

    // The per-position TSV carries only the columns whose module ran — the same rule
    // the adapter columns already follow.
    let text = std::fs::read_to_string(&tsv).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "pos\ta\tc\tg\tt\tn\tdepth\tg_frac", "only the measured columns");
    assert_eq!(lines.len(), 11, "header + 10 positions");
    for l in &lines[1..] {
        assert_eq!(l.split('\t').count(), 8, "{l}");
        assert_eq!(l.split('\t').nth(6).unwrap(), "6", "depth in {l}");
    }
}

#[test]
fn no_quality_skips_the_quality_accumulators_but_keeps_the_framing() {
    let d = tmpdir("no_quality");
    let fq = write_fixture(&d);
    let json = d.join("q.json");
    let (_, _, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--no-quality", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0);
    let v = read_json(&json);
    // Length and composition are unaffected; quality is absent rather than zero.
    assert_eq!(v["metrics"]["reads"], 6);
    assert_eq!(v["metrics"]["bases"], 60);
    assert!(v["metrics"].get("mean_read_mean_q").is_none());
    assert!(v.get("base_position").is_some(), "composition is not a quality module");
    assert_eq!(v["metrics"]["phred_offset"], 33, "a stated default, not a measurement");
    // The report must say the offset was not measured.
    assert!(v["assumptions"].as_array().unwrap().iter().any(|a| a["what"] == "quality"));
    // The quality line is still read for framing, so a record with a short quality
    // string is still an error rather than silent nonsense.
    let bad = d.join("bad.fastq");
    std::fs::write(&bad, b"@r1\nACGTACGT\n+\nIIII\n").unwrap();
    let (_, stderr, code) = run(&d, &[bad.to_str().unwrap(), "--no-quality", "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("quality length"), "{stderr}");
}

#[test]
fn no_heatmap_never_builds_the_grid() {
    let d = tmpdir("no_heatmap");
    let fq = write_fixture(&d);
    let json = d.join("h.json");
    let (_, _, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--no-heatmap", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0);
    let v = read_json(&json);
    assert!(v.get("heatmap").is_none());
    assert!(v["metrics"].get("heatmap_peak").is_none());
    assert_eq!(v["metrics"]["reads"], 6, "everything else still runs");
    assert!(v.get("polyx").is_some(), "one module off does not gate the others");
}

#[test]
fn only_rejects_an_unknown_or_empty_module_list() {
    let d = tmpdir("only_bad");
    let fq = write_fixture(&d);
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--only", "bogus", "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("bogus"), "{stderr}");
    assert!(stderr.contains("dup") && stderr.contains("adapter"), "{stderr}");

    // An empty list is refused, not read as "no modules, empty report".
    for list in ["", ",", " , "] {
        let (stdout, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--only", list, "-q"]);
        assert_ne!(code, 0, "--only {list:?} must be refused");
        assert!(stdout.is_empty(), "no report may be printed on refusal: {stdout}");
        assert!(stderr.contains("no module named"), "{list:?}: {stderr}");
    }

    // A selection that --no-quality empties is refused too...
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--only", "quality", "--no-quality", "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("nothing left to run"), "{stderr}");
    // ... as is --adapter-content left out of an --only list, rather than silently ignored.
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--only", "dup", "--adapter-content", "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("--adapter-content") && stderr.contains("adapter"), "{stderr}");
}

#[test]
fn gating_skips_the_work_it_disables() {
    let d = tmpdir("only_work");
    let fq = write_fixture(&d);

    // No gating flag: the counters are not even printed, because nothing was skipped.
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "-q"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(!stderr.contains("measurement work"), "opt-in only: {stderr}");
    assert!(stderr.contains("duplication fingerprint table projected"), "{stderr}");

    // --only dup fingerprints every read and never touches a position row, the quality
    // bytes, the grid or the k-mer counter. These are the operations themselves, counted
    // where they happen: a zero here is a loop that did not run.
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--only", "dup", "-q"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(stderr.contains("modules: dup"), "{stderr}");
    assert!(
        stderr.contains(
            "6 reads measured: 0 base cells, 0 quality cells, 0 quality bytes, 6 fingerprints, \
             0 polyX reads, 0 heatmap cells, 0 12-mer windows"
        ),
        "{stderr}"
    );

    // --only bases measures the per-position composition, and nothing else — and with
    // no dup module there is no fingerprint table at all: not projected, not built.
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--only", "bases", "-q"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        stderr.contains(
            "6 reads measured: 60 base cells, 0 quality cells, 0 quality bytes, 0 fingerprints"
        ),
        "{stderr}"
    );
    assert!(stderr.contains("the dup module is off"), "{stderr}");
    assert!(!stderr.contains("fingerprint table projected"), "no table is projected: {stderr}");

    // --no-quality keeps every other measurement and drops exactly the quality work.
    let (_, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--no-quality", "-q"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        stderr.contains("60 base cells, 0 quality cells, 0 quality bytes, 6 fingerprints"),
        "{stderr}"
    );
}

#[test]
fn gating_survives_the_multi_file_thread_split() {
    let d = tmpdir("only_threads");
    let a = write_fixture(&d);
    let b = d.join("second.fastq");
    std::fs::write(&b, &FASTQ[..FASTQ.len() / 2]).unwrap();
    let j1 = d.join("g1.json");
    let j8 = d.join("g8.json");
    run(&d, &[a.to_str().unwrap(), b.to_str().unwrap(), "--only", "length,gc",
              "--json", j1.to_str().unwrap(), "-j", "1", "-q"]);
    run(&d, &[a.to_str().unwrap(), b.to_str().unwrap(), "--only", "length,gc",
              "--json", j8.to_str().unwrap(), "-j", "8", "-q"]);
    let mut v1 = read_json(&j1);
    let mut v8 = read_json(&j8);
    blank_run_variance(&mut v1);
    blank_run_variance(&mut v8);
    assert_eq!(v1, v8, "-j must not change a gated report either");
    assert_eq!(v8["modules"].as_array().unwrap().len(), 2);
    assert!(v8["metrics"]["min_len"].is_number(), "length ran");
    assert!(v8["metrics"]["gc_fraction"].is_number(), "gc ran");
    assert!(v8["metrics"].get("mean_read_mean_q").is_none(), "quality did not");
    assert!(v8["metrics"].get("distinct_seqs").is_none(), "dup did not");
}

#[test]
fn gating_refuses_an_output_whose_modules_are_off() {
    let d = tmpdir("only_outputs");
    let fq = write_fixture(&d);
    let fq = fq.to_str().unwrap();

    // A rendered page, a fitted profile and a per-position table cannot say "not
    // measured" field by field, so they are refused while it is still cheap.
    let html = d.join("h.html");
    let (_, stderr, code) = run(&d, &[fq, "--only", "dup", "--html", html.to_str().unwrap(), "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("--html") && stderr.contains("quality"), "{stderr}");
    assert!(!html.exists(), "nothing may be rendered from empty accumulators");

    let prof = d.join("p.json");
    let (_, stderr, code) = run(&d, &[fq, "--only", "dup", "--profile", prof.to_str().unwrap(), "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("--profile") && stderr.contains("quality"), "{stderr}");

    let tsv = d.join("x.tsv");
    let (_, stderr, code) = run(&d, &[fq, "--only", "length", "--positions-tsv", tsv.to_str().unwrap(), "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("--positions-tsv"), "{stderr}");

    // --overrep-mode builds on the fingerprint table, so it needs the dup module.
    let (_, stderr, code) = run(&d, &[fq, "--only", "length", "--overrep-mode", "table", "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("dup"), "{stderr}");

    // The refusals are specific, and they name what is missing: --no-heatmap costs the
    // heatmap section of the page, so the page is refused and says why...
    let (_, stderr, code) = run(&d, &[fq, "--no-heatmap", "--html", html.to_str().unwrap(), "-q"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("heatmap"), "{stderr}");
    assert!(!html.exists());

    // ... while a selection that keeps every module the page renders is allowed, and the
    // page says which modules produced it.
    let all = d.join("all.html");
    let (_, stderr, code) = run(&d, &[fq, "--only", EVERY_MODULE, "--html", all.to_str().unwrap(), "-q"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(all.exists());
}

// ---------------------------------------------------------------------------
// Task 11: --tiles per-tile quality, off by default.
//
// The tile module rides on its own opt-in flag: off (the default) nothing is
// computed and the report is byte for byte what it was; on, quality is separated
// by the tile parsed from the 7-field Illumina header, and headers without
// coordinates are counted, never dropped.
// ---------------------------------------------------------------------------

/// 7-field Illumina headers: `@instrument:run:flowcell:lane:tile:x:y`.
/// Tile 1101 gets good quality, tile 1102 gets bad, on the same read content.
const TILE_FASTQ: &[u8] = b"@LH00916:1:SIMFLOW1:1:1101:1:1 1:N:0:ATCGATCG\nACGTACGTAC\n+\nIIIIIIIIII\n\
@LH00916:1:SIMFLOW1:1:1101:2:2 1:N:0:ATCGATCG\nACGTACGTAC\n+\nIIIIIIIIII\n\
@LH00916:1:SIMFLOW1:1:1102:3:3 1:N:0:ATCGATCG\nACGTACGTAC\n+\n#####AAAAA\n\
@LH00916:1:SIMFLOW1:1:1102:4:4 1:N:0:ATCGATCG\nACGTACGTAC\n+\n#####AAAAA\n\
@simple_read_with_no_coords\nACGTACGTAC\n+\nIIIIIIIIII\n";

#[test]
fn tiles_module_separates_quality_by_tile_and_is_off_by_default() {
    let d = tmpdir("tiles");
    let p = d.join("t.fastq");
    std::fs::write(&p, TILE_FASTQ).unwrap();

    // Off by default: no block, and headers without coords are simply not counted.
    let j0 = d.join("off.json");
    run(&d, &[p.to_str().unwrap(), "--json", j0.to_str().unwrap(), "-q"]);
    assert!(read_json(&j0).get("per_tile").is_none(), "must be opt-in");

    let json = d.join("on.json");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--tiles", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let blk = &v["per_tile"];
    assert_eq!(blk["tiles"].as_array().unwrap().len(), 2, "1101 and 1102");
    assert_eq!(blk["reads_with_coordinates"], 4);
    assert_eq!(blk["reads_without_coordinates"], 1, "the simple header is counted, not dropped");
    // The projection line names the cost.
    assert!(stderr.contains("per-tile"), "{stderr}");

    let by_tile: std::collections::HashMap<String, &serde_json::Value> = blk["tiles"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| (t["tile"].as_str().unwrap().to_string(), t))
        .collect();
    // Tile 1101: every base is 'I' = Q40. Tile 1102: five '#' (Q2) then five 'A' (Q32).
    let good = &by_tile["1101"];
    assert_eq!(good["mean_q"], 40.0);
    let bad = &by_tile["1102"];
    // (5*2 + 5*32) / 10 = 17
    assert!((bad["mean_q"].as_f64().unwrap() - 17.0).abs() < 1e-9, "{bad}");
    // Per-position detail exists for both.
    assert!(good["per_position_mean_q"].as_array().unwrap().len() == 10);
}

#[test]
fn tiles_module_reports_the_memory_projection_when_on() {
    let d = tmpdir("tiles_mem");
    let p = d.join("t.fastq");
    std::fs::write(&p, TILE_FASTQ).unwrap();
    let (_, stderr, code) = run(&d, &[p.to_str().unwrap(), "--tiles", "-q"]);
    assert_eq!(code, 0);
    assert!(stderr.contains("per-tile quality"), "{stderr}");
    assert!(stderr.contains("MB"), "{stderr}");
}

/// Deterministic pseudo-random DNA (the generator insert.rs's tests use), so a
/// fragment cannot accidentally be low-complexity.
fn dna(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed | 1;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Paired FASTQ: `n` pairs of 50 bp reads off a 70 bp fragment, so every pair's
/// insert is 70 and recoverable from the mates' 30-base overlap.
fn write_pairs(dir: &Path, n: usize) -> (PathBuf, PathBuf) {
    let (mut a, mut b) = (String::new(), String::new());
    for i in 0..n {
        let frag = dna(70, 1000 + i as u64);
        let r1 = frag[..50].to_vec();
        let r2 = revcomp(&frag[20..]);
        a.push_str(&format!(
            "@p{i}/1\n{}\n+\n{}\n",
            String::from_utf8(r1).unwrap(),
            "I".repeat(50)
        ));
        b.push_str(&format!(
            "@p{i}/2\n{}\n+\n{}\n",
            String::from_utf8(r2).unwrap(),
            "I".repeat(50)
        ));
    }
    let p1 = dir.join("r1.fastq");
    let p2 = dir.join("r2.fastq");
    std::fs::write(&p1, a).unwrap();
    std::fs::write(&p2, b).unwrap();
    (p1, p2)
}

/// The `<td>` cell of a key/value row (`<tr><th>label</th><td>…</td></tr>`).
fn kv_row<'a>(html: &'a str, label: &str) -> &'a str {
    let marker = format!("<tr><th>{label}</th><td>");
    let start = html
        .find(&marker)
        .unwrap_or_else(|| panic!("no {label:?} row in the page"))
        + marker.len();
    let rest = &html[start..];
    &rest[..rest.find("</td>").unwrap()]
}

#[test]
fn the_html_insert_table_divides_by_the_pairs_it_actually_sampled() {
    let d = tmpdir("insert_html");
    let (r1, r2) = write_pairs(&d, 200);
    let html = d.join("i.html");
    let json = d.join("i.json");
    let (_, stderr, code) = run(
        &d,
        &[
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--insert-hist",
            "--json",
            json.to_str().unwrap(),
            "--html",
            html.to_str().unwrap(),
            "-q",
        ],
    );
    assert_eq!(code, 0, "{stderr}");
    let j = read_json(&json);
    let sampled = j["insert_hist"]["pairs_sampled"].as_u64().unwrap();
    assert_eq!(sampled, 200, "all 200 pairs were read");
    assert!(
        j["insert_hist"]["measured"].as_u64().unwrap() >= 190,
        "the fixture must actually measure inserts: {}",
        j["insert_hist"]["measured"]
    );

    let h = std::fs::read_to_string(&html).unwrap();
    // The denominator is the sample, not the --insert-pairs cap (20000 by
    // default). Dividing by the cap printed "pairs sampled 20000" for a 200-pair
    // run and made the four outcome buckets sum to ~1% instead of ~100%.
    assert_eq!(
        kv_row(&h, "pairs sampled"),
        "200",
        "the page must state the real sample, not the cap"
    );
    let mut sum = 0.0;
    for label in [
        "measured from overlap",
        "measured from the adapter",
        "no overlap (insert above the resolution)",
        "complete overlap, no adapter found",
    ] {
        let cell = kv_row(&h, label);
        let pct: f64 = cell
            .rsplit('(')
            .next()
            .unwrap()
            .trim_end_matches("%)")
            .parse()
            .unwrap();
        sum += pct;
    }
    assert!(
        (sum - 100.0).abs() < 0.5,
        "the outcome buckets partition the sample, so they must sum to ~100%, got {sum}%"
    );
}

#[test]
fn a_partial_run_is_never_called_exact_in_the_html() {
    let d = tmpdir("partial_html");
    let fq = write_fixture(&d);
    let html = d.join("p.html");
    let json = d.join("p.json");
    let (_, stderr, code) = run(
        &d,
        &[
            fq.to_str().unwrap(),
            "--reads",
            "2",
            "--html",
            html.to_str().unwrap(),
            "--json",
            json.to_str().unwrap(),
            "-q",
        ],
    );
    assert_eq!(code, 0, "{stderr}");
    // The JSON already says so; the page next to it used to say the opposite.
    assert_eq!(read_json(&json)["partial"]["is_partial"], true);

    let h = std::fs::read_to_string(&html).unwrap();
    assert!(
        h.contains("reads budget"),
        "a truncated run must name the budget it stopped at"
    );
    assert!(
        !h.contains("no read was skipped"),
        "a run stopped by --reads may not claim every count is exact"
    );
}

#[test]
fn insert_hist_refuses_two_files_that_are_not_in_the_same_pair_order() {
    let d = tmpdir("insert_order");
    let (r1, r2) = write_pairs(&d, 3);
    // Reversed R2: every "pair" is two unrelated reads, which is exactly the
    // case that used to yield a smooth, entirely fictional insert distribution.
    let text = std::fs::read_to_string(&r2).unwrap();
    let mut swapped = String::new();
    for rec in text.trim_end().split('@').filter(|s| !s.is_empty()).rev() {
        swapped.push('@');
        swapped.push_str(rec);
        if !rec.ends_with('\n') {
            swapped.push('\n');
        }
    }
    std::fs::write(&r2, swapped).unwrap();

    let html = d.join("x.html");
    let json = d.join("x.json");
    let (_, stderr, code) = run(
        &d,
        &[
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--insert-hist",
            "--json",
            json.to_str().unwrap(),
            "--html",
            html.to_str().unwrap(),
        ],
    );
    assert_ne!(code, 0, "mis-paired files must be refused, not measured: {stderr}");
    assert!(
        stderr.contains("p0"),
        "the first mismatched pair must be named: {stderr}"
    );
    assert!(
        stderr.contains("p2"),
        "both records of the first mismatch must be named: {stderr}"
    );
    assert!(
        !json.exists() && !html.exists(),
        "no report may be written from mis-paired reads"
    );
}

// ---------------------------------------------------------------------------
// Reproduced-defect regressions (C1-C10) and the output-path contract.
// Each fixture is the one in the defect report, so the numbers below are the
// numbers that were wrong, not a proxy for them.
// ---------------------------------------------------------------------------

/// C1: `base_composition` must sum to `metrics.bases`, including the bases past
/// --max-read-len. Reported: 10 reads x 8 bp at --max-read-len 4 summed to 40 with
/// bases 80, because the block was summed from the truncated per-position rows.
#[test]
fn base_composition_covers_every_base_even_past_max_read_len() {
    let d = tmpdir("composition_all_bases");
    let p = d.join("ten.fastq");
    let mut body = String::new();
    for i in 0..10 {
        body.push_str(&format!("@r{i}\nACGTACGT\n+\nIIIIIIII\n"));
    }
    std::fs::write(&p, body).unwrap();
    let json = d.join("c.json");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--max-read-len", "4", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let comp = &v["base_composition"];
    let sum: u64 = ["a", "c", "g", "t", "other"]
        .iter()
        .map(|k| comp[*k].as_u64().unwrap())
        .sum();
    assert_eq!(v["metrics"]["bases"], 80);
    assert_eq!(sum, 80, "the composition must account for every base: {comp}");
    assert_eq!(comp["a"], 20, "ACGT x 10 reads, 4 As each");
    assert_eq!(v["metrics"]["reads_over_max_read_len"], 10);
    // The per-position matrix still stops at the covered window: the fix widens the
    // totals, it does not attribute untracked positions.
    assert_eq!(v["base_position"].as_object().unwrap().len(), 4);
}

/// C2: `kmer_top_frac` is a share of the k-mer WINDOWS counted, not of the reads.
/// Reported: 2 000 reads of 40 bp gave kmer_top_count 16000 and kmer_top_frac 8.0
/// (29 windows per read -> 58 000 windows, so 16000/58000 = 0.2759).
#[test]
fn kmer_top_fraction_divides_by_windows_and_never_exceeds_one() {
    let d = tmpdir("kmer_windows");
    let p = d.join("k.fastq");
    let read = "ACGT".repeat(10); // 40 bp -> 29 windows
    let mut body = String::new();
    for i in 0..2000 {
        body.push_str(&format!("@r{i}\n{read}\n+\n{}\n", "I".repeat(40)));
    }
    std::fs::write(&p, body).unwrap();
    let json = d.join("k.json");
    let (_, stderr, code) = run(&d, &[p.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&json);
    let m = &v["metrics"];
    assert_eq!(m["kmer_top_count"], 16_000);
    let frac = m["kmer_top_frac"].as_f64().unwrap();
    assert!(
        (frac - 16_000.0 / 58_000.0).abs() < 1e-4,
        "kmer_top_frac must be top/windows = 0.2759, got {frac}"
    );
    assert!(frac <= 1.0, "a share of windows cannot exceed 1, got {frac}");
    let top0 = &v["kmer_spectrum"]["top"][0];
    assert_eq!(top0["count"], 16_000);
    assert!(
        top0["frac"].as_f64().unwrap() <= 1.0,
        "the per-k-mer share is the same denominator: {top0}"
    );
}

/// C3: `--profile` must never publish gc 0.0 for a library with GC. Reported: 800 reads
/// at ~51% GC, every R1 starting with the TruSeq prefix, gave metrics.gc_fraction 0.5138
/// and profile gc 0.0 with no warning — the fragment GC has no denominator when every
/// read is readthrough. It must fall back to the whole-read GC and say so.
#[test]
fn profile_falls_back_to_whole_read_gc_when_every_read_is_readthrough() {
    let d = tmpdir("profile_gc_fallback");
    let p = d.join("adapter.fastq");
    let mut body = String::new();
    for i in 0..800u64 {
        // Adapter prefix on every read (so no fragment base is measurable), body from
        // deterministic DNA so the library has a real, non-zero GC.
        let body_seq = String::from_utf8(dna(37, 7_000 + i)).unwrap();
        body.push_str(&format!("@r{i}\nAGATCGGAAGAGC{body_seq}\n+\n{}\n", "I".repeat(50)));
    }
    std::fs::write(&p, body).unwrap();
    let json = d.join("p.json");
    let prof = d.join("fitted.json");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--profile", prof.to_str().unwrap(),
          "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&prof);
    let report_gc = read_json(&json)["metrics"]["gc_fraction"].as_f64().unwrap();
    let gc = v["gc"].as_f64().unwrap();
    assert!(report_gc > 0.3, "fixture GC: {report_gc}");
    assert!(gc > 0.3, "the profile must not report 0.0 for a {report_gc:.3}-GC library: {gc}");
    assert!(
        (gc - report_gc).abs() < 1e-4,
        "the fallback is the whole-read GC: profile {gc} vs report {report_gc}"
    );
    let fr: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(d.join("fitted.fit_report.json")).unwrap())
            .unwrap();
    assert_eq!(fr["measured"]["gc_is_fragment"], false, "{fr}");
    let warnings = fr["fit_warnings"].as_array().unwrap();
    assert!(
        warnings.iter().any(|w| w.as_str().unwrap_or("").contains("WHOLE-READ GC")),
        "the fallback must be recorded as a fit warning: {warnings:?}"
    );
}

/// C4: `--profile --insert-hist` may claim "MEASURED by --insert-hist" only when the
/// pass resolved an insert. Reported: 30 pairs with no overlap at all (100% censored,
/// every pair on the saturation bin) produced a comment claiming the distribution was
/// measured.
#[test]
fn profile_labels_an_all_censored_insert_distribution_as_censored() {
    let d = tmpdir("profile_insert_censored");
    // 400 bp fragments read from both ends with 50 bp reads: no overlap anywhere.
    let (mut a, mut b) = (String::new(), String::new());
    for i in 0..30u64 {
        let frag = dna(400, 5_000 + i);
        a.push_str(&format!("@q{i}/1\n{}\n+\n{}\n",
            String::from_utf8(frag[..50].to_vec()).unwrap(), "I".repeat(50)));
        b.push_str(&format!("@q{i}/2\n{}\n+\n{}\n",
            String::from_utf8(revcomp(&frag[350..])).unwrap(), "I".repeat(50)));
    }
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    std::fs::write(&r1, a).unwrap();
    std::fs::write(&r2, b).unwrap();
    let json = d.join("c.json");
    let prof = d.join("c.json.prof");
    let (_, stderr, code) = run(
        &d,
        &[r1.to_str().unwrap(), r2.to_str().unwrap(), "--insert-hist",
          "--profile", prof.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let ih = &read_json(&json)["insert_hist"];
    assert_eq!(ih["measured"], 0, "the fixture must measure nothing: {ih}");
    assert_eq!(ih["censored_no_overlap"], 30);
    let v = read_json(&prof);
    let comment = v["_comment"].as_str().unwrap();
    assert!(comment.contains("CENSORED"), "an all-censored CDF is not measured: {comment}");
    assert!(
        !comment.contains("MEASURED by --insert-hist"),
        "the label must not claim a measurement that never happened: {comment}"
    );
    let fr: serde_json::Value = serde_json::from_str(
        // `with_extension` replaces the LAST extension: c.json.prof -> c.json.fit_report.json
        &std::fs::read_to_string(d.join("c.json.fit_report.json")).unwrap(),
    )
    .unwrap();
    assert!(
        fr["invented"].as_array().unwrap().iter().any(|v| v == "insert distribution shape"),
        "the fit report must name the invented shape: {fr}"
    );

    // The other branch: an overlapping pair set IS measured, and says so.
    let (o1, o2) = write_pairs(&d, 100);
    let json2 = d.join("m.json");
    let prof2 = d.join("m.prof");
    let (_, stderr, code) = run(
        &d,
        &[o1.to_str().unwrap(), o2.to_str().unwrap(), "--insert-hist",
          "--profile", prof2.to_str().unwrap(), "--json", json2.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    assert!(read_json(&json2)["insert_hist"]["measured"].as_u64().unwrap() > 90);
    let c2 = read_json(&prof2)["_comment"].as_str().unwrap().to_string();
    assert!(c2.contains("MEASURED by --insert-hist"), "{c2}");
}

/// C5a: mates with different record counts are refused. Reported: R2 truncated to 10 of
/// 50 records gave pairs_sampled 10, pairs_in_input 30, scale_factor 3.0 and exit 0.
#[test]
fn insert_hist_refuses_mates_with_different_record_counts() {
    let d = tmpdir("insert_unequal");
    let (r1, r2) = write_pairs(&d, 20);
    let text = std::fs::read_to_string(&r2).unwrap();
    let mut kept = String::new();
    for (i, rec) in text.trim_end().split('@').filter(|s| !s.is_empty()).enumerate() {
        if i >= 4 {
            break;
        }
        kept.push('@');
        kept.push_str(rec);
        if !rec.ends_with('\n') {
            kept.push('\n');
        }
    }
    std::fs::write(&r2, kept).unwrap();

    let json = d.join("u.json");
    let (stdout, stderr, code) = run(
        &d,
        &[r1.to_str().unwrap(), r2.to_str().unwrap(), "--insert-hist",
          "--json", json.to_str().unwrap(), "-q"],
    );
    assert_ne!(code, 0, "unequal mate counts must be refused, not measured: {stdout}{stderr}");
    assert!(
        stderr.contains("different numbers of records"),
        "the refusal must name the cause: {stderr}"
    );
    assert!(!json.exists(), "no report may be written from a truncated mate file");
}

/// C5b: the reported pair counters are the ones the insert pass actually scanned. Reported:
/// `--reads 4` sampled 50 pairs, reported pairs_in_input 4 and scale_factor 0.08 — a factor
/// derived from the QC scan's budget applied to a different pass's sample.
#[test]
fn insert_hist_reports_the_pairs_it_scanned_and_a_known_input_size() {
    let d = tmpdir("insert_counters");
    let (r1, r2) = write_pairs(&d, 50);

    // A complete run: the pass reaches EOF, so it knows both the sample and the input.
    let full = d.join("full.json");
    let (_, stderr, code) = run(
        &d,
        &[r1.to_str().unwrap(), r2.to_str().unwrap(), "--insert-hist",
          "--json", full.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let ih = &read_json(&full)["insert_hist"];
    assert_eq!(ih["pairs_sampled"], 50);
    assert_eq!(ih["pairs_in_input"], 50, "both files were read to the end: {ih}");
    assert_eq!(ih["scale_factor"], 1.0);

    // --reads 4: the insert pass honours the same budget, and the input's pair count is
    // then genuinely unknown to it (it did not reach EOF), so both are null rather than a
    // rescaled lie.
    let budget = d.join("budget.json");
    let (_, stderr, code) = run(
        &d,
        &[r1.to_str().unwrap(), r2.to_str().unwrap(), "--insert-hist", "--reads", "4",
          "--json", budget.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let ih = &read_json(&budget)["insert_hist"];
    assert_eq!(ih["pairs_sampled"], 4, "the sample stops at the run's own budget: {ih}");
    assert!(
        ih["pairs_in_input"].is_null() && ih["scale_factor"].is_null(),
        "an unknown input size is null, never derived from another pass: {ih}"
    );
}

/// C6: `--overrep-min-count` is THE threshold for over-represented sequences, in the
/// metrics, the `overrep` block and the HTML alike. Reported: at 600 on a 500x sequence
/// the block said [] while repeated_seq_count said 1 and the page listed it.
#[test]
fn overrep_min_count_is_the_same_threshold_everywhere() {
    let d = tmpdir("overrep_threshold");
    let p = d.join("rep.fastq");
    let mut body = String::new();
    for i in 0..500u64 {
        body.push_str(&format!("@d{i}\nACGTACGTACGTACGTACGT\n+\n{}\n", "I".repeat(20)));
    }
    for i in 0..500u64 {
        let seq = String::from_utf8(dna(20, 90_000 + i)).unwrap();
        body.push_str(&format!("@u{i}\n{seq}\n+\n{}\n", "I".repeat(20)));
    }
    std::fs::write(&p, body).unwrap();

    let strict_json = d.join("strict.json");
    let strict_html = d.join("strict.html");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--overrep-mode", "fingerprint", "--overrep-min-count", "600",
          "--json", strict_json.to_str().unwrap(), "--html", strict_html.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&strict_json);
    assert_eq!(v["overrep"]["min_count"], 600, "the block states the flag's threshold");
    assert!(v["overrep"]["sequences"].as_array().unwrap().is_empty());
    assert_eq!(v["metrics"]["repeated_seq_count"], 0, "the metric must honour the flag");
    assert_eq!(v["metrics"]["repeated_seq_frac"], 0.0);
    assert!(v["repeated_sequences"]["top"].as_array().unwrap().is_empty());
    let h = std::fs::read_to_string(&strict_html).unwrap();
    assert!(
        h.contains("at least 600 times"),
        "the page must state the effective threshold and list nothing"
    );
    assert!(!h.contains("<table class=\"rep\">"), "no repeated sequence may be listed");

    // The default threshold (0.1% of 1000 reads = 10) does find the 500x sequence.
    let default_json = d.join("default.json");
    let default_html = d.join("default.html");
    let (_, stderr, code) = run(
        &d,
        &[p.to_str().unwrap(), "--overrep-mode", "fingerprint",
          "--json", default_json.to_str().unwrap(),
          "--html", default_html.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    let v = read_json(&default_json);
    assert_eq!(v["metrics"]["repeated_seq_count"], 1);
    assert_eq!(v["metrics"]["repeated_seq_frac"], 0.5);
    assert_eq!(v["overrep"]["sequences"].as_array().unwrap().len(), 1);
    let h = std::fs::read_to_string(&default_html).unwrap();
    assert!(h.contains("<table class=\"rep\">"), "the default threshold lists it");
}

/// C7: a `--max-mem` below the old fixed 512 MB reserve must admit a run that fits.
/// Reported: `--max-mem 568` refused every run ("needs ~0.0 MB ... limited to 568.0 MB")
/// while 569 proceeded.
///
/// REVERSAL, named deliberately: this test used to end by asserting that `--max-mem 64`
/// REFUSES, which was the fixture pinning the second half of the same bug — the reserve's
/// 64 MiB FLOOR made every ceiling below ~72 MB refuse every run whatever the projection
/// (measured on an 1800-base file whose table projects 0.0 MB: 71 refused, 72 proceeded).
/// The floor is now 8 MiB, so 32 and 64 MB admit a run that needs ~0 bytes, and only a
/// ceiling at the floor itself (8 MB) has nothing left to admit.
#[test]
fn a_max_mem_below_the_old_reserve_admits_a_small_run() {
    let d = tmpdir("maxmem_small");
    let fq = write_fixture(&d);
    let big = write_120bp_fixture(&d);
    for mb in ["32", "64", "256", "512", "568"] {
        let (_, stderr, code) =
            run(&d, &[fq.to_str().unwrap(), "--max-mem", mb, "-q"]);
        assert_eq!(code, 0, "--max-mem {mb} must run: {stderr}");
        assert!(!stderr.contains("limited to"), "{stderr}");
        // ... and the 1800-base file whose table projects ~0.0 MB, which the 71/72 cliff
        // was measured on, stays admitted at the two smallest ceilings.
        let (_, stderr, code) =
            run(&d, &[big.to_str().unwrap(), "--max-mem", mb, "-q"]);
        assert_eq!(code, 0, "--max-mem {mb} on a 1800-base file must run: {stderr}");
    }
    // A ceiling at the reserve floor itself has nothing left for the table: refused, and
    // the refusal says what the limit was.
    let (stdout, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--max-mem", "8", "-q"]);
    assert_ne!(code, 0, "an 8 MB ceiling is the 8 MiB reserve: {stdout}{stderr}");
    assert!(stderr.contains("limited to"), "{stderr}");
}

/// C8: an absurd `--max-read-len` is refused before anything is allocated. Reported:
/// `--insert-hist --max-read-len 1099511627776` aborted with
/// "memory allocation of 17592186044424 bytes failed" AFTER the whole scan.
#[test]
fn an_absurd_max_read_len_is_refused_before_it_is_allocated() {
    let d = tmpdir("max_read_len");
    let (r1, r2) = write_pairs(&d, 2);
    let json = d.join("x.json");
    let (stdout, stderr, code) = run(
        &d,
        &[r1.to_str().unwrap(), r2.to_str().unwrap(), "--insert-hist",
          "--max-read-len", "1099511627776", "--json", json.to_str().unwrap(), "-q"],
    );
    assert_ne!(code, 0, "the value must be refused: {stdout}{stderr}");
    assert!(stderr.contains("--max-read-len"), "name the flag: {stderr}");
    assert!(stderr.contains("maximum"), "say what the maximum is: {stderr}");
    assert!(
        !stderr.contains("memory allocation"),
        "it must refuse, not reach the allocator: {stderr}"
    );
    assert!(!stderr.contains("accumulator state"), "refused before the scan: {stderr}");
    assert!(!json.exists());
}

/// C9: two reports cannot share stdout. Reported: `--json - --positions-tsv -` wrote
/// 596 KB of JSON with TSV threaded through it, and `--json - --html -` appended a page.
#[test]
fn two_outputs_on_one_stdout_are_refused() {
    let d = tmpdir("stdout_conflict");
    let fq = write_fixture(&d);
    for extra in ["--positions-tsv", "--html"] {
        let (stdout, stderr, code) = run(
            &d,
            &[fq.to_str().unwrap(), "--json", "-", extra, "-", "-q"],
        );
        assert_ne!(code, 0, "--json - {extra} - must be refused: {stdout}{stderr}");
        assert!(stderr.contains("stdout"), "name the stream: {stderr}");
        assert!(stdout.is_empty(), "nothing may be written to a refused stream: {stdout}");
    }
    // One stdout output stays legal.
    let (stdout, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--json", "-", "-q"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(stdout.starts_with('{'), "a lone --json - still writes the report: {stdout}");
}

/// C10: an empty input is a refusal, not an all-zero 593 KB report. Reported: a 0-byte
/// file exited 0 and wrote a report of zeros; the other FASTQ tools refuse.
#[test]
fn an_empty_input_is_refused_instead_of_reported_as_zero() {
    let d = tmpdir("empty_input");
    let empty = d.join("empty.fastq");
    std::fs::write(&empty, b"").unwrap();
    let blank = d.join("blank.fastq");
    std::fs::write(&blank, b"\n\n\n\n").unwrap();
    for p in [&empty, &blank] {
        let json = d.join("e.json");
        let (stdout, stderr, code) =
            run(&d, &[p.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"]);
        assert_ne!(code, 0, "an empty input must be refused: {stdout}{stderr}");
        assert!(
            stderr.contains("no FASTQ records"),
            "the refusal must say what was missing: {stderr}"
        );
        assert!(!json.exists(), "no zero report may be written");
    }
    // --profile keeps its own refusal for the same input.
    let prof = d.join("p.json");
    let (_, stderr, code) =
        run(&d, &[empty.to_str().unwrap(), "--profile", prof.to_str().unwrap(), "-q"]);
    assert_ne!(code, 0, "{stderr}");
}

/// An output path's parent directory is created, the way fastmaps does: the report path
/// is a promise to write the report, not a place the user has to mkdir first.
#[test]
fn output_parent_directories_are_created() {
    let d = tmpdir("mkdir_outputs");
    let fq = write_fixture(&d);
    let json = d.join("out/deep/r.json");
    let html = d.join("out/deep/r.html");
    let tsv = d.join("out/deep/r.tsv");
    let prof = d.join("out/deep/r.profile.json");
    let (_, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(),
          "--json", json.to_str().unwrap(),
          "--html", html.to_str().unwrap(),
          "--positions-tsv", tsv.to_str().unwrap(),
          "--profile", prof.to_str().unwrap(), "-q"],
    );
    assert_eq!(code, 0, "{stderr}");
    for p in [&json, &html, &tsv, &prof] {
        assert!(p.is_file(), "{} was not written", p.display());
    }
    assert!(
        d.join("out/deep/r.profile.fit_report.json").is_file(),
        "the fit report travels next to the profile it describes"
    );
}

// ---------------------------------------------------------------------------
// Round-3 defects: `-` on every report output, the summary guard, the --max-mem
// reserve floor, $FASTSUITE_THREADS + busy_cores, the head-dropout denominator,
// --processed-2 on its own, and the summary line's arithmetic.
// ---------------------------------------------------------------------------

/// 120 bp reads whose G fraction is 0.3 in the first 20 positions, in the first 40 AND in
/// the whole read, with an A at position 1. The head-dropout baseline is a G FRACTION, so
/// it must read 0.3 whatever `--max-read-len` is: the per-position planes cover a prefix
/// of the read, and a fraction over that prefix is cap-invariant by construction.
fn write_120bp_fixture(dir: &Path) -> PathBuf {
    let read = format!(
        "A{}{}{}{}{}{}",
        "G".repeat(6),
        "C".repeat(13),
        "G".repeat(6),
        "C".repeat(14),
        "G".repeat(24),
        "C".repeat(56)
    );
    assert_eq!(read.len(), 120);
    assert_eq!(read.matches('G').count(), 36, "36/120 = 0.3");
    let p = dir.join("r120.fastq");
    let mut body = String::new();
    for i in 0..8 {
        body.push_str(&format!("@r{i}\n{read}\n+\n{}\n", "I".repeat(120)));
    }
    std::fs::write(&p, body).unwrap();
    p
}

/// `--html -` was a file named `-` (46 KB) with an empty stdout.
#[test]
fn html_dash_streams_to_stdout_and_writes_no_file_named_dash() {
    let d = tmpdir("html_dash");
    let fq = write_fixture(&d);
    let (stdout, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--html", "-", "-q"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        stdout.trim_start().starts_with("<!DOCTYPE html>"),
        "the page itself must reach stdout: {}",
        &stdout[..stdout.len().min(120)]
    );
    assert!(stdout.trim_end().ends_with("</html>"));
    assert!(!d.join("-").exists(), "no file may be named `-`");
}

/// `--profile -` wrote a file named `-` plus `-.fit_report.json`, and left stdout empty.
#[test]
fn profile_dash_streams_to_stdout_and_writes_no_file_named_dash() {
    let d = tmpdir("profile_dash");
    let fq = write_fixture(&d);
    let (stdout, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--profile", "-", "-q"]);
    assert_eq!(code, 0, "{stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout)
        .unwrap_or_else(|e| panic!("stdout must be the profile JSON: {e}\n{stdout}"));
    assert!(v.is_object() && v.get("read_len").is_some(), "{stdout}");
    assert!(!d.join("-").exists(), "no file may be named `-`");
    assert!(
        !d.join("-.fit_report.json").exists(),
        "the fit report must not be named after the stream"
    );
}

/// Both really stream now, so their collision is genuine — and the refusal must name both.
#[test]
fn html_and_json_on_stdout_are_refused_naming_both_flags() {
    let d = tmpdir("html_json_stream");
    let fq = write_fixture(&d);
    let (stdout, stderr, code) =
        run(&d, &[fq.to_str().unwrap(), "--html", "-", "--json", "-", "-q"]);
    assert_ne!(code, 0, "two reports on one stream must be refused: {stdout}{stderr}");
    assert!(stderr.contains("--json") && stderr.contains("--html"), "{stderr}");
    assert!(stdout.is_empty(), "nothing may be written to a refused stream: {stdout}");
    assert!(!d.join("-").exists());
}

/// `--positions-tsv -` used to get the prose summary appended to the TSV stream: the
/// guard keyed on `--json` only.
#[test]
fn the_summary_line_never_lands_on_a_positions_tsv_stream() {
    let d = tmpdir("tsv_stream");
    let fq = write_fixture(&d);
    let (stdout, stderr, code) = run(&d, &[fq.to_str().unwrap(), "--positions-tsv", "-"]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        stdout.lines().count(),
        11,
        "header + 10 positions and nothing else: {stdout}"
    );
    for line in stdout.lines() {
        assert!(
            !line.starts_with("fastcheck:"),
            "the summary must not be appended to the TSV: {line}"
        );
    }
    assert!(
        stderr.lines().any(|l| l.starts_with("fastcheck: ") && l.contains("reads")),
        "the summary still has to be printed somewhere: {stderr}"
    );
}

/// The head-dropout baseline divided a per-position G sum by EVERY base in the file, so
/// `--max-read-len` scaled it: -0.4639 / -0.4181 / -0.2667 on one 1800-base fixture.
#[test]
fn head_dropout_is_invariant_to_max_read_len() {
    let d = tmpdir("head_dropout_cap");
    let fq = write_120bp_fixture(&d);
    let mut seen = Vec::new();
    for cap in ["20", "40", "1000"] {
        let json = d.join(format!("hd{cap}.json"));
        let (_, stderr, code) = run(
            &d,
            &[fq.to_str().unwrap(), "--max-read-len", cap, "--json", json.to_str().unwrap(), "-q"],
        );
        assert_eq!(code, 0, "{stderr}");
        let m = read_json(&json)["metrics"].clone();
        assert_eq!(m["head_g_frac"], 0.0, "position 1 is A at every cap: {m}");
        assert_eq!(
            m["head_dropout"], 0.3,
            "cap {cap}: baseline - g1 is a G fraction over the bases the planes hold: {m}"
        );
        seen.push(m["head_dropout"].as_f64().unwrap());
    }
    assert!(seen.windows(2).all(|w| w[0] == w[1]), "{seen:?}");
}

/// The profile's `g_mean` had the same denominator bug (0.0361 against a true 0.2333) and
/// is published in the fit report.
#[test]
fn profile_g_mean_is_invariant_to_max_read_len() {
    let d = tmpdir("g_mean_cap");
    let fq = write_120bp_fixture(&d);
    for cap in ["20", "40", "1000"] {
        let prof = d.join(format!("p{cap}.json"));
        let (_, stderr, code) = run(
            &d,
            &[fq.to_str().unwrap(), "--max-read-len", cap, "--profile", prof.to_str().unwrap(), "-q"],
        );
        assert_eq!(code, 0, "{stderr}");
        let fit = prof.with_extension("fit_report.json");
        let v = read_json(&fit);
        assert_eq!(
            v["measured"]["g_mean"], 0.3,
            "cap {cap}: g_mean is a G fraction over the scanned bases: {v}"
        );
    }
}

/// `--processed-2` without `--processed-1` was accepted, then silently ignored
/// (`processed: null`, exit 0).
#[test]
fn processed_2_without_processed_1_is_refused() {
    let d = tmpdir("processed2_alone");
    let fq = write_fixture(&d);
    let p2 = d.join("p2.fastq");
    std::fs::write(&p2, FASTQ).unwrap();
    let json = d.join("r.json");
    let (stdout, stderr, code) = run(
        &d,
        &[fq.to_str().unwrap(), "--processed-2", p2.to_str().unwrap(),
          "--json", json.to_str().unwrap(), "-q"],
    );
    assert_ne!(code, 0, "--processed-2 alone must not be silently ignored: {stdout}{stderr}");
    assert!(stderr.contains("--processed-1"), "name the flag it needs: {stderr}");
    assert!(!json.exists(), "nothing may be measured under a refused flag combination");
}

/// `$FASTSUITE_THREADS` supplies the default (the suite's convention); `-j` wins outright,
/// and unset / unparseable / 0 all mean "not set".
#[test]
fn fastsuite_threads_supplies_the_default_and_explicit_j_wins() {
    let d = tmpdir("fastsuite_threads");
    let fq = write_fixture(&d);
    let default = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) as u64;
    for (env, j, want) in [
        ("4", Some("2"), 2u64),
        ("4", None, 4),
        ("0", None, default),
        ("not-a-number", None, default),
    ] {
        let json = d.join(format!("t_{env}_{}.json", j.unwrap_or("x")));
        let mut args = vec![fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-q"];
        if let Some(j) = j {
            args.push("-j");
            args.push(j);
        }
        let (_, stderr, code) = run_env(&d, &args, &[("FASTSUITE_THREADS", env)]);
        assert_eq!(code, 0, "{stderr}");
        assert_eq!(
            read_json(&json)["metrics"]["threads"], want,
            "FASTSUITE_THREADS={env}, -j={j:?}"
        );
    }
}

/// `metrics.threads` is the `-j` that was ASKED FOR; `busy_cores` is what was used. Without
/// it a report could claim `threads: 8` while one core did the work (fastprep ships it).
#[test]
fn the_report_publishes_the_busy_cores_it_measured() {
    let d = tmpdir("busy_cores");
    let fq = write_fixture(&d);
    let json = d.join("b.json");
    let (_, stderr, code) =
        run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap(), "-j", "8", "-q"]);
    assert_eq!(code, 0, "{stderr}");
    let m = read_json(&json)["metrics"].clone();
    let busy = m["busy_cores"]
        .as_f64()
        .unwrap_or_else(|| panic!("busy_cores must be a number: {m}"));
    assert!(busy >= 0.0, "{m}");
    assert!(
        busy < 8.0,
        "a 6-read plain file cannot keep 8 cores busy, whatever -j says: {m}"
    );
}

/// The summary line and the report now share one implementation (they used to be two
/// copies): the two numbers a user sees in the summary are the report's own.
#[test]
fn the_summary_line_quotes_the_reports_own_gc_and_mean_q() {
    let d = tmpdir("summary_agrees");
    let fq = write_fixture(&d);
    let json = d.join("r.json");
    let (stdout, stderr, code) =
        run(&d, &[fq.to_str().unwrap(), "--json", json.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let m = read_json(&json)["metrics"].clone();
    let line = stdout.lines().next().unwrap();
    let gc = m["gc_fraction"].as_f64().unwrap();
    let q = m["mean_read_mean_q"].as_f64().unwrap();
    assert!(line.contains(&format!("GC {:.2}%", 100.0 * gc)), "{line} vs {m}");
    assert!(line.contains(&format!("mean read Q {:.2}", q)), "{line} vs {m}");
}
