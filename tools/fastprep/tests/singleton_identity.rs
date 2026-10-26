//! Regression: the paired report must satisfy
//! `reads_in == reads_out + reads_dropped` even when a mate survives but has
//! nowhere to go.
//!
//! At 0.5.0 a discarded singleton (a pair whose R1 was dropped by --mean-q and
//! whose R2 was fine, with no `--singletons2`) was counted in `legacy` only:
//! 100 pairs reported `reads_in 200`, `reads_out 0`, `reads_dropped 100`,
//! `r2_reads_dropped 0` while `legacy.discarded_singletons_r2 = 100`.

use std::io::Write;
use std::process::Command;

#[test]
fn discarded_singletons_are_counted_as_dropped() {
    let d = std::env::temp_dir().join(format!(
        "fastprep-singletons-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&d).unwrap();
    let (r1, r2, o1, o2, rep) = (
        d.join("r1.fq"),
        d.join("r2.fq"),
        d.join("o1.fq"),
        d.join("o2.fq"),
        d.join("rep.json"),
    );
    let mut f1 = std::fs::File::create(&r1).unwrap();
    let mut f2 = std::fs::File::create(&r2).unwrap();
    for i in 0..100u32 {
        // Unique per pair (dedup is on by default and would otherwise eat 99
        // of the 100 records before --mean-q ever ran).
        let idx: String = (0..8)
            .map(|k| ['A', 'C', 'G', 'T'][((i >> (2 * k)) & 3) as usize])
            .collect();
        let seq = format!("{idx}ACGTACGTACGTACGTACGTACGT");
        // R1 is Q2 everywhere (below --mean-q 25), R2 is Q40.
        writeln!(f1, "@p{i}\n{seq}\n+\n{}", "#".repeat(32)).unwrap();
        writeln!(f2, "@p{i}\n{seq}\n+\n{}", "I".repeat(32)).unwrap();
    }
    drop(f1);
    drop(f2);

    let out = Command::new(env!("CARGO_BIN_EXE_fastprep"))
        .args([
            "-1",
            r1.to_str().unwrap(),
            "-2",
            r2.to_str().unwrap(),
            "-o",
            o1.to_str().unwrap(),
            "-p",
            o2.to_str().unwrap(),
            "--json",
            rep.to_str().unwrap(),
            "--no-adapter",
            "--allow-total-drop",
            "-q",
        ])
        .output()
        .expect("fastprep ran");
    assert!(
        out.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&rep).unwrap()).unwrap();
    let m = &doc["metrics"];
    let n = |k: &str| m[k].as_u64().unwrap_or(0);
    assert_eq!(n("reads_in"), 200, "reads_in is the sum over mates");
    assert_eq!(
        n("reads_in"),
        n("reads_out") + n("reads_dropped"),
        "the record identity must hold: reads_in {n_rin}, reads_out {n_rout}, reads_dropped {n_rdrop}",
        n_rin = n("reads_in"),
        n_rout = n("reads_out"),
        n_rdrop = n("reads_dropped")
    );
    assert_eq!(n("discarded_singletons_r2"), 100);
    assert_eq!(
        n("r2_reads_dropped"),
        100,
        "a discarded R2 is dropped, not vanished"
    );
    assert_eq!(
        n("r2_reads_in"),
        n("r2_reads_out") + n("r2_reads_dropped"),
        "the identity must also hold per mate"
    );
    let _ = std::fs::remove_dir_all(&d);
}
