//! End-to-end test of `fasthelper compare`.

use std::path::PathBuf;
use std::process::Command;

fn fasthelper() -> &'static str {
    env!("CARGO_BIN_EXE_fasthelper")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fh_compare_{}_{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_report(dir: &PathBuf, name: &str, reads: u64, gc: f64) -> PathBuf {
    let p = dir.join(name);
    let j = format!(
        r#"{{"schema":"fastsuite.report.v1","tool":"fastcheck","metrics":{{"reads":{reads},"gc_fraction":{gc},"mean_read_mean_q":35.0}}}}"#
    );
    std::fs::write(&p, j).unwrap();
    p
}

#[test]
fn compare_produces_before_after_delta() {
    let d = tmpdir("cmp");
    let before = write_report(&d, "before.json", 1000, 0.45);
    let after = write_report(&d, "after.json", 800, 0.40);
    let out = d.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let r = Command::new(fasthelper())
        .args(["compare", "--before", before.to_str().unwrap(), "--after", after.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let tsv = std::fs::read_to_string(out.join("compare.tsv")).unwrap();
    assert!(tsv.contains("reads\t1000\t800\t-200"), "reads row: {tsv}");
    assert!(tsv.contains("gc_fraction\t0.45\t0.4\t-0.0"), "gc row: {tsv}");
    assert!(std::fs::read_to_string(out.join("compare.html")).unwrap().contains("fasthelper compare"));
}
