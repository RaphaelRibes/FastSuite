//! FM-10 — the FASTQ loop reads bytes, not `String`s.
//!
//! `BufRead::lines()` allocates a `String` per line and validates UTF-8;
//! `read_until` into reused `Vec<u8>` does neither. The two must agree byte for
//! byte, including on CRLF input: `lines()` strips a trailing `\r`, so the
//! byte loop has to strip it too or every `\r` read would gain a base of length
//! and move the (x, y) plane.
//!
//! The same path is overwritten between runs so the file path embedded in
//! `heatmap.json` is identical and *all five* outputs can be compared.

use std::path::Path;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn body(eol: &str) -> String {
    let reads = [
        ("@r0", "A".repeat(150)),
        ("@r1", "ACGT".repeat(37) + "AC"),
        ("@r2", "T".repeat(140) + &"G".repeat(10)),
        ("@r3", "N".repeat(150)),
        ("@r4", "ACGTACGTAC".repeat(15)),
    ];
    let mut s = String::new();
    for (h, seq) in reads {
        s.push_str(&format!("{h}{eol}{seq}{eol}+{eol}{}{eol}", "I".repeat(seq.len())));
    }
    s
}

fn run(fq: &Path, out: &Path) {
    let o = Command::new(bin())
        .arg("-o")
        .arg(out)
        .arg("--force")
        .arg(fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
}

#[test]
fn crlf_input_matches_lf_input_byte_for_byte() {
    let dir = std::env::temp_dir().join(format!("maps_crlf_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fq = dir.join("x_R1.fastq");
    let out_lf = dir.join("out_lf");
    let out_crlf = dir.join("out_crlf");

    std::fs::write(&fq, body("\n")).unwrap();
    run(&fq, &out_lf);

    std::fs::write(&fq, body("\r\n")).unwrap();
    run(&fq, &out_crlf);

    for f in [
        "heatmap.json",
        "heatmap.png",
        "heatmap.data",
        "heatmap.preview",
        "heatmap_3d.html",
    ] {
        let a = std::fs::read(out_lf.join(f)).unwrap();
        let b = std::fs::read(out_crlf.join(f)).unwrap();
        assert_eq!(
            a.len(),
            b.len(),
            "{f}: CRLF input must not change the output size ({} vs {})",
            a.len(),
            b.len()
        );
        assert!(a == b, "{f}: CRLF input must be byte-identical to LF input");
    }

    // And the read length really is 150 for every read in both runs.
    let text = std::fs::read_to_string(out_crlf.join("heatmap.json")).unwrap();
    let v: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(v["size"], 151, "a stray \\r would make size 152");
    assert_eq!(v["total_reads"], 5);

    let _ = std::fs::remove_dir_all(&dir);
}
