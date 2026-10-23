//! Golden-stream pin: the sha256 of the three artefacts `fastsim` emits for
//! `--seed 42 --pairs 1000` on the committed `mandrill_E358` profile.
//!
//! fastsim's output is the corpus generator every downstream bench and paper
//! is validated against, so the emitted FASTQ/GT bytes are a contract, not an
//! implementation detail. FS-1's `--gz-level` flag must not move them: the
//! level only changes the DEFLATE encoding, never the payload, so the
//! `.gt.tsv` and the decompressed FASTQ stay identical at every level while
//! the `.gz` bytes legitimately differ.
//!
//! The hashes below were produced by the pre-optimisation binary (0.1.0) and
//! are therefore the 0.1.0 stream, not a snapshot of whatever the current
//! code happens to print. A failure here means the generator moved, which is
//! out of scope for any optimisation batch and needs its own decision.
//!
//! `sha256` is implemented locally (no new dependency): the test also
//! self-checks it against the two NIST vectors so a broken hash cannot be
//! mistaken for a broken stream.
//!
//! Portability caveat: the stream is a pure function of `(seed, profile, code)`
//! plus the platform's `libm` (`exp`/`powf` feed the quality model). The repo
//! pins the compiler in `pixi.toml` for exactly this reason; an unrelated libm
//! change could in principle flip a rounded Phred at a `.5` boundary and move
//! these hashes without any model change. If that ever happens, bisect against
//! the 0.1.0 binary before touching the constants.

use std::path::{Path, PathBuf};
use std::process::Command;

/// Where the dataset-derived profile would live if it were in the repo.
///
/// It deliberately is NOT: `.gitignore` keeps every dataset-derived file out
/// (the repo must stay shippable as a bioconda package / container image), so a
/// clean clone and CI have no `mandrill_E358.json`. These golden tests then SKIP
/// with a visible line instead of failing, and `FASTSIM_GOLDEN_PROFILE` points
/// at the local copy for a full run:
///
/// ```text
/// FASTSIM_GOLDEN_PROFILE=~/hermes_scratch/fastsuite-local-artifacts/mandrill_E358.json \
///     cargo test --test golden_stream
/// ```
const PROFILE_DEFAULT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/profiles/mandrill_E358.json");

/// The path that was searched: `FASTSIM_GOLDEN_PROFILE` when set, else the
/// in-tree location the data policy keeps empty.
fn profile_path() -> PathBuf {
    std::env::var_os("FASTSIM_GOLDEN_PROFILE")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(PROFILE_DEFAULT))
}

fn profile() -> Option<PathBuf> {
    let p = profile_path();
    p.is_file().then_some(p)
}

/// `None` means "the profile is not on this machine": the test reports why it
/// did not run and passes, so the suite count stays honest.
macro_rules! profile_or_skip {
    () => {
        match profile() {
            Some(p) => p,
            None => {
                eprintln!(
                    "SKIP: no profile at {} — the repo ships no dataset-derived file (data \
                     policy), so the pinned stream was not checked on this machine. Point \
                     FASTSIM_GOLDEN_PROFILE at a local copy to run it.",
                    profile_path().display()
                );
                return;
            }
        }
    };
}

/// `--seed 42 --pairs 1000`, gzip level 6 (the default).
const GOLDEN_GZ_R1: &str = "693b6616dab8736f842abb2a4863bdf8a9bfdc6bd921eccc799f90e41d88bd13";
const GOLDEN_GZ_R2: &str = "af6247b21a3003257c016e3b14452176e9d8395a9e157c5d6dbe78b683c0af49";
/// The GT pin was deliberately re-based once, in the 2026-10-21 GT-contract
/// fix: the row gained per-mate `r1_full`/`r2_full` columns and `full_flag`
/// became "BOTH mates are artefact" (it used to be "either"). The FASTQ pins
/// did NOT move, which is the point of pinning them separately: the payload
/// is unchanged, only the labels changed.
const GOLDEN_GT: &str = "56da7e652c0f66b06bab65a7455a97e8b17dbf1c2cf18471aa2c6052cec423bc";
/// The same stream with the compression removed (what any level must inflate
/// back to).
const GOLDEN_PLAIN_R1: &str = "43738c1a013cd9bf40561b9fc7618c90a909799d81ce312fe271f8e855e3f225";
const GOLDEN_PLAIN_R2: &str = "6419bf52f4f65328dd4bc20cafa32d53ab19bccdf8151fd016be9d343f374b12";

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastsim")
}

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fastsim_golden_{}_{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Run one 1000-pair simulation and return the prefix written.
fn simulate(profile: &Path, tag: &str, extra: &[&str]) -> (PathBuf, PathBuf) {
    let dir = scratch(tag);
    let out = dir.join("g");
    let mut cmd = Command::new(bin());
    cmd.args([
        "--pairs",
        "1000",
        "--seed",
        "42",
        "--profile",
        profile.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--force",
        "--quiet",
    ])
    .args(extra);
    let res = cmd.output().unwrap();
    assert!(
        res.status.success(),
        "fastsim {extra:?} failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    (dir, out)
}

fn sha_of(path: &Path) -> String {
    sha256_hex(&std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display())))
}

fn gunzip(path: &Path) -> Vec<u8> {
    use std::io::Read;
    let f = std::fs::File::open(path).unwrap();
    let mut gz = flate2::read::GzDecoder::new(f);
    let mut buf = Vec::new();
    gz.read_to_end(&mut buf).unwrap();
    buf
}

#[test]
fn default_stream_is_the_pinned_0_1_0_stream() {
    let profile = profile_or_skip!();
    let (dir, out) = simulate(&profile, "default", &[]);
    let r1 = PathBuf::from(format!("{}_R1.fastq.gz", out.display()));
    let r2 = PathBuf::from(format!("{}_R2.fastq.gz", out.display()));
    let gt = PathBuf::from(format!("{}.gt.tsv", out.display()));
    assert_eq!(sha_of(&r1), GOLDEN_GZ_R1, "R1 gzip stream moved");
    assert_eq!(sha_of(&r2), GOLDEN_GZ_R2, "R2 gzip stream moved");
    assert_eq!(sha_of(&gt), GOLDEN_GT, "GT rows moved");
    // The default level is 6, so the payload must be the pinned plain stream.
    assert_eq!(sha256_hex(&gunzip(&r1)), GOLDEN_PLAIN_R1);
    assert_eq!(sha256_hex(&gunzip(&r2)), GOLDEN_PLAIN_R2);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn plain_output_is_the_same_stream_without_deflate() {
    let profile = profile_or_skip!();
    let (dir, out) = simulate(&profile, "plain", &["--no-gzip"]);
    let r1 = PathBuf::from(format!("{}_R1.fastq", out.display()));
    let r2 = PathBuf::from(format!("{}_R2.fastq", out.display()));
    let gt = PathBuf::from(format!("{}.gt.tsv", out.display()));
    assert_eq!(sha_of(&r1), GOLDEN_PLAIN_R1, "plain R1 moved");
    assert_eq!(sha_of(&r2), GOLDEN_PLAIN_R2, "plain R2 moved");
    assert_eq!(sha_of(&gt), GOLDEN_GT);
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `--gz-level` changes the DEFLATE encoding only: the payload and the GT stay
/// pinned at every level, and the levels really do differ in size (the
/// documented 17.9% trade is what makes the flag worth having).
#[test]
fn gz_level_changes_bytes_but_not_the_stream() {
    let profile = profile_or_skip!();
    let (dir6, out6) = simulate(&profile, "lvl6", &["--gz-level", "6"]);
    let base6 = std::fs::metadata(format!("{}_R1.fastq.gz", out6.display()))
        .unwrap()
        .len();
    // An explicit --gz-level 6 is the same bytes as the default.
    assert_eq!(
        sha_of(Path::new(&format!("{}_R1.fastq.gz", out6.display()))),
        GOLDEN_GZ_R1,
        "--gz-level 6 must reproduce the default bytes"
    );

    for level in ["1", "4", "9"] {
        let (dir, out) = simulate(&profile, &format!("lvl{level}"), &["--gz-level", level]);
        let r1 = PathBuf::from(format!("{}_R1.fastq.gz", out.display()));
        let r2 = PathBuf::from(format!("{}_R2.fastq.gz", out.display()));
        let gt = PathBuf::from(format!("{}.gt.tsv", out.display()));
        assert_eq!(sha_of(&gt), GOLDEN_GT, "level {level} moved the GT");
        assert_eq!(sha256_hex(&gunzip(&r1)), GOLDEN_PLAIN_R1, "level {level} R1 payload");
        assert_eq!(sha256_hex(&gunzip(&r2)), GOLDEN_PLAIN_R2, "level {level} R2 payload");
        assert_ne!(sha_of(&r1), GOLDEN_GZ_R1, "level {level} must re-encode");
        if level == "1" {
            let s1 = std::fs::metadata(&r1).unwrap().len();
            assert!(
                s1 > base6,
                "level 1 is the fast/larger trade: {s1} vs level 6 {base6}"
            );
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
    std::fs::remove_dir_all(&dir6).unwrap();
}

// ---------------------------------------------------------------------------
// SHA-256 (FIPS 180-4), local so the gate needs no extra crate.
// ---------------------------------------------------------------------------

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

fn sha256_hex(data: &[u8]) -> String {
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bitlen = (data.len() as u64).wrapping_mul(8);
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bitlen.to_be_bytes());

    for chunk in msg.chunks_exact(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[4 * i],
                chunk[4 * i + 1],
                chunk[4 * i + 2],
                chunk[4 * i + 3],
            ]);
        }
        for i in 16..64 {
            let a = w[i - 15];
            let b = w[i - 2];
            let s0 = a.rotate_right(7) ^ a.rotate_right(18) ^ (a >> 3);
            let s1 = b.rotate_right(17) ^ b.rotate_right(19) ^ (b >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let mut v = h;
        for i in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7]
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v[7] = v[6];
            v[6] = v[5];
            v[5] = v[4];
            v[4] = v[3].wrapping_add(t1);
            v[3] = v[2];
            v[2] = v[1];
            v[1] = v[0];
            v[0] = t1.wrapping_add(t2);
        }
        for i in 0..8 {
            h[i] = h[i].wrapping_add(v[i]);
        }
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

#[test]
fn local_sha256_matches_the_nist_vectors() {
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}
