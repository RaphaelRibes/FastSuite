//! Biological insert generation: random DNA (GC-controlled) or sampled from
//! a reference FASTA.
//!
//! An insert is a double-stranded DNA fragment of length I. R1 reads the 5'
//! end of the top strand; R2 reads the 3' end of the top strand (i.e. the
//! reverse complement of the fragment's bottom strand start). The artefact
//! stage (artefact.rs) truncates the insert to the read window.

use anyhow::{bail, Context, Result};
use rand::Rng;

use crate::randutil::sample_gauss;

pub const BASES: [u8; 4] = *b"ACGT";

/// Draw an insert length I ~ N(mean, sd), floored at `floor`.
pub fn sample_insert_len<R: Rng + ?Sized>(mean: f64, sd: f64, floor: u32, rng: &mut R) -> u32 {
    sample_gauss(mean, sd, floor as f64, 1e6, rng).round() as u32
}

/// Generate `len` random DNA bases with target GC fraction.
pub fn random_dna<R: Rng + ?Sized>(len: usize, gc: f64, rng: &mut R) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        let u: f64 = rng.gen();
        // A/T share (1-gc)/2, C/G share gc/2 — same marginal for both pairs
        let b = if u < gc / 2.0 {
            b'C'
        } else if u < gc {
            b'G'
        } else if u < gc + (1.0 - gc) / 2.0 {
            b'A'
        } else {
            b'T'
        };
        out.push(b);
    }
    out
}

/// DNA where the 5′-most `head.len()` bases are sampled from a per-position
/// `[a, c, g, t]` composition and the 3′-most `tail.len()` bases from `tail`
/// (mirrored from the 3′ end), instead of uniform GC. This is how *located*
/// base biases — a 5′ dropout that starts mid-read, a polyG tail, an enriched
/// position — are reproduced: a run-length model cannot, a per-position
/// composition can. `tail` is in the fragment's own 5′→3′ frame, so the caller
/// reverse-complements R2's measured head before passing it.
pub fn random_dna_bias<R: Rng + ?Sized>(
    len: usize,
    gc: f64,
    head: Option<&[[f32; 4]]>,
    tail: Option<&[[f32; 4]]>,
    rng: &mut R,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        // 3′ end: position i is `len-1-i` from the 3′ edge; the tail window
                // covers only the last `tail.len()` positions.
                let pr = head
                    .and_then(|v| v.get(i))
                    .or_else(|| {
                        let d = len - 1 - i;
                        if d < tail.map_or(0, |v| v.len()) {
                            tail.and_then(|v| v.get(d))
                        } else {
                            None
                        }
                    });
        if let Some(pr) = pr {
            let u: f32 = rng.gen();
            let a = pr[0];
            if u < a {
                out.push(b'A');
            } else if u < a + pr[1] {
                out.push(b'C');
            } else if u < a + pr[1] + pr[2] {
                out.push(b'G');
            } else {
                out.push(b'T');
            }
        } else {
            let u: f64 = rng.gen();
            out.push(if u < gc / 2.0 {
                b'C'
            } else if u < gc {
                b'G'
            } else if u < gc + (1.0 - gc) / 2.0 {
                b'A'
            } else {
                b'T'
            });
        }
    }
    out
}

/// Complement a single base (A<->T, C<->G, N stays N).
#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other, // N or lowercase pass through
    }
}

/// Reverse complement a DNA slice.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Ring buffer of recently emitted inserts, the template pool for PCR
/// duplicates.
///
/// A duplicate re-sequences a stored fragment: same template, but fresh
/// quality, sequencing errors and polyG artefacts are drawn independently,
/// because dropout is a per-cluster failure. Emitting byte-identical copies
/// would make dedup-before-trim look correct when it is not — a real
/// duplicate diverges from its original in the same places any two reads
/// of one fragment diverge.
///
/// Only *fresh* (non-duplicate) inserts enter the pool, matching the
/// measured semantics: a duplicate copies an original, and the original was
/// itself sequenced once.
pub struct History {
    buf: Vec<Vec<u8>>,
    cap: usize,
    n: usize,
    w: usize,
}

impl History {
    pub fn new(cap: usize) -> Self {
        History {
            buf: Vec::with_capacity(cap.min(1 << 20)),
            cap,
            n: 0,
            w: 0,
        }
    }

    /// Store a fresh insert (overwrites the oldest when full).
    pub fn push(&mut self, insert: Vec<u8>) {
        if self.n < self.cap {
            self.buf.push(insert);
            self.n += 1;
        } else {
            self.buf[self.w] = insert;
        }
        self.w = (self.w + 1) % self.cap;
    }

    /// Sample a stored insert uniformly; None when the pool is empty.
    pub fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Option<&[u8]> {
        if self.n == 0 {
            return None;
        }
        let i = rng.gen_range(0..self.n);
        Some(&self.buf[i])
    }
}

/// A single contig loaded from a reference FASTA.
#[derive(Debug, Clone)]
pub struct Contig {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Minimal FASTA reader: >header line, then sequence lines (no wrapping
/// enforcement beyond "assume plain ACGTN"). Lowercase bases are uppercased.
pub fn read_fasta(path: &std::path::Path) -> Result<Vec<Contig>> {
    let data = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
    let text = String::from_utf8_lossy(&data);
    let mut contigs = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(h) = line.strip_prefix('>') {
            if !name.is_empty() {
                contigs.push(Contig { name, seq });
            }
            name = h.split_whitespace().next().unwrap_or("").to_string();
            seq = Vec::new();
        } else {
            seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !name.is_empty() {
        contigs.push(Contig { name, seq });
    }
    if contigs.is_empty() {
        bail!("no contigs found in {}", path.display());
    }
    Ok(contigs)
}

/// Sample an insert of `len` bases from a random contig/start position.
/// Rejects positions that would run past a contig end (resamples).
pub fn insert_from_reference<R: Rng + ?Sized>(
    contigs: &[Contig],
    len: u32,
    rng: &mut R,
) -> Result<(String, Vec<u8>)> {
    let len = len as usize;
    let mut attempts = 0;
    loop {
        attempts += 1;
        if attempts > 100_000 {
            bail!("could not place insert of length {len} in reference (contigs too short?)");
        }
        let c = &contigs[rng.gen_range(0..contigs.len())];
        if c.seq.len() < len {
            continue;
        }
        let start = rng.gen_range(0..=c.seq.len() - len);
        return Ok((c.name.clone(), c.seq[start..start + len].to_vec()));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(7)
    }

    #[test]
    fn random_dna_gc_matches_request() {
        let mut rng = rng();
        let seq = random_dna(1_000_000, 0.49, &mut rng);
        let gc = seq.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64 / seq.len() as f64;
        assert!((gc - 0.49).abs() < 0.005, "gc {gc}");
    }

    #[test]
    fn random_dna_only_acgt() {
        let mut rng = rng();
        let seq = random_dna(100_000, 0.5, &mut rng);
        assert!(seq.iter().all(|b| BASES.contains(b)));
    }

    #[test]
    fn random_dna_bias_reproduces_located_composition() {
        let mut rng = rng();
        // positions 0,1 nearly never G; position 2 ~100% G (the ERR10016230 shape)
        let head = [
            [0.33, 0.33, 0.01, 0.33],
            [0.33, 0.33, 0.01, 0.33],
            [0.005, 0.005, 0.985, 0.005],
        ];
        let n_reads = 10_000usize;
        let read_len = 8usize;
        let seq: Vec<u8> = (0..n_reads).flat_map(|_| random_dna_bias(read_len, 0.30, Some(&head), None, &mut rng)).collect();
        let g = |i: usize| seq.iter().skip(i).step_by(read_len).filter(|&&b| b == b'G').count() as f64 / n_reads as f64;
        assert!(g(2) > 0.95, "position 2 should be ~G, got {:.3}", g(2));
        assert!(g(0) < 0.05 && g(1) < 0.05, "positions 0/1 should avoid G: {} {}",
                g(0), g(1));
        // position 3 (beyond the head window) returns to uniform GC
        assert!((g(3) - 0.15).abs() < 0.05, "position 3 should be ~gc/2, got {:.3}", g(3));
    }

    #[test]
    fn revcomp_is_involution() {
        let s = b"ACGTNNAACCGGTT";
        let rc = revcomp(s);
        assert_eq!(revcomp(&rc), s);
        assert_eq!(rc, b"AACCGGTTNNACGT");
    }

    #[test]
    fn sample_insert_len_respects_floor() {
        let mut rng = rng();
        for _ in 0..10_000 {
            let l = sample_insert_len(300.0, 50.0, 150, &mut rng);
            assert!(l >= 150, "insert {l} below floor");
        }
    }

    #[test]
    fn fasta_reader_handles_wrapped_lines() {
        let p = std::env::temp_dir().join("fastsim_test.fa");
        std::fs::write(&p, ">chr1 some description\nACGT\nACGT\n>chr2\nNNN\n").unwrap();
        let cs = read_fasta(&p).unwrap();
        std::fs::remove_file(&p).unwrap();
        assert_eq!(cs.len(), 2);
        assert_eq!(cs[0].name, "chr1");
        assert_eq!(cs[0].seq, b"ACGTACGT");
        assert_eq!(cs[1].name, "chr2");
        assert_eq!(cs[1].seq, b"NNN");
    }

    #[test]
    fn insert_from_reference_returns_valid_slice() {
        let mut rng = rng();
        let contigs = vec![
            Contig {
                name: "c1".into(),
                seq: b"AAAAAAAAAACCCCCCCCCC".to_vec(),
            },
            Contig {
                name: "c2".into(),
                seq: b"GGGGGGGGGGTTTTTTTTTT".to_vec(),
            },
        ];
        let (name, seq) = insert_from_reference(&contigs, 10, &mut rng).unwrap();
        assert!(name == "c1" || name == "c2");
        assert_eq!(seq.len(), 10);
        // the 10-mer must be a slice of one of the contigs
        let c = if name == "c1" {
            &contigs[0]
        } else {
            &contigs[1]
        };
        assert!(c.seq.windows(10).any(|w| w == seq.as_slice()));
    }

    #[test]
    fn insert_from_reference_fails_when_contigs_too_short() {
        let mut rng = rng();
        let contigs = vec![Contig {
            name: "tiny".into(),
            seq: b"ACGT".to_vec(),
        }];
        let err = insert_from_reference(&contigs, 10, &mut rng).unwrap_err();
        assert!(err.to_string().contains("could not place"));
    }
}
