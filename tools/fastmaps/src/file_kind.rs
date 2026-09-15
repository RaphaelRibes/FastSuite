//! Classify input FASTQ files as R1 / R2 / Unknown by glob-pattern matching
//! against their basenames.

use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FileKind {
    Unknown,
    R1,
    R2,
}

impl FileKind {
    /// JSON-friendly tag.
    pub fn as_str(self) -> &'static str {
        match self {
            FileKind::Unknown => "unknown",
            FileKind::R1 => "R1",
            FileKind::R2 => "R2",
        }
    }
}

/// Default glob patterns matching Illumina-style paired-end naming.
/// Comparisons are case-insensitive.
pub fn default_r1_patterns() -> Vec<String> {
    vec!["*R1*".into()]
}

pub fn default_r2_patterns() -> Vec<String> {
    vec!["*R2*".into()]
}

/// Match `name` against a glob pattern using `*` (any sequence) and `?` (one
/// character) wildcards. Other characters are literal. Both inputs are
/// expected to be already lowercased.
fn glob_match(pattern: &[u8], name: &[u8]) -> bool {
    // Iterative two-pointer matcher with backtracking on `*`.
    // O(len(name) × number_of_stars) worst case — fine for filenames.
    let (mut pi, mut ni) = (0, 0);
    let (mut star_p, mut star_n): (Option<usize>, usize) = (None, 0);

    while ni < name.len() {
        if pi < pattern.len() && (pattern[pi] == b'?' || pattern[pi] == name[ni]) {
            pi += 1;
            ni += 1;
        } else if pi < pattern.len() && pattern[pi] == b'*' {
            star_p = Some(pi);
            star_n = ni;
            pi += 1;
        } else if let Some(sp) = star_p {
            // Backtrack: consume one more char with the last `*`.
            pi = sp + 1;
            star_n += 1;
            ni = star_n;
        } else {
            return false;
        }
    }
    // Trailing pattern characters must all be `*`.
    while pi < pattern.len() && pattern[pi] == b'*' {
        pi += 1;
    }
    pi == pattern.len()
}

/// Classify one filename. R1 patterns are checked first; ties go to R1.
pub fn classify(path: &str, r1: &[String], r2: &[String]) -> FileKind {
    let name = Path::new(path)
        .file_name()
        .map(|s| s.to_string_lossy().to_lowercase())
        .unwrap_or_default();
    let name_bytes = name.as_bytes();

    let any_match = |pats: &[String]| -> bool {
        pats.iter()
            .any(|p| glob_match(p.to_lowercase().as_bytes(), name_bytes))
    };

    if any_match(r1) {
        return FileKind::R1;
    }
    if any_match(r2) {
        return FileKind::R2;
    }
    FileKind::Unknown
}

/// Classify every file and report a summary to stderr.
///
/// `customised` indicates whether the user supplied their own patterns.
/// When the autodetect (default patterns) fails for everything, a clear
/// warning is printed.
pub fn classify_all(
    files: &[String],
    r1: &[String],
    r2: &[String],
    customised: bool,
) -> Vec<FileKind> {
    let kinds: Vec<FileKind> = files.iter().map(|f| classify(f, r1, r2)).collect();
    let n_r1 = kinds.iter().filter(|k| **k == FileKind::R1).count();
    let n_r2 = kinds.iter().filter(|k| **k == FileKind::R2).count();
    let n_un = kinds.iter().filter(|k| **k == FileKind::Unknown).count();
    eprintln!("Classified files: {n_r1} R1, {n_r2} R2, {n_un} unknown");

    if n_r1 == 0 && n_r2 == 0 && !customised {
        eprintln!(
            "\n⚠  Autodetection of R1/R2 found NO matches in any of the {} input files.\n   \
             All reads will be aggregated into the combined count only;\n   \
             the R1/R2 split view in the notebook will be empty.\n   \
             If your filenames use a different convention, override with\n   \
             --r1-pattern '<glob>' and --r2-pattern '<glob>' (e.g. '*forward*').\n",
            files.len()
        );
    } else if n_un > 0 {
        eprintln!("Files classified as unknown (counted only in combined view):");
        for (i, k) in kinds.iter().enumerate() {
            if *k == FileKind::Unknown {
                eprintln!("  {}", files[i]);
            }
        }
    }

    kinds
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r1() -> Vec<String> { default_r1_patterns() }
    fn r2() -> Vec<String> { default_r2_patterns() }

    #[test]
    fn glob_basic() {
        assert!(glob_match(b"*r1*", b"e375_seq2_r1_trimmed.fastq"));
        assert!(glob_match(b"*r2*", b"e375_seq2_r2_trimmed.fastq"));
        assert!(!glob_match(b"*r1*", b"e375_seq2_r2_trimmed.fastq"));
    }

    #[test]
    fn glob_question_mark() {
        assert!(glob_match(b"sample?.fastq", b"sample1.fastq"));
        assert!(!glob_match(b"sample?.fastq", b"sample10.fastq"));
    }

    #[test]
    fn glob_anchored() {
        assert!(glob_match(b"e*r1*.fastq", b"e375_seq2_r1_trimmed.fastq"));
        assert!(!glob_match(b"e*r1*.fastq", b"e375_seq2_r1_trimmed.fq"));
    }

    #[test]
    fn glob_trailing_star() {
        assert!(glob_match(b"prefix*", b"prefix_anything"));
        assert!(glob_match(b"prefix*", b"prefix"));
    }

    #[test]
    fn matches_r1_underscore() {
        assert_eq!(classify("E375_seq2_R1_trimmed.fastq", &r1(), &r2()), FileKind::R1);
    }

    #[test]
    fn matches_r2_underscore() {
        assert_eq!(classify("E375_seq2_R2_trimmed.fastq", &r1(), &r2()), FileKind::R2);
    }

    #[test]
    fn matches_r1_dot_suffix() {
        assert_eq!(classify("sample_R1.fastq.gz", &r1(), &r2()), FileKind::R1);
    }

    #[test]
    fn case_insensitive() {
        assert_eq!(classify("sample_r1_001.fastq", &r1(), &r2()), FileKind::R1);
    }

    #[test]
    fn unknown_when_no_match() {
        assert_eq!(classify("just_a_file.fastq", &r1(), &r2()), FileKind::Unknown);
    }

    #[test]
    fn prefers_r1_on_collision() {
        let r1 = vec!["*weird*".into()];
        let r2 = vec!["*weird*".into()];
        assert_eq!(classify("weird_file.fastq", &r1, &r2), FileKind::R1);
    }

    #[test]
    fn custom_patterns_override_defaults() {
        let r1 = vec!["*forward*".into()];
        let r2 = vec!["*reverse*".into()];
        assert_eq!(classify("sample_forward.fastq", &r1, &r2), FileKind::R1);
        assert_eq!(classify("sample_reverse.fastq", &r1, &r2), FileKind::R2);
        assert_eq!(classify("sample_R1.fastq", &r1, &r2), FileKind::Unknown);
    }
}