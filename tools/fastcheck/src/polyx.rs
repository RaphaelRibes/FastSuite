//! polyX run-length measurement: the leading (5′) and trailing (3′) homopolymer
//! run of each read, for every base. This is the generalisation of the polyG
//! artefact (NovaSeq 2-colour 3′ polyG, and the first-cycle 5′ G dropout) to
//! all four bases, so a polyA/polyC/polyT tail or head is caught too.

/// Base code: 0=A, 1=C, 2=G, 3=T. `usize::MAX` for anything else (N, IUPAC).
pub fn base_code(b: u8) -> usize {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => usize::MAX,
    }
}

/// The leading and trailing homopolymer run of a read.
///
/// Returns `(leading_base, leading_len, trailing_base, trailing_len)`. A run
/// stops at the first non-matching base; an N or other IUPAC code ends the run
/// (it cannot be part of a homopolymer). A read with no ACGT at an end has a
/// run of length 0 with `base = usize::MAX`.
pub fn measure_polyx(seq: &[u8]) -> (usize, usize, usize, usize) {
    let mut lb = usize::MAX;
    let mut ll = 0usize;
    for &b in seq {
        let c = base_code(b);
        if c == usize::MAX {
            break;
        }
        if lb == usize::MAX {
            lb = c;
            ll = 1;
        } else if c == lb {
            ll += 1;
        } else {
            break;
        }
    }

    let mut tb = usize::MAX;
    let mut tl = 0usize;
    for &b in seq.iter().rev() {
        let c = base_code(b);
        if c == usize::MAX {
            break;
        }
        if tb == usize::MAX {
            tb = c;
            tl = 1;
        } else if c == tb {
            tl += 1;
        } else {
            break;
        }
    }
    (lb, ll, tb, tl)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trailing_g_run() {
        let (lb, ll, tb, tl) = measure_polyx(b"ACGTGGG");
        assert_eq!(lb, 0); // A
        assert_eq!(ll, 1);
        assert_eq!(tb, 2); // G
        assert_eq!(tl, 3);
    }

    #[test]
    fn leading_g_run() {
        let (lb, ll, tb, tl) = measure_polyx(b"GGGACGT");
        assert_eq!(lb, 2);
        assert_eq!(ll, 3);
        assert_eq!(tb, 3); // T
        assert_eq!(tl, 1);
    }

    #[test]
    fn leading_a_run() {
        let (lb, ll, _, _) = measure_polyx(b"AAAACGT");
        assert_eq!(lb, 0);
        assert_eq!(ll, 4);
    }

    #[test]
    fn no_run() {
        let (lb, ll, tb, tl) = measure_polyx(b"ACGT");
        assert_eq!(lb, 0);
        assert_eq!(ll, 1);
        assert_eq!(tb, 3);
        assert_eq!(tl, 1);
    }

    #[test]
    fn all_same_base_is_full_length() {
        let (lb, ll, tb, tl) = measure_polyx(b"GGGGGG");
        assert_eq!(lb, 2);
        assert_eq!(ll, 6);
        assert_eq!(tb, 2);
        assert_eq!(tl, 6);
    }

    #[test]
    fn n_ends_the_run() {
        let (lb, ll, _, _) = measure_polyx(b"GGGNACGT");
        assert_eq!(lb, 2);
        assert_eq!(ll, 3, "an N cannot continue a homopolymer run");
    }

    #[test]
    fn empty_read_has_no_run() {
        let (lb, ll, tb, tl) = measure_polyx(b"");
        assert_eq!(lb, usize::MAX);
        assert_eq!(ll, 0);
        assert_eq!(tb, usize::MAX);
        assert_eq!(tl, 0);
    }
}
