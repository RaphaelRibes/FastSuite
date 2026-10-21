// Probe: does boundary_anchored remove the mercy-overshoot on clean runs?
use fastprep::algorithm::*;

fn main() {
    // 3' clean run: [non-G flank][run of t G's].
    for t in [6usize, 12, 16, 20] {
        let l = 150usize;
        let mut seq: Vec<u8> = Vec::new();
        while seq.len() < l - t { seq.extend_from_slice(b"ACGT"); }
        seq.truncate(l - t);
        seq.extend(std::iter::repeat_n(b'G', t));
        let q = vec![40u8; l];
        for ba in [false, true] {
            let mut p = Params::default();
            p.min_run = 12;
            p.mercy = 3.0;
            p.scan_interior = false;
            p.boundary_anchored = ba;
            p.prepare();
            let r = polyx_trim(&seq, &q, &p);
            let t3 = match &r { TrimResult::Emit { t3_poly, .. } => *t3_poly, _ => 999 };
            println!("3' clean run t={t} boundary_anchored={ba}: cut_t3={t3} err={}", t3 as i64 - t as i64);
        }
    }
    // 5' clean run: [run of h G's][non-G flank]
    for h in [6usize, 12, 16, 20] {
        let l = 150usize;
        let mut seq: Vec<u8> = vec![b'G'; h];
        while seq.len() < l { seq.extend_from_slice(b"ACGT"); }
        seq.truncate(l);
        let q = vec![40u8; l];
        for ba in [false, true] {
            let mut p = Params::default();
            p.min_run = 12;
            p.mercy = 3.0;
            p.scan_interior = false;
            p.boundary_anchored = ba;
            p.prepare();
            let r = polyx_trim(&seq, &q, &p);
            let t5 = match &r { TrimResult::Emit { t5_poly, .. } => *t5_poly, _ => 999 };
            println!("5' clean run h={h} boundary_anchored={ba}: cut_t5={t5} err={}", t5 as i64 - h as i64);
        }
    }
    // 3' run with an internal substitution error (G-run with one A in it),
    // then clean flank. Boundary-anchored must still include the error base.
    {
        let l = 150usize;
        let t = 14usize;
        let mut seq: Vec<u8> = Vec::new();
        while seq.len() < l - t { seq.extend_from_slice(b"ACGT"); }
        seq.truncate(l - t);
        for i in 0..t { seq.push(if i == 5 { b'A' } else { b'G' }); }
        let q = vec![40u8; l];
        for ba in [false, true] {
            let mut p = Params::default();
            p.min_run = 12; p.mercy = 3.0; p.scan_interior = false; p.boundary_anchored = ba; p.prepare();
            let r = polyx_trim(&seq, &q, &p);
            let t3 = match &r { TrimResult::Emit { t3_poly, .. } => *t3_poly, _ => 999 };
            println!("3' run t={t} w/ internal error boundary_anchored={ba}: cut_t3={t3} err={}", t3 as i64 - t as i64);
        }
    }
}
