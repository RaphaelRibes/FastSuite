// Scratch probe: run polyx_trim on hand-built reads and print the cut.
use fastprep::algorithm::*;

fn main() {
    // read: G + clean alternating biology (G is a biological first base, h=0)
    // Exact bytes from the CLI repro: 158 bp, Q = 'I' = Phred 40.
    let l = 158usize;
    let mut seq: Vec<u8> = Vec::new();
    while seq.len() < l { seq.extend_from_slice(b"ACGT"); }
    seq.truncate(l);
    seq[0] = b'G';  // biological G at position 1, h = 0
    let q = vec![40u8; l];

    for (mr, mercy, hd) in [(6usize, 3.0, None), (12, 1.0, None), (12, 1.0, Some(1)), (6, 3.0, Some(4))] {
        let mut p = Params::default();
        p.min_run = mr;
        p.mercy = mercy;
        p.head_dropout = hd;
        p.scan_interior = false;
        p.prepare();
        let r = polyx_trim(&seq, &q, &p);
        println!("min_run={mr} mercy={mercy} head_dropout={hd:?} -> {:?}", r);
    }
}
