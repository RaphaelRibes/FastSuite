//! Every occupied pixel must be visibly distinct from an empty one.
//!
//! `log_norm` special-cased `max_count <= 1` (the `ln(1) == 0` divide) but not
//! `count == 1` on a grid whose busiest pixel holds more than one read: there
//! `ln(1) / ln(max)` is 0, `heat_color(0)` is the black background, and a pixel
//! holding exactly one read was painted exactly like an empty cell — 8 of 9
//! cell rects `#000000` on a 3×3 grid with three occupied pixels, two of them
//! count 1. The PNG and the SVG both lost those reads.
//!
//! Both are decoded here (the PNG through its filters) and compared against the
//! background and against each other.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_heat_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const BACKGROUND: [u8; 3] = [0, 0, 0];

/// Minimal PNG reader for what the writer emits: 8-bit RGB, no interlace, all
/// five filter types.
fn decode_png(path: &Path) -> (usize, usize, Vec<[u8; 3]>) {
    let raw = std::fs::read(path).unwrap();
    assert_eq!(&raw[..8], b"\x89PNG\r\n\x1a\n", "PNG signature");
    let mut pos = 8usize;
    let (mut w, mut h) = (0usize, 0usize);
    let mut idat: Vec<u8> = Vec::new();
    while pos + 8 <= raw.len() {
        let len = u32::from_be_bytes(raw[pos..pos + 4].try_into().unwrap()) as usize;
        let typ = &raw[pos + 4..pos + 8];
        let data = &raw[pos + 8..pos + 8 + len];
        match typ {
            b"IHDR" => {
                w = u32::from_be_bytes(data[0..4].try_into().unwrap()) as usize;
                h = u32::from_be_bytes(data[4..8].try_into().unwrap()) as usize;
                assert_eq!(data[8], 8, "bit depth");
                assert_eq!(data[9], 2, "colour type: RGB");
                assert_eq!(data[12], 0, "no interlace");
            }
            b"IDAT" => idat.extend_from_slice(data),
            _ => {}
        }
        pos += 12 + len;
    }

    let mut px = Vec::new();
    flate2::read::ZlibDecoder::new(&idat[..])
        .read_to_end(&mut px)
        .expect("IDAT must inflate");

    let bpp = 3usize;
    let stride = w * bpp;
    let mut rows: Vec<Vec<u8>> = Vec::with_capacity(h);
    let mut prev = vec![0u8; stride];
    let mut i = 0usize;
    for _ in 0..h {
        let filter = px[i];
        i += 1;
        let mut line = px[i..i + stride].to_vec();
        i += stride;
        match filter {
            0 => {}
            1 => {
                for x in bpp..stride {
                    line[x] = line[x].wrapping_add(line[x - bpp]);
                }
            }
            2 => {
                for x in 0..stride {
                    line[x] = line[x].wrapping_add(prev[x]);
                }
            }
            3 => {
                for x in 0..stride {
                    let left = if x >= bpp { line[x - bpp] as u16 } else { 0 };
                    line[x] = line[x].wrapping_add(((left + prev[x] as u16) / 2) as u8);
                }
            }
            4 => {
                for x in 0..stride {
                    let a = if x >= bpp { line[x - bpp] as i16 } else { 0 };
                    let b = prev[x] as i16;
                    let c = if x >= bpp { prev[x - bpp] as i16 } else { 0 };
                    let p = a + b - c;
                    let (pa, pb, pc) = ((p - a).abs(), (p - b).abs(), (p - c).abs());
                    let pred = if pa <= pb && pa <= pc { a } else if pb <= pc { b } else { c };
                    line[x] = line[x].wrapping_add(pred as u8);
                }
            }
            other => panic!("unexpected PNG filter type {other}"),
        }
        rows.push(line.clone());
        prev = line;
    }

    // The writer emits grid row y at PNG scanline h-1-y.
    let mut out = vec![BACKGROUND; w * h];
    for y in 0..h {
        let row = &rows[h - 1 - y];
        for x in 0..w {
            out[y * w + x] = [row[x * bpp], row[x * bpp + 1], row[x * bpp + 2]];
        }
    }
    (w, h, out)
}

/// Cell fills from `heatmap.svg`, indexed by `y * w + x`.
fn decode_svg(path: &Path) -> Vec<[u8; 3]> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut cells = Vec::new();
    for line in text.lines() {
        // Cell rects are the only ones 12x12; the colour bar uses height 6.
        if !line.starts_with("<rect ") || !line.contains("width=\"12\" height=\"12\"") {
            continue;
        }
        let fill = line
            .split("fill=\"#")
            .nth(1)
            .and_then(|s| s.strip_suffix("\"/>"))
            .unwrap_or_else(|| panic!("cell rect without a fill: {line}"));
        let byte = |i: usize| u8::from_str_radix(&fill[i..i + 2], 16).unwrap();
        cells.push([byte(0), byte(2), byte(4)]);
    }
    cells
}

/// One run of `n` reads of `seq`; returns (grid width, height, PNG, SVG).
fn run_case(tag: &str, reads: &[(&str, usize)]) -> (usize, usize, Vec<[u8; 3]>, Vec<[u8; 3]>) {
    let dir = tmpdir(tag);
    let fq = dir.join("p_R1.fastq");
    let mut body = String::new();
    for (seq, n) in reads {
        for i in 0..*n {
            body.push_str(&format!("@{seq}{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
        }
    }
    std::fs::write(&fq, body).unwrap();
    let out = dir.join("out");
    let o = Command::new(bin())
        .args(["-o", out.to_str().unwrap(), "-j", "1", "--force", "--svg"])
        .arg(&fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "{tag}: run failed\nstderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let (w, h, png) = decode_png(&out.join("heatmap.png"));
    let svg = decode_svg(&out.join("heatmap.svg"));
    assert_eq!(svg.len(), w * h, "one 12x12 rect per cell");
    (w, h, png, svg)
}

fn visible(tag: &str, rgb: [u8; 3]) {
    assert_ne!(rgb, BACKGROUND, "{tag}: painted exactly like an empty cell");
    assert!(
        rgb.iter().copied().max().unwrap() >= 32,
        "{tag}: {rgb:?} is too dark to read as occupied"
    );
}

/// Three occupied pixels on a 3×3 plane, two of them holding one read and one
/// holding two: every occupied cell must be visible, in the PNG and the SVG.
#[test]
fn count_one_is_not_the_background() {
    // "AA" -> (0,0), "AT" -> (1,1), "TT" -> (2,2).
    let (w, h, png, svg) = run_case("one", &[("AA", 1), ("AT", 1), ("TT", 2)]);
    assert_eq!((w, h), (3, 3));

    let at = |cells: &[[u8; 3]], x: usize, y: usize| cells[y * w + x];
    for (x, y) in [(0, 0), (1, 1), (2, 2)] {
        visible("png", at(&png, x, y));
        visible("svg", at(&svg, x, y));
        assert_eq!(
            at(&png, x, y),
            at(&svg, x, y),
            "the PNG and the SVG must paint ({x},{y}) the same"
        );
    }
    // The two count-1 pixels share a colour, and the count-2 pixel differs.
    assert_eq!(at(&png, 0, 0), at(&png, 1, 1), "equal counts, equal colour");
    assert_ne!(
        at(&png, 0, 0),
        at(&png, 2, 2),
        "a busier pixel must still be distinguishable"
    );
    // Empty cells stay the background: this is not a "paint everything" fix.
    for (x, y) in [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)] {
        assert_eq!(at(&png, x, y), BACKGROUND, "empty ({x},{y}) must stay black");
        assert_eq!(at(&svg, x, y), BACKGROUND, "empty ({x},{y}) must stay black");
    }
}

/// The floor also holds far from the top of the scale: `ln(1)/ln(1000)` is 0
/// too, so a single read on a busy plane used to be black.
#[test]
fn count_one_stays_visible_on_a_wide_scale() {
    let (w, _h, png, svg) = run_case("wide", &[("AA", 1), ("TT", 1000)]);
    visible("png", png[0]);
    visible("svg", svg[0]);
    assert_eq!(png[0], svg[0]);
    // The busiest pixel is the top of the scale, and still distinct.
    let max = png[2 * w + 2];
    visible("png", max);
    assert_ne!(png[0], max, "1 read and 1000 reads must not share a colour");
}
