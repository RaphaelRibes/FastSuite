use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::worker::PixelAcc;

/// Lowest `t` a *non-zero* count may map to.
///
/// `ln(1) / ln(max)` is 0 for every grid whose busiest pixel holds more than one
/// read, and `heat_color(0)` is `[0, 0, 0]` — the background. A pixel holding
/// exactly one read was therefore painted exactly like an empty cell: on a 3×3
/// plane with three occupied pixels, two of them count 1, the PNG and the SVG
/// both showed 8 of 9 cells black and the reads were invisible.
///
/// A non-zero count is never "no reads", so it is lifted to the first clearly
/// visible step of the scale — dark blue, `#000099` — instead of to the
/// background. Counts keep their `ln(c) / ln(max)` order above the floor, and
/// count 0 stays the background.
const MIN_VISIBLE_T: f64 = 0.1;

/// Normalised log-count `t` in `[0, 1]` for the heat scale: `ln(c) / ln(max)`.
///
/// Two degenerate cases, both of which used to black out real pixels:
///
/// * `max_count <= 1` — every occupied pixel holds exactly one read, so
///   `ln(max) == ln(1) == 0`: the `0/0` NaN blacked out the PNG/SVG and emitted
///   an invalid `color: [NaN]` array into the Plotly page. Every occupied pixel
///   *is* the maximum there, so it maps to the top of the scale.
/// * `count == 1` with `max_count > 1` — `ln(1) == 0`, i.e. the background
///   colour for a pixel that is not empty. Lifted to [`MIN_VISIBLE_T`].
fn log_norm(count: u64, max_count: u64) -> f64 {
    if count == 0 {
        return 0.0;
    }
    if max_count <= 1 {
        return 1.0;
    }
    let t = (count as f64).ln() / (max_count as f64).ln();
    t.max(MIN_VISIBLE_T)
}

/// The PNG encoder's error type is not `io::Error`; every caller here is
/// reporting a failed write either way, so map it once (M8).
fn png_error(e: png::EncodingError) -> std::io::Error {
    std::io::Error::other(e.to_string())
}

fn heat_color(t: f64) -> [u8; 3] {
    const STOPS: &[[f64; 3]] = &[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
    ];
    let n = (STOPS.len() - 1) as f64;
    let pos = t.clamp(0.0, 1.0) * n;
    let i = pos.floor() as usize;
    let f = pos.fract();
    let a = STOPS[i.min(STOPS.len() - 1)];
    let b = STOPS[(i + 1).min(STOPS.len() - 1)];
    [
        ((a[0] + f * (b[0] - a[0])) * 255.0) as u8,
        ((a[1] + f * (b[1] - a[1])) * 255.0) as u8,
        ((a[2] + f * (b[2] - a[2])) * 255.0) as u8,
    ]
}

pub fn write_svg(
    counts: &HashMap<(u32, u32), PixelAcc>,
    output_dir: &Path,
    quiet: bool,
) -> std::io::Result<()> {
    let max_x = counts.keys().map(|(x, _)| *x).max().unwrap() as usize;
    let max_y = counts.keys().map(|(_, y)| *y).max().unwrap() as usize;
    let w = max_x + 1;
    let h = max_y + 1;
    let cell: u32 = 12;
    let margin: u32 = 40;
    let y_axis = margin + h as u32 * cell + 28;
    let view_w = margin + w as u32 * cell + 10 + 50;
    let view_h = y_axis;

    let max_count = counts.values().map(|a| a.count).max().unwrap();

    let f = BufWriter::new(File::create(output_dir.join("heatmap.svg"))?);
    let mut f = f;

    let hdr = format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" width=\"{}\" height=\"{}\">\n\
         <rect width=\"{}\" height=\"{}\" fill=\"#000000\"/>\n\
         <text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"10\" fill=\"#ccc\">T+C x T+G log-count heatmap</text>\n",
        view_w, view_h, view_w, view_h,
        view_w, view_h,
        margin, margin
    );
    f.write_all(hdr.as_bytes())?;

    for y in 0..h as u32 {
        for x in 0..w as u32 {
            let px = margin + x * cell;
            let py = view_h - (margin + (y + 1) * cell);
            let (r, g, b) = match counts.get(&(x, y)) {
                Some(acc) => {
                    let t = log_norm(acc.count, max_count);
                    let [r, g, b] = heat_color(t);
                    (r, g, b)
                }
                None => (0, 0, 0),
            };
            let line = format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#{:02X}{:02X}{:02X}\"/>\n",
                px, py, cell, cell, r, g, b
            );
            f.write_all(line.as_bytes())?;
        }
    }

    let cb_x = margin + w as u32 * cell + 10;
    let cb_y = margin;
    let cb_w = 12;
    let cb_h = h as u32 * cell;
    f.write_all(
        format!(
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#808080\" opacity=\"0.25\"/>\n",
            cb_x, cb_y, cb_w, cb_h
        )
        .as_bytes(),
    )?;

    let stops = [("#000", 0), ("#00F", 20), ("#0FF", 40), ("#0F0", 60), ("#FF0", 80), ("#F00", 100)];
    for &(color, pct) in &stops {
        let cy = cb_y + cb_h * pct / 100;
        f.write_all(
            format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"6\" fill=\"{}\"/>\n",
                cb_x, cy, cb_w, color
            )
            .as_bytes(),
        )?;
    }

    for x in (0..=max_x).step_by(5.max(max_x / 10)) {
        let px = margin + x as u32 * cell + cell / 2;
        f.write_all(
            format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"7\" fill=\"#aaa\" text-anchor=\"middle\">{}</text>\n",
                px, y_axis + 8, x
            )
            .as_bytes(),
        )?;
    }
    for y in (0..=max_y).step_by(5.max(max_y / 10)) {
        let py = view_h - (margin + (y as u32 + 1) * cell);
        f.write_all(
            format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"7\" fill=\"#aaa\" text-anchor=\"end\">{}</text>\n",
                margin - 4, py, y
            )
            .as_bytes(),
        )?;
    }

    f.write_all(
        format!(
            "<text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"9\" fill=\"#ccc\" text-anchor=\"middle\">T + C</text>\n",
            margin + w as u32 * cell / 2,
            y_axis + 18
        )
        .as_bytes(),
    )?;
    let y_mid = margin + h as u32 * cell / 2;
    f.write_all(
        format!(
            "<text x=\"8\" y=\"{}\" font-family=\"monospace\" font-size=\"9\" fill=\"#ccc\" text-anchor=\"middle\" transform=\"rotate(-90,8,{})\">T + G</text>\n",
            y_mid, y_mid
        )
        .as_bytes(),
    )?;

    f.write_all(b"</svg>\n")?;
    f.flush()?;

    if !quiet {
        println!("Saved heatmap.svg  ({}x{} cells, {} unique)", w, h, counts.len());
    }
    Ok(())
}
pub fn write_png(
    counts: &HashMap<(u32, u32), PixelAcc>,
    output_dir: &Path,
    quiet: bool,
) -> std::io::Result<()> {
    let max_x = counts.keys().map(|(x, _)| *x).max().unwrap() as usize;
    let max_y = counts.keys().map(|(_, y)| *y).max().unwrap() as usize;
    let max_count = counts.values().map(|a| a.count).max().unwrap();

    let w = max_x + 1;
    let h = max_y + 1;
    let mut pixels = vec![[0u8; 3]; w * h];
    for (&(x, y), acc) in counts {
        let t = log_norm(acc.count, max_count);
        pixels[y as usize * w + x as usize] = heat_color(t);
    }

    let mut rgb: Vec<u8> = Vec::with_capacity(w * h * 3);
    for row in (0..h).rev() {
        for col in 0..w {
            rgb.extend_from_slice(&pixels[row * w + col]);
        }
    }

    let f = BufWriter::new(File::create(output_dir.join("heatmap.png"))?);
    let mut enc = png::Encoder::new(f, w as u32, h as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header().map_err(png_error)?;
    writer.write_image_data(&rgb).map_err(png_error)?;
    if !quiet {
        println!("Saved heatmap.png  ({w}×{h} px)");
    }
    Ok(())
}

/// Write an interactive Plotly 3-D scatter HTML to `output_dir/heatmap_3d.html`.
/// z = mean Phred quality, colour = log(count).
///
/// The series is sorted by `(x, y)` (FM-12). It used to be emitted in `HashMap`
/// iteration order, which std seeds per process — so two runs of the identical
/// command produced different files even at `-j 1`. Together with the integer
/// quality accumulators this makes the page reproducible.
pub fn write_html(
    counts: &HashMap<(u32, u32), PixelAcc>,
    output_dir: &Path,
    quiet: bool,
) -> std::io::Result<()> {
    let max_count = counts.values().map(|a| a.count).max().unwrap();
    let total: u64 = counts.values().map(|a| a.count).sum();

    let mut series: Vec<((u32, u32), PixelAcc)> =
        counts.iter().map(|(&xy, &acc)| (xy, acc)).collect();
    series.sort_unstable_by_key(|&(xy, _)| xy);

    let mut xs = Vec::with_capacity(series.len());
    let mut ys = Vec::with_capacity(series.len());
    let mut zs = Vec::with_capacity(series.len());
    let mut cs = Vec::with_capacity(series.len());

    for ((x, y), acc) in series {
        xs.push(x);
        ys.push(y);
        zs.push(acc.mean_quality());
        cs.push(log_norm(acc.count, max_count));
    }

    let fmt_vec = |v: &[f64]| v.iter().map(|n| format!("{n:.4}")).collect::<Vec<_>>().join(",");
    let fmt_u32 = |v: &[u32]| v.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(",");

    let html = format!(
        r#"<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>FASTQ heatmap 3D</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body {{ margin:0; background:#0e1117; font-family:sans-serif; color:#eee; }}
  #info {{ padding:12px 20px; font-size:13px; color:#aaa; }}
  #plot {{ width:100vw; height:calc(100vh - 40px); }}
</style>
</head>
<body>
<div id="info">{} unique coordinates &nbsp;|&nbsp; {} total reads &nbsp;|&nbsp; z = mean Phred quality &nbsp;|&nbsp; colour = log(count)</div>
<div id="plot"></div>
<script>
const trace = {{
  type: 'scatter3d', mode: 'markers',
  x: [{}], y: [{}], z: [{}],
  marker: {{
    size: 4, color: [{}], colorscale: 'Jet', showscale: true,
    colorbar: {{ title: 'log count', thickness: 14, len: 0.6 }}, opacity: 0.85
  }},
  hovertemplate: 'x: %{{x}}<br>y: %{{y}}<br>mean qual: %{{z:.1f}}<extra></extra>'
}};
const layout = {{
  paper_bgcolor: '#0e1117',
  scene: {{
    bgcolor: '#0e1117',
    xaxis: {{ title: 'x (T+C)', color: '#aaa', gridcolor: '#333', dtick: 5, tick0: 0 }},
    yaxis: {{ title: 'y (T+G)', color: '#aaa', gridcolor: '#333', dtick: 5, tick0: 0 }},
    zaxis: {{ title: 'mean Phred quality', color: '#aaa', gridcolor: '#333' }},
  }},
  margin: {{ l:0, r:0, t:0, b:0 }}
}};
Plotly.newPlot('plot', [trace], layout, {{responsive: true}});
</script>
</body>
</html>
"#,
        counts.len(), total,
        fmt_u32(&xs), fmt_u32(&ys),
        fmt_vec(&zs), fmt_vec(&cs),
    );

    let mut f = File::create(output_dir.join("heatmap_3d.html"))?;
    f.write_all(html.as_bytes())?;
    if !quiet {
        println!("Saved heatmap_3d.html");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression: a grid whose busiest pixel holds exactly one read used to
    /// divide by `ln(1) == 0`, yielding `0/0 == NaN` — a uniformly black PNG /
    /// SVG and an invalid `color: [NaN]` array in the Plotly page.
    ///
    /// `log_norm` special-cases `max_count <= 1` to `1.0` (the occupied pixels
    /// are all the maximum), so the value is finite and maps to a visible
    /// colour instead of black.
    #[test]
    fn single_read_grid_normalises_finite() {
        let max_count = 1u64; // every occupied pixel holds exactly one read
        let t = log_norm(1, max_count);
        assert!(t.is_finite(), "normalised log-count must be finite, got {t}");
        assert_eq!(t, 1.0, "single-read pixel is the maximum → top of scale");
        assert_ne!(
            heat_color(t),
            [0, 0, 0],
            "the single-read pixel must not render black"
        );
    }

    /// The fix must not disturb the ordinary case: with a real `max > 1` the
    /// normalisation stays `ln(c) / ln(max)` (max → 1) for every count whose
    /// value is above the visibility floor.
    #[test]
    fn log_norm_unchanged_for_normal_grids() {
        assert!((log_norm(10, 10) - 1.0).abs() < 1e-12);
        assert!((log_norm(5, 100) - (5f64.ln() / 100f64.ln())).abs() < 1e-12);
        assert!((log_norm(50, 100) - (50f64.ln() / 100f64.ln())).abs() < 1e-12);
        let mid = log_norm(5, 100);
        assert!((mid - (5f64.ln() / 100f64.ln())).abs() < 1e-12);
    }

    /// `ln(1)` is 0, and 0 is the background: a pixel holding one read used to
    /// be painted exactly like an empty one whenever the busiest pixel held more
    /// than one. Count 0 is still the background.
    #[test]
    fn count_one_is_lifted_off_the_background() {
        for max in [2u64, 3, 10, 1_000, 1_000_000] {
            let t = log_norm(1, max);
            assert!(t >= MIN_VISIBLE_T, "count 1 on max {max} gave {t}");
            assert_ne!(heat_color(t), [0, 0, 0], "count 1 on max {max} is black");
            // The scale still orders the counts around it.
            assert!(t < log_norm(max, max), "count 1 must stay below the maximum");
        }
        assert_eq!(log_norm(0, 1_000), 0.0, "no reads is the background");
    }

    /// FM-12: the Plotly series is sorted by `(x, y)`, not emitted in `HashMap`
    /// iteration order (which is randomised per process, so the file used to
    /// change between two runs of the same command).
    #[test]
    fn html_series_is_sorted_by_xy() {
        let dir = std::env::temp_dir().join(format!("maps_html_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let mut counts: HashMap<(u32, u32), PixelAcc> = HashMap::new();
        // Inserted out of order on purpose.
        for (x, y, c) in [(9u32, 1u32, 4u64), (2, 7, 2), (2, 3, 7), (5, 5, 1)] {
            counts.insert(
                (x, y),
                PixelAcc { count: c, sum_q: 30 * c, sum_len: c },
            );
        }
        write_html(&counts, &dir, true).unwrap();
        let text = std::fs::read_to_string(dir.join("heatmap_3d.html")).unwrap();

        let grab = |open: &str, close: &str| -> Vec<u32> {
            let start = text.find(open).unwrap_or_else(|| panic!("no {open}")) + open.len();
            let end = start + text[start..].find(close).unwrap();
            text[start..end]
                .split(',')
                .map(|s| s.trim().parse::<u32>().unwrap())
                .collect()
        };
        let xs = grab("x: [", "], y: [");
        let ys = grab("y: [", "], z: [");
        assert_eq!(xs, vec![2, 2, 5, 9]);
        assert_eq!(ys, vec![3, 7, 5, 1]);
        assert!(
            (1..xs.len()).all(|i| (xs[i - 1], ys[i - 1]) < (xs[i], ys[i])),
            "series must be strictly ordered by (x, y): {xs:?} {ys:?}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }
}
