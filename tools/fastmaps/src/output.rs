use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

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

pub fn write_svg(counts: &HashMap<(u32, u32), (u64, f64)>, output_dir: &Path) {
    let max_x = counts.keys().map(|(x, _)| *x).max().unwrap() as usize;
    let max_y = counts.keys().map(|(_, y)| *y).max().unwrap() as usize;
    let w = max_x + 1;
    let h = max_y + 1;
    let cell: u32 = 12;
    let margin: u32 = 40;
    let y_axis = margin + h as u32 * cell + 28;
    let view_w = margin + w as u32 * cell + 10 + 50;
    let view_h = y_axis;

    let max_count = counts.values().map(|(c, _)| *c).max().unwrap();

    let f = BufWriter::new(File::create(output_dir.join("heatmap.svg")).unwrap());
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
    f.write_all(hdr.as_bytes()).unwrap();

    for y in 0..h as u32 {
        for x in 0..w as u32 {
            let px = margin + x * cell;
            let py = view_h - (margin + (y + 1) * cell);
            let (r, g, b) = match counts.get(&(x, y)) {
                Some((c, _)) => {
                    let t = (*c as f64).ln() / (max_count as f64).ln();
                    let [r, g, b] = heat_color(t);
                    (r, g, b)
                }
                None => (0, 0, 0),
            };
            let line = format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#{:02X}{:02X}{:02X}\"/>\n",
                px, py, cell, cell, r, g, b
            );
            f.write_all(line.as_bytes()).unwrap();
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
    )
    .unwrap();

    let stops = [("#000", 0), ("#00F", 20), ("#0FF", 40), ("#0F0", 60), ("#FF0", 80), ("#F00", 100)];
    for &(color, pct) in &stops {
        let cy = cb_y + cb_h * pct / 100;
        f.write_all(
            format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"6\" fill=\"{}\"/>\n",
                cb_x, cy, cb_w, color
            )
            .as_bytes(),
        )
        .unwrap();
    }

    for x in (0..=max_x).step_by(5.max(max_x / 10)) {
        let px = margin + x as u32 * cell + cell / 2;
        f.write_all(
            format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"7\" fill=\"#aaa\" text-anchor=\"middle\">{}</text>\n",
                px, y_axis + 8, x
            )
            .as_bytes(),
        )
        .unwrap();
    }
    for y in (0..=max_y).step_by(5.max(max_y / 10)) {
        let py = view_h - (margin + (y as u32 + 1) * cell);
        f.write_all(
            format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"7\" fill=\"#aaa\" text-anchor=\"end\">{}</text>\n",
                margin - 4, py, y
            )
            .as_bytes(),
        )
        .unwrap();
    }

    f.write_all(
        format!(
            "<text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"9\" fill=\"#ccc\" text-anchor=\"middle\">T + C</text>\n",
            margin + w as u32 * cell / 2,
            y_axis + 18
        )
        .as_bytes(),
    )
    .unwrap();
    let y_mid = margin + h as u32 * cell / 2;
    f.write_all(
        format!(
            "<text x=\"8\" y=\"{}\" font-family=\"monospace\" font-size=\"9\" fill=\"#ccc\" text-anchor=\"middle\" transform=\"rotate(-90,8,{})\">T + G</text>\n",
            y_mid, y_mid
        )
        .as_bytes(),
    )
    .unwrap();

    f.write_all(b"</svg>\n").unwrap();
    f.flush().unwrap();

    println!("Saved heatmap.svg  ({}x{} cells, {} unique)", w, h, counts.len());
}
pub fn write_png(counts: &HashMap<(u32, u32), (u64, f64)>, output_dir: &Path) {
    let max_x = counts.keys().map(|(x, _)| *x).max().unwrap() as usize;
    let max_y = counts.keys().map(|(_, y)| *y).max().unwrap() as usize;
    let max_count = counts.values().map(|(c, _)| *c).max().unwrap();

    let w = max_x + 1;
    let h = max_y + 1;
    let mut pixels = vec![[0u8; 3]; w * h];
    for (&(x, y), &(c, _)) in counts {
        let t = (c as f64).ln() / (max_count as f64).ln();
        pixels[y as usize * w + x as usize] = heat_color(t);
    }

    let mut rgb: Vec<u8> = Vec::with_capacity(w * h * 3);
    for row in (0..h).rev() {
        for col in 0..w {
            rgb.extend_from_slice(&pixels[row * w + col]);
        }
    }

    let f = BufWriter::new(File::create(output_dir.join("heatmap.png")).unwrap());
    let mut enc = png::Encoder::new(f, w as u32, h as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header().unwrap().write_image_data(&rgb).unwrap();
    println!("Saved heatmap.png  ({w}×{h} px)");
}

/// Write an interactive Plotly 3-D scatter HTML to `output_dir/heatmap_3d.html`.
/// z = mean Phred quality, colour = log(count).
pub fn write_html(counts: &HashMap<(u32, u32), (u64, f64)>, output_dir: &Path) {
    let max_count = counts.values().map(|(c, _)| *c).max().unwrap() as f64;
    let total: u64 = counts.values().map(|(c, _)| c).sum();

    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let mut zs = Vec::new();
    let mut cs = Vec::new();

    for (&(x, y), &(c, cq)) in counts {
        xs.push(x);
        ys.push(y);
        zs.push(cq / c as f64);
        cs.push((c as f64).ln() / max_count.ln());
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

    let mut f = File::create(output_dir.join("heatmap_3d.html")).unwrap();
    f.write_all(html.as_bytes()).unwrap();
    println!("Saved heatmap_3d.html");
}
