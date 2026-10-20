// Chart helpers drawn with cetz primitives — no plotting package needed.
//
// Verified on Typst 0.15.1 + @preview/cetz:0.5.2 by rendering to PNG and measuring
// ink against the page margins (A4, 1.5 cm margins -> 0 ink outside 1.5..19.5 cm).
//
// Rules baked in here, each of which cost real debugging time to find:
//   * `draw.` qualifies every primitive — a star import shadows Typst's own
//     rotate/text/rect and makes rotate() panic with cetz's coordinate rotator;
//   * bar width follows the DATA's x spacing, so a histogram of bases and one of
//     scores work with the same helper (a fixed width draws one giant bar);
//   * the value axis gets numbers, otherwise the chart is unreadable and reads
//     as a clipping bug;
//   * canvas units are centimetres (~1.02 cm per unit); size against the text
//     width, not the page width.
//
// Copy next to the document and `#import "cetz-charts.typ": barchart, histogram,
// groupedbars, scatter`.

#import "@preview/cetz:0.5.2": canvas, draw

#let axis() = 0.5pt
#let label-col = rgb("#555")
#let tick-col = rgb("#666")

/// data: array of (label, value) or (label, value, colour). `ymax` fixes the
/// top of the value axis; without it the axis is padded 18% above the data,
/// which for a percentage chart labels the top tick 118.
#let barchart(data, w: 9.0, h: 3.6, ylabel: none, ymax: none,
              fmt: v => str(calc.round(v, digits: 2))) = canvas({
  let (pad-l, pad-b) = (1.1, 0.9)
  let vals = data.map(d => d.at(1))
  let hi = calc.max(..vals, 0.0)
  let lo = calc.min(..vals, 0.0)
  let span = if hi - lo == 0 { 1.0 } else { hi - lo }
  let hi2 = if ymax != none { ymax } else { hi + span * 0.18 }
  // pad below only when the data really goes negative: padding an all-positive
  // chart pushes the baseline under the axis and the bars float above it
  let lo2 = if lo < 0 { lo - span * 0.12 } else { 0 }
  let (x0, y0) = (pad-l, pad-b)
  let (x1, y1) = (pad-l + w, pad-b + h)
  let n = data.len()
  let bw = w / n * 0.62
  draw.line((x0, y0), (x1, y0), stroke: axis())
  draw.line((x0, y0), (x0, y1), stroke: axis())
  let yz = if lo < 0 { y0 + (0 - lo2) / (hi2 - lo2) * h } else { y0 }
  if lo < 0 {
    draw.line((x0, yz), (x1, yz), stroke: (paint: rgb("#9aa4ad"), thickness: 0.4pt, dash: "dashed"))
  }
  for (i, d) in data.enumerate() {
    let (lab, val) = (d.at(0), d.at(1))
    let col = if d.len() > 2 { d.at(2) } else { rgb("#2f6f9f") }
    let cx = x0 + w / n * (i + 0.5)
    let yv = y0 + (val - lo2) / (hi2 - lo2) * h
    draw.rect((cx - bw / 2, calc.min(yv, yz)), (cx + bw / 2, calc.max(yv, yz)),
              fill: col, stroke: none, radius: 0.03)
    draw.content((cx, yv + 0.16), text(size: 6pt, fill: rgb("#333"), fmt(val)))
    draw.content((cx, y0 - 0.28), text(size: 6pt, fill: label-col)[#lab])
  }
  draw.content((x0 - 0.14, y1), anchor: "east", text(size: 5.5pt, fill: tick-col)[#fmt(hi2)])
  draw.content((x0 - 0.14, y0), anchor: "east", text(size: 5.5pt, fill: tick-col)[#fmt(lo2)])
  if lo < 0 {
    draw.content((x0 - 0.14, yz), anchor: "east", text(size: 5.5pt, fill: tick-col)[0])
  }
  if ylabel != none {
    draw.content((x0 - 0.42, y1 + 0.30), anchor: "west", text(size: 6pt, fill: label-col)[#ylabel])
  }
})

/// groups: array of labels; series: array of (name, colour, array-of-values).
/// `ymax` fixes the top of the value axis; without it the axis is padded 18%
/// above the data, which for a percentage chart labels the top tick 118.
#let groupedbars(groups, series, w: 9.0, h: 3.6, ymax: none,
                 fmt: v => str(calc.round(v, digits: 3))) = canvas({
  let all = series.map(s => s.at(2)).flatten()
  let hi = if ymax != none { ymax } else { calc.max(..all, 0.0) * 1.18 }
  let lo = calc.min(..all, 0.0)
  let (x0, y0) = (1.1, 0.9)
  let (x1, y1) = (x0 + w, y0 + h)
  draw.line((x0, y0), (x1, y0), stroke: axis())
  draw.line((x0, y0), (x0, y1), stroke: axis())
  let ng = groups.len()
  let ns = series.len()
  let gw = w / ng
  let bw = gw * 0.8 / ns
  for (g, gl) in groups.enumerate() {
    for (s, sr) in series.enumerate() {
      let val = sr.at(2).at(g)
      let bx = x0 + gw * g + gw * 0.1 + bw * s
      let yv = y0 + (val - lo) / (hi - lo) * h
      draw.rect((bx, y0), (bx + bw * 0.9, yv), fill: sr.at(1), stroke: none, radius: 0.02)
      draw.content((bx + bw * 0.45, yv + 0.14), text(size: 5.2pt, fill: rgb("#333"), fmt(val)))
    }
    draw.content((x0 + gw * (g + 0.5), y0 - 0.3), text(size: 6pt)[#gl])
  }
  draw.content((x0 - 0.12, y1), anchor: "east", text(size: 5.5pt, fill: tick-col)[#fmt(hi)])
  draw.content((x0 - 0.12, y0), anchor: "east", text(size: 5.5pt, fill: tick-col)[#fmt(lo)])
  for (s, sr) in series.enumerate() {
    let lx = x0 + s * (w / ns)
    draw.rect((lx, y1 + 0.22), (lx + 0.22, y1 + 0.38), fill: sr.at(1), stroke: none, radius: 0.02)
    draw.content((lx + 0.3, y1 + 0.3), anchor: "west", text(size: 6pt)[#sr.at(0)])
  }
})

/// points: array of (x, y). shade-from marks a censored region; marks draws
/// dashed verticals (deciles, thresholds).
#let histogram(points, w: 9.0, h: 3.6, colour: rgb("#2f6f9f"), shade-from: none,
              xlabel: none, ylabel: none, marks: (),
              fmt-y: v => str(calc.round(v / 1000, digits: 0)) + "k") = canvas({
  let (pad-l, pad-b) = (1.1, 0.9)
  let shade-colour = rgb("#e8d6d6")
  let mark-colour = rgb("#b23a48")
  let xs = points.map(p => p.at(0))
  let ys = points.map(p => p.at(1))
  let xmin = calc.min(..xs)
  let xmax = calc.max(..xs)
  let ymax = calc.max(..ys) * 1.15
  let (x0, y0) = (pad-l, pad-b)
  let (x1, y1) = (pad-l + w, pad-b + h)
  let sx = v => x0 + (v - xmin) / (xmax - xmin) * w
  let sy = v => y0 + v / ymax * h
  draw.line((x0, y0), (x1, y0), stroke: axis())
  draw.line((x0, y0), (x0, y1), stroke: axis())
  if shade-from != none {
    draw.rect((sx(shade-from), y0), (x1, y1), fill: shade-colour, stroke: none)
    draw.content(((sx(shade-from) + x1) / 2, y1 - 0.25),
                 text(size: 5.5pt, fill: rgb("#8a4a4a"))[censored above])
  }
  // bar width from the data's own spacing: a bases histogram and a score
  // histogram have ranges ~50x apart and need proportionally different widths
  let bw = (xmax - xmin) / calc.max(points.len(), 1) * 1.1
  for p in points {
    draw.rect((sx(p.at(0)) - bw / 2, y0), (sx(p.at(0)) + bw / 2, sy(p.at(1))),
              fill: colour, stroke: none)
  }
  for m in marks {
    draw.line((sx(m), y0), (sx(m), y1), stroke: (paint: mark-colour, thickness: 0.5pt, dash: "dashed"))
    draw.content((sx(m), y0 - 0.3), text(size: 5.5pt, fill: mark-colour)[#str(calc.round(m, digits: 0))])
  }
  draw.content((x0 - 0.14, y0), anchor: "east", text(size: 5.5pt, fill: tick-col)[0])
  draw.content((x0 - 0.14, sy(calc.max(..ys))), anchor: "east",
               text(size: 5.5pt, fill: tick-col)[#fmt-y(calc.max(..ys))])
  if xlabel != none {
    draw.content(((x0 + x1) / 2, y0 - 0.62), text(size: 6pt)[#xlabel])
  }
  if ylabel != none {
    draw.content((x0 - 0.42, y1 + 0.30), anchor: "west", text(size: 6pt, fill: label-col)[#ylabel])
  }
})

/// points: array of (x, y); `hl` picks the points to draw in `highlight`.
#let scatter(points, w: 9.0, h: 3.6, colour: rgb("#2f6f9f"), highlight: rgb("#b23a48"),
            hl: (p => false), xlabel: none, ylabel: none, labels: ()) = canvas({
  let xs = points.map(p => p.at(0))
  let ys = points.map(p => p.at(1))
  let xmin = calc.min(..xs) * 0.9
  let xmax = calc.max(..xs) * 1.08
  let dy = calc.max(..ys) - calc.min(..ys)
  let dyy = if dy == 0 { 1.0 } else { dy }
  let ymin = calc.min(..ys) - dyy * 0.2
  let ymax = calc.max(..ys) + dyy * 0.2
  let (x0, y0) = (1.2, 0.9)
  let (x1, y1) = (x0 + w, y0 + h)
  let sx = v => x0 + (v - xmin) / (xmax - xmin) * w
  let sy = v => y0 + (v - ymin) / (ymax - ymin) * h
  draw.line((x0, y0), (x1, y0), stroke: axis())
  draw.line((x0, y0), (x0, y1), stroke: axis())
  for (i, p) in points.enumerate() {
    draw.circle((sx(p.at(0)), sy(p.at(1))), radius: 0.10,
                fill: if hl(p) { highlight } else { colour }, stroke: 0.3pt + white)
    if i < labels.len() and labels.at(i) != none {
      draw.content((sx(p.at(0)), sy(p.at(1)) + 0.26), text(size: 5pt, fill: rgb("#444"))[#labels.at(i)])
    }
  }
  if xlabel != none { draw.content(((x0 + x1) / 2, y0 - 0.55), text(size: 6pt)[#xlabel]) }
  if ylabel != none {
    draw.content((x0 - 0.42, y1 + 0.30), anchor: "west", text(size: 6pt, fill: label-col)[#ylabel])
  }
})
