//! PXMT metadata viewer — inspect binary per-read trim metadata.
//!
//! Reads the PXMT binary format produced by `fastprep --metadata` and
//! displays it as a human-readable table.  Supports gzip/zstd-compressed
//! input (auto-detected) and a wide range of filter flags.
//!
//! Memory-efficient: reads records in a single streaming pass, never
//! retaining more than `--limit` records at once.  The summary counters
//! are always O(1) memory.
//!
//! # Examples
//!
//! ```sh
//! # View all records
//! pxmt meta.bin
//!
//! # Show only dropped reads
//! pxmt meta.bin --dropped
//!
//! # Filter by read ID (substring match)
//! pxmt meta.bin --id SRR1234567
//!
//! # Show only adapter-trimmed drops
//! pxmt meta.bin --status adapter
//!
//! # Summary statistics only (constant memory, good for millions of reads)
//! pxmt meta.bin --summary
//!
//! # JSON output for programmatic use
//! pxmt meta.bin --json
//! ```

use clap::Parser;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use needletail::parse_fastx_file;
use fastprep::metadata;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(
    name = "fastprep-meta",
    about = "Inspect PXMT binary metadata files from fastprep",
)]
struct Args {
    /// Path to the PXMT metadata file (.pxmt, .gz, .zst)
    #[arg(required = true)]
    path: PathBuf,

    // Filters
    /// Only show records matching this read-ID substring
    #[arg(long, short = 'i', value_name = "PATTERN")]
    id: Option<String>,

    /// Filter by mate: R1 or R2
    #[arg(long)]
    mate: Option<MateFilter>,

    /// Filter by drop/keep status
    #[arg(long)]
    status: Option<StatusFilter>,

    /// Shorthand: only kept reads (status = kept)
    #[arg(long, conflicts_with = "status")]
    kept: bool,

    /// Shorthand: only dropped reads (status != kept)
    #[arg(long, conflicts_with = "status")]
    dropped: bool,

    /// Only reads with non-zero 5′ adapter trim
    #[arg(long)]
    adapter_5: bool,

    /// Only reads with non-zero 3′ adapter trim
    #[arg(long)]
    adapter_3: bool,

    /// Only reads with any adapter trim
    #[arg(long)]
    adapter: bool,

    /// Only reads with non-zero 5′ polyX trim
    #[arg(long)]
    polyx_5: bool,

    /// Only reads with non-zero 3′ polyX trim
    #[arg(long)]
    polyx_3: bool,

    /// Only reads with any polyX trim
    #[arg(long)]
    polyx: bool,

    /// Only reads with interior runs
    #[arg(long)]
    interior: bool,

    /// Only reads that were untouched (kept with zero trim, no runs)
    #[arg(long)]
    clean: bool,

    // Output modes
    /// Print summary statistics only.  O(1) memory regardless of input size.
    #[arg(long, conflicts_with_all = ["json"])]
    summary: bool,

    /// Print records as JSON array
    #[arg(long)]
    json: bool,

    /// Maximum records to display (0 = unlimited).
    /// Summary counters always reflect the full file.
    #[arg(long, short = 'n', default_value_t = 0)]
    limit: usize,

    /// Sequence view mode: overlay input sequences with trim annotations
    #[arg(long, conflicts_with_all = &["summary", "json", "export_fasta", "out_folder"])]
    view: bool,

    /// Input R1 FASTQ (required with --view, --export-fasta, or --out-folder)
    #[arg(long, short = '1')]
    in_r1: Option<PathBuf>,

    /// Input R2 FASTQ (optional, enables PE mode with --view/--export-fasta/--out-folder)
    #[arg(long, short = '2')]
    in_r2: Option<PathBuf>,

    /// Max sequence characters per line in --view mode (default: 100)
    #[arg(long, default_value_t = 100)]
    view_width: usize,

    /// Hard trim applied to 5′ end (subtracted from P labels in --view)
    #[arg(long, default_value_t = 0)]
    hard_trim_5: usize,

    /// Hard trim applied to 3′ end (subtracted from P labels in --view)
    #[arg(long, default_value_t = 0)]
    hard_trim_3: usize,

    /// Export matching reads to a FASTA file (requires --in-r1)
    #[arg(long, conflicts_with_all = &["summary", "json", "view", "out_folder"])]
    export_fasta: Option<PathBuf>,

    /// Export all categories to separate FASTA files in a directory (requires --in-r1).
    /// Creates: hardtrimmed5.fa, polyx5.fa, polyx3.fa, polyx_interior.fa,
    /// meanq.fa, adapter.fa, duplicate.fa, kept.fa, clean.fa
    #[arg(long, conflicts_with_all = &["summary", "json", "view", "export_fasta"])]
    out_folder: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy)]
enum MateFilter {
    R1,
    R2,
}

impl std::str::FromStr for MateFilter {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "r1" | "1" => Ok(MateFilter::R1),
            "r2" | "2" => Ok(MateFilter::R2),
            _ => Err(format!("invalid mate {:?}, use R1 or R2", s)),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum StatusFilter {
    Kept,
    Dropped,
    Hard,
    Poly5,
    Poly3,
    PolyInt,
    MeanQ,
    Adapter,
    Duplicate,
}

impl std::str::FromStr for StatusFilter {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "kept" | "keep" => Ok(StatusFilter::Kept),
            "dropped" | "drop" => Ok(StatusFilter::Dropped),
            "hard" => Ok(StatusFilter::Hard),
            "poly5" | "poly-5" => Ok(StatusFilter::Poly5),
            "poly3" | "poly-3" => Ok(StatusFilter::Poly3),
            "poly_int" | "polyint" | "interior" => Ok(StatusFilter::PolyInt),
            "meanq" | "mean-q" | "mean_q" => Ok(StatusFilter::MeanQ),
            "adapter" => Ok(StatusFilter::Adapter),
            "duplicate" | "dup" => Ok(StatusFilter::Duplicate),
            _ => Err(format!(
                "invalid status {:?}, use: kept|dropped|hard|poly5|poly3|poly_int|meanq|adapter|duplicate",
                s
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// Record
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct Record {
    id: String,
    mate: u8,   // 0 = R1, 1 = R2
    status: u8, // 0 = kept, 1-6 = drop reasons
    t5_adapter: u32,
    t3_adapter: u32,
    t5_poly: u32,
    t3_poly: u32,
    read_len: u32,
    runs: Vec<(u32, u32)>,
}

impl Record {
    fn status_label(&self) -> &'static str {
        match self.status {
            0 => "kept",
            1 => "hard",
            2 => "poly5",
            3 => "poly3",
            4 => "poly_int",
            5 => "meanq",
            6 => "adapter",
            7 => "duplicate",
            _ => "unknown",
        }
    }

    fn mate_label(&self) -> &'static str {
        match self.mate {
            0 => "R1",
            _ => "R2",
        }
    }

    fn is_dropped(&self) -> bool {
        self.status != 0
    }

    fn has_adapter_trim(&self) -> bool {
        self.t5_adapter > 0 || self.t3_adapter > 0
    }

    fn has_polyx_trim(&self) -> bool {
        self.t5_poly > 0 || self.t3_poly > 0
    }

    fn has_interior(&self) -> bool {
        !self.runs.is_empty()
    }

    fn is_clean(&self) -> bool {
        !self.is_dropped()
            && self.t5_adapter == 0
            && self.t3_adapter == 0
            && self.t5_poly == 0
            && self.t3_poly == 0
            && self.runs.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// Open a file, auto-detecting gzip/zstd compression from magic bytes.
fn open_reader(path: &PathBuf) -> std::io::Result<Box<dyn Read>> {
    // Peek at magic bytes
    let magic = {
        let mut file = std::fs::File::open(path)?;
        let mut m = [0u8; 4];
        file.read_exact(&mut m)?;
        m
    };
    let file = std::fs::File::open(path)?;
    Ok(if magic[..2] == [0x1f, 0x8b] {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if magic[..4] == [0x28, 0xb5, 0x2f, 0xfd] {
        Box::new(zstd::stream::read::Decoder::new(file)?)
    } else {
        Box::new(file)
    })
}

#[derive(Debug, Default)]
struct Summary {
    total: u64,
    kept: u64,
    dropped: u64,
    reason_hard: u64,
    reason_poly5: u64,
    reason_poly3: u64,
    reason_poly_int: u64,
    reason_meanq: u64,
    reason_adapter: u64,
    reason_duplicate: u64,
    r1: u64,
    r2: u64,
    reads_with_adapter: u64,
    reads_with_polyx: u64,
    reads_with_interior: u64,
    clean: u64,
    bases_adapter_5: u64,
    bases_adapter_3: u64,
    bases_polyx_5: u64,
    bases_polyx_3: u64,
    bases_interior: u64,
    total_bases: u64,
    /// Total records that passed the filter (full count, not capped by limit).
    filtered_total: u64,
}

/// Streaming PXMT reader.
///
/// Calls `on_record` for every record in the file.  The caller decides
/// what to store (e.g. collect matching records, increment counters,
/// or both).  Sums are guaranteed complete even when `limit` records
/// have been collected — the read always runs to EOF.
struct PxmtReader<R> {
    buf: std::io::BufReader<R>,
}

impl<R: Read> PxmtReader<R> {
    fn new(reader: R) -> std::io::Result<Self> {
        let mut buf = std::io::BufReader::new(reader);
        let mut magic = [0u8; 4];
        buf.read_exact(&mut magic)?;
        if &magic != b"PXMT" {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "bad magic: expected PXMT, got {:?}",
                    std::str::from_utf8(&magic).unwrap_or("?")
                ),
            ));
        }
        let mut version = [0u8; 1];
        buf.read_exact(&mut version)?;
        if version[0] != metadata::VERSION {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "unsupported PXMT version {}, expected {}",
                    version[0],
                    metadata::VERSION
                ),
            ));
        }
        Ok(Self { buf })
    }

    /// Read the next record, or `None` on EOF.
    fn read_one(&mut self) -> std::io::Result<Option<Record>> {
        let mut id_len_bytes = [0u8; 2];
        if self.buf.read_exact(&mut id_len_bytes).is_err() {
            return Ok(None);
        }
        let id_len = u16::from_le_bytes(id_len_bytes) as usize;

        let mut id_bytes = vec![0u8; id_len];
        if id_len > 0 {
            self.buf.read_exact(&mut id_bytes)?;
        }
        let id = String::from_utf8_lossy(&id_bytes).to_string();

        let mut meta = [0u8; 2];
        self.buf.read_exact(&mut meta)?;
        let mate = meta[0];
        let status = meta[1];

        let mut trim = [0u8; 16];
        self.buf.read_exact(&mut trim)?;
        let t5_adapter = u32::from_le_bytes([trim[0], trim[1], trim[2], trim[3]]);
        let t3_adapter = u32::from_le_bytes([trim[4], trim[5], trim[6], trim[7]]);
        let t5_poly = u32::from_le_bytes([trim[8], trim[9], trim[10], trim[11]]);
        let t3_poly = u32::from_le_bytes([trim[12], trim[13], trim[14], trim[15]]);

        let mut len_bytes = [0u8; 4];
        self.buf.read_exact(&mut len_bytes)?;
        let read_len = u32::from_le_bytes(len_bytes);

        let mut n_runs_bytes = [0u8; 1];
        self.buf.read_exact(&mut n_runs_bytes)?;
        let n_runs = n_runs_bytes[0] as usize;

        let mut runs = Vec::with_capacity(n_runs);
        for _ in 0..n_runs {
            let mut run = [0u8; 8];
            self.buf.read_exact(&mut run)?;
            let start = u32::from_le_bytes([run[0], run[1], run[2], run[3]]);
            let len = u32::from_le_bytes([run[4], run[5], run[6], run[7]]);
            runs.push((start, len));
        }

        Ok(Some(Record {
            id,
            mate,
            status,
            t5_adapter,
            t3_adapter,
            t5_poly,
            t3_poly,
            read_len,
            runs,
        }))
    }
}

/// Stream through all records, collecting summary + matched records.
///
/// `summary` is always fully computed. `matched` is filled (up to `limit`)
/// only for records that pass the filter.  Zero‑copy is achieved by
/// returning records by value — the caller either stores or drops them.
fn stream_records<R: Read>(
    reader: R,
    limit: usize,
    fits: &dyn Fn(&Record) -> bool,
) -> std::io::Result<(Vec<Record>, Summary)> {
    let mut pxmt = PxmtReader::new(reader)?;
    let mut summary = Summary::default();
    let mut matched = Vec::new();

    while let Some(rec) = pxmt.read_one()? {
        // Always update summary counters.
        update_summary(&mut summary, &rec);

        // Collect if it matches the filter (and we still have room).
        if fits(&rec) {
            summary.filtered_total += 1;
            matched.push(rec);
            if limit > 0 && matched.len() >= limit {
                break;
            }
        }
    }

    Ok((matched, summary))
}

#[inline]
fn update_summary(s: &mut Summary, rec: &Record) {
    s.total += 1;
    if rec.is_dropped() {
        s.dropped += 1;
        match rec.status {
            1 => s.reason_hard += 1,
            2 => s.reason_poly5 += 1,
            3 => s.reason_poly3 += 1,
            4 => s.reason_poly_int += 1,
            5 => s.reason_meanq += 1,
            6 => s.reason_adapter += 1,
            7 => s.reason_duplicate += 1,
            _ => {}
        }
    } else {
        s.kept += 1;
    }
    match rec.mate {
        0 => s.r1 += 1,
        _ => s.r2 += 1,
    }
    if rec.has_adapter_trim() {
        s.reads_with_adapter += 1;
    }
    if rec.has_polyx_trim() {
        s.reads_with_polyx += 1;
    }
    if rec.has_interior() {
        s.reads_with_interior += 1;
    }
    if rec.is_clean() {
        s.clean += 1;
    }
    s.bases_adapter_5 += rec.t5_adapter as u64;
    s.bases_adapter_3 += rec.t3_adapter as u64;
    s.bases_polyx_5 += rec.t5_poly as u64;
    s.bases_polyx_3 += rec.t3_poly as u64;
    for (_, l) in &rec.runs {
        s.bases_interior += *l as u64;
    }
    s.total_bases += rec.read_len as u64;
}

// ---------------------------------------------------------------------------
// Output helpers
// ---------------------------------------------------------------------------

fn fmt_num(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.2}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1}K", n as f64 / 1_000.0)
    } else {
        n.to_string()
    }
}

fn pct_str(num: u64, den: u64) -> String {
    if den == 0 {
        "—".to_string()
    } else {
        format!("{:.2}%", 100.0 * num as f64 / den as f64)
    }
}

fn print_summary(s: &Summary) {
    let w = 58;
    eprintln!("{}", "=".repeat(w));
    eprintln!("{:^width$}", "PXMT Metadata Summary", width = w);
    eprintln!("{}", "=".repeat(w));
    eprintln!(
        "  Total records:        {:>10}  {}",
        fmt_num(s.total),
        pct_str(s.total, s.total)
    );
    eprintln!(
        "  Kept:                 {:>10}  {}",
        fmt_num(s.kept),
        pct_str(s.kept, s.total)
    );
    eprintln!(
        "  Dropped:              {:>10}  {}",
        fmt_num(s.dropped),
        pct_str(s.dropped, s.total)
    );
    if s.dropped > 0 {
        eprintln!(
            "    ├─ hard trim:       {:>10}  {}",
            fmt_num(s.reason_hard),
            pct_str(s.reason_hard, s.total)
        );
        eprintln!(
            "    ├─ polyX 5′:        {:>10}  {}",
            fmt_num(s.reason_poly5),
            pct_str(s.reason_poly5, s.total)
        );
        eprintln!(
            "    ├─ polyX 3′:        {:>10}  {}",
            fmt_num(s.reason_poly3),
            pct_str(s.reason_poly3, s.total)
        );
        eprintln!(
            "    ├─ polyX interior:  {:>10}  {}",
            fmt_num(s.reason_poly_int),
            pct_str(s.reason_poly_int, s.total)
        );
        eprintln!(
            "    ├─ mean quality:    {:>10}  {}",
            fmt_num(s.reason_meanq),
            pct_str(s.reason_meanq, s.total)
        );
        eprintln!(
            "    ├─ adapter:         {:>10}  {}",
            fmt_num(s.reason_adapter),
            pct_str(s.reason_adapter, s.total)
        );
        eprintln!(
            "    └─ duplicate:       {:>10}  {}",
            fmt_num(s.reason_duplicate),
            pct_str(s.reason_duplicate, s.total)
        );
    }
    eprintln!(
        "  R1 / single-end:      {:>10}  {}",
        fmt_num(s.r1),
        pct_str(s.r1, s.total)
    );
    eprintln!(
        "  R2:                   {:>10}  {}",
        fmt_num(s.r2),
        pct_str(s.r2, s.total)
    );
    eprintln!(
        "  With adapter trim:    {:>10}  {}",
        fmt_num(s.reads_with_adapter),
        pct_str(s.reads_with_adapter, s.total)
    );
    eprintln!(
        "  With polyX trim:      {:>10}  {}",
        fmt_num(s.reads_with_polyx),
        pct_str(s.reads_with_polyx, s.total)
    );
    eprintln!(
        "  With interior runs:   {:>10}  {}",
        fmt_num(s.reads_with_interior),
        pct_str(s.reads_with_interior, s.total)
    );
    eprintln!(
        "  Clean (untouched):    {:>10}  {}",
        fmt_num(s.clean),
        pct_str(s.clean, s.total)
    );
    eprintln!("{}", "-".repeat(w));
    eprintln!("  Bases in:             {:>10}", fmt_num(s.total_bases));
    eprintln!("  Bases adapter 5′:     {:>10}", fmt_num(s.bases_adapter_5));
    eprintln!("  Bases adapter 3′:     {:>10}", fmt_num(s.bases_adapter_3));
    eprintln!("  Bases polyX 5′:       {:>10}", fmt_num(s.bases_polyx_5));
    eprintln!("  Bases polyX 3′:       {:>10}", fmt_num(s.bases_polyx_3));
    eprintln!("  Bases interior:       {:>10}", fmt_num(s.bases_interior));
    let tot_trim = s.bases_adapter_5
        + s.bases_adapter_3
        + s.bases_polyx_5
        + s.bases_polyx_3
        + s.bases_interior;
    eprintln!(
        "  Bases trimmed (total):{:>10}  {}",
        fmt_num(tot_trim),
        pct_str(tot_trim, s.total_bases)
    );
    eprintln!("{}", "=".repeat(w));
}

fn print_table(records: &[Record]) {
    if records.is_empty() {
        println!("(no matching records)");
        return;
    }

    let mut max_id_len = 4;
    for r in records {
        max_id_len = max_id_len.max(r.id.len().min(60));
    }
    let id_col = max_id_len.max(12).min(60);

    let sep = "─".repeat(id_col + 64);
    println!(
        "  {:<id_col$} {:>4} {:>8} {:>6} {:>6} {:>6} {:>6} {:>7} {:>5} {}",
        "ID",
        "MATE",
        "STATUS",
        "A5",
        "A3",
        "P5",
        "P3",
        "LEN",
        "RUNS",
        "RUNS_DETAIL",
        id_col = id_col,
    );
    println!("{}", sep);

    for r in records {
        let id_display = if r.id.len() > id_col {
            format!("{}…", &r.id[..id_col - 1])
        } else {
            r.id.clone()
        };
        let runs_str = if r.runs.is_empty() {
            "—".to_string()
        } else {
            let total_interior: u32 = r.runs.iter().map(|(_, l)| l).sum();
            format!("{}", total_interior)
        };
        let runs_detail = if r.runs.len() <= 3 {
            r.runs
                .iter()
                .map(|(s, l)| format!("{}:{}", s, l))
                .collect::<Vec<_>>()
                .join(" ")
        } else {
            format!(
                "{} … (+{} more)",
                r.runs[..3]
                    .iter()
                    .map(|(s, l)| format!("{}:{}", s, l))
                    .collect::<Vec<_>>()
                    .join(" "),
                r.runs.len() - 3
            )
        };

        println!(
            "  {:<id_col$} {:>4} {:>8} {:>6} {:>6} {:>6} {:>6} {:>7} {:>5} {}",
            id_display,
            r.mate_label(),
            r.status_label(),
            r.t5_adapter,
            r.t3_adapter,
            r.t5_poly,
            r.t3_poly,
            r.read_len,
            runs_str,
            runs_detail,
            id_col = id_col,
        );
    }
    println!("{}", sep);
}

// ---------------------------------------------------------------------------
// Filter
// ---------------------------------------------------------------------------

fn build_filter(args: &Args, hard_trim_5: usize, hard_trim_3: usize) -> impl Fn(&Record) -> bool + '_ {
    move |r: &Record| -> bool {
        // ID filter
        if let Some(pat) = &args.id {
            if !r.id.contains(pat.as_str()) {
                return false;
            }
        }

        // Mate filter
        if let Some(mf) = args.mate {
            let wanted = match mf {
                MateFilter::R1 => 0,
                MateFilter::R2 => 1,
            };
            if r.mate != wanted {
                return false;
            }
        }

        // Status filter
        if let Some(sf) = args.status {
            let ok = match sf {
                StatusFilter::Kept => !r.is_dropped(),
                StatusFilter::Dropped => r.is_dropped(),
                StatusFilter::Hard => r.status == 1,
                StatusFilter::Poly5 => r.status == 2,
                StatusFilter::Poly3 => r.status == 3,
                StatusFilter::PolyInt => r.status == 4,
                StatusFilter::MeanQ => r.status == 5,
                StatusFilter::Adapter => r.status == 6,
                StatusFilter::Duplicate => r.status == 7,
            };
            if !ok {
                return false;
            }
        }

        if args.kept && r.is_dropped() {
            return false;
        }
        if args.dropped && !r.is_dropped() {
            return false;
        }

        if args.adapter_5 && r.t5_adapter == 0 {
            return false;
        }
        if args.adapter_3 && r.t3_adapter == 0 {
            return false;
        }
        if args.adapter && !r.has_adapter_trim() {
            return false;
        }
        if args.polyx_5 && r.t5_poly <= hard_trim_5 as u32 {
            return false;
        }
        if args.polyx_3 && r.t3_poly <= hard_trim_3 as u32 {
            return false;
        }
        if args.polyx && r.t5_poly <= hard_trim_5 as u32 && r.t3_poly <= hard_trim_3 as u32 {
            return false;
        }
        if args.interior && !r.has_interior() {
            return false;
        }
        if args.clean && !r.is_clean() {
            return false;
        }

        true
    }
}

// ---------------------------------------------------------------------------
// Sequence view mode
// ---------------------------------------------------------------------------

/// Display one record aligned with its FASTQ sequence.
fn display_view(rec: &Record, seq: &[u8], idx: u64, view_width: usize, hard_trim_5: usize, hard_trim_3: usize) {
    let id_display = if rec.id.len() > 32 {
        format!("{}…", &rec.id[..31])
    } else {
        rec.id.clone()
    };

    let seq_len = seq.len();
    let t5 = (rec.t5_adapter + rec.t5_poly) as usize;
    let t3 = (rec.t3_poly + rec.t3_adapter) as usize;
    let out_len = if t5 + t3 > seq_len { 0 } else { seq_len - t5 - t3 };

    // ── Header ──
    if rec.is_dropped() {
        println!(
            "\n── [{}] {} · {} · {}  {}bp ──",
            idx, id_display, rec.mate_label(), rec.status_label(), seq_len
        );
    } else {
        println!(
            "\n── [{}] {} · {} · {}  {}→{}bp ──",
            idx, id_display, rec.mate_label(), rec.status_label(), seq_len, out_len
        );
    }

    // Truncate display to view_width
    let show_len = seq_len.min(view_width);
    let truncated = seq_len > view_width;
    let disp_seq: String = seq[..show_len].iter().map(|&b| b as char).collect();
    let suffix = if truncated { "…" } else { "" };

    // ── Input sequence ──
    println!("  IN:   {}{}", disp_seq, suffix);

    // ── Single annotation bar: █ = hard trim, ▒ = polyX, ▓ = interior run, · = kept ──
    let t5a = (rec.t5_adapter as usize).min(show_len);
    let t5p = (rec.t5_poly as usize).min(show_len.saturating_sub(t5a));
    let t3p = (rec.t3_poly as usize).min(show_len);
    let t3a = (rec.t3_adapter as usize).min(show_len);

    let mut ann_bar: Vec<char> = vec!['·'; show_len];
    let mut label_line: Vec<char> = vec![' '; show_len];

    // Split t5p into hard vs detected
    let t5_hard = hard_trim_5.min(t5p);
    let t5_detected = t5p.saturating_sub(t5_hard);
    // Split t3p into hard vs detected  
    let t3_hard = hard_trim_3.min(t3p);
    let t3_detected = t3p.saturating_sub(t3_hard);

    // 5' trim: █ for hard, ▒ for detected polyX
    for i in 0..t5a.min(show_len) {
        ann_bar[i] = '█';
    }
    for i in 0..t5_hard {
        let pos = t5a + i;
        if pos < show_len { ann_bar[pos] = '█'; }
    }
    for i in 0..t5_detected {
        let pos = t5a + t5_hard + i;
        if pos < show_len { ann_bar[pos] = '▒'; }
    }
    // 3' trim: █ for hard, ▒ for detected polyX
    for i in 0..t3_detected {
        let pos = show_len.saturating_sub(t3p + t3a) + i;
        if pos < show_len { ann_bar[pos] = '▒'; }
    }
    for i in 0..t3_hard {
        let pos = show_len.saturating_sub(t3_hard + t3a) + i;
        if pos < show_len { ann_bar[pos] = '█'; }
    }
    // 3' adapter
    for i in 0..t3a {
        let pos = show_len.saturating_sub(t3a) + i;
        if pos < show_len { ann_bar[pos] = '█'; }
    }

    // Overlay interior runs: ▓ on kept region
    if !rec.runs.is_empty() {
        for &(start, len) in &rec.runs {
            let s = start as usize;
            let e = (s + len as usize).min(show_len);
            for i in s..e {
                if ann_bar[i] == '·' {
                    ann_bar[i] = '▓';
                }
            }
        }
        // Run labels (only for visible runs)
        let mut run_labels = Vec::new();
        for &(start, len) in &rec.runs {
            if start < show_len as u32 {
                run_labels.push(format!("INT {start}:{len}"));
            }
        }
        // Place first run label on the label line
        if let Some(&(first_start, _)) = rec.runs.iter().find(|(s, _)| *s < show_len as u32) {
            let mut current_x = first_start as usize + 1;
            for lbl in &run_labels {
                for (j, c) in lbl.chars().enumerate() {
                    let pos = current_x + j;
                    if pos < show_len { label_line[pos] = c; }
                }
                current_x += lbl.len() + 2; // +2 for separator
            }
        }
    }

    // Trim labels on label line (overwrites run labels where they overlap)
    if t5a > 0 {
        let lbl = format!("A{}", t5a);
        let pos = (t5a / 2).saturating_sub(lbl.len() / 2);
        for (j, c) in lbl.chars().enumerate() {
            if pos + j < show_len { label_line[pos + j] = c; }
        }
    }
    if t5p > 0 {
        let hard = hard_trim_5.min(t5p);
        let detected = t5p.saturating_sub(hard);
        let offset = t5a;
        if hard > 0 {
            let lbl = format!("H{}", hard);
            let pos = offset + (hard / 2).saturating_sub(lbl.len() / 2);
            for (j, c) in lbl.chars().enumerate() {
                if pos + j < show_len { label_line[pos + j] = c; }
            }
        }
        if detected > 0 {
            let lbl = format!("P{}", detected);
            let pos = offset + hard + (detected / 2).saturating_sub(lbl.len() / 2);
            for (j, c) in lbl.chars().enumerate() {
                if pos + j < show_len { label_line[pos + j] = c; }
            }
        }
    }
    if t3p > 0 {
        let hard = hard_trim_3.min(t3p);
        let detected = t3p.saturating_sub(hard);
        let offset = show_len.saturating_sub(t3p + t3a);
        if detected > 0 {
            let lbl = format!("P{}", detected);
            let pos = offset + (detected / 2).saturating_sub(lbl.len() / 2);
            for (j, c) in lbl.chars().enumerate() {
                if pos + j < show_len { label_line[pos + j] = c; }
            }
        }
        if hard > 0 {
            let lbl = format!("H{}", hard);
            let pos = offset + detected + (hard / 2).saturating_sub(lbl.len() / 2);
            for (j, c) in lbl.chars().enumerate() {
                if pos + j < show_len { label_line[pos + j] = c; }
            }
        }
    }
    if t3a > 0 {
        let offset = show_len.saturating_sub(t3a);
        let lbl = format!("A{}", t3a);
        let pos = offset + (t3a / 2).saturating_sub(lbl.len() / 2);
        for (j, c) in lbl.chars().enumerate() {
            if pos + j < show_len { label_line[pos + j] = c; }
        }
    }

    let ann_str: String = ann_bar.iter().collect();
    let label_str: String = label_line.iter().collect();
    println!("  ANN:  {}{}", ann_str, suffix);
    println!("        {}{}", label_str, suffix);

    if rec.is_dropped() {
        println!("        ✗ dropped — {}", rec.status_label());
        return;
    }

    // ── Output sequence (aligned with IN: padded for 5′ trim) ──
    if out_len > 0 {
        let out_start = t5;
        let out_end = seq_len.saturating_sub(t3);
        let pad = out_start.min(view_width);
        let space_pad = " ".repeat(pad);
        let max_out = view_width.saturating_sub(pad);
        let out_seq: String = seq[out_start..out_end]
            .iter()
            .take(max_out)
            .map(|&b| b as char)
            .collect();
        let out_suffix = if out_len > max_out { "…" } else { "" };
        println!("  OUT:  {}{}{}", space_pad, out_seq, out_suffix);
    }
}

/// Stream pxmt + FASTQ files in lockstep, displaying each record.
/// Skips records that don't match the filter (if any).
fn view_records(
    pxmt_reader: Box<dyn Read>,
    in_r1_path: &Path,
    in_r2_path: Option<&Path>,
    limit: usize,
    view_width: usize,
    hard_trim_5: usize,
    hard_trim_3: usize,
    filter: &dyn Fn(&Record) -> bool,
) -> anyhow::Result<()> {
    let mut pxmt = PxmtReader::new(pxmt_reader)?;
    let mut reader_r1 = parse_fastx_file(in_r1_path)
        .map_err(|e| anyhow::anyhow!("cannot open R1: {e}"))?;
    let mut reader_r2 = match in_r2_path {
        Some(p) => Some(parse_fastx_file(p)
            .map_err(|e| anyhow::anyhow!("cannot open R2: {e}"))?),
        None => None,
    };

    let mut shown = 0u64;
    let mut total = 0u64;

    loop {
        let Some(rec) = pxmt.read_one()? else { break };
        total += 1;

        // Read matching FASTQ record (advance both readers in lockstep)
        let fq_rec = if rec.mate == 0 {
            match reader_r1.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => return Err(anyhow::anyhow!("R1 parse error: {e}")),
                None => break,
            }
        } else {
            match &mut reader_r2 {
                Some(r2) => match r2.next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => return Err(anyhow::anyhow!("R2 parse error: {e}")),
                    None => break,
                },
                None => match reader_r1.next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => return Err(anyhow::anyhow!("R1 parse error: {e}")),
                    None => break,
                },
            }
        };

        // Apply filter
        if !filter(&rec) {
            continue;
        }

        display_view(&rec, fq_rec.seq().as_ref(), shown + 1, view_width, hard_trim_5, hard_trim_3);

        shown += 1;
        if limit > 0 && shown >= limit as u64 {
            break;
        }
    }

    if shown == 0 && total > 0 {
        eprintln!("  (no records matched the filter)");
    } else if total > shown {
        eprintln!("\n  Showing {shown} of {total} records");
    }

    Ok(())
}

/// Export matching reads to a FASTA file.
fn export_fasta_records(
    pxmt_reader: Box<dyn Read>,
    in_r1_path: &Path,
    in_r2_path: Option<&Path>,
    out_path: &Path,
    limit: usize,
    filter: &dyn Fn(&Record) -> bool,
) -> anyhow::Result<()> {
    let mut pxmt = PxmtReader::new(pxmt_reader)?;
    let mut reader_r1 =
        parse_fastx_file(in_r1_path).map_err(|e| anyhow::anyhow!("cannot open R1: {e}"))?;
    let mut reader_r2 = match in_r2_path {
        Some(p) => Some(parse_fastx_file(p).map_err(|e| anyhow::anyhow!("cannot open R2: {e}"))?),
        None => None,
    };

    let mut out = std::fs::File::create(out_path)
        .map_err(|e| anyhow::anyhow!("cannot create {out_path:?}: {e}"))?;

    let mut written = 0u64;
    let mut total = 0u64;

    loop {
        let Some(rec) = pxmt.read_one()? else { break };
        total += 1;

        // Read matching FASTQ record in lockstep
        let fq_rec = if rec.mate == 0 {
            match reader_r1.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => return Err(anyhow::anyhow!("R1 parse error: {e}")),
                None => break,
            }
        } else {
            match &mut reader_r2 {
                Some(r2) => match r2.next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => return Err(anyhow::anyhow!("R2 parse error: {e}")),
                    None => break,
                },
                None => match reader_r1.next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => return Err(anyhow::anyhow!("R1 parse error: {e}")),
                    None => break,
                },
            }
        };

        if !filter(&rec) {
            continue;
        }

        // Write FASTA record
        let id = std::str::from_utf8(fq_rec.id()).unwrap_or("?");
        write!(out, ">{id}\n")?;
        out.write_all(fq_rec.seq().as_ref())?;
        out.write_all(b"\n")?;

        written += 1;
        if limit > 0 && written >= limit as u64 {
            break;
        }
    }

    if written == 0 && total > 0 {
        eprintln!("  (no records matched the filter)");
    } else {
        eprintln!("  Exported {written} of {total} records to {out_path:?}");
    }

    Ok(())
}

/// Export all categories to separate FASTA files in a directory.
fn export_categories(
    pxmt_reader: Box<dyn Read>,
    in_r1_path: &Path,
    in_r2_path: Option<&Path>,
    out_dir: &Path,
) -> anyhow::Result<()> {
    use std::io::BufWriter;

    std::fs::create_dir_all(out_dir)
        .map_err(|e| anyhow::anyhow!("cannot create {out_dir:?}: {e}"))?;

    let mut pxmt = PxmtReader::new(pxmt_reader)?;
    let mut reader_r1 =
        parse_fastx_file(in_r1_path).map_err(|e| anyhow::anyhow!("cannot open R1: {e}"))?;
    let mut reader_r2 = match in_r2_path {
        Some(p) => Some(parse_fastx_file(p).map_err(|e| anyhow::anyhow!("cannot open R2: {e}"))?),
        None => None,
    };

    let names = ["hardtrimmed5", "polyx5", "polyx3", "polyx_interior", "meanq", "adapter", "duplicate", "kept", "clean"];
    let mut files: Vec<Option<BufWriter<std::fs::File>>> = names
        .iter()
        .map(|name| {
            let path = out_dir.join(format!("{name}.fa"));
            Some(BufWriter::new(
                std::fs::File::create(&path).unwrap(),
            ))
        })
        .collect();

    const H: usize = 0;
    const P5: usize = 1;
    const P3: usize = 2;
    const PI: usize = 3;
    const MQ: usize = 4;
    const AD: usize = 5;
    const DP: usize = 6;
    const KP: usize = 7;
    const CL: usize = 8;

    let mut total = 0u64;
    let mut counts: [u64; 9] = [0; 9];

    loop {
        let Some(rec) = pxmt.read_one()? else { break };
        total += 1;

        let fq_rec = if rec.mate == 0 {
            match reader_r1.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => return Err(anyhow::anyhow!("R1 parse error: {e}")),
                None => break,
            }
        } else {
            match &mut reader_r2 {
                Some(r2) => match r2.next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => return Err(anyhow::anyhow!("R2 parse error: {e}")),
                    None => break,
                },
                None => match reader_r1.next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => return Err(anyhow::anyhow!("R1 parse error: {e}")),
                    None => break,
                },
            }
        };

        let id = std::str::from_utf8(fq_rec.id()).unwrap_or("?");
        let seq = fq_rec.seq();

        let idx = if rec.is_dropped() {
            match rec.status {
                1 => H,
                2 => P5,
                3 => P3,
                4 => PI,
                5 => MQ,
                6 => AD,
                7 => DP,
                _ => continue,
            }
        } else if rec.t5_adapter == 0 && rec.t3_adapter == 0
            && rec.t5_poly == 0 && rec.t3_poly == 0
            && rec.runs.is_empty()
        {
            CL
        } else {
            KP
        };

        if let Some(ref mut f) = files[idx] {
            write!(f, ">{id}\n")?;
            f.write_all(seq.as_ref())?;
            f.write_all(b"\n")?;
            counts[idx] += 1;
        }
    }

    for (i, f) in files.iter_mut().enumerate() {
        if let Some(fw) = f.take() {
            fw.into_inner().map_err(|e| anyhow::anyhow!("flush error for {}: {e}", names[i]))?;
        }
    }

    eprintln!("  Exported {total} records to {out_dir:?}:");
    for (i, name) in names.iter().enumerate() {
        eprintln!("    {name}.fa — {} reads", counts[i]);
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let reader = open_reader(&args.path)?;

    if args.summary {
        // Pure streaming: no record storage, O(1) memory.
        let mut pxmt = PxmtReader::new(reader)?;
        let mut summary = Summary::default();
        while let Some(rec) = pxmt.read_one()? {
            update_summary(&mut summary, &rec);
        }
        print_summary(&summary);
        return Ok(());
    }

    if args.out_folder.is_some() {
        let in_r1 = args.in_r1.as_ref()
            .ok_or_else(|| anyhow::anyhow!("--in-r1 is required with --out-folder"))?;
        let in_r2 = args.in_r2.as_ref().map(|v| v.as_path());
        let out_dir = args.out_folder.as_ref().unwrap();
        return export_categories(reader, in_r1, in_r2, out_dir);
    }

    if args.export_fasta.is_some() {
        let in_r1 = args.in_r1.as_ref()
            .ok_or_else(|| anyhow::anyhow!("--in-r1 is required with --export-fasta"))?;
        let in_r2 = args.in_r2.as_ref().map(|v| v.as_path());
        let out_path = args.export_fasta.as_ref().unwrap();
        let filter = build_filter(&args, args.hard_trim_5, args.hard_trim_3);
        return export_fasta_records(reader, in_r1, in_r2, out_path, args.limit, &filter);
    }

    if args.view {
        let in_r1 = args.in_r1.as_ref().expect("--in-r1 required with --view");
        let in_r2 = args.in_r2.as_ref().map(|v| v.as_path());
        let filter = build_filter(&args, args.hard_trim_5, args.hard_trim_3);
        return view_records(reader, in_r1, in_r2, args.limit, args.view_width, args.hard_trim_5, args.hard_trim_3, &filter);
    }

    // Normal mode: stream + collect matches.
    let filter = build_filter(&args, 0, 0);
    let (matched, summary) = stream_records(reader, args.limit, &filter)?;

    print_summary(&summary);
    println!();

    if args.json {
        let json_records: Vec<serde_json::Value> = matched
            .iter()
            .map(|r| {
                let mut runs = Vec::with_capacity(r.runs.len());
                for &(s, l) in &r.runs {
                    use serde_json::json;
                    runs.push(json!({"start": s, "len": l}));
                }
                serde_json::json!({
                    "id": r.id,
                    "mate": r.mate_label(),
                    "status": r.status_label(),
                    "t5_adapter": r.t5_adapter,
                    "t3_adapter": r.t3_adapter,
                    "t5_poly": r.t5_poly,
                    "t3_poly": r.t3_poly,
                    "read_len": r.read_len,
                    "runs": runs,
                })
            })
            .collect();

        let output = if matched.len() <= 50 {
            serde_json::to_string_pretty(&serde_json::json!({
                "summary": {
                    "total": summary.total,
                    "kept": summary.kept,
                    "dropped": summary.dropped,
                    "records_shown": matched.len(),
                },
                "records": json_records,
            }))?
        } else {
            serde_json::to_string(&serde_json::json!({
                "summary": {
                    "total": summary.total,
                    "kept": summary.kept,
                    "dropped": summary.dropped,
                    "records_shown": matched.len(),
                },
                "records": json_records,
            }))?
        };
        println!("{}", output);
    } else {
        print_table(&matched);
        if summary.filtered_total > matched.len() as u64 {
            println!(
                "  … and {} more matching records (use -n to show more, --summary for stats only)",
                summary.filtered_total - matched.len() as u64
            );
        }
    }

    Ok(())
}
