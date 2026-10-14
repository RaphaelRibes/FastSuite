use crate::hasher::{AbundanceTracker, HashVerifier, SequenceHasher};
use crate::utils::{base_read_id, compute_write_buffer_size, validate_format_compatibility};
use crate::writer::{
    OutputFormat, preload_existing_hashes, preload_existing_paired_hashes,
    write_fasta_record, write_fasta_record_with_abundance, write_fastq_record_with_abundance,
};
use anyhow::{Context, Result, bail};
use needletail::parse_fastx_file;
use std::fs::{self, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Shared configuration for a deduplication run.
pub(crate) struct DedupConfig {
    /// Whether to overwrite existing output files.
    pub force: bool,
    /// Whether to print diagnostic messages.
    pub verbose: bool,
    /// When `true`, only compute stats without writing output.
    pub dry_run: bool,
    /// GZIP compression level (1–9).
    pub compression_level: u32,
    /// Pre-allocation hint for the internal hash set.
    pub estimated_capacity: usize,
    /// Expected read length in base pairs, used to size I/O buffers.
    pub read_length: usize,
    /// Separator for abundance annotation in output headers (e.g., ";size=").
    pub abundance_sep: String,
    /// Number of parallel workers for multi-threaded mode (1 = single-threaded).
    pub threads: usize,
    /// When `Some`, abundance mode uses a temp file in this directory to avoid
    /// re-decompressing the input gzip in the second pass.
    pub temp_dir: Option<String>,
    /// Output record format, used when the output is a stream; named paths
    /// infer it from their extension.
    pub out_format: OutputFormat,
    /// Output compression, used when the output is a stream; named paths
    /// infer it from their extension.
    pub out_compression: crate::writer::StreamCompression,
}

/// Owned snapshot of a sequence record, used to break needletail's borrow
/// chain when the same reader needs to be advanced multiple times.
struct OwnedRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

// ==========================================
// SINGLE-END MODE
// ==========================================

/// Open the output for a run: stdout when the path is `-` (D-4), the file
/// otherwise. The stream branch carries the explicit format/compression
/// decision from `--out-format` / `--out-compression`.
fn open_output(path: &str, cfg: &DedupConfig) -> Result<(Box<dyn Write>, OutputFormat)> {
    if crate::stream::is_stream(path) {
        Ok((
            crate::writer::prepare_stream_writer(cfg.out_compression, cfg.compression_level),
            cfg.out_format,
        ))
    } else {
        crate::writer::prepare_writer(Path::new(path), cfg.force, cfg.compression_level)
    }
}

/// Executes deduplication for single-end reads.
pub(crate) fn execute_deduplication<T: SequenceHasher + 'static>(
    input_path: &str,
    output_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    let mut verifier = HashVerifier::<T>::new(cfg.estimated_capacity)?;

    if cfg.dry_run {
        if cfg.verbose {
            println!("--dry-run enabled: calculating duplication rate without writing files.");
        }
        let mut reader = crate::stream::open_fastx(input_path)?;
        let mut processed_sequences = 0;
        let mut duplicates = 0;

        while let Some(record) = reader.next() {
            let seq_record = record.context("Invalid sequence data")?;
            let seq = seq_record.seq();
            let hash = T::hash_sequence(&seq);
            let is_unique = verifier.verify(hash);

            if !is_unique {
                duplicates += 1;
            }
            processed_sequences += 1;
        }

        return Ok((processed_sequences, duplicates));
    }

    // --- PRELOADING ---
    let output_path_bound = Path::new(output_path);
    let is_gz = OutputFormat::is_gz(output_path_bound);

    if crate::stream::is_stream(output_path) {
        // D-4: a stream is written fresh — nothing exists to preload hashes
        // from and nothing to truncate. Resume is a file-only feature.
        if cfg.verbose {
            println!("Output is a stream: written fresh (no preload, no resume).");
        }
    } else if cfg.force {
        if cfg.verbose {
            println!("--force enabled: overwriting output.");
        }
    } else {
        let (preloads, valid_bytes) =
            preload_existing_hashes::<T>(output_path, &mut verifier, cfg.verbose)?;
        if preloads > 0 && cfg.verbose {
            println!("{} sequences preloaded.", preloads);
        }

        if output_path_bound.exists() {
            let current_size = fs::metadata(output_path_bound)?.len();

            if valid_bytes < current_size {
                if is_gz {
                    bail!(
                        "The output file (.gz) is corrupted and cannot be truncated. Use --force to restart."
                    );
                }
                if cfg.verbose {
                    println!(
                        "Truncating corrupted file from {} to {} bytes.",
                        current_size, valid_bytes
                    );
                }
                let file = OpenOptions::new().write(true).open(output_path_bound)?;
                file.set_len(valid_bytes)?;
            }
        }
    }

    // --- WRITER PREPARATION ---
    let (writer, output_format) = open_output(output_path, cfg)?;
    let buf_size = compute_write_buffer_size(cfg.read_length);
    let mut buffered_writer = BufWriter::with_capacity(buf_size, writer);

    if cfg.verbose {
        println!("Write buffer size: {} KB", buf_size / 1024);
    }

    // --- READING AND PARSING ---
    let mut reader = crate::stream::open_fastx(input_path)?;
    let mut processed_sequences = 0;
    let mut duplicates = 0;
    let mut format_validated = false;

    while let Some(record) = reader.next() {
        let seq_record = record.context("Invalid sequence data")?;

        if !format_validated {
            let has_quality = seq_record.qual().is_some_and(|q| !q.is_empty());
            validate_format_compatibility(has_quality, output_format)?;
            format_validated = true;
        }

        let seq = seq_record.seq();
        let hash = T::hash_sequence(&seq);
        let is_unique = verifier.verify(hash);

        if is_unique {
            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buffered_writer, seq_record.id(), &seq)?;
                }
                OutputFormat::Fastq => {
                    seq_record
                        .write(&mut buffered_writer, None)
                        .context("Error writing FASTQ")?;
                }
            }
        } else {
            duplicates += 1;
        }
        processed_sequences += 1;
    }

    buffered_writer
        .flush()
        .context("Failed to flush output buffer")?;

    Ok((processed_sequences, duplicates))
}

// ==========================================
// PAIRED-END MODE
// ==========================================

/// Executes deduplication for paired-end reads.
pub(crate) fn execute_paired_deduplication<T: SequenceHasher + 'static>(
    input_r1_path: &str,
    input_r2_path: &str,
    output_r1_path: &str,
    output_r2_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    let mut verifier = HashVerifier::<T>::new(cfg.estimated_capacity)?;

    let mut reader_r1 = crate::stream::open_fastx(input_r1_path)?;
    let mut reader_r2 = crate::stream::open_fastx(input_r2_path)?;

    let mut processed_sequences = 0;
    let mut duplicates = 0;

    if cfg.dry_run {
        if cfg.verbose {
            println!("--dry-run enabled: calculating Paired-End duplication rate.");
        }

        while let (Some(record_r1_res), Some(record_r2_res)) = (reader_r1.next(), reader_r2.next())
        {
            let seq_r1 = record_r1_res.context("Invalid sequence in R1")?;
            let seq_r2 = record_r2_res.context("Invalid sequence in R2")?;

            let base_id_r1 = base_read_id(seq_r1.id());
            let base_id_r2 = base_read_id(seq_r2.id());

            if base_id_r1 != base_id_r2 {
                bail!(
                    "Critical desynchronization detected! R1: {}, R2: {}",
                    String::from_utf8_lossy(base_id_r1),
                    String::from_utf8_lossy(base_id_r2)
                );
            }

            let r1_seq = seq_r1.seq();
            let r2_seq = seq_r2.seq();
            let combined_hash = T::hash_pair(&r1_seq, &r2_seq);
            if !verifier.verify(combined_hash) {
                duplicates += 1;
            }
            processed_sequences += 1;
        }

        // Count remaining records to detect length mismatch
        let mut r1_remaining = 0;
        while reader_r1.next().is_some() {
            r1_remaining += 1;
        }
        let mut r2_remaining = 0;
        while reader_r2.next().is_some() {
            r2_remaining += 1;
        }
        if r1_remaining != r2_remaining {
            bail!(
                "Desynchronization detected at the end: R1 has {} unmatched reads, R2 has {} unmatched reads",
                r1_remaining,
                r2_remaining
            );
        }

        return Ok((processed_sequences, duplicates));
    }

    // --- PRELOADING AND SYNCHRONIZED TRUNCATION ---
    let path_r1_bound = Path::new(output_r1_path);
    let path_r2_bound = Path::new(output_r2_path);

    let is_gz_r1 = OutputFormat::is_gz(path_r1_bound);
    let is_gz_r2 = OutputFormat::is_gz(path_r2_bound);

    if crate::stream::is_stream(output_r1_path) || crate::stream::is_stream(output_r2_path) {
        // D-4: a stream is written fresh — nothing exists to preload hashes
        // from and nothing to truncate. Resume is a file-only feature.
        if cfg.verbose {
            println!("Output is a stream: written fresh (no preload, no resume).");
        }
    } else if cfg.force {
        if cfg.verbose {
            println!("--force enabled: overwriting Paired-End outputs.");
        }
    } else {
        let (preloads, bytes_r1, bytes_r2) = preload_existing_paired_hashes::<T>(
            output_r1_path,
            output_r2_path,
            &mut verifier,
            cfg.verbose,
        )?;

        if preloads > 0 && cfg.verbose {
            println!("{} pairs preloaded and synchronized.", preloads);
        }

        if path_r1_bound.exists() {
            let current_size_r1 = fs::metadata(path_r1_bound)?.len();
            if bytes_r1 < current_size_r1 {
                if is_gz_r1 {
                    bail!("R1 file (.gz) is desynchronized and cannot be truncated. Use --force.");
                }
                if cfg.verbose {
                    println!(
                        "Truncating R1 for resynchronization ({} to {} bytes).",
                        current_size_r1, bytes_r1
                    );
                }
                let file = OpenOptions::new().write(true).open(path_r1_bound)?;
                file.set_len(bytes_r1)?;
            }
        }

        if path_r2_bound.exists() {
            let current_size_r2 = fs::metadata(path_r2_bound)?.len();
            if bytes_r2 < current_size_r2 {
                if is_gz_r2 {
                    bail!("R2 file (.gz) is desynchronized and cannot be truncated. Use --force.");
                }
                if cfg.verbose {
                    println!(
                        "Truncating R2 for resynchronization ({} to {} bytes).",
                        current_size_r2, bytes_r2
                    );
                }
                let file = OpenOptions::new().write(true).open(path_r2_bound)?;
                file.set_len(bytes_r2)?;
            }
        }
    }

    // --- WRITER PREPARATION ---
    let (writer_r1, output_format_r1) = open_output(output_r1_path, cfg)?;
    let (writer_r2, output_format_r2) = open_output(output_r2_path, cfg)?;

    if output_format_r1 != output_format_r2 {
        bail!("Output files R1 and R2 must have the same format (FASTA or FASTQ)");
    }

    let buf_size = compute_write_buffer_size(cfg.read_length);

    if cfg.verbose {
        println!("Write buffer size: {} KB", buf_size / 1024);
    }

    let mut buf_writer_r1 = BufWriter::with_capacity(buf_size, writer_r1);
    let mut buf_writer_r2 = BufWriter::with_capacity(buf_size, writer_r2);

    // --- SYNCHRONIZED MAIN WRITE LOOP ---
    let mut format_validated = false;

    while let (Some(record_r1_res), Some(record_r2_res)) = (reader_r1.next(), reader_r2.next()) {
        let seq_r1 = record_r1_res.context("Invalid sequence in R1")?;
        let seq_r2 = record_r2_res.context("Invalid sequence in R2")?;

        if !format_validated {
            let has_quality = seq_r1.qual().is_some_and(|q| !q.is_empty());
            validate_format_compatibility(has_quality, output_format_r1)?;
            format_validated = true;
        }

        let base_id_r1 = seq_r1
            .id()
            .split(|&b| b == b' ')
            .next()
            .unwrap_or(seq_r1.id());
        let base_id_r2 = seq_r2
            .id()
            .split(|&b| b == b' ')
            .next()
            .unwrap_or(seq_r2.id());

        if base_id_r1 != base_id_r2 {
            bail!(
                "Critical desynchronization detected at pair #{}! R1: {}, R2: {}",
                processed_sequences + 1,
                String::from_utf8_lossy(base_id_r1),
                String::from_utf8_lossy(base_id_r2)
            );
        }

        let r1_seq = seq_r1.seq();
        let r2_seq = seq_r2.seq();
        let combined_hash = T::hash_pair(&r1_seq, &r2_seq);
        let is_unique = verifier.verify(combined_hash);

        if is_unique {
            match output_format_r1 {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buf_writer_r1, seq_r1.id(), &r1_seq)?;
                    write_fasta_record(&mut buf_writer_r2, seq_r2.id(), &r2_seq)?;
                }
                OutputFormat::Fastq => {
                    seq_r1
                        .write(&mut buf_writer_r1, None)
                        .context("Error writing R1 FASTQ")?;
                    seq_r2
                        .write(&mut buf_writer_r2, None)
                        .context("Error writing R2 FASTQ")?;
                }
            }
        } else {
            duplicates += 1;
        }

        processed_sequences += 1;
    }

    if reader_r1.next().is_some() || reader_r2.next().is_some() {
        bail!(
            "Desynchronization detected at the end: one file contains more reads than the other!"
        );
    }

    buf_writer_r1
        .flush()
        .context("Failed to flush R1 output buffer")?;
    buf_writer_r2
        .flush()
        .context("Failed to flush R2 output buffer")?;

    Ok((processed_sequences, duplicates))
}

// ==========================================
// ABUNDANCE TRACKING MODE (SINGLE-END)
// ==========================================

/// Executes abundance-tracking deduplication for single-end reads.
///
/// Two-pass strategy. When `cfg.temp_dir` is set, the first pass writes unique
/// records to an uncompressed temp file so the second pass avoids re-decompressing
/// the input gzip — trading disk I/O for CPU time.
pub(crate) fn execute_deduplication_with_abundance<T: SequenceHasher + 'static>(
    input_path: &str,
    output_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    let mut tracker = AbundanceTracker::<T>::new(cfg.estimated_capacity)?;
    let mut total_reads = 0;

    if let Some(ref temp_dir) = cfg.temp_dir {
        // === TEMP FILE PATH (first pass writes to disk, second pass reads raw) ===

        // --- PASS 1: count abundances + write unique records to temp file ---
        let temp_path = std::path::Path::new(temp_dir)
            .join(format!("fastdedup_abundance_{}.fastq", std::process::id()));
        let temp_file = fs::File::create(&temp_path)
            .with_context(|| format!("Failed to create temp: {temp_path:?}"))?;
        let mut temp_writer = BufWriter::with_capacity(256 * 1024, temp_file);

        {
            let mut reader = crate::stream::open_fastx(input_path)?;
            let mut is_fastq = false;
            let mut format_detected = false;

            while let Some(record) = reader.next() {
                let seq_record = record.context("Invalid sequence data")?;

                if !format_detected {
                    is_fastq = seq_record.qual().is_some_and(|q| !q.is_empty());
                    format_detected = true;
                }

                let seq = seq_record.seq();
                let hash = T::hash_sequence(&seq);
                let (is_first, _) = tracker.track(hash);
                total_reads += 1;

                if is_first {
                    if is_fastq {
                        let qual = seq_record.qual().unwrap_or(&[]);
                        temp_writer.write_all(b"@")?;
                        temp_writer.write_all(seq_record.id())?;
                        temp_writer.write_all(b"\n")?;
                        temp_writer.write_all(&seq)?;
                        temp_writer.write_all(b"\n+\n")?;
                        temp_writer.write_all(qual)?;
                        temp_writer.write_all(b"\n")?;
                    } else {
                        temp_writer.write_all(b">")?;
                        temp_writer.write_all(seq_record.id())?;
                        temp_writer.write_all(b"\n")?;
                        temp_writer.write_all(&seq)?;
                        temp_writer.write_all(b"\n")?;
                    }
                }
            }
        }

        temp_writer.flush().context("Failed to flush temp file")?;
        drop(temp_writer);

        if cfg.verbose {
            println!(
                "Abundance tracking: {} unique / {} total reads (temp: {:?})",
                tracker.len(),
                total_reads,
                temp_path
            );
        }

        // --- PASS 2: read temp file (raw, no decompression) ---
        let (writer, output_format) = open_output(output_path, cfg)?;
        let buf_size = compute_write_buffer_size(cfg.read_length);
        let mut buffered_writer = BufWriter::with_capacity(buf_size, writer);

        let sep = cfg.abundance_sep.as_bytes();
        let temp_path_str = temp_path.to_string_lossy().to_string();
        let mut reader = parse_fastx_file(&temp_path_str).context("Failed to read temp file")?;
        let mut format_validated = false;

        while let Some(record) = reader.next() {
            let seq_record = record.context("Invalid sequence data in temp file")?;

            if !format_validated {
                let has_quality = seq_record.qual().is_some_and(|q| !q.is_empty());
                validate_format_compatibility(has_quality, output_format)?;
                format_validated = true;
            }

            let seq = seq_record.seq();
            let hash = T::hash_sequence(&seq);
            let count = tracker.get_count(hash);

            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record_with_abundance(
                        &mut buffered_writer,
                        seq_record.id(),
                        &seq,
                        count,
                        sep,
                    )?;
                }
                OutputFormat::Fastq => {
                    let qual = seq_record.qual().unwrap_or(&[]);
                    write_fastq_record_with_abundance(
                        &mut buffered_writer,
                        seq_record.id(),
                        &seq,
                        qual,
                        count,
                        sep,
                    )?;
                }
            }
        }

        buffered_writer
            .flush()
            .context("Failed to flush output buffer")?;

        let _ = fs::remove_file(&temp_path);
        let duplicates = total_reads - tracker.len();
        Ok((total_reads, duplicates))

    } else {
        // === CLASSIC PATH (re-read input gzip for second pass) ===

        // --- PASS 1: count abundances only ---
        {
            let mut reader = crate::stream::open_fastx(input_path)?;
            while let Some(record) = reader.next() {
                let seq_record = record.context("Invalid sequence data")?;
                let seq = seq_record.seq();
                let hash = T::hash_sequence(&seq);
                tracker.track(hash);
                total_reads += 1;
            }
        }

        // --- PASS 2: re-read gzip, annotate with abundance ---
        let (writer, output_format) = open_output(output_path, cfg)?;
        let buf_size = compute_write_buffer_size(cfg.read_length);
        let mut buffered_writer = BufWriter::with_capacity(buf_size, writer);

        if cfg.verbose {
            println!(
                "Abundance tracking: {} unique / {} total reads (re-reading input)",
                tracker.len(),
                total_reads
            );
        }

        let sep = cfg.abundance_sep.as_bytes();
        let mut reader = crate::stream::open_fastx(input_path)?;
        let mut written = rustc_hash::FxHashSet::default();
        let mut format_validated = false;

        while let Some(record) = reader.next() {
            let seq_record = record.context("Invalid sequence data")?;

            if !format_validated {
                let has_quality = seq_record.qual().is_some_and(|q| !q.is_empty());
                validate_format_compatibility(has_quality, output_format)?;
                format_validated = true;
            }

            let seq = seq_record.seq();
            let hash = T::hash_sequence(&seq);
            if !written.insert(hash) {
                continue;
            }
            let count = tracker.get_count(hash);

            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record_with_abundance(
                        &mut buffered_writer,
                        seq_record.id(),
                        &seq,
                        count,
                        sep,
                    )?;
                }
                OutputFormat::Fastq => {
                    let qual = seq_record.qual().unwrap_or(&[]);
                    write_fastq_record_with_abundance(
                        &mut buffered_writer,
                        seq_record.id(),
                        &seq,
                        qual,
                        count,
                        sep,
                    )?;
                }
            }
        }

        buffered_writer
            .flush()
            .context("Failed to flush output buffer")?;

        let duplicates = total_reads - tracker.len();
        Ok((total_reads, duplicates))
    }
}

// ==========================================
// ABUNDANCE TRACKING MODE (PAIRED-END)
// ==========================================

/// Executes abundance-tracking deduplication for paired-end reads.
///
/// Two-pass strategy. When `cfg.temp_dir` is set, the first pass writes unique
/// pairs to an uncompressed temp file so the second pass avoids re-decompressing
/// the input gzip — trading disk I/O for CPU time.
pub(crate) fn execute_paired_deduplication_with_abundance<T: SequenceHasher + 'static>(
    input_r1_path: &str,
    input_r2_path: &str,
    output_r1_path: &str,
    output_r2_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    let mut tracker = AbundanceTracker::<T>::new(cfg.estimated_capacity)?;
    let mut total_reads = 0;

    if let Some(ref temp_dir) = cfg.temp_dir {
        // === TEMP FILE PATH ===

        // --- PASS 1: count abundances + write unique pairs to temp file ---
        let temp_path = std::path::Path::new(temp_dir)
            .join(format!("fastdedup_abundance_{}.fastq", std::process::id()));
        let temp_file = fs::File::create(&temp_path)
            .with_context(|| format!("Failed to create temp: {temp_path:?}"))?;
        let mut temp_writer = BufWriter::with_capacity(256 * 1024, temp_file);

        {
            let mut reader_r1 = crate::stream::open_fastx(input_r1_path)?;
            let mut reader_r2 = crate::stream::open_fastx(input_r2_path)?;
            let mut is_fastq = false;
            let mut format_detected = false;

            while let (Some(r1_res), Some(r2_res)) = (reader_r1.next(), reader_r2.next()) {
                let seq_r1 = r1_res.context("Invalid sequence in R1")?;

                if !format_detected {
                    is_fastq = seq_r1.qual().is_some_and(|q| !q.is_empty());
                    format_detected = true;
                }

                let seq_r2 = r2_res.context("Invalid sequence in R2")?;

                let base_id_r1 = base_read_id(seq_r1.id());
                let base_id_r2 = base_read_id(seq_r2.id());
                if base_id_r1 != base_id_r2 {
                    bail!(
                        "Critical desynchronization detected at pair #{}! R1: {}, R2: {}",
                        total_reads + 1,
                        String::from_utf8_lossy(base_id_r1),
                        String::from_utf8_lossy(base_id_r2)
                    );
                }

                let combined_hash = T::hash_pair(&seq_r1.seq(), &seq_r2.seq());
                let (is_first, _) = tracker.track(combined_hash);
                total_reads += 1;

                if is_first {
                    write_paired_record(&mut temp_writer, &seq_r1, &seq_r2, is_fastq)?;
                }
            }

            // Check for trailing records (desynchronization at end)
            if reader_r1.next().is_some() || reader_r2.next().is_some() {
                bail!(
                    "Desynchronization detected at the end: one file contains more reads than the other!"
                );
            }
        }

        temp_writer.flush().context("Failed to flush temp file")?;
        drop(temp_writer);

        if cfg.verbose {
            println!(
                "Abundance tracking: {} unique / {total_reads} pairs (temp: {:?})",
                tracker.len(),
                temp_path
            );
        }

        // --- PASS 2: read temp file in pairs, annotate with abundance ---
        let (writer_r1, output_format_r1) = open_output(output_r1_path, cfg)?;
        let (writer_r2, output_format_r2) = open_output(output_r2_path, cfg)?;

        if output_format_r1 != output_format_r2 {
            bail!("Output files R1 and R2 must have the same format (FASTA or FASTQ)");
        }

        let buf_size = compute_write_buffer_size(cfg.read_length);
        let mut buf_writer_r1 = BufWriter::with_capacity(buf_size, writer_r1);
        let mut buf_writer_r2 = BufWriter::with_capacity(buf_size, writer_r2);

        let sep = cfg.abundance_sep.as_bytes();
        let temp_path_str = temp_path.to_string_lossy().to_string();
        let mut reader = parse_fastx_file(&temp_path_str).context("Failed to read temp file")?;
        let mut format_validated = false;

        // Temp file is interleaved: R1,R2,R1,R2...
        // Extract owned data to avoid borrow conflicts.
        loop {
            let r1_data = match reader.next() {
                Some(r) => {
                    let rec = r.context("Invalid R1 in temp file")?;
                    OwnedRecord {
                        id: rec.id().to_vec(),
                        seq: rec.seq().to_vec(),
                        qual: rec.qual().map(|q| q.to_vec()),
                    }
                }
                None => break,
            };
            let r2_data = match reader.next() {
                Some(r) => {
                    let rec = r.context("Invalid R2 in temp file")?;
                    OwnedRecord {
                        id: rec.id().to_vec(),
                        seq: rec.seq().to_vec(),
                        qual: rec.qual().map(|q| q.to_vec()),
                    }
                }
                None => bail!("Desynchronization in temp file: odd number of records"),
            };

            if !format_validated {
                let has_quality = r1_data.qual.is_some();
                validate_format_compatibility(has_quality, output_format_r1)?;
                format_validated = true;
            }

            let combined_hash = T::hash_pair(&r1_data.seq, &r2_data.seq);
            let count = tracker.get_count(combined_hash);

            match output_format_r1 {
                OutputFormat::Fasta => {
                    write_fasta_record_with_abundance(
                        &mut buf_writer_r1,
                        &r1_data.id,
                        &r1_data.seq,
                        count,
                        sep,
                    )?;
                    write_fasta_record_with_abundance(
                        &mut buf_writer_r2,
                        &r2_data.id,
                        &r2_data.seq,
                        count,
                        sep,
                    )?;
                }
                OutputFormat::Fastq => {
                    let qual_r1 = r1_data.qual.as_deref().unwrap_or(&[]);
                    let qual_r2 = r2_data.qual.as_deref().unwrap_or(&[]);
                    write_fastq_record_with_abundance(
                        &mut buf_writer_r1,
                        &r1_data.id,
                        &r1_data.seq,
                        qual_r1,
                        count,
                        sep,
                    )?;
                    write_fastq_record_with_abundance(
                        &mut buf_writer_r2,
                        &r2_data.id,
                        &r2_data.seq,
                        qual_r2,
                        count,
                        sep,
                    )?;
                }
            }
        }

        buf_writer_r1
            .flush()
            .context("Failed to flush R1 output buffer")?;
        buf_writer_r2
            .flush()
            .context("Failed to flush R2 output buffer")?;

        let _ = fs::remove_file(&temp_path);
        let duplicates = total_reads - tracker.len();
        Ok((total_reads, duplicates))

    } else {
        // === CLASSIC PATH (re-read input gzip for second pass) ===

        // --- PASS 1: count abundances only ---
        {
            let mut reader_r1 = crate::stream::open_fastx(input_r1_path)?;
            let mut reader_r2 = crate::stream::open_fastx(input_r2_path)?;

            while let (Some(r1_res), Some(r2_res)) = (reader_r1.next(), reader_r2.next()) {
                let seq_r1 = r1_res.context("Invalid sequence in R1")?;
                let seq_r2 = r2_res.context("Invalid sequence in R2")?;

                let base_id_r1 = base_read_id(seq_r1.id());
                let base_id_r2 = base_read_id(seq_r2.id());
                if base_id_r1 != base_id_r2 {
                    bail!(
                        "Critical desynchronization detected at pair #{}! R1: {}, R2: {}",
                        total_reads + 1,
                        String::from_utf8_lossy(base_id_r1),
                        String::from_utf8_lossy(base_id_r2)
                    );
                }

                let combined_hash = T::hash_pair(&seq_r1.seq(), &seq_r2.seq());
                tracker.track(combined_hash);
                total_reads += 1;
            }

            // Check for trailing records (desynchronization at end)
            if reader_r1.next().is_some() || reader_r2.next().is_some() {
                bail!(
                    "Desynchronization detected at the end: one file contains more reads than the other!"
                );
            }
        }

        // --- PASS 2: re-read gzip, annotate with abundance ---
        let (writer_r1, output_format_r1) = open_output(output_r1_path, cfg)?;
        let (writer_r2, output_format_r2) = open_output(output_r2_path, cfg)?;

        if output_format_r1 != output_format_r2 {
            bail!("Output files R1 and R2 must have the same format (FASTA or FASTQ)");
        }

        let buf_size = compute_write_buffer_size(cfg.read_length);
        let mut buf_writer_r1 = BufWriter::with_capacity(buf_size, writer_r1);
        let mut buf_writer_r2 = BufWriter::with_capacity(buf_size, writer_r2);

        if cfg.verbose {
            println!(
                "Abundance tracking: {} unique / {total_reads} pairs (re-reading input)",
                tracker.len()
            );
        }

        let sep = cfg.abundance_sep.as_bytes();
        let mut reader_r1 = crate::stream::open_fastx(input_r1_path)?;
        let mut reader_r2 = crate::stream::open_fastx(input_r2_path)?;
        let mut written = rustc_hash::FxHashSet::default();
        let mut format_validated = false;

        while let (Some(r1_res), Some(r2_res)) = (reader_r1.next(), reader_r2.next()) {
            let seq_r1 = r1_res.context("Invalid sequence in R1")?;
            let seq_r2 = r2_res.context("Invalid sequence in R2")?;

            let base_id_r1 = base_read_id(seq_r1.id());
            let base_id_r2 = base_read_id(seq_r2.id());
            if base_id_r1 != base_id_r2 {
                bail!(
                    "Critical desynchronization detected at pair #{}! R1: {}, R2: {}",
                    total_reads + 1,
                    String::from_utf8_lossy(base_id_r1),
                    String::from_utf8_lossy(base_id_r2)
                );
            }

            if !format_validated {
                let has_quality = seq_r1.qual().is_some_and(|q| !q.is_empty());
                validate_format_compatibility(has_quality, output_format_r1)?;
                format_validated = true;
            }

            let combined_hash = T::hash_pair(&seq_r1.seq(), &seq_r2.seq());
            if !written.insert(combined_hash) {
                continue;
            }
            let count = tracker.get_count(combined_hash);

            match output_format_r1 {
                OutputFormat::Fasta => {
                    write_fasta_record_with_abundance(
                        &mut buf_writer_r1,
                        seq_r1.id(),
                        &seq_r1.seq(),
                        count,
                        sep,
                    )?;
                    write_fasta_record_with_abundance(
                        &mut buf_writer_r2,
                        seq_r2.id(),
                        &seq_r2.seq(),
                        count,
                        sep,
                    )?;
                }
                OutputFormat::Fastq => {
                    let qual_r1 = seq_r1.qual().unwrap_or(&[]);
                    let qual_r2 = seq_r2.qual().unwrap_or(&[]);
                    write_fastq_record_with_abundance(
                        &mut buf_writer_r1,
                        seq_r1.id(),
                        &seq_r1.seq(),
                        qual_r1,
                        count,
                        sep,
                    )?;
                    write_fastq_record_with_abundance(
                        &mut buf_writer_r2,
                        seq_r2.id(),
                        &seq_r2.seq(),
                        qual_r2,
                        count,
                        sep,
                    )?;
                }
                }
        }

        buf_writer_r1
            .flush()
            .context("Failed to flush R1 output buffer")?;
        buf_writer_r2
            .flush()
            .context("Failed to flush R2 output buffer")?;

        Ok((total_reads, total_reads - tracker.len()))
    }
}

/// Helper: write a pair of records (R1 then R2) to a temp file in raw FASTQ/FASTA.
fn write_paired_record(
    writer: &mut impl Write,
    r1: &needletail::parser::SequenceRecord,
    r2: &needletail::parser::SequenceRecord,
    is_fastq: bool,
) -> Result<()> {
    if is_fastq {
        let qual_r1 = r1.qual().unwrap_or(&[]);
        let qual_r2 = r2.qual().unwrap_or(&[]);
        writer.write_all(b"@")?;
        writer.write_all(r1.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(&r1.seq())?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(qual_r1)?;
        writer.write_all(b"\n")?;
        writer.write_all(b"@")?;
        writer.write_all(r2.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(&r2.seq())?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(qual_r2)?;
        writer.write_all(b"\n")?;
    } else {
        writer.write_all(b">")?;
        writer.write_all(r1.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(&r1.seq())?;
        writer.write_all(b"\n")?;
        writer.write_all(b">")?;
        writer.write_all(r2.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(&r2.seq())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
