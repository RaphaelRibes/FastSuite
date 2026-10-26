use crate::hasher::{HashVerifier, SequenceHasher};
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use needletail::parse_fastx_file;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::Path;

/// Supported output formats for sequences.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OutputFormat {
    Fasta,
    Fastq,
}

impl OutputFormat {
    /// Determines the output format based on the file extension.
    pub fn from_extension(path: &Path) -> Self {
        let path_str = path.to_string_lossy().to_lowercase();

        let without_gz = path_str.strip_suffix(".gz").unwrap_or(&path_str);

        if without_gz.ends_with(".fasta")
            || without_gz.ends_with(".fa")
            || without_gz.ends_with(".fna")
        {
            OutputFormat::Fasta
        } else {
            OutputFormat::Fastq
        }
    }

    /// Checks if the file path indicates a gzipped file.
    pub fn is_gz(path: &Path) -> bool {
        path.to_string_lossy().to_lowercase().ends_with(".gz")
    }
}

/// Compression decision for stream output (D-4): a stream has no
/// filename, so `-o -` cannot infer gzip from an extension.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StreamCompression {
    None,
    Gz,
}

/// Open stdout as an output stream. The record format is chosen per
/// record by the caller from `--out-format`; this only decides the bytes'
/// compression. Gzip finishes on drop, exactly like the file path.
pub fn prepare_stream_writer(
    compression: StreamCompression,
    compression_level: u32,
) -> Box<dyn Write> {
    const STREAM_BUFFER: usize = 1 << 20; // 1 MiB, same as the file path
    match compression {
        StreamCompression::None => {
            Box::new(io::BufWriter::with_capacity(STREAM_BUFFER, io::stdout()))
        }
        StreamCompression::Gz => Box::new(GzEncoder::new(
            io::BufWriter::with_capacity(STREAM_BUFFER, io::stdout()),
            Compression::new(compression_level),
        )),
    }
}

/// Prepares a writer for the output file based on the format and whether it should
/// be forced (overwritten) or appended.
///
/// # Errors
///
/// Returns an error if:
/// - Append mode is requested on a `.gz` file (multi-member gzip is unreliable
///   for downstream bioinformatics tools).
/// - The file cannot be opened or created.
pub fn prepare_writer(
    path: &Path,
    force: bool,
    compression_level: u32,
) -> Result<(Box<dyn Write>, OutputFormat)> {
    // Create the parent directories a named output needs. Every other tool in
    // the suite does this (fastprep, fastcheck, fastmaps, fastsim), and
    // `-o results/run/o.fq` here failed with "Failed to create file" on a
    // directory the caller had not made by hand: same flag, same meaning, one
    // tool refusing. A stream never reaches this function.
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating {} for {}", parent.display(), path.display()))?;
        }
    }

    let format = OutputFormat::from_extension(path);
    let is_gz = OutputFormat::is_gz(path);

    let file = if path.exists() && !force {
        if is_gz {
            bail!(
                "Cannot append to gzipped output {:?}. Multi-member gzip is unreliable for \
                 downstream tools. Use --force to overwrite, or use an uncompressed output \
                 format for incremental/resumable runs.",
                path
            );
        }
        OpenOptions::new()
            .append(true)
            .open(path)
            .with_context(|| format!("Failed to open file in append mode: {:?}", path))?
    } else {
        File::create(path).with_context(|| format!("Failed to create file: {:?}", path))?
    };

    let writer: Box<dyn Write> = if is_gz {
        // FD-7: deflate emits ~32 KiB units, so an unbuffered `File` here costs
        // one write syscall (and a possible short-write memmove) per block.
        // The plain branch is buffered by the caller; this branch was not.
        Box::new(GzEncoder::new(
            io::BufWriter::with_capacity(1 << 20, file),
            Compression::new(compression_level),
        ))
    } else {
        Box::new(file)
    };

    Ok((writer, format))
}

/// True when an existing output is too short to hold a FASTX record.
///
/// The output file is created before the first record is written, so a run that
/// crashes — or that is refused after `File::create` — leaves a 0-byte file
/// behind. That is the common resume artifact, and there is nothing in it to
/// preload. Two bytes is needletail's minimum: it sniffs the first two.
pub(crate) fn preload_is_empty(path: &Path) -> Result<bool> {
    Ok(fs::metadata(path)?.len() < 2)
}

/// Opens an existing output for the resume preload pass.
///
/// `what` names the file in the message ("preload" for a single output, "R1
/// preload" for a mate). A file that is not FASTX at all gets a message naming
/// the fix instead of needletail's bare parse error.
pub(crate) fn open_preload_reader(
    path: &str,
    what: &str,
) -> Result<Box<dyn needletail::FastxReader>> {
    parse_fastx_file(path).with_context(|| {
        format!("Error opening {what} file {path}: not a FASTX file; use --force to overwrite it")
    })
}

/// Preloads existing hashes into memory from an existing output file (Single-End).
pub fn preload_existing_hashes<T: SequenceHasher>(
    path: &str,
    verifier: &mut HashVerifier<T>,
    verbose: bool,
) -> Result<(usize, u64)> {
    if !Path::new(path).exists() {
        return Ok((0, 0));
    }

    // A 0- or 1-byte output holds no record: preload nothing and let the caller
    // truncate it (it does when `valid_bytes < current_size`), so the append
    // starts on a record boundary.
    if preload_is_empty(Path::new(path))? {
        return Ok((0, 0));
    }

    if verbose {
        println!("Preloading sequences from existing output...");
    }

    let mut reader = open_preload_reader(path, "preload")?;
    let mut count = 0;
    let mut valid_bytes = ByteCounter(0);

    while let Some(record_result) = reader.next() {
        match record_result {
            Ok(record) => {
                let hash = T::hash_sequence(&record.seq());
                verifier.verify(hash);

                count += 1;
                let _ = record.write(&mut valid_bytes, None);
            }
            Err(e) => {
                if verbose {
                    eprintln!("Incomplete sequence detected at end of file ({}).", e);
                    eprintln!("Calculating fail-safe truncation point...");
                }
                break;
            }
        }
    }

    Ok((count, valid_bytes.0))
}

/// Preloads existing hashes into memory from existing output files (Paired-End)
/// and computes valid byte sizes to synchronize truncation.
pub fn preload_existing_paired_hashes<T: SequenceHasher>(
    path_r1: &str,
    path_r2: &str,
    verifier: &mut HashVerifier<T>,
    verbose: bool,
) -> Result<(usize, u64, u64)> {
    let path_r1_ref = Path::new(path_r1);
    let path_r2_ref = Path::new(path_r2);

    if !path_r1_ref.exists() || !path_r2_ref.exists() {
        return Ok((0, 0, 0));
    }

    // Either mate being empty means no pair can be preloaded: return nothing,
    // and the caller truncates both to 0 — the same state as a fresh run.
    if preload_is_empty(path_r1_ref)? || preload_is_empty(path_r2_ref)? {
        return Ok((0, 0, 0));
    }

    if verbose {
        println!("Preloading and synchronizing pairs from existing outputs...");
    }

    let mut reader_r1 = open_preload_reader(path_r1, "R1 preload")?;
    let mut reader_r2 = open_preload_reader(path_r2, "R2 preload")?;

    let mut count = 0;
    let mut valid_bytes_r1 = ByteCounter(0);
    let mut valid_bytes_r2 = ByteCounter(0);

    while let (Some(record_r1_res), Some(record_r2_res)) = (reader_r1.next(), reader_r2.next()) {
        match (record_r1_res, record_r2_res) {
            (Ok(record_r1), Ok(record_r2)) => {
                let combined_hash = T::hash_pair(&record_r1.seq(), &record_r2.seq());
                verifier.verify(combined_hash);

                count += 1;
                let _ = record_r1.write(&mut valid_bytes_r1, None);
                let _ = record_r2.write(&mut valid_bytes_r2, None);
            }
            _ => {
                if verbose {
                    eprintln!("Incomplete sequence or desynchronization detected at end of file.");
                    eprintln!("Calculating fail-safe truncation points for R1 and R2...");
                }
                break;
            }
        }
    }

    Ok((count, valid_bytes_r1.0, valid_bytes_r2.0))
}

/// A no-op writer that counts bytes written, used to compute valid byte offsets
/// without actually performing I/O.
pub(crate) struct ByteCounter(pub u64);

impl Write for ByteCounter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0 += buf.len() as u64;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Writes a single FASTA record using raw bytes.
#[inline]
pub fn write_fasta_record(writer: &mut impl Write, id: &[u8], seq: &[u8]) -> Result<()> {
    writer
        .write_all(b">")
        .context("Error writing FASTA '>' prefix")?;
    writer.write_all(id).context("Error writing FASTA ID")?;
    writer
        .write_all(b"\n")
        .context("Error writing FASTA ID newline")?;
    writer
        .write_all(seq)
        .context("Error writing FASTA sequence")?;
    writer
        .write_all(b"\n")
        .context("Error writing FASTA sequence newline")?;
    Ok(())
}

/// Writes a single FASTA record with an abundance annotation suffix on the header.
#[inline]
pub fn write_fasta_record_with_abundance(
    writer: &mut impl Write,
    id: &[u8],
    seq: &[u8],
    abundance: u32,
    sep: &[u8],
) -> Result<()> {
    writer
        .write_all(b">")
        .context("Error writing FASTA '>' prefix")?;
    writer.write_all(id).context("Error writing FASTA ID")?;
    writer
        .write_all(sep)
        .context("Error writing abundance separator")?;
    let size = abundance.to_string();
    writer
        .write_all(size.as_bytes())
        .context("Error writing abundance value")?;
    writer
        .write_all(b"\n")
        .context("Error writing FASTA ID newline")?;
    writer
        .write_all(seq)
        .context("Error writing FASTA sequence")?;
    writer
        .write_all(b"\n")
        .context("Error writing FASTA sequence newline")?;
    Ok(())
}

/// Writes a single FASTQ record with an abundance annotation suffix on the header.
#[inline]
pub fn write_fastq_record_with_abundance(
    writer: &mut impl Write,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    abundance: u32,
    sep: &[u8],
) -> Result<()> {
    writer
        .write_all(b"@")
        .context("Error writing FASTQ '@' prefix")?;
    writer.write_all(id).context("Error writing FASTQ ID")?;
    writer
        .write_all(sep)
        .context("Error writing abundance separator")?;
    let size = abundance.to_string();
    writer
        .write_all(size.as_bytes())
        .context("Error writing abundance value")?;
    writer
        .write_all(b"\n")
        .context("Error writing FASTQ ID newline")?;
    writer
        .write_all(seq)
        .context("Error writing FASTQ sequence")?;
    writer
        .write_all(b"\n+\n")
        .context("Error writing FASTQ '+' line")?;
    writer
        .write_all(qual)
        .context("Error writing FASTQ quality")?;
    writer
        .write_all(b"\n")
        .context("Error writing FASTQ quality newline")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- OutputFormat ---

    #[test]
    fn test_output_format_fastq_default() {
        assert_eq!(
            OutputFormat::from_extension(Path::new("out.fastq")),
            OutputFormat::Fastq
        );
        assert_eq!(
            OutputFormat::from_extension(Path::new("out.fq.gz")),
            OutputFormat::Fastq
        );
        assert_eq!(
            OutputFormat::from_extension(Path::new("out.unknown")),
            OutputFormat::Fastq
        );
    }

    #[test]
    fn test_output_format_fasta() {
        assert_eq!(
            OutputFormat::from_extension(Path::new("out.fasta")),
            OutputFormat::Fasta
        );
        assert_eq!(
            OutputFormat::from_extension(Path::new("out.fa")),
            OutputFormat::Fasta
        );
        assert_eq!(
            OutputFormat::from_extension(Path::new("out.fna.gz")),
            OutputFormat::Fasta
        );
    }

    #[test]
    fn test_is_gz() {
        assert!(OutputFormat::is_gz(Path::new("file.fastq.gz")));
        assert!(!OutputFormat::is_gz(Path::new("file.fastq")));
    }
}
