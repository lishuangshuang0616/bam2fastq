//! writer.rs — FASTQ output writing
//!
//! Contains `GenWriter` (gzip or raw output), `FastqWriter` (per-lane file management),
//! and `FastqManager` (dispatches writes to the appropriate writer by read-group).

use crate::formatter::{FormatBamRecords, FqRecord, Rg};
use anyhow::{anyhow, Error};
use gzp::{
    deflate::Gzip,
    par::compress::{ParCompress, ParCompressBuilder},
};
use once_cell::sync::Lazy;
use shardio::helper::ThreadProxyWriter;
use std::collections::HashMap;
use std::fs::{create_dir, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

// ─── Available system memory (cached once at startup) ─────────────────────

/// Available system memory at startup (bytes).
pub static AVAILABLE_MEMORY: Lazy<usize> = Lazy::new(get_available_memory);

/// Public accessor for available system memory (used by go() to size caches).
pub fn get_memory() -> usize {
    get_available_memory()
}

#[cfg(target_os = "macos")]
fn get_available_memory() -> usize {
    use std::process::Command;
    let output = Command::new("vm_stat")
        .output()
        .expect("Failed to execute vm_stat");
    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        for line in stdout.lines() {
            if line.contains("Pages free") || line.contains("Pages inactive") {
                let parts: Vec<&str> = line.split(':').collect();
                if parts.len() > 1 {
                    let pages: usize = parts[1].trim().trim_end_matches('.').parse().unwrap_or(0);
                    return pages * 4096;
                }
            }
        }
    }
    1024 * 1024 * 1024
}

#[cfg(target_os = "linux")]
fn get_available_memory() -> usize {
    use std::fs;
    let meminfo = fs::read_to_string("/proc/meminfo").unwrap_or_default();
    for line in meminfo.lines() {
        if line.starts_with("MemAvailable:") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() > 1 {
                let kb: usize = parts[1].parse().unwrap_or(0);
                return kb * 1024;
            }
        }
    }
    1024 * 1024 * 1024
}

#[cfg(not(any(target_os = "macos", target_os = "linux")))]
fn get_available_memory() -> usize {
    1024 * 1024 * 1024
}

// ─── GenWriter ────────────────────────────────────────────────────────────

/// Output writer that is either parallel-gzip compressed or raw uncompressed.
///
/// The `Gz` variant wraps a `BufWriter<ParCompress<Gzip>>`:
///   1MB BufWriter → ParCompress (N parallel 64KB blocks) → output file
///
/// The `Raw` variant uses `ThreadProxyWriter<BufWriter<File>>` for uncompressed output.
pub enum GenWriter {
    Gz(BufWriter<ParCompress<Gzip>>),
    Raw(ThreadProxyWriter<BufWriter<File>>),
}

impl Write for GenWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            GenWriter::Gz(bw) => bw.write(buf),
            GenWriter::Raw(w) => w.write(buf),
        }
    }

    fn write_vectored(&mut self, bufs: &[std::io::IoSlice<'_>]) -> std::io::Result<usize> {
        match self {
            GenWriter::Gz(bw) => bw.write_vectored(bufs),
            GenWriter::Raw(w) => w.write_vectored(bufs),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            GenWriter::Gz(bw) => bw.flush(),
            GenWriter::Raw(w) => w.flush(),
        }
    }
}

/// Alias for `GenWriter` used at call sites.
pub type Bgw = GenWriter;

// ─── FastqManager ─────────────────────────────────────────────────────────

/// Manages a collection of `FastqWriter` instances, one per read-group.
/// Dispatches `write()` calls to the correct writer based on the `Rg` tag.
pub struct FastqManager {
    pub writers: HashMap<Rg, FastqWriter>,
    pub out_dir: PathBuf,
}

impl FastqManager {
    pub fn new(
        out_dir: &Path,
        formatter: FormatBamRecords,
        _sample_name: String,
        reads_per_fastq: Option<usize>,
        no_compress: bool,
    ) -> Self {
        let mut sample_def_paths: HashMap<String, PathBuf> = HashMap::new();
        let mut writers = HashMap::new();

        for (_, &(ref _sample, lane)) in formatter.rg_spec.iter() {
            let sample = _sample.clone();
            let path = sample_def_paths
                .entry(sample)
                .or_insert_with(|| out_dir.join(_sample));

            let writer = FastqWriter::new(
                path,
                formatter.clone(),
                "bam2fastq".to_string(),
                1,
                reads_per_fastq,
                no_compress,
            );
            writers.insert((_sample.clone(), lane), writer);
        }

        FastqManager {
            writers,
            out_dir: out_dir.to_path_buf(),
        }
    }

    /// Write a read pair directly to the appropriate `FastqWriter`.
    pub fn write(
        &mut self,
        rg: &Option<Rg>,
        r1: &FqRecord,
        r2: &FqRecord,
        i1: &Option<FqRecord>,
        i2: &Option<FqRecord>,
    ) {
        match rg {
            Some(ref rg) => {
                if let Some(w) = self.writers.get_mut(rg) {
                    w.write(r1, r2, i1, i2).expect("Failed to write records");
                }
            }
            None => {
                if let Some(w) = self.writers.values_mut().next() {
                    w.write(r1, r2, i1, i2).expect("Failed to write records");
                }
            }
        }
    }

    /// No-op kept for call-site compatibility; writes are now synchronous.
    #[inline]
    pub fn flush_buffer(&mut self) {}

    pub fn total_written(&self) -> usize {
        self.writers.values().map(|w| w.total_written).sum()
    }

    pub fn paths(&self) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {
        self.writers
            .values()
            .flat_map(|w| w.path_sets.clone())
            .collect()
    }
}

// ─── FastqWriter ──────────────────────────────────────────────────────────

/// Manages the actual file handles for one read-group's FASTQ output.
/// Supports optional chunking via `reads_per_fastq`.
pub struct FastqWriter {
    formatter: FormatBamRecords,
    out_dir: PathBuf,
    sample_name: String,
    lane: u32,

    pub r1: Option<Bgw>,
    pub r2: Option<Bgw>,
    pub i1: Option<Bgw>,
    pub i2: Option<Bgw>,

    chunk_written: usize,
    pub total_written: usize,
    n_chunks: usize,
    reads_per_fastq: Option<usize>,
    pub path_sets: Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>,
    no_compress: bool,
}

impl FastqWriter {
    pub fn new(
        out_dir: &Path,
        formatter: FormatBamRecords,
        sample_name: String,
        lane: u32,
        reads_per_fastq: Option<usize>,
        no_compress: bool,
    ) -> Self {
        Self {
            formatter,
            out_dir: out_dir.to_path_buf(),
            sample_name,
            lane,
            r1: None,
            r2: None,
            i1: None,
            i2: None,
            chunk_written: 0,
            total_written: 0,
            n_chunks: 0,
            reads_per_fastq,
            path_sets: vec![],
            no_compress,
        }
    }

    fn get_paths(
        out_dir: &Path,
        sample_name: &str,
        lane: u32,
        n_files: usize,
        formatter: &FormatBamRecords,
        no_compress: bool,
    ) -> (PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>) {
        let ext = if no_compress { "fastq" } else { "fastq.gz" };
        let mk = |suffix: &str| {
            out_dir.join(format!(
                "{}_L{:02}_{:01}_{}.{}",
                sample_name,
                lane,
                n_files + 1,
                suffix,
                ext
            ))
        };

        if let Some(ref names) = formatter.rename {
            (
                mk(&names[0]),
                mk(&names[1]),
                if !formatter.i1_spec.is_empty() {
                    Some(mk(&names[2]))
                } else {
                    None
                },
                if !formatter.i2_spec.is_empty() {
                    Some(mk(&names[3]))
                } else {
                    None
                },
            )
        } else {
            let mk_default = |suffix: &str| {
                out_dir.join(format!(
                    "{}_L{:02}_{:01}_{}.{}",
                    sample_name,
                    lane,
                    n_files + 1,
                    suffix,
                    ext
                ))
            };
            (
                mk_default("1"),
                mk_default("2"),
                if !formatter.i1_spec.is_empty() {
                    Some(mk_default("I1"))
                } else {
                    None
                },
                if !formatter.i2_spec.is_empty() {
                    Some(mk_default("I2"))
                } else {
                    None
                },
            )
        }
    }

    /// Write a FASTQ record using `write_vectored` for maximum IO throughput.
    /// Falls back to `write_all` on partial writes (rare with BufWriter backing).
    pub fn write_rec(w: &mut Bgw, rec: &FqRecord) -> Result<(), Error> {
        use std::io::IoSlice;
        let mut slices: [IoSlice<'_>; 7] = [
            IoSlice::new(b"@"),
            IoSlice::new(&rec.head),
            IoSlice::new(b"\n"),
            IoSlice::new(&rec.seq),
            IoSlice::new(b"\n+\n"),
            IoSlice::new(&rec.qual),
            IoSlice::new(b"\n"),
        ];
        let mut bufs = &mut slices[..];
        while !bufs.is_empty() {
            let n = w.write_vectored(bufs)?;
            let mut remaining = n;
            while remaining > 0 && !bufs.is_empty() {
                if remaining >= bufs[0].len() {
                    remaining -= bufs[0].len();
                    bufs = &mut bufs[1..];
                } else {
                    let offset = remaining;
                    w.write_all(&bufs[0][offset..])?;
                    bufs = &mut bufs[1..];
                    remaining = 0;
                }
            }
        }
        Ok(())
    }

    pub fn try_write_rec(w: &mut Option<Bgw>, rec: &Option<FqRecord>) -> Result<(), Error> {
        if let Some(ref mut w) = w {
            if let Some(r) = rec {
                FastqWriter::write_rec(w, r)?;
            } else {
                return Err(anyhow!("No record to write"));
            }
        }
        Ok(())
    }

    pub fn try_write_recs(w: &mut Option<Bgw>, rec: &FqRecord) -> Result<(), Error> {
        if let Some(ref mut w) = w {
            FastqWriter::write_rec(w, rec)?;
        }
        Ok(())
    }

    pub fn write(
        &mut self,
        r1: &FqRecord,
        r2: &FqRecord,
        i1: &Option<FqRecord>,
        i2: &Option<FqRecord>,
    ) -> Result<(), Error> {
        if self.total_written == 0 {
            let _ = create_dir(&self.out_dir);
            self.cycle_writers();
        }

        FastqWriter::try_write_recs(&mut self.r1, r1)?;
        FastqWriter::try_write_recs(&mut self.r2, r2)?;
        FastqWriter::try_write_rec(&mut self.i1, i1)?;
        FastqWriter::try_write_rec(&mut self.i2, i2)?;

        self.chunk_written += 1;
        self.total_written += 1;

        if let Some(max_reads) = self.reads_per_fastq {
            if self.chunk_written == max_reads {
                self.cycle_writers();
            }
        }
        Ok(())
    }

    fn cycle_writers(&mut self) {
        let paths = Self::get_paths(
            &self.out_dir,
            &self.sample_name,
            self.lane,
            self.n_chunks,
            &self.formatter,
            self.no_compress,
        );
        self.r1 = Some(Self::open_writer(&paths.0, self.no_compress));
        self.r2 = Some(Self::open_writer(&paths.1, self.no_compress));
        self.i1 = paths
            .2
            .as_ref()
            .map(|p| Self::open_writer(p, self.no_compress));
        self.i2 = paths
            .3
            .as_ref()
            .map(|p| Self::open_writer(p, self.no_compress));
        self.n_chunks += 1;
        self.chunk_written = 0;
        self.path_sets.push(paths);
    }

    #[allow(dead_code)]
    pub fn close_current_writers(&mut self) {
        Self::close_writer(self.r1.take());
        Self::close_writer(self.r2.take());
        Self::close_writer(self.i1.take());
        Self::close_writer(self.i2.take());
    }

    #[allow(dead_code)]
    fn close_writer(w: Option<GenWriter>) {
        if let Some(mut writer) = w {
            let _ = writer.flush();
            drop(writer);
        }
    }

    fn open_writer<P: AsRef<Path>>(path: P, no_compress: bool) -> GenWriter {
        let file = File::create(path).expect("Failed to create output file");

        let available_memory = *AVAILABLE_MEMORY;
        let out_buf_size = if available_memory < 2 * 1024 * 1024 * 1024 {
            1 << 21 // 2MB
        } else {
            1 << 23 // 8MB
        };

        if no_compress {
            GenWriter::Raw(ThreadProxyWriter::new(
                BufWriter::with_capacity(out_buf_size, file),
                out_buf_size / 4,
            ))
        } else {
            let in_buf_size = 1 << 20; // 1MB input buffer to gzp
                                       // gzp is IO-bound: use threads/4 to keep CPU budget for Rayon format work
            let gzp_threads = (rayon::current_num_threads() / 4).max(1);
            let gz: ParCompress<Gzip> = ParCompressBuilder::new()
                .num_threads(gzp_threads)
                .expect("Failed to initialize parallel gzip compressor")
                .from_writer(BufWriter::with_capacity(out_buf_size, file));
            GenWriter::Gz(BufWriter::with_capacity(in_buf_size, gz))
        }
    }
}
