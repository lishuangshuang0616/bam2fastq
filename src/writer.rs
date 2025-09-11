//! FASTQ file writing functionality

use std::collections::HashMap;
use std::fs::{create_dir, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use anyhow::{Error, Context};
use flate2::write::GzEncoder;
use shardio::helper::ThreadProxyWriter;

use crate::types::*;
use crate::formatter::FormatBamRecords;
use crate::memory_pool::ThreadSafeByteBufferPool;

/// Type alias for buffered gzip writer
type Bgw = ThreadProxyWriter<BufWriter<GzEncoder<File>>>;

/// Manager for multiple FASTQ writers
pub struct FastqManager {
    writers: HashMap<Rg, FastqWriter>,
    out_dir: PathBuf,
    buffer_pool: ThreadSafeByteBufferPool,
}

impl FastqManager {
    /// Create a new FastqManager
    pub fn new(
        out_dir: &Path,
        formatter: FormatBamRecords,
        _sample_name: String,
        reads_per_fastq: Option<usize>,
    ) -> Self {
        let buffer_pool = ThreadSafeByteBufferPool::new(100, 8192); // Pool for write buffers
        let mut sample_def_paths = HashMap::new();
        let mut writers = HashMap::new();

        for (_, &(ref _sample, lane)) in formatter.rg_spec.iter() {
            let sample = _sample.clone();
            let path = sample_def_paths.entry(sample).or_insert_with(|| {
                let suffix = _sample;
                out_dir.join(suffix)
            });

            let writer = FastqWriter::new(
                path,
                formatter.clone(),
                "bam2fastq".to_string(),
                1,
                reads_per_fastq,
            );

            writers.insert((_sample.clone(), lane), writer);
        }

        FastqManager {
            writers,
            out_dir: out_dir.to_path_buf(),
            buffer_pool,
        }
    }

    /// Write records to appropriate writer
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
            },
            None => {
                // If no RG, use the first available writer
                if let Some(w) = self.writers.values_mut().next() {
                    w.write(r1, r2, i1, i2).expect("Failed to write records");
                }
            }
        }
    }

    /// Get total number of records written
    pub fn total_written(&self) -> usize {
        self.writers.values().map(|w| w.total_written).sum()
    }

    /// Get all output file paths
    pub fn paths(&self) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {
        self.writers
            .iter()
            .flat_map(|(_, w)| w.path_sets.clone())
            .collect()
    }
}

/// FASTQ file writer for a single read group
pub struct FastqWriter {
    formatter: FormatBamRecords,
    out_dir: PathBuf,
    _sample_name: String,
    lane: u32,
    
    r1: Option<Bgw>,
    r2: Option<Bgw>,
    i1: Option<Bgw>,
    i2: Option<Bgw>,

    chunk_written: usize,
    pub total_written: usize,
    n_chunks: usize,
    reads_per_fastq: Option<usize>,
    pub path_sets: Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>,
    // Pre-allocated write buffer for better performance
    write_buffer: Vec<u8>,
}

impl FastqWriter {
    /// Create a new FastqWriter
    pub fn new(
        out_dir: &Path,
        formatter: FormatBamRecords,
        sample_name: String,
        lane: u32,
        reads_per_fastq: Option<usize>,
    ) -> Self {
        Self {
            formatter,
            out_dir: out_dir.to_path_buf(),
            _sample_name: sample_name,
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
            write_buffer: Vec::with_capacity(8192),
        }
    }

    /// Get file paths for current chunk
    fn get_paths(
        out_dir: &Path,
        sample_name: &str,
        lane: u32,
        n_files: usize,
        formatter: &FormatBamRecords,
    ) -> (PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>) {
        if formatter.rename.is_none() {
            let r1 = out_dir.join(format!(
                "{}_L{:02}_{:01}_1.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1
            ));
            let r2 = out_dir.join(format!(
                "{}_L{:02}_{:01}_2.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1
            ));
            let i1 = out_dir.join(format!(
                "{}_L{:02}_{:01}_I1.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1
            ));
            let i2 = out_dir.join(format!(
                "{}_L{:02}_{:01}_I2.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1
            ));
            (
                r1, 
                r2,
                if !formatter.i1_spec.is_empty() {
                    Some(i1)
                } else {
                    None
                },
                if !formatter.i2_spec.is_empty() {
                    Some(i2)
                } else {
                    None
                },
            )
        } else {
            let rename = formatter.rename.as_ref().unwrap();
            let r1 = out_dir.join(format!(
                "{}_L{:02}_{:01}_{}.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1,
                rename[0]
            ));
            let r2 = out_dir.join(format!(
                "{}_L{:02}_{:01}_{}.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1,
                rename[1]
            ));
            let i1 = if rename.len() > 2 && !formatter.i1_spec.is_empty() {
                Some(out_dir.join(format!(
                    "{}_L{:02}_{:01}_{}.fastq.gz", 
                    sample_name, 
                    lane,
                    n_files + 1,
                    rename[2]
                )))
            } else {
                None
            };
            let i2 = if rename.len() > 3 && !formatter.i2_spec.is_empty() {
                Some(out_dir.join(format!(
                    "{}_L{:02}_{:01}_{}.fastq.gz", 
                    sample_name, 
                    lane,
                    n_files + 1,
                    rename[3]
                )))
            } else {
                None
            };
            (r1, r2, i1, i2)
        }
    }

    /// Write a single FASTQ record
    pub fn write_rec(
        w: &mut Bgw,
        rec: &FqRecord,
    ) -> Result<(), Error> {
        // Pre-calculate total size to minimize allocations
        let total_size = 1 + rec.head.len() + rec.seq.len() + rec.qual.len() + 5; // +1 for @, +5 for newlines and "+"
        let mut buffer = Vec::with_capacity(total_size);
        
        buffer.push(b'@');
        buffer.extend_from_slice(&rec.head);
        buffer.push(b'\n');
        buffer.extend_from_slice(&rec.seq);
        buffer.extend_from_slice(b"\n+\n");
        buffer.extend_from_slice(&rec.qual);
        buffer.push(b'\n');
        
        w.write_all(&buffer)?;
        Ok(())
    }

    /// Write optional FASTQ record
    pub fn write_fq_record(
        w: &mut Option<Bgw>,
        rec: &Option<FqRecord>,
    ) -> Result<(), Error> {
        match (w, rec) {
            (Some(ref mut writer), Some(ref record)) => {
                Self::write_rec(writer, record)?;
            }
            (None, None) => {}
            _ => {
                return Err(anyhow::anyhow!("Mismatch between writer and record availability"));
            }
        }
        Ok(())
    }

    /// Write mandatory FASTQ record
    pub fn write_mandatory_fq_record(
        w: &mut Option<Bgw>,
        rec: &FqRecord,
    ) -> Result<(), Error> {
        match w {
            Some(ref mut writer) => Self::write_rec(writer, rec),
            None => Err(anyhow::anyhow!("Writer not available for mandatory record")),
        }
    }

    /// Write records to files
    fn write(
        &mut self,
        r1: &FqRecord,
        r2: &FqRecord,
        i1: &Option<FqRecord>,
        i2: &Option<FqRecord>,
    ) -> Result<(), Error> {
        if self.reads_per_fastq.is_some() && self.chunk_written >= self.reads_per_fastq.unwrap() {
            self.cycle_writers()?;
        }

        if self.r1.is_none() {
            self.cycle_writers()?;
        }

        // Use optimized batch writing when possible
        Self::write_mandatory_fq_record(&mut self.r1, r1)?;
        Self::write_mandatory_fq_record(&mut self.r2, r2)?;
        Self::write_fq_record(&mut self.i1, i1)?;
        Self::write_fq_record(&mut self.i2, i2)?;

        self.chunk_written += 1;
        self.total_written += 1;

        Ok(())
    }

    /// Cycle to new output files
    fn cycle_writers(&mut self) -> Result<(), Error> {
        let paths = Self::get_paths(
            &self.out_dir,
            &self._sample_name,
            self.lane,
            self.n_chunks,
            &self.formatter,
        );

        if let Some(parent) = paths.0.parent() {
            if !parent.exists() {
                create_dir(parent).with_context(|| format!("Failed to create directory: {:?}", parent))?;
            }
        }

        self.r1 = Some(Self::open_gzip_writer(&paths.0)?);
        self.r2 = Some(Self::open_gzip_writer(&paths.1)?);
        self.i1 = paths.2.as_ref().map(Self::open_gzip_writer).transpose()?;
        self.i2 = paths.3.as_ref().map(Self::open_gzip_writer).transpose()?;

        self.path_sets.push(paths);
        self.chunk_written = 0;
        self.n_chunks += 1;

        Ok(())
    }

    /// Open a gzip writer for the given path
    fn open_gzip_writer<P: AsRef<Path>>(
        path: P
    ) -> Result<ThreadProxyWriter<BufWriter<GzEncoder<File>>>, Error> {
        let file = File::create(&path)
            .with_context(|| format!("Failed to create file: {:?}", path.as_ref()))?;
        let gz = GzEncoder::new(file, flate2::Compression::default());
        let buf = BufWriter::new(gz);
        Ok(ThreadProxyWriter::new(buf, 1 << 21))
    }
}