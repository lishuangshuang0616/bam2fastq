use std::backtrace;
use std::io::{BufWriter, Write};
use std::{borrow::Cow, collections::HashMap, path::PathBuf, str::FromStr};
use rpcache::RpCache;
use serde::{Serialize, Deserialize};
use shardio::helper::ThreadProxyWriter;
use shardio::SortKey;
use shardio::{ShardReader, ShardWriter};
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam::{self, Read};
use regex::Regex;
use itertools::Itertools;
use anyhow::{Error, anyhow, Context};
use flate2::write::GzEncoder;
use std::fs::{create_dir, File};
use std::panic;
use std::path::Path;
use std::str;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use tempfile;
use once_cell::sync::Lazy;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

// Function to get available system memory in bytes
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
                    // Page size is typically 4096 bytes on macOS
                    return pages * 4096;
                }
            }
        }
    }
    1024 * 1024 * 1024 // Default to 1GB if we can't determine
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
                return kb * 1024; // Convert KB to bytes
            }
        }
    }
    1024 * 1024 * 1024 // Default to 1GB if we can't determine
}

// Fallback for other platforms
#[cfg(not(any(target_os = "macos", target_os = "linux")))]
fn get_available_memory() -> usize {
    1024 * 1024 * 1024 // Default to 1GB
}

static RG_REGEX: Lazy<Regex> = Lazy::new(|| Regex::new(r"^([0-9]+)-[0-9A-F]+$").unwrap());

mod bx_index;
mod locus;
mod rpcache;
use bx_index::BxListIter;

type OutPaths = (
    PathBuf,
    PathBuf,
    Option<PathBuf>,
    Option<PathBuf>,
);

type Rg = (String, u32);

type FormattedReadPair = (
    Option<Rg>,
    FqRecord,
    FqRecord,
    Option<FqRecord>,
    Option<FqRecord>,
);


#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct FqRecord {
    #[serde(with = "serde_bytes")]
    head: Vec<u8>,
    #[serde(with = "serde_bytes")]
    seq: Vec<u8>,
    #[serde(with = "serde_bytes")]
    qual: Vec<u8>,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
enum ReadNum {
    R1,
    R2
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SerFq {
    read_group: Option<Rg>,
    #[serde(with = "serde_bytes")]
    header_key: Vec<u8>,
    rec: FqRecord,
    read_num: ReadNum,
    i1: Option<FqRecord>,
    i2: Option<FqRecord>,
}

struct SerFqSort;

impl SortKey<SerFq> for SerFqSort {
    type Key = Vec<u8>;
    fn sort_key(t: &SerFq) -> Cow<Vec<u8>> {
        Cow::Borrowed(&t.header_key)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SpecEntry {
    Tags(String, String),
    #[allow(dead_code)]
    Ns(usize),
    Read,
}

#[derive(Debug, Clone)]
struct FormatBamRecords {
    rg_spec: HashMap<String, Rg>,
    r1_spec: Vec<SpecEntry>,
    r2_spec: Vec<SpecEntry>,
    i1_spec: Vec<SpecEntry>,
    i2_spec: Vec<SpecEntry>,
    rename: Option<Vec<String>>,
    order: [u32; 4],
}

pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        _ => panic!("invalid nucleotide: {}", b as char),
    }
}

impl FormatBamRecords {

    /// Automatically detect if BAM contains paired-end reads by sampling records
    pub fn detect_paired_end<R: bam::Read>(reader: &mut R) -> Result<bool, Error> {
        let mut paired_count = 0;
        let mut total_count = 0;
        let max_samples = 100_000; // Sample up to 1000 reads for detection
        
        for result in reader.records() {
            let record = result.context("Failed to read BAM record")?;
            total_count += 1;
            
            // Check if read is paired (flag 0x1) and properly paired (flag 0x2)
            if record.is_paired() {
                paired_count += 1;
            }
            
            // Stop sampling after enough records
            if total_count >= max_samples {
                break;
            }
        }
        
        // Consider paired-end if more than 50% of sampled reads are paired
        let paired_ratio = if total_count > 0 {
            paired_count as f64 / total_count as f64
        } else {
            0.0
        };
        
        Ok(paired_ratio > 0.5)
    }

    /// Create C4 configuration with automatic paired-end detection
    pub fn c4head_auto<R: bam::Read>(reader: &mut R) -> Result<FormatBamRecords, Error> {
        let is_paired = Self::detect_paired_end(reader)?;
        
        if is_paired {
            // For true paired-end data, both R1 and R2 contain read sequences
            Ok(FormatBamRecords {
                rg_spec: Self::parse_rgs(reader),
                r1_spec: vec![
                    SpecEntry::Tags("CR".to_string(), "CY".to_string()),
                    SpecEntry::Tags("UR".to_string(), "UY".to_string()),
                    SpecEntry::Read  // Include read sequence for R1
                ],
                r2_spec: vec![SpecEntry::Read],
                i1_spec: vec![],
                i2_spec: vec![],
                rename: None,
                order: [1, 2, 0, 0],
            })
        } else {
            // For single-end data, only R2 contains read sequence
            Ok(FormatBamRecords {
                rg_spec: Self::parse_rgs(reader),
                r1_spec: vec![
                    SpecEntry::Tags("CR".to_string(), "CY".to_string()),
                    SpecEntry::Tags("UR".to_string(), "UY".to_string())
                ],
                r2_spec: vec![SpecEntry::Read],
                i1_spec: vec![],
                i2_spec: vec![],
                rename: None,
                order: [1, 2, 0, 0],
            })
        }
    }

    pub fn c4head<R: bam::Read>(reader: &R) -> FormatBamRecords {
        FormatBamRecords {
            rg_spec: Self::parse_rgs(reader),
            r1_spec: vec![
                SpecEntry::Tags("CR".to_string(), "CY".to_string()),
                SpecEntry::Tags("UR".to_string(), "UY".to_string())
            ],
            r2_spec: vec![SpecEntry::Read],
            i1_spec: vec![],
            i2_spec: vec![],
            rename: None,
            order: [1, 2, 0, 0],
        }
    }

    fn parse_rgs<R: bam::Read>(reader: &R) -> HashMap<String, Rg> {
        let text = match std::str::from_utf8(reader.header().as_bytes()) {
            Ok(t) => t,
            Err(_) => return HashMap::new(), // Return empty map on UTF-8 error
        };

        let mut rg_items = text
            .lines()
            .filter(|l| l.starts_with("@RG"))
            .filter_map(Self::parse_rg_line)
            .collect::<HashMap<_, _>>();

        if rg_items.is_empty() {
            for _i in 1..2 {
                let name = format!("bam2fastq_output");
                rg_items.insert(name.clone(), (name, 0));
            }
        }
        //println!("{:?}",rg_items);

        rg_items
    }

    fn parse_rg_line(line: &str) -> Option<(String, (String, u32))> {
        let mut entries = line.split('\t');
        entries.next()?; // consume @RG entry

        let mut tags = entries
            .map(|entry| entry.split_once(':').ok_or_else(|| anyhow!("Invalid RG entry format")))
            .collect::<Result<HashMap<_, _>, _>>()
            .ok()?;

        let v = tags.remove("ID")?;
        let (rg, lane) = v.rsplit_once(':')?;

        match u32::from_str(lane) {
            Ok(n) => Some((v.to_string(), (rg.to_string(), n))),
            Err(_) => {
                let cap = RG_REGEX.captures(lane)?;
                let lane_u32 = u32::from_str(cap.get(1).unwrap().as_str())
                    .map_err(|e| anyhow!("Failed to parse lane number: {}", e))
                    .ok()?;
                Some((v.to_string(), (rg.to_string(), lane_u32)))
            }
        }
    }

    fn try_get_rg(&self, rec: &Record) -> Option<Rg> {
        let rg = rec.aux(b"RG");
        match rg {
            Ok(Aux::String(s)) => {
                let key = match String::from_utf8(Vec::from(s)) {
                    Ok(k) => k,
                    Err(_) => return None, // Return None on UTF-8 error
                };
                self.rg_spec.get(&key).cloned()
            }
            Ok(..) => {
                eprintln!(
                    "invalid type of RG header. record: {}",
                    match str::from_utf8(rec.qname()) {
                        Ok(s) => s,
                        Err(_) => "invalid UTF-8",
                    }
                );
                None
            },
            Err(_) => None,
        }
    }

    pub fn find_rg(&self, rec: &Record) -> Option<Rg> {
        let main_rg_tag = self.try_get_rg(rec);

        if main_rg_tag.is_some() {
            main_rg_tag
        } else {
            let emit = |tag| {
                let corrected_bc = match String::from_utf8(Vec::from(tag)) {
                    Ok(s) => s,
                    Err(_) => return None, // Return None on UTF-8 error
                };
                let mut parts = corrected_bc.split('-');
                let _ = parts.next();
                match parts.next() {
                    Some(v) => {
                        match u32::from_str(v) {
                            Ok(_v) => {
                                //println!("got gg: {}", v);
                                let name = format!("bam2fastq_output");
                                self.rg_spec.get(&name).cloned()
                            }
                            _ => None,
                        }
                    }
                    _ => None,
                }
            };

            if let Ok(Aux::String(s)) = rec.aux(b"CB") {
                return emit(s);
            }

            if let Ok(Aux::String(s)) = rec.aux(b"BX") {
                return emit(s);
            }

            None
        }
    }

    /// Convert a BAM record to a Fq record, for internal caching
    pub fn bam_rec_to_ser(&self, rec: &Record) -> Result<SerFq, Error> {
        Ok(
            match (rec.is_first_in_template(), rec.is_last_in_template()) {
                (true, false) => SerFq {
                    header_key: rec.qname().to_vec(),
                    read_group: self.find_rg(rec),
                    read_num: ReadNum::R1,
                    rec: self
                        .bam_rec_to_fq(rec, &self.r1_spec, self.order[0])
                        .unwrap(),
                    i1: if !self.i1_spec.is_empty() {
                        Some(self.bam_rec_to_fq(rec, &self.i1_spec, self.order[2])?)
                    } else {
                        None
                    },
                    i2: if !self.i2_spec.is_empty() {
                        Some(self.bam_rec_to_fq(rec, &self.i2_spec, self.order[3])?)
                    } else {
                        None
                    },
                },
                (false, true) => SerFq {
                    header_key: rec.qname().to_vec(),
                    read_group: self.find_rg(rec),
                    read_num: ReadNum::R2,
                    rec: self
                        .bam_rec_to_fq(rec, &self.r2_spec, self.order[1])
                        .unwrap(),
                    i1: if !self.i1_spec.is_empty() {
                        Some(self.bam_rec_to_fq(rec, &self.i1_spec, self.order[2])?)
                    } else {
                        None
                    },
                    i2: if !self.i2_spec.is_empty() {
                        Some(self.bam_rec_to_fq(rec, &self.i2_spec, self.order[3])?)
                    } else {
                        None
                    },
                },
                _ => {
                    let e = anyhow!(
                        "Not a valid read pair: {}, {}",
                        rec.is_first_in_template(),
                        rec.is_last_in_template()
                    );
                    return Err(e);
                }
            },
        )
    }

    fn fetch_tag(rec: &Record, tag: &str, last_tag: bool, dest: &mut Vec<u8>) -> Result<(), Error> {
        match rec.aux(tag.as_bytes()) {
            Ok(Aux::String(s)) => dest.extend_from_slice(s.as_bytes()),
            // old BAM files have single-char strings as Char
            Ok(Aux::Char(c)) => dest.push(c),
            Err(_) => {
                if last_tag {
                    return Ok(());
                }
                let e = anyhow!(
                    "BAM record missing tag: {:?} on read {:?}. You do not appear to have an original C4 BAM file.",
                    tag,
                    str::from_utf8(rec.qname()).map_err(|e| anyhow!("Invalid read name UTF-8: {}", e))?
                );
                return Err(e);
            }
            Ok(tag_val) => {
                let e = anyhow!("Invalid BAM record: read: {:?} unexpected tag type. Expected string for {:?}, got {:?}.\n ", str::from_utf8(rec.qname()).unwrap(), tag, tag_val);
                return Err(e);
            }
        }

        Ok(())
    }

    /// Convert a BAM record to Fq record ready to be written
    pub fn bam_rec_to_fq(
        &self,
        rec: &Record,
        spec: &[SpecEntry],
        read_number: u32,
    ) -> Result<FqRecord, Error> {
        let mut head = Vec::new();
        let qname = rec.qname();
                // Find the position of the slash (if it exists)
        let base_name = if let Some(pos) = qname.iter().position(|&x| x == b'/') {
            &qname[..pos]
        } else {
            qname
        };
        // Construct the new header
        head.extend_from_slice(base_name);
        let head_suffix = format!("/{}", read_number);
        head.extend(head_suffix.as_bytes());

        // Reconstitute read and QVs
        let mut read = Vec::new();
        let mut qv = Vec::new();

        for (idx, item) in spec.iter().enumerate() {
            // It OK for the final tag in the spec to be missing from the read
            let last_item = idx == spec.len() - 1;

            match *item {
                // Data from a tag
                SpecEntry::Tags(ref read_tag, ref qv_tag) => {
                    Self::fetch_tag(rec, read_tag, last_item, &mut read)?;
                    Self::fetch_tag(rec, qv_tag, last_item, &mut qv)?;
                }

                // Just hardcode some Ns -- for cases where we didn't retain the required data
                SpecEntry::Ns(len) => {
                    for _ in 0..len {
                        read.push(b'N');
                        qv.push(b'J');
                    }
                }

                SpecEntry::Read => {
                    // The underlying read
                    let mut seq = rec.seq().as_bytes();
                    let mut qual: Vec<u8> = rec.qual().iter().map(|x| x + 33).collect();

                    if rec.is_reverse() {
                        seq.reverse();
                        for b in seq.iter_mut() {
                            *b = complement(*b);
                        }

                        qual.reverse();
                    }

                    read.extend(seq);
                    qv.extend(qual);
                }
            }
        }

        let fq_rec = FqRecord {
            head,
            seq: read,
            qual: qv,
        };

        Ok(fq_rec)
    }

    pub fn format_read_pair(
        &self,
        r1_rec: &Record,
        r2_rec: &Record,
    ) -> Result<FormattedReadPair, Error> {
        let r1 = self.bam_rec_to_fq(r1_rec, &self.r1_spec, self.order[0])?;
        let r2 = self.bam_rec_to_fq(r2_rec, &self.r2_spec, self.order[1])?;

        let i1 = if !self.i1_spec.is_empty() {
            Some(self.bam_rec_to_fq(r1_rec, &self.i1_spec, self.order[2])?)
        } else {
            None
        };

        let i2 = if !self.i2_spec.is_empty() {
            Some(self.bam_rec_to_fq(r1_rec, &self.i2_spec, self.order[3])?)
        } else {
            None
        };

        let rg = self.find_rg(r1_rec);
        Ok((rg, r1, r2, i1, i2))
    }

    pub fn format_read(&self, rec: &Record) -> Result<FormattedReadPair, Error> {
        let r1 = self.bam_rec_to_fq(rec, &self.r1_spec, self.order[0])?;
        let r2 = self.bam_rec_to_fq(rec, &self.r2_spec, self.order[1])?;

        let i1 = if !self.i1_spec.is_empty() {
            Some(self.bam_rec_to_fq(rec, &self.i1_spec, self.order[2])?)
        } else {
            None
        };

        let i2 = if !self.i2_spec.is_empty() {
            Some(self.bam_rec_to_fq(rec, &self.i2_spec, self.order[3])?)
        } else {
            None
        };

        let rg = self.find_rg(rec);
        Ok((rg, r1, r2, i1, i2))
    }

    /// A spec implies double-ended reads if both the R1 and R2 reads generate different BAM records.
    /// If not the R1 and R2 sequences can be derived from a single BAM entry.
    pub fn is_double_ended(&self) -> bool {
        self.r1_spec.contains(&SpecEntry::Read) && self.r2_spec.contains(&SpecEntry::Read)
    }
}



enum GenWriter {
    Gz(ThreadProxyWriter<BufWriter<GzEncoder<File>>>),
    Raw(ThreadProxyWriter<BufWriter<File>>),
}

impl Write for GenWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            GenWriter::Gz(w) => w.write(buf),
            GenWriter::Raw(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            GenWriter::Gz(w) => w.flush(),
            GenWriter::Raw(w) => w.flush(),
        }
    }
}

type Bgw = GenWriter;

struct FastqManager {
    writers: HashMap<Rg, FastqWriter>,
    out_dir: PathBuf,
    // Add batch writing buffer
    write_buffer: Vec<(Option<Rg>, FqRecord, FqRecord, Option<FqRecord>, Option<FqRecord>)>,
    buffer_size: usize,
}

impl FastqManager {
    pub fn new (
        out_dir: &Path,
        formatter: FormatBamRecords,
        _sample_name: String,
        reads_per_fastq: Option<usize>,
        no_compress: bool,
    ) -> Self {
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
                no_compress,
            );

            writers.insert((_sample.clone(), lane), writer);
        }

        // Determine optimal buffer size based on available memory
        let available_memory = get_available_memory();
        let buffer_size = if available_memory < 2 * 1024 * 1024 * 1024 { // Less than 2GB
            100 // Smaller buffer for low memory systems
        } else if available_memory < 8 * 1024 * 1024 * 1024 { // Less than 8GB
            500 // Medium buffer
        } else {
            1000 // Large buffer for high memory systems
        };

        FastqManager {
            writers,
            out_dir: out_dir.to_path_buf(),
            write_buffer: Vec::with_capacity(buffer_size),
            buffer_size,
        }
    }

    pub fn write(
        &mut self,
        rg: &Option<Rg>,
        r1: &FqRecord,
        r2: &FqRecord,
        i1: &Option<FqRecord>,
        i2: &Option<FqRecord>,
    ) {
        // Add to buffer instead of writing immediately
        self.write_buffer.push((
            rg.clone(),
            r1.clone(),
            r2.clone(),
            i1.clone(),
            i2.clone(),
        ));

        // Flush buffer when it reaches capacity
        if self.write_buffer.len() >= self.buffer_size {
            self.flush_buffer();
        }
    }

    pub fn flush_buffer(&mut self) {
        // Process all buffered writes in batch
        for (rg, r1, r2, i1, i2) in self.write_buffer.drain(..) {
            match rg {
                Some(ref rg) => {
                    if let Some(w) = self.writers.get_mut(&rg) {
                        w.write(&r1, &r2, &i1, &i2).expect("Failed to write records");
                    }
                },
                None => {
                    // If there is no RG, use the first available writer
                    if let Some(w) = self.writers.values_mut().next() {
                        w.write(&r1, &r2, &i1, &i2).expect("Failed to write records");
                    }
                }
            }
        }
    }

    pub fn total_written(&self) -> usize {
        self.writers.iter().map(|(_, w)| w.total_written).sum()
    }

    pub fn paths(&self) -> Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)> {
        self.writers
            .iter()
            .flat_map(|(_, w)| w.path_sets.clone())
            .collect()
    }

    pub fn flush_all_writers(&mut self) {
        // First flush the buffer
        self.flush_buffer();
        
        // Then flush all writers
        for writer in self.writers.values_mut() {
            writer.close_current_writers();
        }
    }
}



struct FastqWriter {
    formatter: FormatBamRecords,
    out_dir: PathBuf,
    sample_name: String,
    lane: u32,
    
    r1: Option<Bgw>,
    r2: Option<Bgw>,
    i1: Option<Bgw>,
    i2: Option<Bgw>,

    chunk_written: usize,
    total_written: usize,
    n_chunks: usize,
    reads_per_fastq: Option<usize>,
    path_sets: Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>, 
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
        let extension = if no_compress { "fastq" } else { "fastq.gz" };
        if formatter.rename.is_none() {
            let r1 = out_dir.join(format!(
                "{}_L{:02}_{:01}_1.{}", 
                sample_name, 
                lane,
                n_files + 1,
                extension
            ));
            let r2 = out_dir.join(format!(
                "{}_L{:02}_{:01}_2.{}", 
                sample_name, 
                lane,
                n_files + 1,
                extension
            ));
            let i1 = out_dir.join(format!(
                "{}_L{:02}_{:01}_I1.{}", 
                sample_name, 
                lane,
                n_files + 1,
                extension
            ));
            let i2 = out_dir.join(format!(
                "{}_L{:02}_{:01}_I2.{}", 
                sample_name, 
                lane,
                n_files + 1,
                extension
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
            let new_read_names = formatter.rename.as_ref().unwrap();

            let r1 = out_dir.join(format!(
                "{}_L{:02}_{:01}_{}.{}", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[0],
                extension
            ));
            let r2 = out_dir.join(format!(
                "{}_L{:02}_{:01}_{}.{}", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[1],
                extension
            ));
            let i1 = out_dir.join(format!(
                "{}_L{:02}_{:01}_{}.{}", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[2],
                extension
            ));
            let i2 = out_dir.join(format!(
                "{}_L{:02}_{:01}_{}.{}", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[3],
                extension
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
        }
    }

    pub fn write_rec(
        w: &mut Bgw,
        rec: &FqRecord,
    ) -> Result<(), Error> {
        w.write_all(b"@")?;
        w.write_all(&rec.head)?;
        w.write_all(b"\n")?;
        w.write_all(&rec.seq)?;
        w.write_all(b"\n")?;
        w.write_all(b"+\n")?;
        w.write_all(&rec.qual)?;
        w.write_all(b"\n")?;
        Ok(())
    }

    pub fn try_write_rec(
        w: &mut Option<Bgw>,
        rec: &Option<FqRecord>,
    ) -> Result<(), Error> {
        if let Some(ref mut w) = w {
            if let Some(r) = rec {
                FastqWriter::write_rec(w, r)?;
            } else {
                return Err(anyhow!("No record to write"));
            }
        }
        Ok(())
    }

    pub fn try_write_recs(
        w: &mut Option<Bgw>,
        recs: &FqRecord
    ) -> Result<(), Error> {
        if let Some(ref mut w) = w {
            FastqWriter::write_rec(w, recs)?;
        }
        Ok(())
    }

    fn write(
        &mut self,
        r1: &FqRecord,
        r2: &FqRecord,
        i1: & Option<FqRecord>,
        i2: & Option<FqRecord>,
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

                // Only split if reads_per_fastq is specified
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
        self.i1 = paths.2.as_ref().map(|p| Self::open_writer(p, self.no_compress));
        self.i2 = paths.3.as_ref().map(|p| Self::open_writer(p, self.no_compress));

        self.n_chunks += 1;
        self.chunk_written = 0;
        self.path_sets.push(paths);
    }

    pub fn close_current_writers(&mut self) {
        // Explicitly close and flush current writers
        if let Some(mut writer) = self.r1.take() {
            let _ = writer.flush();
        }
        if let Some(mut writer) = self.r2.take() {
            let _ = writer.flush();
        }
        if let Some(mut writer) = self.i1.take() {
            let _ = writer.flush();
        }
        if let Some(mut writer) = self.i2.take() {
            let _ = writer.flush();
        }
    }

    fn open_writer<P: AsRef<Path>>(
        path: P,
        no_compress: bool,
    ) -> GenWriter {
        let file = File::create(path)
            .map_err(|e| anyhow!("Failed to create output file: {}", e))
            .expect("Failed to create output file");
        // Use adaptive compression settings based on available memory
        let available_memory = get_available_memory();
        
        // Adjust buffer size based on available memory
        let buffer_size = if available_memory < 2 * 1024 * 1024 * 1024 { // Less than 2GB
            1 << 20 // 1MB buffer
        } else {
            1 << 24 // 16MB buffer (original size)
        };

        if no_compress {
            GenWriter::Raw(ThreadProxyWriter::new(
                BufWriter::with_capacity(buffer_size, file), 
                buffer_size / 4
            ))
        } else {
            let compression_level = if available_memory < 4 * 1024 * 1024 * 1024 { // Less than 4GB
                flate2::Compression::fast() // Faster, less memory
            } else {
                flate2::Compression::new(6) // Balanced compression
            };
            let gz = GzEncoder::new(file, compression_level);
            GenWriter::Gz(ThreadProxyWriter::new(
                BufWriter::with_capacity(buffer_size, gz), 
                buffer_size / 4
            ))
        }
    }
}

#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "BAM to FASTQ Converter for C4 Single Cell RNA seq Data", long_about = None)]
pub struct Args {
    /// Input BAM file path
    #[arg(value_name = "BAM", help = "Path to the input BAM file")]
    bam: String,

    /// Output directory for FASTQ files
    #[arg(value_name = "OUTPUT", help = "Directory where FASTQ files will be written")] 
    outputpath: String,

    /// Number of CPU threads to use
    #[arg(
        short = 't',
        long,
        value_name = "THREADS",
        default_value = "4",
        help = "Number of CPU threads for parallel processing"
    )]
    threads: usize,

    /// Process specific genomic region
    #[arg(
        short = 'r',
        long,
        value_name = "REGION",
        help = "Process reads from a specific genomic region (format: chr1:1000-2000)"
    )]
    locus: Option<String>,

    /// BX tag list file (hidden option)
    #[arg(long, hide = true)]
    bx_list: Option<String>,

    /// Number of reads per FASTQ file
    #[arg(
        short = 'n',
        long,
        value_name = "READS",
        help = "Maximum number of reads per FASTQ file. All reads go to a single file if not specified."
    )]
    reads_per_fastq: Option<usize>,

    /// Maximum memory to use in MB (default: automatically determined)
    #[arg(
        long,
        value_name = "MEMORY",
        help = "Maximum memory to use in MB. If not specified, will be automatically determined based on system resources."
    )]
    max_memory: Option<usize>,

    /// Show detailed error traceback
    #[arg(long, hide = true)]
    traceback: bool,

    /// Relaxed mode for unpaired reads
    #[arg(
        long,
        hide = true,
        default_value_t = true,
        help = "Skip unpaired reads instead of throwing an error"
    )]
    relaxed: bool,

    /// Automatically detect if BAM contains paired-end reads
    #[arg(
        long,
        hide = true,
        default_value_t = true,
        help = "Automatically detect if BAM contains paired-end reads based on BAM flags"
    )]
    auto_detect: bool,

    /// Disable output file compression
    #[arg(long, help = "Disable gzip compression for output FASTQ files")]
    no_compress: bool,
}

fn set_panic_handler() {
    panic::set_hook(Box::new(move |info| {
        let backtrace = backtrace::Backtrace::capture();

        let msg = match info.payload().downcast_ref::<&'static str>() {
            Some(s) => *s,
            None => match info.payload().downcast_ref::<String>() {
                Some(s) => &**s,
                None => "Box<Any>",
            }
        };

        let msg = match info.location() {
            Some(location) => format!(
                "bam2fastq failed unexpectedly: '{}' at {}:{}\nBacktrace:\n{:?}",
                msg,
                location.file(),
                location.line(),
                backtrace
            ),
            None => format!(
                "bam2fastq failed unexpectedly: '{}'\nBacktrace:\n{:?}", 
                msg,
                backtrace
            ),
        };

        println!("{}", msg);
    }))
}

pub fn go(args: Args, cache_size: Option<usize>) -> Result<Vec<OutPaths>, Error> {
    // Calculate optimal cache size based on available memory or user-specified limit
    let cache_size = match cache_size {
        Some(size) => size,
        None => {
            let available_memory = match args.max_memory {
                Some(mb) => mb * 1024 * 1024, // Convert MB to bytes
                None => get_available_memory(),
            };
            // Use 1/8 of available memory for cache, with minimum of 100k and maximum of 2M entries
            let calculated_size = (available_memory / 8) / 1024; // Rough estimate assuming 1KB per entry
            calculated_size.clamp(100_000, 2_000_000)
        }
    };

    let path = std::path::PathBuf::from(args.bam.clone());
    if !path.exists() {
        return Err(anyhow!("BAM file doesn't exist: {:?}", path));
    }

    match args.locus {
        Some(ref locus) => {
            let loc = locus::Locus::from_str(locus)
                .context("Invalid locus argument. Please use format: 'chr1:123-456'")?;
            let mut bam = bam::IndexedReader::from_path(&args.bam)
                .context(
                    "Error opening BAM file. The BAM file must be indexed when using --locus",
                )?;
            let tid = bam
                .header()
                .tid(loc.chrom.as_bytes())
                .ok_or_else(|| anyhow!("Requested chromosome not present: {}", loc.chrom))?;

            bam.fetch((tid, loc.start, loc.end))?;
            inner(args.clone(), cache_size, bam)
        }
        None => {
            let _bam = bam::Reader::from_path(&args.bam);
            let bam = _bam.context("Error opening BAM file")?;
            inner(args, cache_size, bam)
        }
    }
}

pub fn inner<R: bam::Read>(
    args: Args,
    cache_size: usize,
    mut bam: R,
) -> Result<Vec<OutPaths>, Error> {
    bam.set_threads(args.threads)?;
    
    let formatter = if args.auto_detect {
        // For auto-detection, we need to create a separate reader to sample the BAM
        // Since we can't easily clone the reader, we'll open the file again for detection
        let mut detection_bam = bam::Reader::from_path(&args.bam)
            .context("Failed to open BAM file for auto-detection")?;
        
        println!("Auto-detecting paired-end status from BAM file...");
        let detected_formatter = FormatBamRecords::c4head_auto(&mut detection_bam)?;
        let is_paired = detected_formatter.is_double_ended();
        println!("Detected: {} data", if is_paired { "paired-end" } else { "single-end" });
        
        detected_formatter
    } else {
        FormatBamRecords::c4head(&bam)
    };

    let out_path = Path::new(&args.outputpath);
    if !out_path.exists() {
        create_dir(&args.outputpath).context(anyhow!(
            "error creating output dir"
        ))?;
    }

    let fq = FastqManager::new(
        out_path, 
        formatter.clone(), 
        "bam2fastq".to_string(), 
        args.reads_per_fastq,
        args.no_compress,
    );

    if formatter.is_double_ended() {
        if args.bx_list.is_some() {
            let bxi = bx_index::BxIndex::new(args.bam)?;
            let bx_iter = BxListIter::from_path(
                args.bx_list.unwrap(), 
                bxi, 
                bam
            )?;
            proc_double_ended(bx_iter, formatter, fq, cache_size, false, args.relaxed, args.threads)
        } else {
            proc_double_ended(
                bam.records(),
                formatter,
                fq,
                cache_size,
                args.locus.is_some(),
                args.relaxed,
                args.threads,
            )
        }
    } else if args.bx_list.is_some() {
        let bxi = bx_index::BxIndex::new(args.bam)?;
        let bx_iter = BxListIter::from_path(args.bx_list.unwrap(), bxi, bam)?;
        proc_double_ended(bx_iter, formatter, fq, cache_size, false, args.relaxed, args.threads)
    } else {
        proc_single_ended(bam.records(), formatter, fq)
    }
}


fn proc_double_ended<I, E> (
    records: I,
    formatter: FormatBamRecords,
    fq: FastqManager,
    cache_size: usize,
    retricted_locus: bool,
    relaxed: bool,
    num_threads: usize,
) -> Result<Vec<OutPaths>, Error>
where 
    I: Iterator<Item = Result<Record, E>>,
    Result<Record, E>: Context<Record, E>,
{
    // Set up thread pool for parallel processing
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .map_err(|e| anyhow!("Failed to initialize thread pool: {}", e))?;

    // Create progress bar
    let progress_bar = ProgressBar::new_spinner();
    progress_bar.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {pos} reads processed ({per_sec}/s) {msg}")
        .map_err(|e| anyhow!("Failed to set progress bar style: {}", e))?
        .progress_chars("#>-"));
    progress_bar.enable_steady_tick(std::time::Duration::from_millis(100));
    
    let temp_file = tempfile::NamedTempFile::new_in(&fq.out_dir)?;
    
    // Adjust ShardWriter parameters based on available memory
    let available_memory = get_available_memory();
    let shard_count = if available_memory < 2 * 1024 * 1024 * 1024 { // Less than 2GB
        16 // Fewer shards to reduce memory usage
    } else {
        32 // Original shard count
    };
    
    let buffer_size = if available_memory < 4 * 1024 * 1024 * 1024 { // Less than 4GB
        1024 // Smaller buffer
    } else {
        2048 // Original buffer size
    };
    
    // Phase 1: Parallel BAM record processing with async orphan handling
    let (total_read_pairs, mut fq) = {
        let mut rp_cache = RpCache::new(cache_size, relaxed);
        let w: ShardWriter<SerFq, SerFqSort> = ShardWriter::new(
            temp_file.path(), 
            shard_count, 
            buffer_size, 
            1 << 20 // Reduced from 1<<21 to reduce memory pressure
        )?;
        let sender = w.get_sender();
        let mut total_read_pairs = 0;
        let mut processed_reads = 0;

        // Create channels for async orphan processing
        let (orphan_tx, orphan_rx): (Sender<Record>, Receiver<Record>) = mpsc::channel();
        let formatter_clone = formatter.clone();
        let mut orphan_sender = sender.clone();
        
        // Spawn background thread for orphan processing
        let orphan_handle = thread::spawn(move || -> Result<(), Error> {
            while let Ok(orphan) = orphan_rx.recv() {
                let ser = formatter_clone.bam_rec_to_ser(&orphan)
                    .map_err(|e| anyhow!("Failed to serialize orphaned read: {}", e))?;
                orphan_sender.send(ser)?;
            }
            Ok(())
        });

        // Create shared FastqManager for parallel writing
        let fq_shared = Arc::new(Mutex::new(fq));
        let formatter_shared = Arc::new(formatter.clone());

        // Process records sequentially but with async orphan handling
        for _rec in records {
            let rec = _rec.context("Error when reading BAM")?;
            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            processed_reads += 1;
            
            match (rec.is_first_in_template(), rec.is_last_in_template()) {
                (false, false) => {
                    return Err(anyhow!("Not single-end read {}", 
                        str::from_utf8(rec.qname()).map_err(|e| anyhow!("Invalid read name UTF-8: {}", e))?))
                }
                (true, true) => {
                    return Err(anyhow!("Read has both r1 and r2 flags: {}", 
                        str::from_utf8(rec.qname()).map_err(|e| anyhow!("Invalid read name UTF-8: {}", e))?))
                }
                (true, false) => total_read_pairs += 1, // R1 read
                (false, true) => (), // R2 read
            }

            let tid = rec.tid();
            let pos = rec.pos();

            if let Some((r1, r2)) = rp_cache.cache_rec(rec) {
                let (rg, fq1, fq2, fq_i1, fq_i2) = formatter_shared
                    .format_read_pair(&r1, &r2)
                    .map_err(|e| anyhow!("Failed to format read pair: {}", e))?;
                
                // Write immediately without blocking
                {
                    let mut fq_guard = fq_shared.lock().unwrap();
                    fq_guard.write(&rg, &fq1, &fq2, &fq_i1, &fq_i2);
                }
            }

            // Async orphan handling - send to background thread
            if rp_cache.len() > cache_size * 3 / 4 {
                for orphan in rp_cache.clear_orphans(tid, pos) {
                    orphan_tx.send(orphan).map_err(|e| anyhow!("Failed to send orphan: {}", e))?;
                }
            }
            
            // Update progress bar every 1000 records
            if processed_reads % 1000 == 0 {
                progress_bar.set_position(processed_reads);
                progress_bar.set_message(format!("Processing BAM records... (Parallel)"));
                
                // Report memory usage periodically
                if processed_reads % 100000 == 0 {
                    let current_memory = get_available_memory();
                    progress_bar.set_message(format!("Processing BAM records... (Memory: {} MB free, {} threads)", 
                        current_memory / (1024 * 1024), num_threads));
                }
            }
        }



        progress_bar.set_message("Processing remaining orphaned reads...");
        for (_, orphan) in rp_cache.cache.drain() {
            orphan_tx.send(orphan).map_err(|e| anyhow!("Failed to send orphan: {}", e))?;
        }

        // Close orphan channel and wait for background thread
        drop(orphan_tx);
        orphan_handle.join().map_err(|e| anyhow!("Orphan processing thread failed: {:?}", e))??;

        progress_bar.finish_with_message(format!("Processed {} reads, found {} read pairs (Parallel)", processed_reads, total_read_pairs));
        
        // Extract FastqManager from Arc<Mutex<>>
        let fq = Arc::try_unwrap(fq_shared).map_err(|_| anyhow!("Failed to unwrap FastqManager"))?.into_inner().unwrap();
        
        (total_read_pairs, fq)
    };

    // Phase 2: Parallel orphan read processing
    let write_progress = ProgressBar::new_spinner();
    write_progress.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.blue} [{elapsed_precise}] Writing orphan FASTQ files... {msg}")
        .map_err(|e| anyhow!("Failed to set progress bar style: {}", e))?
        .progress_chars("#>-"));
    write_progress.enable_steady_tick(std::time::Duration::from_millis(100));
    
    let reader = ShardReader::<SerFq, SerFqSort>::open(temp_file.path())?;
    let mut ncached = 0;
    
    // Parallel processing of orphan reads with immediate writing
    let reader_iter = reader.iter()?;
    let chunk_groups = reader_iter.chunk_by(|x| x.as_ref().ok().map(|x| x.header_key.clone()));
    
    // Use the extracted FastqManager
    
    // Process chunks in parallel batches
    let mut chunk_batch = Vec::new();
    const CHUNK_BATCH_SIZE: usize = 100;
    
    for (_, items) in &chunk_groups {
        let item_vec: Result<Vec<SerFq>, _> = items.collect();
        let item_vec = item_vec?;
        
        if item_vec.len() != 2 && !retricted_locus {
            let header = str::from_utf8(&item_vec[0].rec.head)
                .map_err(|e| anyhow!("Invalid UTF-8 in read header: {}", e))?;
            if !relaxed {
                let msg = anyhow!("Didn't find both records for a paired end read. Is your BAM file complete?\nRead name of unpaired record: {}", header);
                return Err(msg);
            } else {
                println!("Didn't find both records for a paired end read. Skipping. Read name of unpaired record: {}", header);
                continue;
            }
        }
        
        if item_vec.len() != 2 && retricted_locus {
            continue;
        }

        chunk_batch.push(item_vec);
        
        // Process batch when it reaches target size
        if chunk_batch.len() >= CHUNK_BATCH_SIZE {
            let batch = std::mem::take(&mut chunk_batch);
            
            // Process batch in parallel
            let processed_pairs: Vec<_> = batch.into_par_iter().map(|mut item_vec| {
                // Remove sorting requirement - just ensure R1/R2 assignment is correct
                let (r1, r2) = if item_vec.len() == 2 {
                    if item_vec[0].read_num == ReadNum::R1 {
                        (item_vec.swap_remove(0), item_vec.swap_remove(0))
                    } else {
                        (item_vec.swap_remove(1), item_vec.swap_remove(0))
                    }
                } else {
                    // Single read case - treat as R1
                    let read = item_vec.swap_remove(0);
                    let dummy_r2 = SerFq {
                        read_group: read.read_group.clone(),
                        header_key: read.header_key.clone(),
                        rec: FqRecord {
                            head: read.rec.head.clone(),
                            seq: Vec::new(),
                            qual: Vec::new(),
                        },
                        read_num: ReadNum::R2,
                        i1: None,
                        i2: None,
                    };
                    (read, dummy_r2)
                };
                (r1, r2)
            }).collect();
            
            // Write all processed pairs immediately
            for (r1, r2) in processed_pairs {
                fq.write(&r1.read_group, &r1.rec, &r2.rec, &r1.i1, &r1.i2);
                ncached += 1;
            }
        }
        
        // Update progress bar every 100 records written
        if ncached % 100 == 0 {
            write_progress.set_message(format!("Written {} read pairs (Parallel)", ncached));
        }
        
        // Report memory usage and trigger cleanup if needed
        if ncached % 50000 == 0 {
            let current_memory = get_available_memory();
            write_progress.set_message(format!("Written {} read pairs (Memory: {} MB free, {} threads)", 
                ncached, current_memory / (1024 * 1024), num_threads));
            
            // Aggressive memory management based on available memory
            if current_memory < 200 * 1024 * 1024 {
                // Critical memory pressure - flush all writers immediately
                fq.flush_all_writers();
                write_progress.set_message(format!("Written {} read pairs (Critical memory pressure - flushed all buffers)", ncached));
            } else if current_memory < 500 * 1024 * 1024 {
                // Moderate memory pressure - trigger garbage collection
                std::hint::black_box(());
                write_progress.set_message(format!("Written {} read pairs (Low memory detected, optimizing...)", ncached));
            }
        }
    }
    
    // Process remaining chunks in the last batch
    if !chunk_batch.is_empty() {
        let batch = chunk_batch;
        
        let processed_pairs: Vec<_> = batch.into_par_iter().map(|mut item_vec| {
            let (r1, r2) = if item_vec.len() == 2 {
                if item_vec[0].read_num == ReadNum::R1 {
                    (item_vec.swap_remove(0), item_vec.swap_remove(0))
                } else {
                    (item_vec.swap_remove(1), item_vec.swap_remove(0))
                }
            } else {
                let read = item_vec.swap_remove(0);
                let dummy_r2 = SerFq {
                    read_group: read.read_group.clone(),
                    header_key: read.header_key.clone(),
                    rec: FqRecord {
                        head: read.rec.head.clone(),
                        seq: Vec::new(),
                        qual: Vec::new(),
                    },
                    read_num: ReadNum::R2,
                    i1: None,
                    i2: None,
                };
                (read, dummy_r2)
            };
            (r1, r2)
        }).collect();
        
        for (r1, r2) in processed_pairs {
            fq.write(&r1.read_group, &r1.rec, &r2.rec, &r1.i1, &r1.i2);
            ncached += 1;
        }
    }
    
    write_progress.finish_with_message(format!("Completed! Written {} read pairs (Parallel)", ncached));
    
    // Ensure all buffered writes are flushed
    fq.flush_buffer();
    
    println!(
        "Writing finished. \nObserved {} unique read ids. \nWrote {} read pairs ({} cached) using {} threads",
        total_read_pairs,
        fq.total_written(),
        ncached,
        num_threads
    );
    Ok(fq.paths())
}

fn proc_single_ended<I>(
    records: I,
    formatter: FormatBamRecords,
    mut fq: FastqManager,
) -> Result<Vec<OutPaths>, Error>
where
    I: Iterator<Item = Result<Record, rust_htslib::errors::Error>>,
{
    // Create progress bar
    let progress_bar = ProgressBar::new_spinner();
    progress_bar.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {pos} reads processed ({per_sec}/s) {msg}")
        .map_err(|e| anyhow!("Failed to set progress bar style: {}", e))?
        .progress_chars("#>-"));
    progress_bar.enable_steady_tick(std::time::Duration::from_millis(100));
    
    let mut processed_count = 0;
    let mut written_count = 0;
    
    // Batch processing with controlled memory usage
    let mut batch = Vec::new();
    let batch_size = 5000; // Process in batches of 5000 records
    
    for rec_result in records {
        let rec = rec_result.context("Error when reading BAM")?;
        
        if rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        
        processed_count += 1;
        batch.push(rec);
        
        // Process batch when it's full or check memory periodically
        if batch.len() >= batch_size || processed_count % 10000 == 0 {
            // Process current batch
            for rec in batch.drain(..) {
                let formatted = formatter.format_read(&rec)
                    .map_err(|e| anyhow!("Failed to format single read: {}", e))?;
                let (rg, r1, r2, i1, i2) = formatted;
                fq.write(&rg, &r1, &r2, &i1, &i2);
                written_count += 1;
            }
            
            // Report progress and memory usage
            if processed_count % 10000 == 0 {
                progress_bar.set_position(processed_count as u64);
                let current_memory = get_available_memory();
                progress_bar.set_message(format!("Processing single-end reads... (Memory: {} MB free, {} written)", 
                    current_memory / (1024 * 1024), written_count));
                
                // Force flush if memory is getting low
                if current_memory < 500 * 1024 * 1024 { // Less than 500MB
                    fq.flush_buffer();
                    progress_bar.set_message(format!("Low memory detected - flushed buffer ({} written)", written_count));
                }
            } else if processed_count % 1000 == 0 {
                progress_bar.set_position(processed_count as u64);
                progress_bar.set_message(format!("Processing single-end reads... ({} written)", written_count));
            }
        }
    }
    
    // Process any remaining records in the final batch
    for rec in batch {
        let formatted = formatter.format_read(&rec)
            .map_err(|e| anyhow!("Failed to format single read: {}", e))?;
        let (rg, r1, r2, i1, i2) = formatted;
        fq.write(&rg, &r1, &r2, &i1, &i2);
        written_count += 1;
    }

    progress_bar.finish_with_message(format!("Processed {} reads, written {} reads", processed_count, written_count));
    
    // Ensure all buffered writes are flushed
    fq.flush_buffer();
    
    // make sure we have the right number of output reads
    println!(
        "Writing finished. \nObserved {} reads. \nWrote {} reads",
        processed_count,
        fq.total_written()
    );
    Ok(fq.paths())
}


fn main() {
    set_panic_handler();
    std::env::set_var("RUST_BACKTRACE", "1");

    //println!("bam2fastq v{}", VERSION);
        // Use clap to parse command line arguments
    let args = Args::parse();

    let traceback = args.traceback;
    let res = go(args, None);

    if let Err(ref e) = res {
        println!("bam2fastq error: {e}\n");

        //println!("If this error is unexpected. Please re-run with --traceback and include stack trace with an error report");

        if traceback {
            println!("see below for more details:");
            println!("==========================");
            println!("{}\n{}", e, e.backtrace());
        };
        ::std::process::exit(1);
    }
}


#[cfg(test)]
mod tests {

    #[test]
    fn test_lr21() {
        // Create a fixed output directory
        let output_dir = "target/fastq_results";

        let args = crate::Args {
            threads: 2,
            bam: "/Volumes/mac_up/anno_decon_sorted.bam".to_string(),
            outputpath: output_dir.to_string(),
            reads_per_fastq: None,
            locus: None,
            bx_list: None,
            traceback: false,
            relaxed: false,
            max_memory: None,
            auto_detect: true,
            no_compress: false,
        };

                // Run conversion
        let out_path_sets = super::go(args, None).unwrap();
        
                // Print result file paths
        println!("\nGenerated FASTQ files:");
        for (r1, r2, i1, i2) in out_path_sets {
            println!("R1: {}", r1.display());
            println!("R2: {}", r2.display());
            if let Some(i1_path) = i1 {
                println!("I1: {}", i1_path.display());
            }
            if let Some(i2_path) = i2 {
                println!("I2: {}", i2_path.display());
            }
            println!("---");
        }

                // Optional: Check if files are generated and print file sizes
        println!("\nFile size information:");
        for entry in std::fs::read_dir(output_dir).unwrap() {
            let entry = entry.unwrap();
            let metadata = entry.metadata().unwrap();
            println!("{}: {} bytes", entry.file_name().to_string_lossy(), metadata.len());
        }
    }

    

    #[test]
    fn bad_bam() {
        let tempdir = tempfile::Builder::new()
            .prefix("bam_to_fq_test")
            .tempdir()
            .expect("create temp dir");
        let tmp_path = tempdir.path().join("outs");

        let args = crate::Args  {
            threads: 2,
            bam: "/Users/lishuangshuang/Documents/scrna/dnbc4tools/target/my_test_3/pos_sortednon_multiplexed.bam".to_string(),
            outputpath: tmp_path.to_str().unwrap().to_string(),
            reads_per_fastq: None,
            locus: None,
            bx_list: None,
            traceback: false,
            relaxed: false,
            max_memory: None,
            auto_detect: true,
            no_compress: false,
        };

        let res = super::go(args, Some(2));

        println!("res: {:?}", res);
    }
}