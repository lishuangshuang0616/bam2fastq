use std::backtrace;
use std::io::{BufWriter, Write};
use std::{borrow::Cow, collections::HashMap, path::PathBuf, str::FromStr};
use rpcache::RpCache;
use serde::{Serialize, Deserialize};
use shardio::helper::ThreadProxyWriter;
use shardio::SortKey;
use shardio::{ShardReader, ShardWriter};
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam::{self, header, Read};
use regex::Regex;
use itertools::Itertools;
use anyhow::{Error, anyhow, Context};
use flate2::write::GzEncoder;
use std::fs::{create_dir, File};
use std::panic;
use std::path::Path;
use std::str;
use docopt::Docopt;

mod bx_index;
mod locus;
mod rpcache;

use bx_index::BxListIter;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const USAGE: &str = "Usage: cargo run --release <input>";

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


#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)]
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

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum SpecEntry {
    Tags(String, String),
    Ns(usize),
    Read
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
    pub fn from_headers<R: bam::Read>(reader: &R) -> Option<Self> {
        let mut spec = Self::parse_spec(reader);
        let seq_names = Self::parse_seq_names(reader);

        if spec.is_empty() {
            None
        } else {
            Some(Self {
                rg_spec: HashMap::new(),
                r1_spec: spec.remove("R1").unwrap(),
                r2_spec: spec.remove("R2").unwrap(),
                i1_spec: spec.remove("I1").unwrap_or_default(),
                i2_spec: spec.remove("I2").unwrap_or_default(),
                rename: seq_names,
                order: [1, 3, 2, 4],
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
        let text = std::str::from_utf8(reader.header().as_bytes()).unwrap();

        let mut rg_items = text
            .lines()
            .filter(|l| l.starts_with("@RG"))
            .filter_map(Self::parse_rg_line)
            .collect::<HashMap<_, _>>();

        if rg_items.is_empty() {
            println!("WARNING: no @RG (read group) headers found in BAM file. Splitting data by the GEM well marked in the corrected barcode tag.");
            println!("Reads without a corrected barcode will not appear in output FASTQs");
            // No RG items in header -- invent a set fixed set of RGs
            // each observed Gem group in the BAM file will get mapped to these.
            for i in 1..100 {
                let name = format!("gemgroup{:03}", i);
                rg_items.insert(name.clone(), (name, 0));
            }
        }

        rg_items
    }

    fn parse_rg_line(line: &str) -> Option<(String, (String, u32))> {
        let mut entries = line.split('\t');
        entries.next()?; // consume @RG entry

        let mut tags = entries
            .map(|entry| entry.split_once(':').unwrap())
            .collect::<HashMap<_, _>>();

        let v = tags.remove("ID")?;
        let (rg, lane) = v.rsplit_once(':')?;

        match u32::from_str(lane) {
            Ok(n) => Some((v.to_string(), (rg.to_string(), n))),
            Err(_) => {
                // Handle case in ALIGNER pipeline prior to 2.1.3 -- samtools merge would append a unique identifier to each RG ID tags
                // Detect this condition and remove from lane
                let re = Regex::new(r"^([0-9]+)-[0-9A-F]+$").unwrap();
                let cap = re.captures(lane)?;
                let lane_u32 = u32::from_str(cap.get(1).unwrap().as_str()).unwrap();
                Some((v.to_string(), (rg.to_string(), lane_u32)))
            }
        }
    }

    /// Parse the specs from BAM headers if available
    fn parse_spec<R: bam::Read>(reader: &R) -> HashMap<String, Vec<SpecEntry>> {
        // Example header line:
        // @CO	10x_bam_to_fastq:R1(RX:QX,TR:TQ,SEQ:QUAL)
        let re = Regex::new(r"@CO\t10x_bam_to_fastq:(\S+)\((\S+)\)").unwrap();
        let text = String::from_utf8(Vec::from(reader.header().as_bytes())).unwrap();

        text.lines()
            .into_iter()
            .filter_map(|l| {
                re.captures(l).map(|c| {
                    let read = c.get(1).unwrap().as_str().to_string();
                    let tag_list = c.get(2).unwrap().as_str();

                    let spec_entries = tag_list
                        .split(',')
                        .into_iter()
                        .map(|el| {
                            if el == "SEQ:QUAL" {
                                SpecEntry::Read
                            } else {
                                let (rtag, qtag) =
                                    el.split(':').map(ToString::to_string).next_tuple().unwrap();
                                SpecEntry::Tags(rtag, qtag)
                            }
                        })
                        .collect();

                    (read, spec_entries)
                })
            })
            .collect()
    }

    // Example header line:
    // @CO	10x_bam_to_fastq_seqnames:R1,R3,I1,R2
    // In this case, the @CO header lines marked R1, R2, I1, I2 will
    // be used to write reads to output files R1, R3, I1, and R2, respectively
    fn parse_seq_names<R: bam::Read>(reader: &R) -> Option<Vec<String>> {
        let text = String::from_utf8(Vec::from(reader.header().as_bytes())).unwrap();
        let re = Regex::new(r"@CO\t10x_bam_to_fastq_seqnames:(\S+)").unwrap();

        for l in text.lines() {
            if let Some(c) = re.captures(l) {
                let names = c.get(1).unwrap().as_str().split(',');
                let seq_names = names
                    .into_iter()
                    .map(std::string::ToString::to_string)
                    .collect();
                return Some(seq_names);
            }
        }
        None
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
                    "BAM record missing tag: {:?} on read {:?}. You do not appear to have an original 10x BAM file.\nIf you downloaded this BAM file from SRA, you likely need to download the 'Original Format' version of the BAM available for most 10x datasets.",
                    tag,
                    str::from_utf8(rec.qname()).unwrap()
                );
                return Err(e);
            }
            Ok(tag_val) => {
                let e = anyhow!("Invalid BAM record: read: {:?} unexpected tag type. Expected string for {:?}, got {:?}.\n You do not appear to have the original 10x BAM file. If you downloaded this BAM file from SRA, you likely need to download the 'Original Format' version of the BAM available for most 10x datasets.", str::from_utf8(rec.qname()).unwrap(), tag, tag_val);
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
        head.extend_from_slice(rec.qname());
        let head_suffix = format!(" {}:N:0:0", read_number);
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



type Bgw = ThreadProxyWriter<BufWriter<GzEncoder<File>>>;

struct FastqManager {
    writers: HashMap<Rg, FastqWriter>,
    out_dir: PathBuf,
}

impl FastqManager {
    pub fn new (
        out_dir: &Path,
        formatter: FormatBamRecords,
        _sample_name: String,
        reads_per_fastq: usize,
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
                lane,
                reads_per_fastq,
            );

            writers.insert((_sample.clone(), lane), writer);
        }

        FastqManager {
            writers,
            out_dir: out_dir.to_path_buf(),
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
        if let &Some(ref rg) = rg {
            self.writers.get_mut(rg).map(|w| w.write(r1, r2, i1, i2));
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
    reads_per_fastq: usize,
    path_sets: Vec<(PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>)>,
}

impl FastqWriter {
    pub fn new(
        out_dir: &Path,
        formatter: FormatBamRecords,
        sample_name: String,
        lane: u32,
        reads_per_fastq: usize,
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
        }
    }

    fn get_paths(
        out_dir: &Path,
        sample_name: &str,
        lane: u32,
        n_files: usize,
        formatter: &FormatBamRecords,
    ) -> (PathBuf, PathBuf, Option<PathBuf>, Option<PathBuf>) {
        if formatter.rename.is_none() {
            let r1 = out_dir.join(format!(
                "{}_L{:02}_{:02}_1.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1
            ));
            let r2 = out_dir.join(format!(
                "{}_L{:02}_{:02}_2.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1
            ));
            let i1 = out_dir.join(format!(
                "{}_L{:02}_{:02}_I1.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1
            ));
            let i2 = out_dir.join(format!(
                "{}_L{:02}_{:02}_I2.fastq.gz", 
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
            let new_read_names = formatter.rename.as_ref().unwrap();

            let r1 = out_dir.join(format!(
                "{}_L{:02}_{:02}_{}.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[0]
            ));
            let r2 = out_dir.join(format!(
                "{}_L{:02}_{:02}_{}.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[1]
            ));
            let i1 = out_dir.join(format!(
                "{}_L{:02}_{:02}_{}.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[2]
            ));
            let i2 = out_dir.join(format!(
                "{}_L{:02}_{:02}_{}.fastq.gz", 
                sample_name, 
                lane,
                n_files + 1,
                new_read_names[3]
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
                panic!("No record to write");
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

        if self.chunk_written == self.reads_per_fastq {
            self.cycle_writers();
        }
        Ok(())
    }

    fn cycle_writers(&mut self) {
        let paths = Self::get_paths(
            &self.out_dir, 
            &self.sample_name, 
            self.lane, 
            self.n_chunks, 
            &self.formatter
        );

        self.r1 = Some(Self::open_gzip_writer(paths.0));
        self.r2 = Some(Self::open_gzip_writer(paths.1));
        self.i1 = paths.2.as_ref().map(Self::open_gzip_writer);
        self.i2 = paths.3.as_ref().map(Self::open_gzip_writer);

        self.n_chunks += 1;
    }

    fn open_gzip_writer<P: AsRef<Path>>(
        path: P
    ) -> ThreadProxyWriter<BufWriter<GzEncoder<File>>> {
        let file = File::create(path).unwrap();
        let gz = GzEncoder::new(file, flate2::Compression::default());
        ThreadProxyWriter::new(
            BufWriter::with_capacity(1 << 22, gz), 1 << 19
        )
    }
}

#[derive(Debug, Deserialize, Clone)]
pub struct Args {
    arg_bam: String,
    arg_output_path: String,
    flag_nthreads: usize,
    flag_locus: Option<String>,
    flag_bx_list: Option<String>,
    flag_reads_per_fastq: usize,
    flag_traceback: bool,
    flag_relaxed: bool,
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
    let cache_size = cache_size.unwrap_or(500_000);

    let path = std::path::PathBuf::from(args.arg_bam.clone());
    if !path.exists() {
        return Err(anyhow!("BAM file doesn't exist: {:?}", path));
    }

    match args.flag_locus {
        Some(ref locus) => {
            let loc = locus::Locus::from_str(locus)
                .context("Invalid locus argument. Please use format: 'chr1:123-456'")?;
            let mut bam = bam::IndexedReader::from_path(&args.arg_bam)
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
            let _bam = bam::Reader::from_path(&args.arg_bam);
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
    bam.set_threads( args.flag_nthreads)?;
    let formatter = {
        FormatBamRecords::c4head(&bam)
    };

    let out_path = Path::new(&args.arg_output_path);
    create_dir(&args.arg_output_path).context(anyhow!(
        "error creating output dir"
    ))?;

    let fq = FastqManager::new(
        out_path, 
        formatter.clone(), 
        "bam2fastq".to_string(), 
        args.flag_reads_per_fastq
    );

    if formatter.is_double_ended() {
        if args.flag_bx_list.is_some() {
            let bxi = bx_index::BxIndex::new(args.arg_bam)?;
            let bx_iter = BxListIter::from_path(
                args.flag_bx_list.unwrap(), 
                bxi, 
                bam
            )?;
            proc_double_ended(bx_iter, formatter, fq, cache_size, false, args.flag_relaxed)
        } else {
            proc_double_ended(
                bam.records(),
                formatter,
                fq,
                cache_size,
                args.flag_locus.is_some(),
                args.flag_relaxed,
            )
        }
    } else if args.flag_bx_list.is_some() {
        let bxi = bx_index::BxIndex::new(args.arg_bam)?;
        let bx_iter = BxListIter::from_path(args.flag_bx_list.unwrap(), bxi, bam)?;
        proc_double_ended(bx_iter, formatter, fq, cache_size, false, args.flag_relaxed)
    } else {
        proc_single_ended(bam.records(), formatter, fq)
    }
}


fn proc_double_ended<I, E> (
    records: I,
    formatter: FormatBamRecords,
    mut fq: FastqManager,
    cache_size: usize,
    retricted_locus: bool,
    relaxed: bool,
) -> Result<Vec<OutPaths>, Error>
where 
    I: Iterator<Item = Result<Record, E>>,
    Result<Record, E>: Context<Record, E>,
{
    let temp_file = tempfile::NamedTempFile::new_in(&fq.out_dir)?;
    let total_read_pairs = {
        let mut rp_cache = RpCache::new(cache_size, relaxed);
        let w: ShardWriter<SerFq, SerFqSort> = ShardWriter::new(
            temp_file.path(), 
            32, 
            2048, 
            1 << 21
        )?;
        let mut sender = w.get_sender();
        let mut totble_read_pairs = 0;

        for _rec in records {
            let rec = _rec.context("Error when reading BAM")?;
            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            match (rec.is_first_in_template(), rec.is_last_in_template()) {
                (false, false) => {
                    return Err(anyhow!("Not single-end read {}", str::from_utf8(rec.qname()).unwrap()))
                }
                (true, true) => {
                    return Err(anyhow!("Read has both r1 and r2 flags: {}", str::from_utf8(rec.qname()).unwrap()))
                }
                (true, false) => totble_read_pairs += 1,
                (false, true) => (),
            }

            let tid = rec.tid();
            let pos = rec.pos();

            if let Some((r1, r2)) = rp_cache.cache_rec(rec) {
                let (rg, fq1, fq2, fq_i1, fq_i2) = formatter
                    .format_read_pair(&r1, &r2)
                    .unwrap();
                fq.write(&rg, &fq1, &fq2, &fq_i1, &fq_i2)
            }

            if rp_cache.len() > cache_size {
                for orphan in rp_cache.clear_orphans(tid, pos) {
                    let ser = formatter.bam_rec_to_ser(&orphan)?;
                    sender.send(ser)?;
                }
            }
        }

        for (_, orphan) in rp_cache.cache.drain() {
            let ser = formatter.bam_rec_to_ser(&orphan)?;
            sender.send(ser)?;
        }

        totble_read_pairs
    };

    let reader = ShardReader::<SerFq, SerFqSort>::open(temp_file.path())?;
    let mut ncached = 0;
    for (_, items) in &reader
        .iter()?
        .chunk_by(|x| x.as_ref().ok().map(|x| x.header_key.clone()))
    {
        let _item_vec: Result<Vec<SerFq>, _> = items.collect();
        let mut item_vec = _item_vec?;
        if item_vec.len() != 2 && !retricted_locus {
            let header = str::from_utf8(&item_vec[0].rec.head).unwrap();
            if !relaxed {
                let msg = anyhow!("Didn't find both records for a paired end read. Is your BAM file complete?\nRead name of unpaired record: {}", header);
                return Err(msg);
            } else {
                println!("Didn't find both records for a paired end read. Skipping. Read name of unpaired record: {}", header);
            }
        }
        if item_vec.len() != 2 && retricted_locus{
            continue;
        }

        item_vec.sort_by_key(|x| x.read_num);
        let r1 = item_vec.swap_remove(0);
        let r2 = item_vec.swap_remove(0);
        fq.write(&r1.read_group, &r1.rec, &r2.rec, &r1.i1, &r1.i2);
        ncached += 1;
    }
    println!(
        "Writing finished.  Observed {} unique read ids. Wrote {} read pairs ({} cached)",
        total_read_pairs,
        fq.total_written(),
        ncached
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
    let total_reads = {
        // Count total R1s observed, so we can make sure we've preserved all read pairs
        let mut total_reads = 0;

        for _rec in records {
            let rec = _rec.context("Error when reading BAM")?;

            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            total_reads += 1;

            let (rg, fq1, fq2, fq_i1, fq_i2) = formatter.format_read(&rec)?;
            fq.write(&rg, &fq1, &fq2, &fq_i1, &fq_i2);
        }

        total_reads
    };

    // make sure we have the right number of output reads
    println!(
        "Writing finished.  Observed {} read pairs. Wrote {} read pairs",
        total_reads,
        fq.total_written()
    );
    Ok(fq.paths())
}


fn main() {
    set_panic_handler();
    std::env::set_var("RUST_BACKTRACE", "1");

    println!("bam2fastq v{}", VERSION);
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    let traceback = args.flag_traceback;
    let res = go(args, None);

    if let Err(ref e) = res {
        println!("bam2fastq error: {e}\n");

        println!("If this error is unexpected. Please re-run with --traceback and include stack trace with an error report");

        if traceback {
            println!("see below for more details:");
            println!("==========================");
            println!("{}\n{}", e, e.backtrace());
        };
        ::std::process::exit(1);
    }
}