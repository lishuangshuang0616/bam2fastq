//! formatter.rs — BAM → FASTQ record conversion
//!
//! Contains all core data types and the `FormatBamRecords` struct that handles
//! reading BAM tags and sequences and producing `FqRecord` / `SerFq` outputs.

use anyhow::{anyhow, Error};
use once_cell::sync::Lazy;
use regex::Regex;
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam::{self};
use serde::{Deserialize, Serialize};
use shardio::SortKey;
use std::{borrow::Cow, collections::HashMap, str, str::FromStr};

// ─── Public type aliases ───────────────────────────────────────────────────

pub type Rg = (String, u32);

pub type OutPaths = (
    std::path::PathBuf,
    std::path::PathBuf,
    Option<std::path::PathBuf>,
    Option<std::path::PathBuf>,
);

pub type FormattedReadPair = (
    Option<Rg>,
    FqRecord,
    FqRecord,
    Option<FqRecord>,
    Option<FqRecord>,
);

// ─── Core data structures ──────────────────────────────────────────────────

/// A single FASTQ record (header, sequence, quality).
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub struct FqRecord {
    #[serde(with = "serde_bytes")]
    pub head: Vec<u8>,
    #[serde(with = "serde_bytes")]
    pub seq: Vec<u8>,
    #[serde(with = "serde_bytes")]
    pub qual: Vec<u8>,
}

/// Whether a read is R1 or R2 in a pair.
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum ReadNum {
    R1,
    R2,
}

/// Serialised FASTQ read, stored in the shardio shard for paired-end matching.
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SerFq {
    pub read_group: Option<Rg>,
    #[serde(with = "serde_bytes")]
    pub header_key: Vec<u8>,
    pub rec: FqRecord,
    pub read_num: ReadNum,
    pub i1: Option<FqRecord>,
    pub i2: Option<FqRecord>,
}

/// Sort key for shardio: sort by QNAME byte string.
pub struct SerFqSort;

impl SortKey<SerFq> for SerFqSort {
    type Key = Vec<u8>;
    fn sort_key(t: &SerFq) -> Cow<'_, Vec<u8>> {
        Cow::Borrowed(&t.header_key)
    }
}

/// One entry in a read spec: either BAM tags, N padding, or the raw read sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpecEntry {
    Tags(String, String),
    #[allow(dead_code)]
    Ns(usize),
    Read,
}

/// Holds all configuration needed to convert BAM records to FASTQ records.
#[derive(Debug, Clone)]
pub struct FormatBamRecords {
    pub rg_spec: HashMap<String, Rg>,
    pub r1_spec: Vec<SpecEntry>,
    pub r2_spec: Vec<SpecEntry>,
    pub i1_spec: Vec<SpecEntry>,
    pub i2_spec: Vec<SpecEntry>,
    pub rename: Option<Vec<String>>,
    pub order: [u32; 4],
}

// ─── Complement lookup table ───────────────────────────────────────────────

/// 256-element lookup table for fast nucleotide complement calculation.
static COMPLEMENT_LUT: [u8; 256] = {
    let mut lut = [b'N'; 256];
    lut[b'A' as usize] = b'T';
    lut[b'a' as usize] = b't';
    lut[b'C' as usize] = b'G';
    lut[b'c' as usize] = b'g';
    lut[b'G' as usize] = b'C';
    lut[b'g' as usize] = b'c';
    lut[b'T' as usize] = b'A';
    lut[b't' as usize] = b'a';
    lut[b'N' as usize] = b'N';
    lut[b'n' as usize] = b'n';
    lut
};

#[inline(always)]
pub fn complement(b: u8) -> u8 {
    COMPLEMENT_LUT[b as usize]
}

/// Regex to parse lane numbers from RG IDs like "SAMPLENAME:1-ABCDEF".
pub static RG_REGEX: Lazy<Regex> = Lazy::new(|| Regex::new(r"^([0-9]+)-[0-9A-F]+$").unwrap());

// ─── FormatBamRecords implementation ──────────────────────────────────────

impl FormatBamRecords {
    /// Automatically detect if BAM contains paired-end reads by sampling records.
    pub fn detect_paired_end<R: bam::Read>(reader: &mut R) -> Result<bool, Error> {
        let mut paired_count = 0;
        let mut total_count = 0;
        let max_samples = 100_000;

        for result in reader.records() {
            let record = result.map_err(|e| anyhow!("Failed to read BAM record: {}", e))?;
            total_count += 1;
            if record.is_paired() {
                paired_count += 1;
            }
            if total_count >= max_samples {
                break;
            }
        }

        let paired_ratio = if total_count > 0 {
            paired_count as f64 / total_count as f64
        } else {
            0.0
        };
        Ok(paired_ratio > 0.5)
    }

    /// Create C4 configuration with automatic paired-end detection.
    pub fn c4head_auto<R: bam::Read>(reader: &mut R) -> Result<FormatBamRecords, Error> {
        let is_paired = Self::detect_paired_end(reader)?;
        if is_paired {
            Ok(FormatBamRecords {
                rg_spec: Self::parse_rgs(reader),
                r1_spec: vec![
                    SpecEntry::Tags("CR".to_string(), "CY".to_string()),
                    SpecEntry::Tags("UR".to_string(), "UY".to_string()),
                    SpecEntry::Read,
                ],
                r2_spec: vec![SpecEntry::Read],
                i1_spec: vec![],
                i2_spec: vec![],
                rename: None,
                order: [1, 2, 0, 0],
            })
        } else {
            Ok(FormatBamRecords {
                rg_spec: Self::parse_rgs(reader),
                r1_spec: vec![
                    SpecEntry::Tags("CR".to_string(), "CY".to_string()),
                    SpecEntry::Tags("UR".to_string(), "UY".to_string()),
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
                SpecEntry::Tags("UR".to_string(), "UY".to_string()),
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
            Err(_) => return HashMap::new(),
        };

        let mut rg_items = text
            .lines()
            .filter(|l| l.starts_with("@RG"))
            .filter_map(Self::parse_rg_line)
            .collect::<HashMap<_, _>>();

        if rg_items.is_empty() {
            let name = "bam2fastq_output".to_string();
            rg_items.insert(name.clone(), (name, 0));
        }
        rg_items
    }

    fn parse_rg_line(line: &str) -> Option<(String, (String, u32))> {
        let mut entries = line.split('\t');
        entries.next()?; // consume @RG

        let tags = entries
            .map(|entry| {
                entry
                    .split_once(':')
                    .ok_or_else(|| anyhow!("Invalid RG entry format"))
            })
            .collect::<Result<HashMap<_, _>, _>>()
            .ok()?;

        let v = tags.get("ID")?;
        let (rg, lane) = v.rsplit_once(':')?;

        match u32::from_str(lane) {
            Ok(n) => Some((v.to_string(), (rg.to_string(), n))),
            Err(_) => {
                let cap = RG_REGEX.captures(lane)?;
                let lane_u32 = u32::from_str(cap.get(1).unwrap().as_str()).ok()?;
                Some((v.to_string(), (rg.to_string(), lane_u32)))
            }
        }
    }

    fn try_get_rg(&self, rec: &Record) -> Option<Rg> {
        match rec.aux(b"RG") {
            Ok(Aux::String(s)) => {
                let key = String::from_utf8(Vec::from(s)).ok()?;
                self.rg_spec.get(&key).cloned()
            }
            Ok(..) => {
                eprintln!(
                    "invalid type of RG header. record: {}",
                    str::from_utf8(rec.qname()).unwrap_or("invalid UTF-8")
                );
                None
            }
            Err(_) => None,
        }
    }

    pub fn find_rg(&self, rec: &Record) -> Option<Rg> {
        let main_rg_tag = self.try_get_rg(rec);
        if main_rg_tag.is_some() {
            return main_rg_tag;
        }

        let emit = |tag: &[u8]| {
            let corrected_bc = String::from_utf8(Vec::from(tag)).ok()?;
            let mut parts = corrected_bc.split('-');
            let _ = parts.next();
            match parts.next() {
                Some(v) if u32::from_str(v).is_ok() => {
                    let name = "bam2fastq_output".to_string();
                    self.rg_spec.get(&name).cloned()
                }
                _ => None,
            }
        };

        if let Ok(Aux::String(s)) = rec.aux(b"CB") {
            return emit(s.as_bytes());
        }
        if let Ok(Aux::String(s)) = rec.aux(b"BX") {
            return emit(s.as_bytes());
        }
        None
    }

    /// Convert a BAM record to a `SerFq` (for internal paired-end caching).
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
                    return Err(anyhow!(
                        "Not a valid read pair: is_first={}, is_last={}",
                        rec.is_first_in_template(),
                        rec.is_last_in_template()
                    ));
                }
            },
        )
    }

    fn fetch_tag(rec: &Record, tag: &str, last_tag: bool, dest: &mut Vec<u8>) -> Result<(), Error> {
        match rec.aux(tag.as_bytes()) {
            Ok(Aux::String(s)) => dest.extend_from_slice(s.as_bytes()),
            Ok(Aux::Char(c)) => dest.push(c),
            Err(_) => {
                if last_tag {
                    return Ok(());
                }
                return Err(anyhow!(
                    "BAM record missing tag: {:?} on read {:?}. \
                     You do not appear to have an original C4 BAM file.",
                    tag,
                    str::from_utf8(rec.qname())
                        .map_err(|e| anyhow!("Invalid read name UTF-8: {}", e))?
                ));
            }
            Ok(tag_val) => {
                return Err(anyhow!(
                    "Invalid BAM record: read: {:?} unexpected tag type. \
                     Expected string for {:?}, got {:?}.",
                    str::from_utf8(rec.qname()).unwrap_or("?"),
                    tag,
                    tag_val
                ));
            }
        }
        Ok(())
    }

    /// Convert a BAM record to an `FqRecord` ready to be written.
    pub fn bam_rec_to_fq(
        &self,
        rec: &Record,
        spec: &[SpecEntry],
        read_number: u32,
    ) -> Result<FqRecord, Error> {
        let qname = rec.qname();
        let base_name = if let Some(pos) = qname.iter().position(|&x| x == b'/') {
            &qname[..pos]
        } else {
            qname
        };

        // Pre-allocate header: qname + "/" + single digit
        let mut head = Vec::with_capacity(base_name.len() + 2);
        head.extend_from_slice(base_name);
        head.push(b'/');
        head.push(b'0' + (read_number as u8));

        let mut read = Vec::new();
        let mut qv = Vec::new();

        for (idx, item) in spec.iter().enumerate() {
            let last_item = idx == spec.len() - 1;

            match *item {
                SpecEntry::Tags(ref read_tag, ref qv_tag) => {
                    Self::fetch_tag(rec, read_tag, last_item, &mut read)?;
                    Self::fetch_tag(rec, qv_tag, last_item, &mut qv)?;
                }
                SpecEntry::Ns(len) => {
                    for _ in 0..len {
                        read.push(b'N');
                        qv.push(b'J');
                    }
                }
                SpecEntry::Read => {
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

        Ok(FqRecord {
            head,
            seq: read,
            qual: qv,
        })
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
        Ok((self.find_rg(r1_rec), r1, r2, i1, i2))
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
        Ok((self.find_rg(rec), r1, r2, i1, i2))
    }

    /// Returns `true` if both R1 and R2 require separate BAM records (paired-end).
    pub fn is_double_ended(&self) -> bool {
        self.r1_spec.contains(&SpecEntry::Read) && self.r2_spec.contains(&SpecEntry::Read)
    }
}
