//! BAM record formatting functionality

use std::{borrow::Cow, collections::HashMap, str::FromStr};
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam;
use regex::Regex;
use anyhow::{Error, anyhow};
use shardio::SortKey;

use crate::types::*;
use crate::utils::complement;
use crate::memory_pool::ByteBufferPool;

/// Sort key implementation for SerFq
pub struct SerFqSort;

impl SortKey<SerFq> for SerFqSort {
    type Key = Vec<u8>;
    fn sort_key(t: &SerFq) -> Cow<Vec<u8>> {
        Cow::Borrowed(&t.header_key)
    }
}

/// BAM record formatter for C4 single-cell data
#[derive(Debug, Clone)]
pub struct FormatBamRecords {
    pub rg_spec: HashMap<String, Rg>,
    pub r1_spec: Vec<SpecEntry>,
    pub r2_spec: Vec<SpecEntry>,
    pub i1_spec: Vec<SpecEntry>,
    pub i2_spec: Vec<SpecEntry>,
    pub rename: Option<Vec<String>>,
    pub order: [u32; 4],
    _buffer_pool: ByteBufferPool,
}

impl FormatBamRecords {
    /// Create formatter for C4 single-cell data
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
            _buffer_pool: ByteBufferPool::new(50, 1024),
        }
    }

    /// Parse read groups from BAM header
    fn parse_rgs<R: bam::Read>(reader: &R) -> HashMap<String, Rg> {
        let text = std::str::from_utf8(reader.header().as_bytes()).unwrap();

        let mut rg_items = text
            .lines()
            .filter(|l| l.starts_with("@RG"))
            .filter_map(Self::parse_rg_line)
            .collect::<HashMap<_, _>>();

        if rg_items.is_empty() {
            for _i in 1..2 {
                let name = "bam2fastq_output".to_string();
                rg_items.insert(name.clone(), (name, 0));
            }
        }

        rg_items
    }

    /// Parse a single read group line
    fn parse_rg_line(line: &str) -> Option<(String, (String, u32))> {
        let mut entries = line.split('\t');
        entries.next()?; // consume @RG entry

        let mut tags = HashMap::new();
        for entry in entries {
            if let Some((key, value)) = entry.split_once(':') {
                tags.insert(key, value);
            } else {
                return None; // Invalid entry format
            }
        }

        let v = tags.remove("ID")?;
        let (rg, lane) = v.rsplit_once(':')?;

        match u32::from_str(lane) {
            Ok(n) => Some((v.to_string(), (rg.to_string(), n))),
            Err(_) => {
                let re = Regex::new(r"^([0-9]+)-[0-9A-F]+$").ok()?;
                let cap = re.captures(lane)?;
                let lane_str = cap.get(1)?.as_str();
                let lane_u32 = u32::from_str(lane_str).ok()?;
                Some((v.to_string(), (rg.to_string(), lane_u32)))
            }
        }
    }

    /// Try to get read group from RG tag
    fn try_get_rg(&self, rec: &Record) -> Option<Rg> {
        let rg = rec.aux(b"RG");
        match rg {
            Ok(Aux::String(s)) => {
                let key = String::from_utf8(Vec::from(s)).unwrap();
                self.rg_spec.get(&key).cloned()
            }
            Ok(..) => panic!(
                "invalid type of RG header. record: {}",
                std::str::from_utf8(rec.qname()).unwrap()
            ),
            Err(_) => None,
        }
    }

    /// Find read group for a record
    pub fn find_rg(&self, rec: &Record) -> Option<Rg> {
        let main_rg_tag = self.try_get_rg(rec);

        if main_rg_tag.is_some() {
            main_rg_tag
        } else {
            let emit = |tag| {
                let corrected_bc = String::from_utf8(Vec::from(tag)).unwrap();
                let mut parts = corrected_bc.split('-');
                let _ = parts.next();
                match parts.next() {
                    Some(v) => {
                        match u32::from_str(v) {
                            Ok(_v) => {
                                let name = "bam2fastq_output".to_string();
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

    /// Convert a BAM record to a SerFq record for internal caching
    pub fn bam_rec_to_ser(&self, rec: &Record) -> Result<SerFq, Error> {
        Ok(
            match (rec.is_first_in_template(), rec.is_last_in_template()) {
                (true, false) => SerFq {
                    header_key: rec.qname().to_vec(),
                    read_group: self.find_rg(rec),
                    read_num: ReadNum::R1,
                    rec: self
                        .bam_rec_to_fq(rec, &self.r1_spec, self.order[0])?,
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
                        .bam_rec_to_fq(rec, &self.r2_spec, self.order[1])?,
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

    /// Fetch tag from BAM record
    fn fetch_tag(&self, rec: &Record, tag: &str, last_tag: bool, dest: &mut Vec<u8>) -> Result<(), Error> {
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
                    std::str::from_utf8(rec.qname()).unwrap()
                );
                return Err(e);
            }
            Ok(tag_val) => {
                let e = anyhow!("Invalid BAM record: read: {:?} unexpected tag type. Expected string for {:?}, got {:?}.\n ", std::str::from_utf8(rec.qname()).unwrap(), tag, tag_val);
                return Err(e);
            }
        }

        Ok(())
    }

    /// Convert a BAM record to FqRecord ready to be written
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
        // Build new header
        head.extend_from_slice(base_name);
        let head_suffix = format!("/{}", read_number);
        head.extend(head_suffix.as_bytes());

        // Reconstitute read and QVs
        let mut read = self._buffer_pool.get();
        let mut qv = self._buffer_pool.get();

        for (idx, item) in spec.iter().enumerate() {
            // It OK for the final tag in the spec to be missing from the read
            let last_item = idx == spec.len() - 1;

            match *item {
                // Data from a tag
                SpecEntry::Tags(ref read_tag, ref qv_tag) => {
                    self.fetch_tag(rec, read_tag, last_item, &mut read)?;
                    self.fetch_tag(rec, qv_tag, last_item, &mut qv)?;
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
                            *b = complement(*b)?;
                        }

                        qual.reverse();
                    }

                    read.extend(seq);
                    qv.extend(qual);
                }
            }
        }

        let result = FqRecord {
            head,
            seq: read.clone(),
            qual: qv.clone(),
        };
        
        // Return buffers to pool
        self._buffer_pool.put(read);
        self._buffer_pool.put(qv);
        
        Ok(result)
    }

    /// Format a read pair
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

    /// Format a single read
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

    /// Check if spec implies double-ended reads
    pub fn is_double_ended(&self) -> bool {
        self.r1_spec.contains(&SpecEntry::Read) && self.r2_spec.contains(&SpecEntry::Read)
    }
}