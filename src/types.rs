//! Type definitions for BAM to FASTQ conversion

use std::path::PathBuf;
use serde::{Serialize, Deserialize};

/// Output paths tuple: (R1, R2, I1, I2)
pub type OutPaths = (
    PathBuf,
    PathBuf,
    Option<PathBuf>,
    Option<PathBuf>,
);

/// Read group identifier: (sample_name, lane)
pub type Rg = (String, u32);

/// Formatted read pair: (read_group, R1, R2, I1, I2)
pub type FormattedReadPair = (
    Option<Rg>,
    FqRecord,
    FqRecord,
    Option<FqRecord>,
    Option<FqRecord>,
);

/// FASTQ record structure
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct FqRecord {
    #[serde(with = "serde_bytes")]
    pub head: Vec<u8>,
    #[serde(with = "serde_bytes")]
    pub seq: Vec<u8>,
    #[serde(with = "serde_bytes")]
    pub qual: Vec<u8>,
}

/// Read number enumeration
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum ReadNum {
    R1,
    R2,
}

/// Serialized FASTQ record for caching
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

/// Specification entry for read formatting
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpecEntry {
    Tags(String, String),
    #[allow(dead_code)]
    Ns(usize),
    Read,
}