use std::{borrow::Cow, collections::HashMap, path::PathBuf};
use serde::{Serialize, Deserialize};
use shardio::SortKey;
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam::{self, Read};



mod bx_index;


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

    pub fn find_rg(&self, _rec: &Record) -> Option<Rg> {
        None
    }
}































fn main() {
    println!("Hello, world!");
}
