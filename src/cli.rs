//! Command line interface definitions

use clap::Parser;

/// Command line arguments for BAM to FASTQ converter
#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "BAM to FASTQ Converter for C4 Single Cell RNA seq Data", long_about = None)]
pub struct Args {
    /// Path to the input BAM file
    #[arg(value_name = "BAM", help = "Path to the input BAM file")]
    pub bam: String,

    /// Directory where FASTQ files will be written
    #[arg(value_name = "OUTPUT", help = "Directory where FASTQ files will be written")] 
    pub outputpath: String,

    /// Number of CPU threads for parallel processing
    #[arg(
        short = 't',
        long,
        value_name = "THREADS",
        default_value = "4",
        help = "Number of CPU threads for parallel processing"
    )]
    pub threads: usize,

    /// Process reads from a specific genomic region
    #[arg(
        short = 'r',
        long,
        value_name = "REGION",
        help = "Process reads from a specific genomic region (format: chr1:1000-2000)"
    )]
    pub locus: Option<String>,

    /// Barcode list file (hidden option)
    #[arg(long, hide = true)]
    pub bx_list: Option<String>,

    /// Maximum number of reads per FASTQ file
    #[arg(
        short = 'n',
        long,
        value_name = "READS",
        help = "Maximum number of reads per FASTQ file. All reads go to a single file if not specified."
    )]
    pub reads_per_fastq: Option<usize>,

    /// Enable traceback on errors (hidden option)
    #[arg(long, hide = true)]
    pub traceback: bool,

    /// Skip unpaired reads instead of throwing an error
    #[arg(
        long,
        hide = true,
        default_value_t = true,
        help = "Skip unpaired reads instead of throwing an error"
    )]
    pub relaxed: bool,
}