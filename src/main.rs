//! main.rs — Entry point and orchestration for bam2fastq
//!
//! This binary converts C4 single-cell RNA-seq BAM files to FASTQ format.
//! The heavy logic lives in the submodules below:
//!   - `formatter` — BAM record → FqRecord conversion
//!   - `writer`    — FASTQ file management (GenWriter, FastqManager)
//!   - `processing`— end-to-end pipelines (single-ended, double-ended)

mod bx_index;
mod formatter;
mod locus;
mod processing;
mod rpcache;
mod writer;

use anyhow::{anyhow, Context, Error};
use bx_index::BxListIter;
use clap::Parser;
use formatter::{FormatBamRecords, OutPaths};
use processing::{proc_double_ended, proc_single_ended};
use rust_htslib::bam::{self, Read};
use std::{backtrace, panic, path::Path, str::FromStr};
use writer::FastqManager;

// ─── CLI arguments ────────────────────────────────────────────────────────

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "BAM to FASTQ Converter for C4 Single Cell RNA seq Data",
    long_about = None
)]
pub struct Args {
    /// Input BAM file path
    #[arg(value_name = "BAM", help = "Path to the input BAM file")]
    bam: String,

    /// Output directory for FASTQ files
    #[arg(
        value_name = "OUTPUT",
        help = "Directory where FASTQ files will be written"
    )]
    outputpath: String,

    /// Number of CPU threads to use
    #[arg(
        short = 't',
        long,
        value_name = "THREADS",
        default_value_t = num_cpus::get(),
        help = "Number of CPU threads for parallel processing (default: all available cores)"
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

    /// Maximum memory to use in MB (default: auto-detected)
    #[arg(
        long,
        value_name = "MEMORY",
        help = "Maximum memory to use in MB. Auto-determined if not specified."
    )]
    max_memory: Option<usize>,

    /// Show detailed error traceback
    #[arg(long, hide = true)]
    traceback: bool,

    /// Relaxed mode: skip unpaired reads instead of erroring
    #[arg(long, hide = true, default_value_t = true)]
    relaxed: bool,

    /// Automatically detect paired-end status from BAM flags
    #[arg(long, hide = true, default_value_t = true)]
    auto_detect: bool,

    /// Disable gzip compression for output FASTQ files
    #[arg(long, help = "Disable gzip compression for output FASTQ files")]
    no_compress: bool,
}

// ─── Panic handler ────────────────────────────────────────────────────────

fn set_panic_handler() {
    panic::set_hook(Box::new(move |info| {
        let backtrace = backtrace::Backtrace::capture();
        let msg = match info.payload().downcast_ref::<&'static str>() {
            Some(s) => *s,
            None => match info.payload().downcast_ref::<String>() {
                Some(s) => &**s,
                None => "Box<Any>",
            },
        };
        let msg = match info.location() {
            Some(loc) => format!(
                "bam2fastq failed unexpectedly: '{}' at {}:{}\nBacktrace:\n{:?}",
                msg,
                loc.file(),
                loc.line(),
                backtrace
            ),
            None => format!(
                "bam2fastq failed unexpectedly: '{}'\nBacktrace:\n{:?}",
                msg, backtrace
            ),
        };
        println!("{}", msg);
    }))
}

// ─── Public entry points ──────────────────────────────────────────────────

/// Top-level entry point used by tests and the CLI. Opens the BAM file,
/// applies any locus filter, and dispatches to [`inner`].
pub fn go(args: Args, cache_size: Option<usize>) -> Result<Vec<OutPaths>, Error> {
    let cache_size = match cache_size {
        Some(size) => size,
        None => {
            let available_memory = match args.max_memory {
                Some(mb) => mb * 1024 * 1024,
                None => writer::get_memory(),
            };
            let calculated_size = (available_memory / 8) / 1024;
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
            let mut bam = bam::IndexedReader::from_path(&args.bam).context(
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
            let bam = bam::Reader::from_path(&args.bam).context("Error opening BAM file")?;
            inner(args, cache_size, bam)
        }
    }
}

/// Configure thread pools, open the formatter, create the `FastqManager`,
/// then dispatch to the appropriate processing pipeline.
pub fn inner<R: bam::Read>(
    args: Args,
    cache_size: usize,
    mut bam: R,
) -> Result<Vec<OutPaths>, Error> {
    // Thread budget allocation:
    //   Rayon par_iter (format): args.threads      — CPU-intensive, gets the full budget
    //   BAM BGZF decompression:  max(1, threads/4) — IO-bound, diminishing returns
    //   gzp compression:         max(1, threads/4) — IO-bound, diminishing returns
    let io_threads = (args.threads / 4).max(1);
    bam.set_threads(io_threads)?;

    // Initialize Rayon global thread pool to the full --threads count.
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap_or_else(|e| eprintln!("Warning: failed to configure Rayon thread pool: {}", e));

    let formatter = if args.auto_detect {
        let mut detection_bam = bam::Reader::from_path(&args.bam)
            .context("Failed to open BAM file for auto-detection")?;
        let _ = detection_bam.set_threads(io_threads);
        println!("Auto-detecting paired-end status from BAM file...");
        let detected = FormatBamRecords::c4head_auto(&mut detection_bam)?;
        println!(
            "Detected: {} data",
            if detected.is_double_ended() {
                "paired-end"
            } else {
                "single-end"
            }
        );
        detected
    } else {
        FormatBamRecords::c4head(&bam)
    };

    let out_path = Path::new(&args.outputpath);
    if !out_path.exists() {
        std::fs::create_dir(&args.outputpath).context("error creating output dir")?;
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
            let bx_iter = BxListIter::from_path(args.bx_list.unwrap(), bxi, bam)?;
            proc_double_ended(
                bx_iter,
                formatter,
                fq,
                cache_size,
                false,
                args.relaxed,
                args.threads,
            )
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
        proc_double_ended(
            bx_iter,
            formatter,
            fq,
            cache_size,
            false,
            args.relaxed,
            args.threads,
        )
    } else {
        proc_single_ended(bam.records(), formatter, fq, args.threads)
    }
}

// ─── Binary entry point ───────────────────────────────────────────────────

fn main() {
    set_panic_handler();
    std::env::set_var("RUST_BACKTRACE", "1");

    let args = Args::parse();
    let traceback = args.traceback;
    let res = go(args, None);

    if let Err(ref e) = res {
        println!("bam2fastq error: {e}\n");
        if traceback {
            println!("see below for more details:");
            println!("==========================");
            println!("{}\n{}", e, e.backtrace());
        }
        ::std::process::exit(1);
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    #[test]
    fn test_lr21() {
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
            no_compress: true,
        };
        let out_path_sets = super::go(args, None).unwrap();
        println!("\nGenerated FASTQ files:");
        for (r1, r2, i1, i2) in out_path_sets {
            println!("R1: {}", r1.display());
            println!("R2: {}", r2.display());
            if let Some(p) = i1 {
                println!("I1: {}", p.display());
            }
            if let Some(p) = i2 {
                println!("I2: {}", p.display());
            }
            println!("---");
        }
    }

    #[test]
    fn bad_bam() {
        let tempdir = tempfile::Builder::new()
            .prefix("bam_to_fq_test")
            .tempdir()
            .expect("create temp dir");
        let tmp_path = tempdir.path().join("outs");
        let args = crate::Args {
            threads: 2,
            bam: "/my_test_3/pos_sortednon_multiplexed.bam".to_string(),
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
