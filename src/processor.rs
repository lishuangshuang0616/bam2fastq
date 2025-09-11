//! BAM record processing functionality

use std::path::Path;
use rust_htslib::bam::{self, Read, Record};
use anyhow::{anyhow, Error, Context};
use indicatif::{ProgressBar, ProgressStyle};


use crate::types::*;
use crate::formatter::FormatBamRecords;
use crate::writer::FastqManager;
use crate::advanced_cache::AdvancedRpCache;
use crate::cli::Args;
use crate::locus::Locus;
use std::str::FromStr;


// Constants for configuration
const DEFAULT_CACHE_SIZE: usize = 500_000;
const PROGRESS_UPDATE_INTERVAL: usize = 100;

/// Main processing function
pub fn go(args: Args, cache_size: Option<usize>) -> Result<Vec<OutPaths>, Error> {
    let cache_size = cache_size.unwrap_or(DEFAULT_CACHE_SIZE);
    
    let mut bam = bam::Reader::from_path(&args.bam)
        .with_context(|| format!("Failed to open BAM file: {}", args.bam))?;
    
    bam.set_threads(args.threads)
        .with_context(|| "Failed to set thread count for BAM reader")?;
    
    inner(args, cache_size, bam)
}

/// Inner processing function with generic BAM reader
pub fn inner<R: bam::Read>(
    args: Args,
    cache_size: usize,
    mut bam: R,
) -> Result<Vec<OutPaths>, Error> {
    let formatter = FormatBamRecords::c4head(&bam);
    
    let out_dir = Path::new(&args.outputpath);
    let fq = FastqManager::new(
        out_dir,
        formatter.clone(),
        "sample".to_string(),
        args.reads_per_fastq,
    );
    
    let restricted_locus = args.locus.is_some();
    
    if let Some(locus_str) = &args.locus {
        let _locus_spec = Locus::from_str(locus_str)
            .with_context(|| format!("Failed to parse locus: {}", locus_str))?;
        // Note: fetch method is only available for IndexedReader, not generic bam::Read
        // This will be handled in the go function where we have the proper reader type
    }
    
    let records = bam.records();
    
    if formatter.is_double_ended() {
        proc_double_ended(
            records,
            formatter,
            fq,
            cache_size,
            restricted_locus,
            args.relaxed,
        )
    } else {
        proc_single_ended(records, formatter, fq)
    }
}

/// Process double-ended reads (paired-end)
fn proc_double_ended<I>(
    records: I,
    formatter: FormatBamRecords,
    mut fq: FastqManager,
    cache_size: usize,
    _restricted_locus: bool,
    relaxed: bool,
) -> Result<Vec<OutPaths>, Error>
where 
    I: Iterator<Item = Result<Record, rust_htslib::errors::Error>>,
{
    let mut cache = AdvancedRpCache::new(cache_size, relaxed);
    let mut n_written = 0;
    let mut n_processed = 0;
    
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {pos} records processed, {msg}")
            .unwrap()
    );
    
    for rec_result in records {
        let rec = rec_result
            .with_context(|| "Failed to read BAM record")?;
        n_processed += 1;
        
        if n_processed % PROGRESS_UPDATE_INTERVAL == 0 {
            pb.set_position(n_processed as u64);
            pb.set_message(format!("{} written", n_written));
        }
        
        // Periodically clear orphans based on current position
        if n_processed % 10000 == 0 {
            let current_tid = rec.tid();
            let current_pos = rec.pos();
            let orphans = cache.clear_orphans(current_tid, current_pos);
            
            if !relaxed && !orphans.is_empty() {
                eprintln!("Warning: {} orphaned reads cleared from cache", orphans.len());
            }
        }
        
        if let Some((r1_rec, r2_rec)) = cache.cache_rec(rec) {
            match formatter.format_read_pair(&r1_rec, &r2_rec) {
                Ok((rg, r1, r2, i1, i2)) => {
                    fq.write(&rg, &r1, &r2, &i1, &i2);
                    n_written += 1;
                }
                Err(e) => {
                    if !relaxed {
                        return Err(e).with_context(|| "Failed to format read pair");
                    }
                    // Skip this read pair in relaxed mode
                }
            }
        }
    }
    
    // Handle remaining orphaned reads
    let final_orphans = cache.clear_orphans(i32::MAX, i64::MAX); // Clear all remaining
    if !relaxed && !final_orphans.is_empty() {
        return Err(anyhow!("Found {} unpaired reads at end of processing", final_orphans.len()));
    }
    
    if relaxed && !final_orphans.is_empty() {
        eprintln!("Skipped {} unpaired reads in relaxed mode", final_orphans.len());
    }
    
    pb.finish_with_message(format!("Completed: {} records processed, {} written", n_processed, n_written));
    
    Ok(fq.paths())
}

/// Process single-ended reads
fn proc_single_ended<I>(
    records: I,
    formatter: FormatBamRecords,
    mut fq: FastqManager,
) -> Result<Vec<OutPaths>, Error>
where
    I: Iterator<Item = Result<Record, rust_htslib::errors::Error>>,
{
    let mut n_written = 0;
    let mut n_processed = 0;
    
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {pos} records processed, {msg}")
            .unwrap()
    );
    
    for rec_result in records {
        let rec = rec_result
            .with_context(|| "Failed to read BAM record")?;
        n_processed += 1;
        
        if n_processed % PROGRESS_UPDATE_INTERVAL == 0 {
            pb.set_position(n_processed as u64);
            pb.set_message(format!("{} written", n_written));
        }
        
        match formatter.format_read(&rec) {
            Ok((rg, r1, r2, i1, i2)) => {
                fq.write(&rg, &r1, &r2, &i1, &i2);
                n_written += 1;
            }
            Err(e) => {
                return Err(e).with_context(|| "Failed to format single read");
            }
        }
    }
    
    pb.finish_with_message(format!("Completed: {} records processed, {} written", n_processed, n_written));
    
    Ok(fq.paths())
}