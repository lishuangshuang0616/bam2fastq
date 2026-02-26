//! processing.rs — BAM → FASTQ conversion pipeline
//!
//! Contains `proc_single_ended` (streaming, Rayon-parallel) and
//! `proc_double_ended` (RpCache-based read-pair matching).

use crate::formatter::{FormatBamRecords, FormattedReadPair, OutPaths, ReadNum, SerFq, SerFqSort};
use crate::rpcache::RpCache;
use crate::writer::{FastqManager, AVAILABLE_MEMORY};
use anyhow::{anyhow, Context, Error};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use rayon::prelude::*;
use rust_htslib::bam::record::Record;
use shardio::{ShardReader, ShardWriter};
use std::str;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

// ─── Single-ended processing ───────────────────────────────────────────────

/// Process a single-ended BAM stream.
///
/// Pipeline stages (all concurrent via channel):
/// 1. Main thread: collect BAM records into batches of `BATCH_SIZE`
/// 2. Rayon `par_iter`: format each batch in parallel across all available cores
/// 3. Dedicated write thread: receives formatted batches and writes to disk
pub fn proc_single_ended<I>(
    records: I,
    formatter: FormatBamRecords,
    fq: FastqManager,
    num_threads: usize,
) -> Result<Vec<OutPaths>, Error>
where
    I: Iterator<Item = Result<Record, rust_htslib::errors::Error>>,
{
    let progress_bar = ProgressBar::new_spinner();
    progress_bar.set_style(
        ProgressStyle::default_spinner()
            .template(
                "{spinner:.green} [{elapsed_precise}] {pos} reads processed ({per_sec}/s) {msg}",
            )
            .map_err(|e| anyhow!("Failed to set progress bar style: {}", e))?
            .progress_chars("#>-"),
    );
    progress_bar.enable_steady_tick(std::time::Duration::from_millis(100));

    // Buffer scales with thread count: more threads → faster format → need bigger buffer
    // Formula: max(4, threads / 2) e.g. 4 threads→4, 16 threads→8, 32 threads→16
    let channel_buffer = (num_threads / 2).max(4);
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<FormattedReadPair>>(channel_buffer);

    // Dedicated write thread — owns FastqManager so writes never block the formatter
    let write_thread = std::thread::spawn(move || -> Result<FastqManager, Error> {
        let mut fq = fq;
        while let Ok(batch) = rx.recv() {
            for (rg, r1, r2, i1, i2) in batch {
                fq.write(&rg, &r1, &r2, &i1, &i2);
            }
        }
        fq.flush_buffer();
        Ok(fq)
    });

    const BATCH_SIZE: usize = 5_000;
    let mut batch: Vec<Record> = Vec::with_capacity(BATCH_SIZE);
    let mut processed_count: u64 = 0;

    for rec_result in records {
        let rec = rec_result.context("Error when reading BAM")?;
        if rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        batch.push(rec);

        if batch.len() >= BATCH_SIZE {
            let formatted: Vec<FormattedReadPair> = batch
                .par_iter()
                .filter_map(|rec| {
                    formatter
                        .format_read(rec)
                        .map_err(|e| eprintln!("Format error: {}", e))
                        .ok()
                })
                .collect();

            processed_count += formatted.len() as u64;
            progress_bar.set_position(processed_count);
            tx.send(formatted)
                .map_err(|_| anyhow!("Write thread died unexpectedly"))?;
            batch.clear();
        }
    }

    // Flush last partial batch
    if !batch.is_empty() {
        let formatted: Vec<FormattedReadPair> = batch
            .par_iter()
            .filter_map(|rec| {
                formatter
                    .format_read(rec)
                    .map_err(|e| eprintln!("Format error: {}", e))
                    .ok()
            })
            .collect();

        processed_count += formatted.len() as u64;
        progress_bar.set_position(processed_count);
        tx.send(formatted)
            .map_err(|_| anyhow!("Write thread died unexpectedly"))?;
    }

    drop(tx); // Signal write thread: no more data
    let fq = write_thread
        .join()
        .map_err(|_| anyhow!("Write thread panicked"))?
        .context("Write thread error")?;

    progress_bar.finish_with_message(format!("Processed {} reads", processed_count));
    println!(
        "Writing finished. \nObserved {} reads. \nWrote {} reads",
        processed_count,
        fq.total_written()
    );
    Ok(fq.paths())
}

// ─── Double-ended processing ───────────────────────────────────────────────

/// Process a paired-end BAM stream.
///
/// Phase 1: Read BAM records sequentially, matching R1+R2 pairs via `RpCache`.
///          Matched pairs are written immediately; orphaned reads go to a shardio temp file.
/// Phase 2: Read the sorted shardio shard, re-pair orphans, and write remaining pairs.
///          Writing is pipelined via a dedicated write thread.
pub fn proc_double_ended<I, E>(
    records: I,
    formatter: FormatBamRecords,
    mut fq: FastqManager,
    cache_size: usize,
    retricted_locus: bool,
    relaxed: bool,
    num_threads: usize,
) -> Result<Vec<OutPaths>, Error>
where
    I: Iterator<Item = Result<Record, E>>,
    Result<Record, E>: Context<Record, E>,
{
    // Thread pool for parallel processing
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .map_err(|e| anyhow!("Failed to initialize thread pool: {}", e))?;

    let progress_bar = ProgressBar::new_spinner();
    progress_bar.set_style(
        ProgressStyle::default_spinner()
            .template(
                "{spinner:.green} [{elapsed_precise}] {pos} reads processed ({per_sec}/s) {msg}",
            )
            .map_err(|e| anyhow!("Failed to set progress bar style: {}", e))?
            .progress_chars("#>-"),
    );
    progress_bar.enable_steady_tick(std::time::Duration::from_millis(100));

    let temp_file = tempfile::NamedTempFile::new_in(&fq.out_dir)?;

    let available_memory = *AVAILABLE_MEMORY;
    let shard_count = if available_memory < 2 * 1024 * 1024 * 1024 {
        16
    } else {
        32
    };
    let buffer_size = if available_memory < 4 * 1024 * 1024 * 1024 {
        1024
    } else {
        2048
    };

    // ── Phase 1: sequential BAM scan with async orphan handling ──────────
    let (total_read_pairs, mut fq) = {
        let mut rp_cache = RpCache::new(cache_size, relaxed);
        let w: ShardWriter<SerFq, SerFqSort> =
            ShardWriter::new(temp_file.path(), shard_count, buffer_size, 1 << 20)?;
        let sender = w.get_sender();
        let mut total_read_pairs = 0u64;
        let mut processed_reads = 0u64;

        let (orphan_tx, orphan_rx): (Sender<Record>, Receiver<Record>) = mpsc::channel();
        let formatter_clone = formatter.clone();
        let mut orphan_sender = sender.clone();

        let orphan_handle = thread::spawn(move || -> Result<(), Error> {
            while let Ok(orphan) = orphan_rx.recv() {
                let ser = formatter_clone
                    .bam_rec_to_ser(&orphan)
                    .map_err(|e| anyhow!("Failed to serialize orphaned read: {}", e))?;
                orphan_sender.send(ser)?;
            }
            Ok(())
        });

        let formatter_shared = std::sync::Arc::new(formatter.clone());

        for _rec in records {
            let rec = _rec.context("Error when reading BAM")?;
            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }
            processed_reads += 1;

            match (rec.is_first_in_template(), rec.is_last_in_template()) {
                (false, false) => {
                    return Err(anyhow!(
                        "Not single-end read {}",
                        str::from_utf8(rec.qname())
                            .map_err(|e| anyhow!("Invalid read name UTF-8: {}", e))?
                    ))
                }
                (true, true) => {
                    return Err(anyhow!(
                        "Read has both r1 and r2 flags: {}",
                        str::from_utf8(rec.qname())
                            .map_err(|e| anyhow!("Invalid read name UTF-8: {}", e))?
                    ))
                }
                (true, false) => total_read_pairs += 1,
                (false, true) => (),
            }

            let tid = rec.tid();
            let pos = rec.pos();

            if let Some((r1, r2)) = rp_cache.cache_rec(rec) {
                let (rg, fq1, fq2, fq_i1, fq_i2) = formatter_shared
                    .format_read_pair(&r1, &r2)
                    .map_err(|e| anyhow!("Failed to format read pair: {}", e))?;
                fq.write(&rg, &fq1, &fq2, &fq_i1, &fq_i2);
            }

            if rp_cache.len() > cache_size * 3 / 4 {
                for orphan in rp_cache.clear_orphans(tid, pos) {
                    orphan_tx
                        .send(orphan)
                        .map_err(|e| anyhow!("Failed to send orphan: {}", e))?;
                }
            }

            if processed_reads % 10_000 == 0 {
                progress_bar.set_position(processed_reads);
            }
        }

        progress_bar.set_message("Processing remaining orphaned reads...");
        for (_, orphan) in rp_cache.cache.drain() {
            orphan_tx
                .send(orphan)
                .map_err(|e| anyhow!("Failed to send orphan: {}", e))?;
        }

        drop(orphan_tx);
        orphan_handle
            .join()
            .map_err(|e| anyhow!("Orphan processing thread failed: {:?}", e))??;

        progress_bar.finish_with_message(format!(
            "Processed {} reads, found {} read pairs",
            processed_reads, total_read_pairs
        ));
        (total_read_pairs, fq)
    };

    // ── Phase 2: read sorted orphan shard, re-pair, write ────────────────
    let write_progress = ProgressBar::new_spinner();
    write_progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.blue} [{elapsed_precise}] Writing orphan FASTQ files... {msg}")
            .map_err(|e| anyhow!("Failed to set progress bar style: {}", e))?
            .progress_chars("#>-"),
    );
    write_progress.enable_steady_tick(std::time::Duration::from_millis(100));

    let reader = ShardReader::<SerFq, SerFqSort>::open(temp_file.path())?;
    let mut ncached = 0usize;

    let channel_buffer = (num_threads / 2).max(4);
    type PairBatch = Vec<(SerFq, SerFq)>;
    let (pair_tx, pair_rx) = std::sync::mpsc::sync_channel::<PairBatch>(channel_buffer);

    let write_thread = std::thread::spawn(move || -> FastqManager {
        while let Ok(pairs) = pair_rx.recv() {
            for (r1, r2) in pairs {
                fq.write(&r1.read_group, &r1.rec, &r2.rec, &r1.i1, &r1.i2);
            }
        }
        fq.flush_buffer();
        fq
    });

    let reader_iter = reader.iter()?;
    let chunk_groups = reader_iter.chunk_by(|x| x.as_ref().ok().map(|x| x.header_key.clone()));

    let mut chunk_batch: Vec<Vec<SerFq>> = Vec::new();
    const CHUNK_BATCH_SIZE: usize = 1000;

    for (_, items) in &chunk_groups {
        let item_vec: Result<Vec<SerFq>, _> = items.collect();
        let item_vec = item_vec?;

        if item_vec.len() != 2 && !retricted_locus {
            let header = str::from_utf8(&item_vec[0].rec.head)
                .map_err(|e| anyhow!("Invalid UTF-8 in read header: {}", e))?;
            if !relaxed {
                return Err(anyhow!(
                    "Didn't find both records for a paired end read. \
                     Is your BAM file complete?\nRead name of unpaired record: {}",
                    header
                ));
            } else {
                println!(
                    "Didn't find both records for a paired end read. \
                     Skipping. Read name: {}",
                    header
                );
                continue;
            }
        }
        if item_vec.len() != 2 && retricted_locus {
            continue;
        }

        chunk_batch.push(item_vec);

        if chunk_batch.len() >= CHUNK_BATCH_SIZE {
            let batch = std::mem::take(&mut chunk_batch);
            let processed_pairs: PairBatch = batch
                .into_par_iter()
                .map(|mut item_vec| {
                    if item_vec[0].read_num == ReadNum::R1 {
                        (item_vec.swap_remove(0), item_vec.swap_remove(0))
                    } else {
                        (item_vec.swap_remove(1), item_vec.swap_remove(0))
                    }
                })
                .collect();

            ncached += processed_pairs.len();
            pair_tx.send(processed_pairs).ok();

            if ncached % 10_000 == 0 {
                write_progress.set_message(format!("Written {} read pairs", ncached));
            }
        }
    }

    // Flush remaining
    if !chunk_batch.is_empty() {
        let processed_pairs: PairBatch = chunk_batch
            .into_par_iter()
            .map(|mut item_vec| {
                if item_vec.len() == 2 {
                    if item_vec[0].read_num == ReadNum::R1 {
                        (item_vec.swap_remove(0), item_vec.swap_remove(0))
                    } else {
                        (item_vec.swap_remove(1), item_vec.swap_remove(0))
                    }
                } else {
                    let read = item_vec.swap_remove(0);
                    use crate::formatter::{FqRecord, ReadNum as Rn};
                    let dummy_r2 = SerFq {
                        read_group: read.read_group.clone(),
                        header_key: read.header_key.clone(),
                        rec: FqRecord {
                            head: read.rec.head.clone(),
                            seq: Vec::new(),
                            qual: Vec::new(),
                        },
                        read_num: Rn::R2,
                        i1: None,
                        i2: None,
                    };
                    (read, dummy_r2)
                }
            })
            .collect();

        ncached += processed_pairs.len();
        pair_tx.send(processed_pairs).ok();
    }

    drop(pair_tx);
    let fq = write_thread
        .join()
        .map_err(|_| anyhow!("Write thread panicked"))?;

    write_progress.finish_with_message(format!(
        "Completed! Written {} read pairs (Parallel)",
        ncached
    ));
    println!(
        "Writing finished. \nObserved {} unique read ids. \nWrote {} read pairs ({} cached) using {} threads",
        total_read_pairs,
        fq.total_written(),
        ncached,
        num_threads
    );
    Ok(fq.paths())
}
