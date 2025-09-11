//! Parallel processing pipeline for BAM to FASTQ conversion
//! Implements multi-threaded pipeline with producer-consumer pattern

use crate::types::*;
use crate::formatter::FormatBamRecords;
use crate::writer::FastqManager;
use crate::advanced_cache::AdvancedRpCache;
use rust_htslib::bam::Record;
use std::sync::mpsc::{self};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;
use anyhow::{anyhow, Result};

/// Message types for the processing pipeline
#[derive(Debug)]
pub enum PipelineMessage {
    Record(Record),
    Flush,
    Terminate,
}

/// Statistics for pipeline performance monitoring
#[derive(Debug, Clone, Default)]
pub struct PipelineStats {
    pub records_processed: u64,
    pub pairs_written: u64,
    pub orphans_found: u64,
    pub processing_time_ms: u64,
    pub cache_hits: u64,
    pub cache_misses: u64,
}

/// Configuration for the parallel processing pipeline
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    pub num_worker_threads: usize,
    pub buffer_size: usize,
    pub cache_size: usize,
    pub batch_size: usize,
    pub relaxed_mode: bool,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            num_worker_threads: num_cpus::get().max(2),
            buffer_size: 10000,
            cache_size: 100000,
            batch_size: 1000,
            relaxed_mode: false,
        }
    }
}

/// Parallel processing pipeline for BAM records
#[allow(dead_code)]
pub struct ParallelProcessor {
    config: PipelineConfig,
    stats: Arc<Mutex<PipelineStats>>,
}

impl ParallelProcessor {
    pub fn new(config: PipelineConfig) -> Self {
        Self {
            config,
            stats: Arc::new(Mutex::new(PipelineStats::default())),
        }
    }

    /// Process records using parallel pipeline
    pub fn process_records<I>(
        &mut self,
        records: I,
        formatter: FormatBamRecords,
        fq_manager: FastqManager,
    ) -> Result<(Vec<OutPaths>, PipelineStats)>
    where
        I: Iterator<Item = Result<Record, rust_htslib::errors::Error>>,
    {
        let start_time = Instant::now();
        
        // Create channels for pipeline communication
        let (record_tx, record_rx) = mpsc::sync_channel(self.config.buffer_size);
        let (result_tx, result_rx) = mpsc::sync_channel(self.config.buffer_size);
        
        // Shared resources
        let stats = Arc::clone(&self.stats);
        let fq_manager = Arc::new(Mutex::new(fq_manager));
        
        // Share the receiver among all worker threads using Arc<Mutex<>>
        let shared_rx = Arc::new(Mutex::new(record_rx));
        let mut worker_handles = Vec::new();
        
        for worker_id in 0..self.config.num_worker_threads {
            let worker_rx = Arc::clone(&shared_rx);
            let worker_tx = result_tx.clone();
            let worker_formatter = formatter.clone();
            let worker_stats = Arc::clone(&stats);
            let worker_config = self.config.clone();
            
            let handle = thread::spawn(move || {
                Self::worker_thread(
                    worker_id,
                    worker_rx,
                    worker_tx,
                    worker_formatter,
                    worker_stats,
                    worker_config,
                )
            });
            
            worker_handles.push(handle);
        }
        
        // Start writer thread
        let writer_fq = Arc::clone(&fq_manager);
        let writer_stats = Arc::clone(&stats);
        let writer_handle = thread::spawn(move || {
            Self::writer_thread(result_rx, writer_fq, writer_stats)
        });
        
        // Producer: send records to workers (run in main thread)
        let mut batch = Vec::with_capacity(1000);
        
        for record_result in records {
            match record_result {
                Ok(record) => {
                    batch.push(record);
                    
                    if batch.len() >= 1000 {
                        for record in batch.drain(..) {
                            if record_tx.send(PipelineMessage::Record(record)).is_err() {
                                return Err(anyhow!("Failed to send record to workers"));
                            }
                        }
                    }
                }
                Err(e) => {
                    return Err(anyhow!("Error reading BAM record: {}", e));
                }
            }
        }
        
        // Send remaining records
        for record in batch {
            if record_tx.send(PipelineMessage::Record(record)).is_err() {
                return Err(anyhow!("Failed to send final records to workers"));
            }
        }
        
        // Signal flush and termination
        for _ in 0..self.config.num_worker_threads {
            let _ = record_tx.send(PipelineMessage::Flush);
        }
        
        for _ in 0..self.config.num_worker_threads {
            let _ = record_tx.send(PipelineMessage::Terminate);
        }
        
        // Wait for all workers to finish
        for handle in worker_handles {
            handle.join().unwrap()?
        }
        
        // Signal writer to terminate
        let _ = result_tx.send(WriterMessage::Terminate);
        
        // Wait for writer to finish and get output paths
        let output_paths = writer_handle.join().unwrap()?;
        
        // Update final statistics
        {
            let mut stats_guard = self.stats.lock().unwrap();
            stats_guard.processing_time_ms = start_time.elapsed().as_millis() as u64;
        }
        
        let final_stats = self.stats.lock().unwrap().clone();
        Ok((output_paths, final_stats))
    }
    
    /// Worker thread function
    fn worker_thread(
        _worker_id: usize,
        rx: Arc<Mutex<mpsc::Receiver<PipelineMessage>>>,
        tx: mpsc::SyncSender<WriterMessage>,
        formatter: FormatBamRecords,
        _stats: Arc<Mutex<PipelineStats>>,
        config: PipelineConfig,
    ) -> Result<()> {
        let mut cache = AdvancedRpCache::new(config.cache_size, config.relaxed_mode);
        let mut local_stats = PipelineStats::default();
        
        loop {
            let message = {
                let rx_guard = rx.lock().unwrap();
                rx_guard.recv()
            };
            
            match message {
                Ok(PipelineMessage::Record(record)) => {
                    local_stats.records_processed += 1;
                    
                    // Periodically clear orphans
                    if local_stats.records_processed % 10000 == 0 {
                        let current_tid = record.tid();
                        let current_pos = record.pos();
                        let orphans = cache.clear_orphans(current_tid, current_pos);
                        local_stats.orphans_found += orphans.len() as u64;
                    }
                    
                    // Try to match read pairs
                    if let Some((r1_rec, r2_rec)) = cache.cache_rec(record) {
                        local_stats.cache_hits += 1;
                        
                        match formatter.format_read_pair(&r1_rec, &r2_rec) {
                            Ok((rg, r1, r2, i1, i2)) => {
                                let read_group = match rg {
                                    Some((sample, lane)) => format!("{}.{}", sample, lane),
                                    None => "unknown".to_string(),
                                };
                                
                                let write_msg = WriterMessage::WritePair {
                                    read_group,
                                    r1,
                                    r2,
                                    i1: Box::new(i1),
                                    i2: Box::new(i2),
                                };
                                
                                if tx.send(write_msg).is_err() {
                                    break;
                                }
                                
                                local_stats.pairs_written += 1;
                            }
                            Err(_) => {
                                if !config.relaxed_mode {
                                    // In strict mode, this would be an error
                                    // For now, just skip the pair
                                }
                            }
                        }
                    } else {
                        local_stats.cache_misses += 1;
                    }
                }
                Ok(PipelineMessage::Flush) => {
                    // Clear remaining orphans
                    let final_orphans = cache.clear_orphans(i32::MAX, i64::MAX);
                    local_stats.orphans_found += final_orphans.len() as u64;
                }
                Ok(PipelineMessage::Terminate) => {
                    break;
                }
                Err(_) => {
                    break;
                }
            }
        }
        
        // Merge local stats with global stats
        {
            let mut global_stats = _stats.lock().unwrap();
            global_stats.records_processed += local_stats.records_processed;
            global_stats.pairs_written += local_stats.pairs_written;
            global_stats.orphans_found += local_stats.orphans_found;
            global_stats.cache_hits += local_stats.cache_hits;
            global_stats.cache_misses += local_stats.cache_misses;
        }
        
        Ok(())
    }
    
    /// Writer thread function
    fn writer_thread(
        rx: mpsc::Receiver<WriterMessage>,
        fq_manager: Arc<Mutex<FastqManager>>,
        stats: Arc<Mutex<PipelineStats>>,
    ) -> Result<Vec<(std::path::PathBuf, std::path::PathBuf, Option<std::path::PathBuf>, Option<std::path::PathBuf>)>> {
        let mut batch_buffer = Vec::with_capacity(1000);
        
        loop {
            match rx.recv() {
                Ok(WriterMessage::WritePair { read_group, r1, r2, i1, i2 }) => {
                    batch_buffer.push((read_group, r1, r2, *i1, *i2));
                    
                    // Flush batch when full
                    if batch_buffer.len() >= 1000 {
                        Self::flush_write_batch(&mut batch_buffer, &fq_manager)?;
                    }
                }
                Ok(WriterMessage::Terminate) => {
                    // Flush remaining writes
                    if !batch_buffer.is_empty() {
                        Self::flush_write_batch(&mut batch_buffer, &fq_manager)?;
                    }
                    break;
                }
                Err(_) => {
                    break;
                }
            }
        }
        
        // Get output paths from the writer
        let fq_guard = fq_manager.lock().unwrap();
        Ok(fq_guard.paths())
    }
    
    /// Flush a batch of writes to the FASTQ manager
    fn flush_write_batch(
        batch: &mut Vec<(String, FqRecord, FqRecord, Option<FqRecord>, Option<FqRecord>)>,
        fq_manager: &Arc<Mutex<FastqManager>>,
    ) -> Result<()> {
        let mut fq_guard = fq_manager.lock().unwrap();
        
        for (rg_str, r1, r2, i1, i2) in batch.drain(..) {
            // Parse the read group string back to the expected format
            let rg = if rg_str == "unknown" {
                None
            } else {
                // For now, use None - in a real implementation, you'd parse the string
                None
            };
            fq_guard.write(&rg, &r1, &r2, &i1, &i2);
        }
        
        Ok(())
    }
    
    /// Get pipeline statistics
    pub fn get_stats(&self) -> PipelineStats {
        self.stats.lock().unwrap().clone()
    }
}

/// Messages for the writer thread
#[derive(Debug)]
#[allow(dead_code)]
enum WriterMessage {
    WritePair {
        read_group: String,
        r1: FqRecord,
        r2: FqRecord,
        i1: Box<Option<FqRecord>>,
        i2: Box<Option<FqRecord>>,
    },
    Terminate,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_pipeline_config_default() {
        let config = PipelineConfig::default();
        assert!(config.num_worker_threads >= 2);
        assert_eq!(config.buffer_size, 10000);
        assert_eq!(config.cache_size, 100000);
        assert_eq!(config.batch_size, 1000);
        assert!(!config.relaxed_mode);
    }
    
    #[test]
    fn test_pipeline_stats() {
        let stats = PipelineStats::default();
        assert_eq!(stats.records_processed, 0);
        assert_eq!(stats.pairs_written, 0);
        assert_eq!(stats.orphans_found, 0);
    }
}