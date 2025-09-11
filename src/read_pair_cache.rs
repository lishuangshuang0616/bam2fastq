
use rust_htslib::bam::record::Record;
use std::collections::HashMap;
use crate::memory_pool::{ByteBufferPool, StringPool};

/// Read-pair cache with optimized memory management
/// Uses pre-allocated buffers and efficient data structures
pub struct RpCache {
    pub cache_size: usize,
    pub cache: HashMap<Vec<u8>, Record>,
    pub relaxed: bool,
    buffer_pool: ByteBufferPool,
    string_pool: StringPool,
    // Pre-allocated vectors for batch operations
    orphan_keys_buffer: Vec<Vec<u8>>,
    orphans_buffer: Vec<Record>,
}

impl RpCache {
    pub fn new(cache_size: usize, relaxed: bool) -> RpCache {
        RpCache {
            cache: HashMap::with_capacity(cache_size),
            cache_size,
            relaxed,
            buffer_pool: ByteBufferPool::new(20, 256), // Pool for qname buffers
            string_pool: StringPool::new(10, 128), // Pool for string operations
            orphan_keys_buffer: Vec::with_capacity(cache_size / 4),
            orphans_buffer: Vec::with_capacity(cache_size / 4),
        }
    }

    pub fn cache_rec(&mut self, rec: Record) -> Option<(Record, Record)> {
        // If cache already has entry, we have a pair! Return both
        match self.cache.remove(rec.qname()) {
            Some(old_rec) => {
                if rec.is_first_in_template() && old_rec.is_last_in_template() {
                    Some((rec, old_rec))
                } else if old_rec.is_first_in_template() && rec.is_last_in_template() {
                    Some((old_rec, rec))
                } else {
                    if self.relaxed {
                        println!(
                            "Found extra BAM record for qname: {}. Skipping due to --relaxed",
                            String::from_utf8_lossy(rec.qname())
                        );
                        return None;
                    }

                    println!(
                        "Found invalid set of BAM record for qname: {}.",
                        String::from_utf8_lossy(rec.qname())
                    );
                    println!("This may be caused by inputting the same FASTQ record to Long Ranger twice");
                    panic!(
                        "invalid BAM record detected: {}",
                        String::from_utf8_lossy(rec.qname())
                    );
                }
            }
            None => {
                self.cache.insert(Vec::from(rec.qname()), rec);
                None
            }
        }
    }

    pub fn clear_orphans(&mut self, current_tid: i32, current_pos: i64) -> Vec<Record> {
        // Clear and reuse pre-allocated buffers
        self.orphans_buffer.clear();
        self.orphan_keys_buffer.clear();

        let mut dist = 5000;
        let target_size = self.cache_size / 2;

        while self.cache.len() > target_size && dist > 100 {
            // Collect keys to remove in pre-allocated buffer
            for (key, rec) in self.cache.iter() {
                // Evict unmapped reads, reads on a previous chromosome, or reads that are >dist behind the current position
                if rec.tid() == -1
                    || (current_pos - rec.pos()).abs() > dist
                    || rec.tid() != current_tid
                {
                    self.orphan_keys_buffer.push(key.clone());
                }
            }

            // Remove records and collect orphans
            self.orphans_buffer.reserve(self.orphan_keys_buffer.len());
            for key in self.orphan_keys_buffer.drain(..) {
                if let Some(rec) = self.cache.remove(&key) {
                    self.orphans_buffer.push(rec);
                }
            }

            dist /= 2;
        }

        // Cache got too full -- clear everything
        if dist <= 100 {
            self.orphans_buffer.reserve(self.cache.len());
            for (_, rec) in self.cache.drain() {
                self.orphans_buffer.push(rec);
            }
            // Recreate HashMap with proper capacity
            self.cache = HashMap::with_capacity(self.cache_size);
        }

        // Return orphans by moving from buffer
        std::mem::take(&mut self.orphans_buffer)
    }

    pub fn len(&self) -> usize {
        self.cache.len()
    }
}