//! Advanced caching implementation with efficient data structures

use rust_htslib::bam::record::Record;
use std::collections::{HashMap, BTreeMap};
use crate::memory_pool::ByteBufferPool;

/// Cache entry with metadata for efficient eviction
#[derive(Debug, Clone)]
struct CacheEntry {
    record: Record,
    tid: i32,
    _pos: i64,
    insert_time: u64,
}

impl CacheEntry {
    fn new(record: Record, insert_time: u64) -> Self {
        let tid = record.tid();
        let _pos = record.pos();
        Self {
            record,
            tid,
            _pos,
            insert_time,
        }
    }
}

/// Position-based key for efficient spatial queries
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct PositionKey {
    tid: i32,
    pos: i64,
    _qname_hash: u64,
}

impl PositionKey {
    fn new(tid: i32, pos: i64, qname: &[u8]) -> Self {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        let mut hasher = DefaultHasher::new();
        qname.hash(&mut hasher);
        let qname_hash = hasher.finish();
        
        Self {
            tid,
            pos,
            _qname_hash: qname_hash,
        }
    }
}

/// Advanced read-pair cache with multiple indexing strategies
pub struct AdvancedRpCache {
    // Primary index by qname for fast pair lookup
    qname_index: HashMap<Vec<u8>, CacheEntry>,
    
    // Secondary index by position for spatial queries
    position_index: BTreeMap<PositionKey, Vec<u8>>,
    
    // Time-based index for LRU eviction
    time_index: BTreeMap<u64, Vec<Vec<u8>>>,
    
    // Configuration
    cache_size: usize,
    relaxed: bool,
    current_time: u64,
    
    // Memory optimization
    _buffer_pool: ByteBufferPool,
    
    // Pre-allocated buffers for batch operations
    eviction_candidates: Vec<Vec<u8>>,
    orphans_buffer: Vec<Record>,
}

impl AdvancedRpCache {
    pub fn new(cache_size: usize, relaxed: bool) -> Self {
        Self {
            qname_index: HashMap::with_capacity(cache_size),
            position_index: BTreeMap::new(),
            time_index: BTreeMap::new(),
            cache_size,
            relaxed,
            current_time: 0,
            _buffer_pool: ByteBufferPool::new(50, 256),
            eviction_candidates: Vec::with_capacity(cache_size / 4),
            orphans_buffer: Vec::with_capacity(cache_size / 4),
        }
    }

    /// Cache a record and return paired records if found
    pub fn cache_rec(&mut self, rec: Record) -> Option<(Record, Record)> {
        self.current_time += 1;
        let qname = Vec::from(rec.qname());
        
        // Check if we already have the mate
        if let Some(cached_entry) = self.qname_index.remove(&qname) {
            // Remove from secondary indices
            let pos_key = PositionKey::new(cached_entry.tid, cached_entry._pos, &qname);
            self.position_index.remove(&pos_key);
            
            // Remove from time index
            if let Some(time_entries) = self.time_index.get_mut(&cached_entry.insert_time) {
                time_entries.retain(|q| q != &qname);
                if time_entries.is_empty() {
                    self.time_index.remove(&cached_entry.insert_time);
                }
            }
            
            // Return the pair in correct order
            if rec.is_first_in_template() && cached_entry.record.is_last_in_template() {
                Some((rec, cached_entry.record))
            } else if cached_entry.record.is_first_in_template() && rec.is_last_in_template() {
                Some((cached_entry.record, rec))
            } else {
                if self.relaxed {
                    eprintln!(
                        "Found extra BAM record for qname: {}. Skipping due to --relaxed",
                        String::from_utf8_lossy(&qname)
                    );
                    return None;
                }
                panic!(
                    "Invalid BAM record pair for qname: {}",
                    String::from_utf8_lossy(&qname)
                );
            }
        } else {
            // Cache the record
            let entry = CacheEntry::new(rec, self.current_time);
            let pos_key = PositionKey::new(entry.tid, entry._pos, &qname);
            
            // Insert into all indices
            self.qname_index.insert(qname.clone(), entry);
            self.position_index.insert(pos_key, qname.clone());
            
            // Add to time index
            self.time_index
                .entry(self.current_time)
                .or_default()
                .push(qname);
            
            // Check if we need to evict
            if self.qname_index.len() > self.cache_size {
                self.evict_old_entries();
            }
            
            None
        }
    }

    /// Advanced cache eviction using multiple strategies
    pub fn clear_orphans(&mut self, current_tid: i32, current_pos: i64) -> Vec<Record> {
        self.orphans_buffer.clear();
        self.eviction_candidates.clear();
        
        // Strategy 1: Evict by distance (spatial locality)
        self.evict_by_distance(current_tid, current_pos);
        
        // Strategy 2: Evict by age if still over capacity
        if self.qname_index.len() > self.cache_size / 2 {
            self.evict_by_age();
        }
        
        // Strategy 3: Emergency eviction if still over capacity
        if self.qname_index.len() > self.cache_size * 3 / 4 {
            self.emergency_evict();
        }
        
        std::mem::take(&mut self.orphans_buffer)
    }
    
    /// Evict entries based on genomic distance
    fn evict_by_distance(&mut self, current_tid: i32, current_pos: i64) {
        let max_distance = 10000; // 10kb
        
        // Find entries to evict based on position
        for (pos_key, qname) in self.position_index.iter() {
            if pos_key.tid != current_tid || (current_pos - pos_key.pos).abs() > max_distance {
                self.eviction_candidates.push(qname.clone());
            }
        }
        
        // Remove evicted entries
        for qname in self.eviction_candidates.drain(..) {
            if let Some(entry) = self.qname_index.remove(&qname) {
                let pos_key = PositionKey::new(entry.tid, entry._pos, &qname);
                self.position_index.remove(&pos_key);
                
                // Remove from time index
                if let Some(time_entries) = self.time_index.get_mut(&entry.insert_time) {
                    time_entries.retain(|q| q != &qname);
                    if time_entries.is_empty() {
                        self.time_index.remove(&entry.insert_time);
                    }
                }
                
                self.orphans_buffer.push(entry.record);
            }
        }
    }
    
    /// Evict oldest entries (LRU strategy)
    fn evict_by_age(&mut self) {
        let target_size = self.cache_size / 2;
        let mut _evicted_count = 0;
        
        // Collect old entries
        let mut old_times: Vec<_> = self.time_index.keys().cloned().collect();
        old_times.sort();
        
        for time in old_times {
            if self.qname_index.len() <= target_size {
                break;
            }
            
            if let Some(qnames) = self.time_index.remove(&time) {
                for qname in qnames {
                    if let Some(entry) = self.qname_index.remove(&qname) {
                        let pos_key = PositionKey::new(entry.tid, entry._pos, &qname);
                        self.position_index.remove(&pos_key);
                        self.orphans_buffer.push(entry.record);
                        _evicted_count += 1;
                        
                        if self.qname_index.len() <= target_size {
                            break;
                        }
                    }
                }
            }
        }
    }
    
    /// Emergency eviction when cache is critically full
    fn emergency_evict(&mut self) {
        let _emergency_target = self.cache_size / 4;
        
        // Clear everything if we're in emergency mode
        self.orphans_buffer.reserve(self.qname_index.len());
        for (_, entry) in self.qname_index.drain() {
            self.orphans_buffer.push(entry.record);
        }
        
        self.position_index.clear();
        self.time_index.clear();
        
        // Recreate with proper capacity
        self.qname_index = HashMap::with_capacity(self.cache_size);
    }
    
    /// Evict old entries when cache is full
    fn evict_old_entries(&mut self) {
        let target_size = self.cache_size * 3 / 4;
        
        while self.qname_index.len() > target_size {
            // Find oldest entry
            if let Some((&oldest_time, _)) = self.time_index.iter().next() {
                if let Some(qnames) = self.time_index.remove(&oldest_time) {
                    for qname in qnames {
                        if let Some(entry) = self.qname_index.remove(&qname) {
                            let pos_key = PositionKey::new(entry.tid, entry._pos, &qname);
                            self.position_index.remove(&pos_key);
                            // Don't add to orphans buffer here as this is just capacity management
                        }
                        
                        if self.qname_index.len() <= target_size {
                            break;
                        }
                    }
                }
            } else {
                break; // No more entries to evict
            }
        }
    }

    /// Get current cache size
    pub fn len(&self) -> usize {
        self.qname_index.len()
    }
    
    /// Check if cache is empty
    pub fn is_empty(&self) -> bool {
        self.qname_index.is_empty()
    }
    
    /// Get cache statistics
    pub fn stats(&self) -> CacheStats {
        CacheStats {
            total_entries: self.qname_index.len(),
            position_index_size: self.position_index.len(),
            time_index_size: self.time_index.len(),
            current_time: self.current_time,
        }
    }
}

/// Cache statistics for monitoring
#[derive(Debug, Clone)]
pub struct CacheStats {
    pub total_entries: usize,
    pub position_index_size: usize,
    pub time_index_size: usize,
    pub current_time: u64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Record;
    
    fn create_test_record(qname: &str, tid: i32, pos: i64, is_first: bool) -> Record {
        let mut record = Record::new();
        record.set_qname(qname.as_bytes());
        record.set_tid(tid);
        record.set_pos(pos);
        if is_first {
            record.set_flags(0x40); // First in template
        } else {
            record.set_flags(0x80); // Last in template
        }
        record
    }
    
    #[test]
    fn test_advanced_cache_pairing() {
        let mut cache = AdvancedRpCache::new(100, false);
        
        let r1 = create_test_record("read1", 0, 1000, true);
        let r2 = create_test_record("read1", 0, 1100, false);
        
        // Cache first read
        assert!(cache.cache_rec(r1).is_none());
        assert_eq!(cache.len(), 1);
        
        // Cache second read - should return pair
        let pair = cache.cache_rec(r2);
        assert!(pair.is_some());
        assert_eq!(cache.len(), 0);
    }
    
    #[test]
    fn test_distance_based_eviction() {
        let mut cache = AdvancedRpCache::new(10, true);
        
        // Add records at different positions
        for i in 0..5 {
            let r = create_test_record(&format!("read{}", i), 0, i * 1000, true);
            cache.cache_rec(r);
        }
        
        assert_eq!(cache.len(), 5);
        
        // Clear orphans with current position far away
        let orphans = cache.clear_orphans(0, 20000);
        assert!(!orphans.is_empty());
    }
}