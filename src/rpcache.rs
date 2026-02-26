use rust_htslib::bam::record::Record;
use rustc_hash::FxHashMap;

/// Position-based key for efficient spatial queries
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct PositionKey {
    tid: i32,
    pos: i64,
}

impl std::hash::Hash for PositionKey {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_i32(self.tid);
        state.write_i64(self.pos);
    }
}

/// Read-pair cache with spatial indexing for efficient eviction
/// Let's us stream through the BAM and find nearby mates so we can write them out immediately
/// Reads whose mate is not found promptly are kept in memory with efficient indexing
pub struct RpCache {
    pub cache_size: usize,
    pub cache: FxHashMap<Vec<u8>, Record>,
    pub relaxed: bool,
    // Spatial index for efficient eviction
    // FxHashMap is O(1) compared to BTreeMap's O(log N). Order isn't strictly needed for eviction grouping.
    position_index: FxHashMap<PositionKey, Vec<Vec<u8>>>,
}

impl RpCache {
    pub fn new(cache_size: usize, relaxed: bool) -> RpCache {
        RpCache {
            cache: FxHashMap::with_capacity_and_hasher(cache_size, Default::default()),
            cache_size,
            relaxed,
            position_index: FxHashMap::default(),
        }
    }

    pub fn cache_rec(&mut self, rec: Record) -> Option<(Record, Record)> {
        let qname = Vec::from(rec.qname());

        // If cache already has entry, we have a pair! Return both
        match self.cache.remove(&qname) {
            Some(old_rec) => {
                // Remove from position index
                let pos_key = PositionKey {
                    tid: old_rec.tid(),
                    pos: old_rec.pos(),
                };
                if let Some(qnames) = self.position_index.get_mut(&pos_key) {
                    qnames.retain(|q| q != &qname);
                    if qnames.is_empty() {
                        self.position_index.remove(&pos_key);
                    }
                }

                if rec.is_first_in_template() && old_rec.is_last_in_template() {
                    Some((rec, old_rec))
                } else if old_rec.is_first_in_template() && rec.is_last_in_template() {
                    Some((old_rec, rec))
                } else {
                    if self.relaxed {
                        println!(
                            "Found extra BAM record for qname: {}. Skipping due to --relaxed",
                            String::from_utf8_lossy(&qname)
                        );
                        return None;
                    }

                    println!(
                        "Found invalid set of BAM record for qname: {}.",
                        String::from_utf8_lossy(&qname)
                    );
                    println!("This may be caused by inputting the same FASTQ record to Long Ranger twice");
                    panic!(
                        "invalid BAM record detected: {}",
                        String::from_utf8_lossy(&qname)
                    );
                }
            }
            None => {
                // Add to cache and position index
                let pos_key = PositionKey {
                    tid: rec.tid(),
                    pos: rec.pos(),
                };
                self.position_index
                    .entry(pos_key)
                    .or_default()
                    .push(qname.clone());
                self.cache.insert(qname, rec);
                None
            }
        }
    }

    pub fn clear_orphans(&mut self, current_tid: i32, current_pos: i64) -> Vec<Record> {
        let mut orphans = Vec::new();
        let max_distance = 10000; // 10kb

        // Single-pass: drain orphans and clean position index simultaneously
        self.position_index.retain(|pos_key, qnames| {
            if pos_key.tid != current_tid || (current_pos - pos_key.pos).abs() > max_distance {
                // Evict all reads at this position
                for key in qnames.iter() {
                    if let Some(rec) = self.cache.remove(key) {
                        orphans.push(rec);
                    }
                }
                false // Remove this position from index
            } else {
                // Still in range — clean up any stale qnames
                qnames.retain(|qname| self.cache.contains_key(qname));
                !qnames.is_empty()
            }
        });

        // Emergency eviction if cache still too large after positional eviction
        if self.cache.len() > self.cache_size * 3 / 4 {
            let excess = self.cache.len() - self.cache_size / 2;
            let mut count = 0;

            while count < excess && !self.cache.is_empty() {
                if let Some(key) = self.cache.keys().next().cloned() {
                    if let Some(rec) = self.cache.remove(&key) {
                        orphans.push(rec);
                        count += 1;
                    }
                }
            }

            // Rebuild position index after emergency eviction
            self.position_index.clear();
            for (key, rec) in &self.cache {
                let pos_key = PositionKey {
                    tid: rec.tid(),
                    pos: rec.pos(),
                };
                self.position_index
                    .entry(pos_key)
                    .or_default()
                    .push(key.clone());
            }
        }

        orphans
    }

    pub fn len(&self) -> usize {
        self.cache.len()
    }
}
