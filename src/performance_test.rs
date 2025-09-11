//! Performance testing utilities for memory optimization validation

use std::time::Instant;
use std::collections::HashMap;
use crate::memory_pool::{ByteBufferPool, ThreadSafeByteBufferPool};

/// Benchmark buffer allocation patterns
pub struct MemoryBenchmark {
    iterations: usize,
}

impl MemoryBenchmark {
    pub fn new(iterations: usize) -> Self {
        Self { iterations }
    }

    /// Test traditional Vec allocation vs pool allocation
    pub fn benchmark_buffer_allocation(&self) -> (u128, u128) {
        // Traditional allocation benchmark
        let start = Instant::now();
        for _ in 0..self.iterations {
            let mut buffer = Vec::with_capacity(1024);
            buffer.extend_from_slice(&[0u8; 512]);
            buffer.clear();
        }
        let traditional_time = start.elapsed().as_nanos();

        // Pool allocation benchmark
        let pool = ByteBufferPool::new(50, 1024);
        let start = Instant::now();
        for _ in 0..self.iterations {
            let mut buffer = pool.get();
            buffer.extend_from_slice(&[0u8; 512]);
            pool.put(buffer);
        }
        let pool_time = start.elapsed().as_nanos();

        (traditional_time, pool_time)
    }

    /// Test HashMap with pre-allocation vs dynamic allocation
    pub fn benchmark_hashmap_allocation(&self) -> (u128, u128) {
        // Dynamic allocation
        let start = Instant::now();
        for _ in 0..self.iterations {
            let mut map: HashMap<String, i32> = HashMap::new();
            for i in 0..100 {
                map.insert(format!("key_{}", i), i);
            }
        }
        let dynamic_time = start.elapsed().as_nanos();

        // Pre-allocated
        let start = Instant::now();
        for _ in 0..self.iterations {
            let mut map: HashMap<String, i32> = HashMap::with_capacity(100);
            for i in 0..100 {
                map.insert(format!("key_{}", i), i);
            }
        }
        let preallocated_time = start.elapsed().as_nanos();

        (dynamic_time, preallocated_time)
    }

    /// Test thread-safe buffer pool performance
    pub fn benchmark_thread_safe_pool(&self) -> u128 {
        let pool = ThreadSafeByteBufferPool::new(50, 1024);
        
        let start = Instant::now();
        for _ in 0..self.iterations {
            let mut buffer = pool.get();
            buffer.extend_from_slice(&[0u8; 512]);
            pool.put(buffer);
        }
        start.elapsed().as_nanos()
    }

    /// Print benchmark results
    pub fn run_all_benchmarks(&self) {
        println!("Running memory optimization benchmarks with {} iterations...", self.iterations);
        
        let (trad_time, pool_time) = self.benchmark_buffer_allocation();
        println!("Buffer Allocation:");
        println!("  Traditional: {} ns", trad_time);
        println!("  Pool-based:  {} ns", pool_time);
        println!("  Improvement: {:.2}%", 
                 ((trad_time as f64 - pool_time as f64) / trad_time as f64) * 100.0);
        
        let (dyn_time, pre_time) = self.benchmark_hashmap_allocation();
        println!("\nHashMap Allocation:");
        println!("  Dynamic:      {} ns", dyn_time);
        println!("  Pre-allocated: {} ns", pre_time);
        println!("  Improvement: {:.2}%", 
                 ((dyn_time as f64 - pre_time as f64) / dyn_time as f64) * 100.0);
        
        let ts_time = self.benchmark_thread_safe_pool();
        println!("\nThread-safe Pool: {} ns", ts_time);
        
        println!("\nMemory optimization benchmarks completed.");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_memory_benchmark() {
        let benchmark = MemoryBenchmark::new(1000);
        let (trad_time, pool_time) = benchmark.benchmark_buffer_allocation();
        
        // Pool should generally be faster or at least not significantly slower
        assert!(pool_time <= trad_time * 2, "Pool allocation is too slow compared to traditional");
    }

    #[test]
    fn test_hashmap_benchmark() {
        let benchmark = MemoryBenchmark::new(100);
        let (dyn_time, pre_time) = benchmark.benchmark_hashmap_allocation();
        
        // Pre-allocated should be faster
        assert!(pre_time <= dyn_time, "Pre-allocated HashMap should be faster");
    }
}