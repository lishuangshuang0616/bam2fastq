//! Memory pool implementation for efficient memory management

use std::collections::VecDeque;
use std::sync::{Arc, Mutex};
use std::cell::RefCell;

/// A simple object pool for reusing Vec<u8> buffers
#[derive(Debug)]
pub struct ByteBufferPool {
    pool: RefCell<VecDeque<Vec<u8>>>,
    max_size: usize,
    buffer_capacity: usize,
}

impl Clone for ByteBufferPool {
    fn clone(&self) -> Self {
        Self::new(self.max_size, self.buffer_capacity)
    }
}

impl ByteBufferPool {
    /// Create a new buffer pool
    pub fn new(max_size: usize, buffer_capacity: usize) -> Self {
        Self {
            pool: RefCell::new(VecDeque::with_capacity(max_size)),
            max_size,
            buffer_capacity,
        }
    }

    /// Get a buffer from the pool or create a new one
    pub fn get(&self) -> Vec<u8> {
        let mut pool = self.pool.borrow_mut();
        if let Some(mut buffer) = pool.pop_front() {
            buffer.clear();
            buffer
        } else {
            Vec::with_capacity(self.buffer_capacity)
        }
    }

    /// Return a buffer to the pool
    pub fn put(&self, mut buffer: Vec<u8>) {
        let mut pool = self.pool.borrow_mut();
        if pool.len() < self.max_size {
            // Only keep buffers that haven't grown too large
            if buffer.capacity() <= self.buffer_capacity * 2 {
                buffer.clear();
                pool.push_back(buffer);
            }
        }
    }

    /// Get current pool size
    pub fn size(&self) -> usize {
        self.pool.borrow().len()
    }
}

/// Thread-safe buffer pool
pub struct ThreadSafeByteBufferPool {
    pool: Arc<Mutex<VecDeque<Vec<u8>>>>,
    max_size: usize,
    buffer_capacity: usize,
}

impl ThreadSafeByteBufferPool {
    /// Create a new thread-safe buffer pool
    pub fn new(max_size: usize, buffer_capacity: usize) -> Self {
        Self {
            pool: Arc::new(Mutex::new(VecDeque::with_capacity(max_size))),
            max_size,
            buffer_capacity,
        }
    }

    /// Get a buffer from the pool or create a new one
    pub fn get(&self) -> Vec<u8> {
        let mut pool = self.pool.lock().unwrap();
        if let Some(mut buffer) = pool.pop_front() {
            buffer.clear();
            buffer
        } else {
            Vec::with_capacity(self.buffer_capacity)
        }
    }

    /// Return a buffer to the pool
    pub fn put(&self, mut buffer: Vec<u8>) {
        let mut pool = self.pool.lock().unwrap();
        if pool.len() < self.max_size {
            // Only keep buffers that haven't grown too large
            if buffer.capacity() <= self.buffer_capacity * 2 {
                buffer.clear();
                pool.push_back(buffer);
            }
        }
    }

    /// Clone the pool for sharing across threads
    pub fn clone(&self) -> Self {
        Self {
            pool: Arc::clone(&self.pool),
            max_size: self.max_size,
            buffer_capacity: self.buffer_capacity,
        }
    }
}

/// RAII wrapper for automatic buffer return to pool
pub struct PooledBuffer<'a> {
    buffer: Option<Vec<u8>>,
    pool: &'a ByteBufferPool,
}

impl<'a> PooledBuffer<'a> {
    pub fn new(pool: &'a ByteBufferPool) -> Self {
        let buffer = pool.get();
        Self {
            buffer: Some(buffer),
            pool,
        }
    }

    /// Get mutable reference to the buffer
    pub fn as_mut(&mut self) -> &mut Vec<u8> {
        self.buffer.as_mut().unwrap()
    }

    /// Get reference to the buffer
    pub fn as_ref(&self) -> &Vec<u8> {
        self.buffer.as_ref().unwrap()
    }

    /// Take ownership of the buffer (prevents automatic return to pool)
    pub fn take(mut self) -> Vec<u8> {
        self.buffer.take().unwrap()
    }
}

impl Drop for PooledBuffer<'_> {
    fn drop(&mut self) {
        if let Some(buffer) = self.buffer.take() {
            self.pool.put(buffer);
        }
    }
}

/// Pre-allocated string pool for common operations
pub struct StringPool {
    pool: RefCell<VecDeque<String>>,
    max_size: usize,
    initial_capacity: usize,
}

impl StringPool {
    pub fn new(max_size: usize, initial_capacity: usize) -> Self {
        Self {
            pool: RefCell::new(VecDeque::with_capacity(max_size)),
            max_size,
            initial_capacity,
        }
    }

    pub fn get(&self) -> String {
        let mut pool = self.pool.borrow_mut();
        if let Some(mut s) = pool.pop_front() {
            s.clear();
            s
        } else {
            String::with_capacity(self.initial_capacity)
        }
    }

    pub fn put(&self, mut s: String) {
        let mut pool = self.pool.borrow_mut();
        if pool.len() < self.max_size && s.capacity() <= self.initial_capacity * 2 {
            s.clear();
            pool.push_back(s);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buffer_pool() {
        let pool = ByteBufferPool::new(5, 1024);
        
        // Get a buffer
        let mut buf1 = pool.get();
        buf1.extend_from_slice(b"test");
        assert_eq!(buf1.len(), 4);
        
        // Return it
        pool.put(buf1);
        assert_eq!(pool.size(), 1);
        
        // Get it back
        let buf2 = pool.get();
        assert_eq!(buf2.len(), 0); // Should be cleared
        assert!(buf2.capacity() >= 1024);
    }

    #[test]
    fn test_pooled_buffer_raii() {
        let pool = ByteBufferPool::new(5, 1024);
        
        {
            let mut pooled = PooledBuffer::new(&pool);
            pooled.as_mut().extend_from_slice(b"test");
        } // Buffer should be automatically returned here
        
        assert_eq!(pool.size(), 1);
    }

    #[test]
    fn test_thread_safe_pool() {
        let pool = ThreadSafeByteBufferPool::new(5, 1024);
        let pool_clone = pool.clone();
        
        let buf = pool.get();
        pool_clone.put(buf);
        
        let buf2 = pool.get();
        assert!(buf2.capacity() >= 1024);
    }
}