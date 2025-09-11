//! Utility functions for BAM to FASTQ conversion

use std::backtrace;
use std::panic;
use anyhow::{Error, anyhow};

/// Complement a nucleotide base
pub fn complement(b: u8) -> Result<u8, Error> {
    match b {
        b'A' => Ok(b'T'),
        b'T' => Ok(b'A'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'N' => Ok(b'N'),
        _ => Err(anyhow!("invalid nucleotide: {}", b as char)),
    }
}

/// Set up panic handler for better error reporting
pub fn set_panic_handler() {
    panic::set_hook(Box::new(|panic_info| {
        let backtrace = backtrace::Backtrace::capture();
        
        eprintln!("\n=== PANIC OCCURRED ===");
        
        if let Some(location) = panic_info.location() {
            eprintln!("Location: {}:{}", location.file(), location.line());
        }
        
        if let Some(message) = panic_info.payload().downcast_ref::<&str>() {
            eprintln!("Message: {}", message);
        } else if let Some(message) = panic_info.payload().downcast_ref::<String>() {
            eprintln!("Message: {}", message);
        }
        
        eprintln!("\nBacktrace:\n{}", backtrace);
        eprintln!("=====================\n");
    }));
}