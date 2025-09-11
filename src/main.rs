use clap::Parser;

// Module declarations
mod bx_index;
mod locus;
mod read_pair_cache;
mod types;
mod utils;
mod cli;
mod formatter;
mod writer;
mod processor;
mod memory_pool;
mod performance_test;
mod advanced_cache;
mod parallel_processor;
mod integration_tests;

// Re-exports
use cli::Args;
use processor::go;
use utils::set_panic_handler;

fn main() {
    set_panic_handler();
    std::env::set_var("RUST_BACKTRACE", "1");

    // Parse command line arguments using clap
    let args = Args::parse();

    let traceback = args.traceback;
    let res = go(args, None);

    if let Err(ref e) = res {
        println!("bam2fastq error: {e}\n");

        if traceback {
            println!("see below for more details:");
            println!("==========================");
            println!("{}\n{}", e, e.backtrace());
        };
        ::std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_lr21() {
        // Create fixed output directory
        let output_dir = "target/fastq_results";

        let args = Args {
            threads: 10,
            bam: "/Users/lishuangshuang/Documents/scrna/dnbc4tools/target/my_test_3/pos_sortednon_multiplexed.bam".to_string(),
            outputpath: output_dir.to_string(),
            reads_per_fastq: None,
            locus: None,
            bx_list: None,
            traceback: false,
            relaxed: false,
        };

        // Run conversion
        let out_path_sets = super::go(args, None).unwrap();
        
        // Print resulting FASTQ file paths
        println!("\nGenerated FASTQ files:");
        for (r1, r2, i1, i2) in out_path_sets {
            println!("R1: {}", r1.display());
            println!("R2: {}", r2.display());
            if let Some(i1_path) = i1 {
                println!("I1: {}", i1_path.display());
            }
            if let Some(i2_path) = i2 {
                println!("I2: {}", i2_path.display());
            }
            println!("---");
        }

        // Optional: Check files are generated and print file sizes
        println!("\nFile size information:");
        for entry in std::fs::read_dir(output_dir).unwrap() {
            let entry = entry.unwrap();
            let metadata = entry.metadata().unwrap();
            println!("{}: {} bytes", entry.file_name().to_string_lossy(), metadata.len());
        }
    }

    #[test]
    fn bad_bam() {
        let tempdir = tempfile::Builder::new()
            .prefix("bam_to_fq_test")
            .tempdir()
            .expect("create temp dir");
        let tmp_path = tempdir.path().join("outs");

        let args = Args {
            threads: 2,
            bam: "/Users/lishuangshuang/Documents/scrna/dnbc4tools/target/my_test_3/pos_sortednon_multiplexed.bam".to_string(),
            outputpath: tmp_path.to_str().unwrap().to_string(),
            reads_per_fastq: None,
            locus: None,
            bx_list: None,
            traceback: false,
            relaxed: false,
        };

        let res = super::go(args, Some(2));

        println!("res: {:?}", res);
    }
}