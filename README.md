# bam2fastq

A high-performance BAM to FASTQ converter for **C4 single-cell RNA-seq data** (DNBelab C4 / DNBC4tools). Designed for speed with parallel processing pipelines, it converts coordinate-sorted or name-sorted BAM files into compressed FASTQ outputs.

---

## Features

- 🚀 **Parallel processing** — Rayon-based batch formatting + gzp parallel gzip compression
- 🔄 **Paired-end support** — RpCache-based read pair matching with orphan handling
- 🧵 **Full thread control** — `--threads` correctly bounds Rayon, htslib, and gzp
- 🔍 **Auto-detection** — Automatically detects single-end vs. paired-end BAM
- 🗺️ **Region filtering** — Process only a specific genomic locus
- 📦 **chunked FASTQ** — Optional splitting into fixed-size FASTQ files

---

## Installation

### Requirements

- Rust ≥ 1.70
- `htslib` (for BAM reading; usually installed via the `rust-htslib` crate automatically)

### Build

```bash
git clone https://github.com/lishuangshuang0616/bam2fastq.git
cd bam2fastq
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

The binary will be at `target/release/bam2fastq`.

---

## Usage

```
bam2fastq [OPTIONS] <BAM> <OUTPUT>
```

### Positional arguments

| Argument | Description |
|----------|-------------|
| `<BAM>`  | Path to the input BAM file |
| `<OUTPUT>` | Output directory for FASTQ files |

### Options

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--threads` | `-t` | `auto` | Number of CPU threads for parallel processing (default: all available cores) |
| `--locus` | `-r` | — | Process a specific region, e.g. `chr1:1000-2000` |
| `--reads-per-fastq` | `-n` | — | Maximum reads per output FASTQ file |
| `--max-memory` | — | auto | Maximum memory in MB (auto-detected if omitted) |
| `--no-compress` | — | off | Disable gzip compression (outputs `.fastq` instead of `.fastq.gz`) |

### Examples

**Basic conversion (paired-end, auto-detected):**
```bash
bam2fastq -t 16 sample.bam ./output/
```

**Specific locus (requires BAM index):**
```bash
bam2fastq -t 8 -r chr1:1000000-2000000 sample.bam ./output_chr1/
```

**Uncompressed output:**
```bash
bam2fastq -t 8 --no-compress sample.bam ./output/
```

**Split into 500M-read chunks:**
```bash
bam2fastq -t 16 -n 500000000 sample.bam ./output/
```

---

## Output

Files are written to `<OUTPUT>/` using the naming convention:

```
<sample>_L<lane>_<chunk>_1.fastq.gz   # R1
<sample>_L<lane>_<chunk>_2.fastq.gz   # R2
<sample>_L<lane>_<chunk>_I1.fastq.gz  # Index 1 (if present)
<sample>_L<lane>_<chunk>_I2.fastq.gz  # Index 2 (if present)
```

---

## Architecture

The codebase is organized into focused modules:

| Module | Responsibility |
|--------|----------------|
| `formatter.rs` | BAM record → `FqRecord` conversion, complement LUT, tag parsing |
| `writer.rs` | `GenWriter` (gzip/raw), `FastqWriter` (file handles), `FastqManager` (dispatch by RG) |
| `processing.rs` | `proc_single_ended` and `proc_double_ended` pipeline implementations |
| `main.rs` | CLI args (`Args`), `go()` / `inner()` orchestration, entry point |

### Processing Pipeline

**Single-ended:**
```
BAM records → [batch] → Rayon par_iter (format) → channel → write thread → disk
```

**Paired-end:**
```
Phase 1: BAM records → RpCache (R1+R2 matching) → matched pairs → disk
                     → orphans → shardio temp file
Phase 2: shardio shard → par_iter (sort/re-pair) → channel → write thread → disk
```

### Thread Budget

With `--threads N`:
- **Rayon** (format work): `N` threads — the CPU-intensive part gets the full budget
- **BAM BGZF decompression**: `max(1, N/4)` threads — IO-bound, diminishing returns
- **gzp compression**: `max(1, N/4)` threads — IO-bound, diminishing returns

Total active threads is approximately `N × 1.5` rather than `N × 3`.

---

## Performance Tips

| Storage Type | Recommended Settings |
|---|---|
| Network/HDD (≤150 MB/s) | `--threads 16` (default batch/buffer settings are tuned for this) |
| NVMe SSD (≥500 MB/s) | `--threads 20+`, increase `BATCH_SIZE` in `processing.rs` to 10000-20000 |

---

## License

MIT
