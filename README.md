# rust-seq2kminmers

Rust crate to convert any DNA sequence to an iterator of its k-min-mers. It is useful for rust-mdbg.

# Intro

For non-bioinformaticians, this library essentially performs a variant of Run Length Encoding + a fast rolling hash of DNA sequences.

For a primer on k-min-mers, see: https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00332-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS240547122100332X%3Fshowall%3Dtrue

# Features

* Uses NtHash1 behind the scenes. 

* Allows fast iteration in HPC-space but reports positions in original sequence space.

* SIMD implementation for Homopolymer Compression 

* SIMD implementation for 32-bit NtHash1 and some 31-bit variant of NtHash2 (modified from https://github.com/bcgsc/ntHash/pull/9)

# Performance

Very good, thanks for asking. Scalar ntHash is around 100-200 MB/sec, SIMD ntHash around 1 GB/sec, HPC ntHash around 2-3 GB/sec

# Testing

    cargo run # to run a small test in main.rs

    cargo run -- [filename.fasta] [nb_threads] # to run a bigger test using an input file (also main.rs)

    cargo test # run some tests 

    cargo bench # run some throughput benchmarks
