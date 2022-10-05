# rust-seq2kminmers

Rust crate to convert any DNA sequence to an iterator of its kminmers. It is useful for rust-mdbg.

* Uses NtHash1 behind the scenes. 

* Allows fast iteration in HPC-space but reports positions in original sequence space.

# Testing

cargo run # to run a small test in main.rs

cargo run -- [filename.fasta] [nb_threads] # to run a bigger test using an input file (also main.rs)

cargo test # run some tests 

cargo bench # run some throughput benchmarks
