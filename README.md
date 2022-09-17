# rust-seq2kminmers

Rust crate to convert any DNA sequence to an iterator of its kminmers. It is useful for rust-mdbg.

* Uses NtHash1 behind the scenes. 

* Allows fast iteration in HPC-space but reports positions in original sequence space.

