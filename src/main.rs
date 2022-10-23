use std::time::Instant;
#[allow(unused_imports)]
use rand::distributions::{Distribution, Uniform};
use rust_seq2kminmers::{KminmersIterator, HashMode};
use rust_parallelfastx::{parallel_fastx};
#[allow(unused_imports)]
use rust_seq2kminmers::{encode_rle, hpc, encode_rle_simd};

fn main() {
    let args: Vec<String> = std::env::args().collect();

    //let mode = HashMode::Regular;
    //let mode = HashMode::Simd;
    //let mode = HashMode::Hpc;
    let mode = HashMode::HpcSimd;

    // A simple example given just a sequence in a string
    if args.len() < 2
    {
        let seq = "AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
        //let seq = "AACTTTTTGGGGGGCAAAAAACCCCCCCTGCCCCCCAAACTTTTTGGGGGGCAAAAAACCCCCCCTGCCCCCCAAACTTTTTGGGGGGCAAAAAACCCCCCCTGCCCCCCA";

        /*
        let range = Uniform::from(0..4);
        let mut rng = rand::thread_rng();
        let seq_len = 10000;
        let randseq = (0..seq_len)
            .map(|_| match range.sample(&mut rng) {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => 'N',
            })
        .collect::<String>();
        let seq = &randseq;
        */

        println!("seq:    {:?}",seq);
        println!("HPC:    {:?}",hpc(seq));
        // no point displaying that, correctness is tested in tests/main.src
        //println!("HPCopt: {:?}",encode_rle_simd(seq));
        //println!("encode_rle:{:?}",encode_rle(seq));
        println!("Demonstrating how to construct k-min-mers (l=31, k=5, d=0.1) out of a test sequence: {}",seq);
        let iter = KminmersIterator::new(seq.as_bytes(), 31, 5, 0.1, mode).unwrap();
        for kminmer in iter
        {
            println!("kminmer: {:?}",kminmer);
        }
    }

    // A complete example with FASTA parsing
    else
    {
        let k = 5;
        let l = 31;
        let d = 0.01;
        let filename = &args[1];
        let nb_threads = args[2].parse::<i32>().unwrap();
        println!("Enumerating k-min-mers for the input file {} in parallel ({} threads)", filename, nb_threads);


        let task = |seq_str: &[u8], _seq_id: &str|  {
            let iter = KminmersIterator::new(seq_str, l, k, d, mode).unwrap();
            //let iter :Vec<u64> = vec![];
            let mut _count = 0;
            for _kminmer in iter
            {
                //println!("seq_id {} kminmer: {:?}",_seq_id,kminmer);
                _count += 1;
            }
            //println!("seq_id {} nb kminmers: {}",_seq_id,_count);

        };

        let start = Instant::now();
        parallel_fastx(&filename, nb_threads.try_into().unwrap(), task);
        let duration = start.elapsed();
        println!("FASTA to kminmers in {:?}.", duration);

    }
}
