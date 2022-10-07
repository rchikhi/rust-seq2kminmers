#![feature(bench_black_box)]
use std::time::Instant;
use rust_seq2kminmers::{KminmersIterator};
use rust_parallelfastx::{parallel_fastx};

fn main() {
    let args: Vec<String> = std::env::args().collect();

    
    // A simple example given just a sequence in a string
    if args.len() < 2
    {
        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
        println!("Demonstrating how to construct k-min-mers (k=10, l=5, d=0.1) out of a test sequence: {}",std::str::from_utf8(seq).unwrap());
        let iter = KminmersIterator::new(seq, 10, 5, 0.1, true).unwrap();
        for kminmer in iter
        {
            println!("kminmer: {:?}",kminmer);
        }
    }


    // A complete example with FASTA parsing
    else
    {
        let k = 7;
        let l = 31;
        let d = 0.01;
        let filename = &args[1];
        let nb_threads = args[2].parse::<i32>().unwrap();
        println!("Enumerating k-min-mers for the input file {} in parallel ({} threads)", filename, nb_threads);


        let task = |seq_str: &[u8], _seq_id: &str|  {
            let iter = KminmersIterator::new(seq_str, k, l, d, true).unwrap();
            //let iter :Vec<u64> = vec![];
            for kminmer in iter
            {
                //println!("seq_id {} kminmer: {:?}",seq_id,kminmer);
                
                // prevent from being optimized away
                std::hint::black_box(&kminmer);
            }
        };

        let start = Instant::now();
        parallel_fastx(&filename, nb_threads.try_into().unwrap(), task);
        let duration = start.elapsed();
        println!("FASTA to kminmers in {:?}.", duration);

    }
}
