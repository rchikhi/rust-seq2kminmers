use std::fs;
use rust_seq2kminmers::{KminmersIterator, HashMode, Kminmer};
use rust_seq2kminmers::{hpc, encode_rle, encode_rle_simd};

// these tests were made for H=u64
// if they fail when H=u32, that's normal
#[test]
fn main() {

    let genome_file = "tests/ecoli.genome.100k.fa";
    let contents = fs::read_to_string(genome_file).expect("Failed to read test genome file").split("\n").collect::<Vec<&str>>()[1].to_string();

    let hashes : Vec<u64> = vec![ 
         6097375827354318,
         5077268723048817,
         17093614815813553,
         13932651659877218,
         2254626575123847,
         4725847317728813,
         10971942364167709,
         1406844240705087,
         15284878278949327,
         13429516156719180,
         10760699289819902,
         11244197813995113,
         6993910349997344,
         22098843726082404,
         4944933674400292,
         14212811059278321,
         9310664830401458,
         11232758307960192,
         9720472733789719,
         13210101786532125,
    ];

    let iter = KminmersIterator::new(contents.as_bytes(), 10, 5, 0.0001, HashMode::Regular).unwrap();
    let mut count = 0;
    for kminmer in iter
    {
        println!("kminmer: {:?}",kminmer.get_hash());
        assert!(kminmer.get_hash() == hashes[count]);
        count += 1;
    }

    // a big HPC SIMD correctness test
    assert!(encode_rle(&contents).0 == hpc(&contents));
    assert!(encode_rle(&contents).0 == encode_rle_simd(&contents).0);
    
    assert!(encode_rle(&contents).1 == encode_rle_simd(&contents).1.into_iter().map(|p| p as usize).collect::<Vec<usize>>());
}
