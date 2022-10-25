use std::fs;
use rust_seq2kminmers::{KminmersIterator, HashMode, Kminmer, KminmerHash, H};
use rust_seq2kminmers::{hpc, encode_rle, encode_rle_simd};
use std::any::type_name;

fn type_of<T>(_: T) -> &'static str {
        type_name::<T>()
}

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

    let hashes32 : Vec<u64> = vec! [
        143479479014703,
        1415094313937202,
        7085699921625713,
        2731023262850893,
        3529660833839258,
        2520689800435504,
        3515165585325381,
        2855190423625803,
        5122855536061684,
        244022361441902,
        2856446528761135,
        906939906227534,
        2115341643533671,
        246274980452770,
        159737436030657,
    ];


    let iter = KminmersIterator::new(contents.as_bytes(), 10, 5, 0.0001, HashMode::Regular).unwrap();
    let mut count = 0;
    let dummy_hash: H = 0;
    for kminmer in iter
    {
        println!("kminmer: {:?}",kminmer.get_hash());
        if type_of(dummy_hash) == "u64" {
            assert!(kminmer.get_hash() == hashes[count]);
        } else {
            assert!(type_of(dummy_hash) == "u32");
            assert!(kminmer.get_hash() == hashes32[count]);
        }
        count += 1;
    }

    // a big HPC SIMD correctness test
    assert!(encode_rle(&contents).0 == hpc(&contents));
    assert!(encode_rle(&contents).0 == encode_rle_simd(&contents.as_bytes()).0);
    assert!(encode_rle(&contents).1 == encode_rle_simd(&contents.as_bytes()).1.into_iter().map(|p| p as usize).collect::<Vec<usize>>());


    // a stresstest over some parameters
    for l in vec![5,7,11,17,25,31] {
        for k in vec![2,5,8] {
            assert!(KminmersIterator::new(contents.as_bytes(), l, k, 0.01, HashMode::Regular).unwrap().collect::<Vec<KminmerHash>>() == 
                    KminmersIterator::new(contents.as_bytes(), l, k, 0.01, HashMode::Simd).unwrap().collect::<Vec<KminmerHash>>());
            assert!(KminmersIterator::new(contents.as_bytes(), l, k, 0.01, HashMode::Hpc).unwrap().collect::<Vec<KminmerHash>>() == 
                    KminmersIterator::new(contents.as_bytes(), l, k, 0.01, HashMode::HpcSimd).unwrap().collect::<Vec<KminmerHash>>());
        }
    }
}
