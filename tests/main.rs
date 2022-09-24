use std::fs;
use rust_seq2kminmers::KminmersIterator;

#[test]
fn main() {

    let genome_file = "tests/ecoli.genome.100k.fa";
    let contents = fs::read_to_string(genome_file).expect("Failed to read test genome file").split("\n").collect::<Vec<&str>>()[1].to_string();

    let hashes : Vec<u64> = vec![ 5238461524779056570,
     14895685229490448351,
     5103682387073363600,
     11643183865033904618,
     5256025917849251646,
     10448852938929078343,
     11987499807632838007,
     15567715264670518846,
     2580426438973794393,
     1838538700886040430,
     2570222892016362111,
     16127441378801417204,
     16009455758200025226,
     8042481517955484350,
     17023110711370061084,
     12295444363597936549,
     3365939716928708239,
     12820581993041321376,
     4080354370492127023,
     11125575058059394899
    ];

    let iter = KminmersIterator::new(contents.as_bytes(), 10, 5, 0.0001, false).unwrap();
    let mut count = 0;
    for kminmer in iter
    {
        println!("kminmer: {:?}",kminmer.get_hash());
        assert!(kminmer.get_hash() == hashes[count]);
        count += 1;
    }
}
