
use rust_seq2kminmers::{KminmersIterator, KminmersHashIterator};

fn main() {
    let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
    let iter = KminmersIterator::new(seq, 10, 5, 0.1, true).unwrap();
    for kminmer in iter
    {
        println!("kminmer: {:?}",kminmer);
    }

    let iter = KminmersHashIterator::new(seq, 10, 5, 0.1, true).unwrap();
    for kminmer in iter
    {
        println!("kminmer: {:?}",kminmer);
    }
}

