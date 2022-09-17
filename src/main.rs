
use rust_seq2kminmers::KminmersIterator;

fn main() {
    let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
    let iter = KminmersIterator::new(seq, 10, 5, 0.1, false).unwrap();
    for kminmer in iter
    {
        println!("kminmer: {:?}",kminmer);
    }
}

