
use rust_seq2kminmers::{KminmersIterator, KminmersWriteIterator, KminmersReadIterator};

fn main() {
    let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
    let iter = KminmersIterator::new(seq, 10, 5, 0.1, true).unwrap();
    for kminmer in iter
    {
        println!("kminmer: {:?}",kminmer);
    }
    println!("Write minimizers while generating k-min-mers:");
    let iter = KminmersWriteIterator::new(seq, 10, 5, 0.1, false, "seq").unwrap();
    for kminmer in iter
    {
        println!("kminmer-write: {:?}",kminmer);
    }
    println!("Read minimizers and generate k-min-mers:");
    let iter2 = KminmersReadIterator::new(10, 4, 0.1, "seq").unwrap();
    for kminmer in iter2
    {
        println!("kminmer-read: {:?}",kminmer);
    }
}

