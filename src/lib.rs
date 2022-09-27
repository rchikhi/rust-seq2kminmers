use nthash::NtHashIterator;
mod kminmer;
pub use kminmer::{Kminmer, KminmerHash};
mod nthash_hpc;
pub use nthash_hpc::NtHashHPCIterator;
use std::io::{Result, Write, Lines, BufReader, BufRead};
use std::fs::{File, OpenOptions};
// An iterator for getting k-min-mers out of a DNA sequence
///
/// Parameters:
/// l: minimizer length
/// k: k-min-mer length
/// density: density of minimizer scheme
/// hpc: true if minimizers are computed as if the sequence was homopolymer compressed but all positions reported
///      in original sequence space, false if everything is done in original sequence space
///
/// Hashing is performed by NtHash1
///
/// Example usage:
/// ```
///     use rust_seq2kminmers::KminmersIterator;
///     
///     fn main() {
///        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
///        let iter = KminmersIterator::new(seq, 10, 5, 0.1, false).unwrap();
///        for kminmer in iter
///        {
///            println!("kminmer: {:?}",kminmer);
///        }
///     }
/// ```

pub struct KminmersIterator<'a> {
    seq_pos : usize, 
    k: usize,
    l: usize,
    hash_bound: u64,
    hpc: bool,
    nthash_hpc_iterator: Option<NtHashHPCIterator<'a>>,
    nthash_iterator: Option<NtHashIterator<'a>>,
    curr_sk : Vec::<u64>,
    curr_pos : Vec::<usize>,
    count : usize,
}

impl<'a> KminmersIterator<'a> {
    pub fn new(seq: &'a [u8], l: usize, k: usize, density: f64, hpc: bool) -> Result<KminmersIterator<'a>> {

        let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;

        let mut nthash_hpc_iterator = None;
        let mut nthash_iterator = None;
        if seq.len() > l {
            if hpc
            {
                nthash_hpc_iterator = Some(NtHashHPCIterator::new(seq, l, hash_bound).unwrap());
            }
            else
            { 
                nthash_iterator = Some(NtHashIterator::new(seq, l).unwrap());
            }
        }

        let curr_sk = Vec::<u64>::new();
        let curr_pos = Vec::<usize>::new();

        Ok(KminmersIterator {
            seq_pos: 0,
            k,
            l,
            hash_bound,
            hpc,
            nthash_hpc_iterator: nthash_hpc_iterator,
            nthash_iterator: nthash_iterator,
            curr_pos,
            curr_sk,
            count : 0
        })
    }
}

impl<'a> Iterator for KminmersIterator<'a> {
    type Item = Kminmer;

    fn next(&mut self) -> Option<Kminmer> {
        let kminmer;
        loop
        {
            let mut j;
            let mut hash;
            if self.hpc
            {
                match self.nthash_hpc_iterator.as_mut().unwrap().next()
                {
                    Some(n) => { (j,hash) = n; } 
                    None => return None
                };
            }
            else
            {
                loop
                {
                    match self.nthash_iterator.as_mut().unwrap().next()
                    {
                        Some(x) => { hash = x;}
                        None => return None
                    };
                    self.seq_pos += 1;
                    j = self.seq_pos;
                    if hash < self.hash_bound { break; }
                }
            }

            self.curr_pos.push(j); // raw sequence position
            self.curr_sk.push(hash);
            if self.curr_sk.len() == self.k { 
                kminmer = Kminmer::new(&self.curr_sk, self.curr_pos[0], self.curr_pos[self.k - 1] + self.l - 1, self.count);
                self.curr_sk = self.curr_sk[1..self.k].to_vec();
                self.curr_pos = self.curr_pos[1..self.k].to_vec();
                self.count += 1;
                break; 
            }
        }
        Some(kminmer)
    }
}


pub struct KminmersWriteIterator<'a> {
    seq_pos : usize, 
    k: usize,
    l: usize,
    hash_bound: u64,
    hpc: bool,
    nthash_hpc_iterator: Option<NtHashHPCIterator<'a>>,
    nthash_iterator: Option<NtHashIterator<'a>>,
    curr_sk : Vec::<u64>,
    curr_pos : Vec::<usize>,
    count : usize,
    file: File,
}

impl<'a> KminmersWriteIterator<'a> {
    pub fn new(seq: &'a [u8], l: usize, k: usize, density: f64, hpc: bool, prefix: &str) -> Result<KminmersWriteIterator<'a>> {
        let buf = format!("{}-{}-{}.mers", prefix, l, density);
        let file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(buf)
        .unwrap();
        let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;

        let mut nthash_hpc_iterator = None;
        let mut nthash_iterator = None;
        if seq.len() > l {
            if hpc
            {
                nthash_hpc_iterator = Some(NtHashHPCIterator::new(seq, l, hash_bound).unwrap());
            }
            else
            { 
                nthash_iterator = Some(NtHashIterator::new(seq, l).unwrap());
            }
        }

        let curr_sk = Vec::<u64>::new();
        let curr_pos = Vec::<usize>::new();

        Ok(KminmersWriteIterator {
            seq_pos: 0,
            k,
            l,
            hash_bound,
            hpc,
            nthash_hpc_iterator: nthash_hpc_iterator,
            nthash_iterator: nthash_iterator,
            curr_pos,
            curr_sk,
            count : 0,
            file,
        })
    }
}

impl<'a> Iterator for KminmersWriteIterator<'a> {
    type Item = Kminmer;

    fn next(&mut self) -> Option<Kminmer> {
        let kminmer;
        loop
        {
            let mut j;
            let mut hash;
            if self.hpc
            {
                match self.nthash_hpc_iterator.as_mut().unwrap().next()
                {
                    Some(n) => { (j,hash) = n; } 
                    None => return None
                };
            }
            else
            {
                loop
                {
                    match self.nthash_iterator.as_mut().unwrap().next()
                    {
                        Some(x) => { hash = x;}
                        None => return None
                    };
                    self.seq_pos += 1;
                    j = self.seq_pos;
                    if hash < self.hash_bound { break; }
                }
            }
            write!(self.file, "{}\t{}\n", j, hash)
            .expect("Unable to write minimizer.");
            self.curr_pos.push(j); // raw sequence position
            self.curr_sk.push(hash);
            if self.curr_sk.len() == self.k { 
                kminmer = Kminmer::new(&self.curr_sk, self.curr_pos[0], self.curr_pos[self.k - 1] + self.l - 1, self.count);
                self.curr_sk = self.curr_sk[1..self.k].to_vec();
                self.curr_pos = self.curr_pos[1..self.k].to_vec();
                self.count += 1;
                break; 
            }
        }
        Some(kminmer)
    }
}

pub struct KminmersReadIterator {
    k: usize,
    l: usize,
    line_iterator: Lines<BufReader<File>>,
    curr_sk : Vec::<u64>,
    curr_pos : Vec::<usize>,
    count : usize,
}

impl KminmersReadIterator {
    pub fn new(l: usize, k: usize, density: f64, prefix: &str) -> Result<KminmersReadIterator> {
        let buf = format!("{}-{}-{}.mers", prefix, l, density);
        let file = OpenOptions::new().read(true).open(buf).expect("Could not open minimizer index.");
        let lines = BufReader::new(file).lines();
        let curr_sk = Vec::<u64>::new();
        let curr_pos = Vec::<usize>::new();
        Ok(KminmersReadIterator {
            k,
            l,
            line_iterator: lines,
            curr_pos,
            curr_sk,
            count : 0,
        })
    }
}

impl Iterator for KminmersReadIterator {
    type Item = Kminmer;

    fn next(&mut self) -> Option<Kminmer> {
        let kminmer;
        loop
        {
            let mut j;
            let mut hash;
            let line = self.line_iterator.next();
            if let Some(Ok(l)) = line {
                let v : Vec<&str> = l.split("\t").collect();
                j = v[0].parse::<usize>().unwrap();
                hash = v[1].parse::<u64>().unwrap();

            }
            else {return None;}
            self.curr_pos.push(j); // raw sequence position
            self.curr_sk.push(hash);
            if self.curr_sk.len() == self.k { 
                kminmer = Kminmer::new(&self.curr_sk, self.curr_pos[0], self.curr_pos[self.k - 1] + self.l - 1, self.count);
                self.curr_sk = self.curr_sk[1..self.k].to_vec();
                self.curr_pos = self.curr_pos[1..self.k].to_vec();
                self.count += 1;
                break; 
            }
        }
        Some(kminmer)
    }
}

pub struct KminmersHashIterator<'a> {
    seq_pos : usize, 
    k: usize,
    l: usize,
    hash_bound: u64,
    hpc: bool,
    nthash_hpc_iterator: Option<NtHashHPCIterator<'a>>,
    nthash_iterator: Option<NtHashIterator<'a>>,
    curr_sk : Vec::<u64>,
    curr_pos : Vec::<usize>,
    count : usize,
}

impl<'a> KminmersHashIterator<'a> {
    pub fn new(seq: &'a [u8], l: usize, k: usize, density: f64, hpc: bool) -> Result<KminmersHashIterator<'a>> {

        let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;

        let mut nthash_hpc_iterator = None;
        let mut nthash_iterator = None;
        if seq.len() > l {
            if hpc
            {
                nthash_hpc_iterator = Some(NtHashHPCIterator::new(seq, l, hash_bound).unwrap());
            }
            else
            { 
                nthash_iterator = Some(NtHashIterator::new(seq, l).unwrap());
            }
        }

        let curr_sk = Vec::<u64>::new();
        let curr_pos = Vec::<usize>::new();

        Ok(KminmersHashIterator {
            seq_pos: 0,
            k,
            l,
            hash_bound,
            hpc,
            nthash_hpc_iterator: nthash_hpc_iterator,
            nthash_iterator: nthash_iterator,
            curr_pos,
            curr_sk,
            count : 0
        })
    }
}

impl<'a> Iterator for KminmersHashIterator<'a> {
    type Item = KminmerHash;

    fn next(&mut self) -> Option<KminmerHash> {
        let kminmer;
        loop
        {
            let mut j;
            let mut hash;
            if self.hpc
            {
                match self.nthash_hpc_iterator.as_mut().unwrap().next()
                {
                    Some(n) => { (j,hash) = n; } 
                    None => return None
                };
            }
            else
            {
                loop
                {
                    match self.nthash_iterator.as_mut().unwrap().next()
                    {
                        Some(x) => { hash = x;}
                        None => return None
                    };
                    self.seq_pos += 1;
                    j = self.seq_pos;
                    if hash < self.hash_bound { break; }
                }
            }

            self.curr_pos.push(j); // raw sequence position
            self.curr_sk.push(hash);
            if self.curr_sk.len() == self.k { 
                kminmer = KminmerHash::new(&self.curr_sk, self.curr_pos[0], self.curr_pos[self.k - 1] + self.l - 1, self.count);
                self.curr_sk = self.curr_sk[1..self.k].to_vec();
                self.curr_pos = self.curr_pos[1..self.k].to_vec();
                self.count += 1;
                break; 
            }
        }
        Some(kminmer)
    }
}
