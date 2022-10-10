#![feature(stdsimd)]
#![feature(core_intrinsics)]
use nthash::NtHashIterator;
mod kminmer;
pub use kminmer::{Kminmer, KminmerVec, KminmerHash};
mod nthash_hpc;
pub use nthash_hpc::NtHashHPCIterator;
mod nthash_simd;
pub use nthash_simd::NtHashSIMDIterator;
mod nthash_hpc_simd;
pub use nthash_hpc_simd::NtHashHPCSIMDIterator;
mod nthash_c;
pub use nthash_c::nthash_c;
mod hpc;
pub use hpc::{hpc, encode_rle, encode_rle_simd};

use std::io::{Result};

#[derive(PartialEq, Clone, Copy)]
pub enum HashMode {
    Regular,
    Hpc,
    Simd,
    HpcSimd
}

//pub type H  = u32; // hash precision
//pub type FH = f32;
pub type H  = u64; 
pub type FH = f64;

// An iterator for getting k-min-mers out of a DNA sequence
///
/// Parameters:
/// l: minimizer length
/// k: k-min-mer length
/// density: density of minimizer scheme
/// mode: 
///     Hpc or Hpcsimd:
///     minimizers are computed as if the sequence was homopolymer compressed but all positions reported
///      in original sequence space,
///     Regular: everything is done in original sequence space
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
    hash_bound: H,
    mode: HashMode,
    nthash_hpc_simd_iterator: Option<NtHashHPCSIMDIterator>,
    nthash_hpc_iterator: Option<NtHashHPCIterator<'a>>,
    nthash_simd_iterator: Option<NtHashSIMDIterator<'a>>,
    nthash_iterator: Option<NtHashIterator<'a>>,
    curr_sk : Vec::<H>,
    curr_pos : Vec::<usize>,
    kminmer_fhash: H,
    kminmer_rhash: H,
    count : usize,
}


impl<'a> KminmersIterator<'a> {
    pub fn new(seq: &'a [u8], l: usize, k: usize, density: FH, mode: HashMode) -> Result<KminmersIterator<'a>> {

        let hash_bound = ((density as FH) * (H::max_value() as FH)) as H;

        let mut nthash_hpc_simd_iterator = None;
        let mut nthash_hpc_iterator = None;
        let mut nthash_simd_iterator = None;
        let mut nthash_iterator = None;
        if seq.len() > l {
            if mode == HashMode::Hpc {
                nthash_hpc_iterator = Some(NtHashHPCIterator::new(seq, l, hash_bound).unwrap());
            }
            else if mode == HashMode::Simd {
                nthash_simd_iterator = Some(NtHashSIMDIterator::new(seq, l, hash_bound));
            }
            else if mode == HashMode::HpcSimd {
                nthash_hpc_simd_iterator = Some(NtHashHPCSIMDIterator::new(seq, l, hash_bound));
            }
            else { 
                nthash_iterator = Some(NtHashIterator::new(seq, l).unwrap());
            }
        }
        
        let curr_sk  = Vec::<H>    ::with_capacity((seq.len() as f64 * density) as usize);
        let curr_pos = Vec::<usize>::with_capacity((seq.len() as f64 * density) as usize);

        Ok(KminmersIterator {
            seq_pos: 0,
            k,
            l,
            hash_bound,
            mode,
            nthash_hpc_simd_iterator,
            nthash_hpc_iterator,
            nthash_simd_iterator,
            nthash_iterator,
            curr_pos,
            curr_sk,
            kminmer_fhash: 0,
            kminmer_rhash: 0,
            count : 0
        })
    }
}

pub type KminmerType = KminmerHash;

impl<'a> Iterator for KminmersIterator<'a> {
    type Item = KminmerType;
    fn next(&mut self) -> Option<KminmerType> {
        let res;
        loop
        {
            let mut j;
            let mut hash: H;
            if self.mode == HashMode::HpcSimd {
                match self.nthash_hpc_simd_iterator.as_mut().unwrap().next()
                {
                    Some(n) => { (j,hash) = n; } 
                    None => return None
                };
            }
            else if self.mode == HashMode::Simd {
                match self.nthash_simd_iterator.as_mut().unwrap().next()
                {
                    Some(n) => { (j,hash) = n; } 
                    None => return None
                };
            }
            else if self.mode == HashMode::Hpc
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
                        Some(x) => { hash = x as H;}
                        None => return None
                    };
                    self.seq_pos += 1;
                    j = self.seq_pos;
                    if hash < self.hash_bound { break; }
                }
            }

            self.curr_pos.push(j); // raw sequence position
            self.curr_sk.push(hash);
            if self.curr_sk.len() >= self.k { 
                self.kminmer_fhash ^= self.kminmer_fhash.rotate_left(1) ^ (hash as H).rotate_left(self.k as u32) ^ (self.curr_pos[self.count] as H);
                self.kminmer_rhash ^= self.kminmer_rhash.rotate_right(1) ^ (hash as H) .rotate_right(1) ^ (self.curr_pos[self.count] as H).rotate_left((self.k-1) as u32);
                let hash = if self.kminmer_fhash < self.kminmer_rhash { self.kminmer_fhash } else { self.kminmer_rhash } as u32;
                //res = Some(Kminmer::new(&self.curr_sk[self.count..self.count+self.k], self.curr_pos[self.count], self.curr_pos[self.count+self.k - 1] + self.l - 1, self.count));
                res = Some(KminmerHash::new_from_hash(hash, self.curr_pos[self.count], self.curr_pos[self.count+self.k - 1] + self.l - 1, self.count, self.kminmer_rhash < self.kminmer_fhash));
                //res= Some(Kminmer::new(&self.curr_sk[0..0], 0,0,0));
                self.count += 1;
                break; 
            }
            else
            {
                self.kminmer_fhash ^= (hash as H).rotate_left((self.k-self.curr_sk.len() + 1) as u32);
                self.kminmer_rhash ^= (hash as H).rotate_left(self.curr_sk.len() as u32)
            }
        }
        res
    }
}


