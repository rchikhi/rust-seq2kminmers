#![feature(stdsimd)]
#![feature(core_intrinsics)]
use nthash32::NtHashIterator;
mod kminmer;
pub use kminmer::{Kminmer, KminmerVec, KminmerHash};
mod nthash_hpc;
pub use nthash_hpc::NtHashHPCIterator;
//mod nthash2_avx512_32;
//pub use nthash2_avx512_32::NtHashSIMDIterator;
mod nthash_avx512_32;
pub use nthash_avx512_32::NtHashSIMDIterator;
mod nthash_hpc_simd;
pub use nthash_hpc_simd::NtHashHPCSIMDIterator;
//mod nthash_c;
//pub use nthash_c::nthash_c;
mod hpc;
pub use hpc::{hpc, encode_rle, encode_rle_simd};

use std::io::{Result};

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum HashMode {
    Regular,
    Hpc,
    Simd,
    HpcSimd
}

// minimizer hash type
//pub type H  = u16; 
pub type H  = u32; // hash precision
//pub type H  = u64; 
//pub type FH = f32;
pub type FH = f64;

// kminmer hash type 
pub type KH = u64; 

pub type KminmerType = KminmerHash;

// An iterator for getting k-min-mers out of a DNA sequence
///
/// Parameters:
///     * l: minimizer length
///     * k: k-min-mer length
///     * density: density of minimizer scheme
///
/// mode: 
///     * Hpc or Hpcsimd:
///          minimizers are computed as if the sequence was homopolymer compressed but all positions reported
///          in original sequence space,
///     * Regular: everything is done in original sequence space
///
/// Hashing is performed by NtHash1
///
/// Example usage:
/// ```
///     use rust_seq2kminmers::{KminmersIterator, HashMode, Kminmer};
///     
///     fn main() {
///        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
///        let iter = KminmersIterator::new(seq, 10, 5, 0.1, HashMode::Hpc).unwrap();
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
    nthash_hpc_simd_iterator: Option<NtHashHPCSIMDIterator<'a>>,
    nthash_hpc_iterator: Option<NtHashHPCIterator<'a>>,
    nthash_simd_iterator: Option<NtHashSIMDIterator<'a>>,
    nthash_iterator: Option<NtHashIterator<'a>>,
    curr_sk : Vec::<KH>,
    curr_pos : Vec::<usize>,
    kminmer_fhash: KH,
    kminmer_rhash: KH,
    count : usize,
}


impl<'a> KminmersIterator<'a> {
    pub fn new(seq: &'a [u8], l: usize, k: usize, density: FH, mode: HashMode) -> Result<KminmersIterator<'a>> {

        let hash_bound = ((density as FH) * (H::MAX as FH)) as H;

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
                nthash_hpc_simd_iterator = Some(NtHashHPCSIMDIterator::new(seq, l, hash_bound).unwrap());
            }
            else { 
                nthash_iterator = Some(NtHashIterator::new(seq, l).unwrap());
            }
        }
        
        let curr_sk  = Vec::<KH>   ::with_capacity((seq.len() as FH * density) as usize);
        let curr_pos = Vec::<usize>::with_capacity((seq.len() as FH * density) as usize);

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


// the madness one needs to go through for function overloading
// (thanks for the help, https://stackoverflow.com/a/56100816)
trait MixHash{
    fn mixhash(&self) -> KH 
        where Self: Sized;
}

impl MixHash for u16
{
    fn mixhash(&self) -> KH 
    {
        let mut x :KH= *self as KH;
        // need a stronger hash this time (murmur64)
        x ^= x.rotate_left(33);
        x*= 0xff51afd7ed558ccd;
        x ^= x.rotate_left(33);
        x *= 0xc4ceb9fe1a85ec53;
        x ^= x.rotate_left(33);
        x
    }   
}

impl MixHash for u32
{
    fn mixhash(&self) -> KH 
    {
        let mut x :KH= *self as KH;
        // this is xorshift
        // but maybe use murmur here too (for less kminmer hashes collisions)
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    }   
}

impl MixHash for u64
{
    fn mixhash(&self) -> KH
    {
        *self
    }
}

impl<'a> Iterator for KminmersIterator<'a> {
    type Item = KminmerType;
    fn next(&mut self) -> Option<KminmerType> {
        let res;
        loop
        {
            let mut j;
            let mut hash: H;
            if self.mode == HashMode::HpcSimd {
                // the unwrap_or_else() magic is just an optimization:
                // https://www.reddit.com/r/rust/comments/dmws17/new_to_rust_is_unwrap_free/f55lalo/
                if self.nthash_hpc_simd_iterator.is_none() { return None; }
                match self.nthash_hpc_simd_iterator.as_mut().unwrap_or_else(|| unsafe { std::hint::unreachable_unchecked() }).next()
                {
                    Some(n) => { (j,hash) = n; } 
                    None => return None
                };
            }
            else if self.mode == HashMode::Simd {
                if self.nthash_simd_iterator.is_none() { return None; }
                match self.nthash_simd_iterator.as_mut().unwrap_or_else(|| unsafe { std::hint::unreachable_unchecked() }).next()
                {
                    Some(n) => { (j,hash) = n; }
                    None => return None
                };
            }
            else if self.mode == HashMode::Hpc
            {
                if self.nthash_hpc_iterator.is_none() { return None; }
                match self.nthash_hpc_iterator.as_mut().unwrap_or_else(|| unsafe { std::hint::unreachable_unchecked() }).next()
                {
                    Some(n) => { (j,hash) = n; } 
                    None => return None
                };
            }
            else 
            {
                loop
                {
                    if self.nthash_iterator.is_none() { return None; }
                    match self.nthash_iterator.as_mut().unwrap().next()
                    {
                        Some(x) => { hash = x as H;}
                        None => return None
                    };
                    j = self.seq_pos;
                    self.seq_pos += 1;
                    if hash <= self.hash_bound { break; }
                }
            }
            self.curr_pos.push(j); // raw sequence position
            let hash : KH = if self.mode == HashMode::Simd { (hash as u32).mixhash() } else { hash.mixhash() }; // only necessary if input hashes are u32
            self.curr_sk.push(hash);

            if self.curr_sk.len() >= self.k { 
                //assert!(self.curr_sk.len() == self.count+self.k) ;  // commented, but it's true
                //assert!(self.curr_pos[self.count+self.k - 1] == j); // same
                if self.curr_sk.len() == self.k
                {
                    self.kminmer_fhash ^= hash.rotate_left((self.k-1-(self.curr_sk.len()-1)) as u32);
                    self.kminmer_rhash ^= hash.rotate_left((self.curr_sk.len()-1) as u32);
                }
                else
                {
                    self.kminmer_fhash = self.kminmer_fhash.rotate_left(1)  ^ hash
                        ^ self.curr_sk[self.count-1].rotate_left(self.k as u32);
                    self.kminmer_rhash = self.kminmer_rhash.rotate_right(1) ^ hash.rotate_left((self.k-1) as u32)  
                        ^ self.curr_sk[self.count-1].rotate_right(1);
                }
                let nthash = if self.kminmer_fhash < self.kminmer_rhash { self.kminmer_fhash } else { self.kminmer_rhash };
                let rev = self.kminmer_rhash < self.kminmer_fhash;
                //let (hash,rev) = nthash1_minimizer_space(&self.curr_sk[self.count..self.count+self.k]);
                //res = Some(Kminmer::new(&self.curr_sk[self.count..self.count+self.k], self.curr_pos[self.count], j + self.l - 1, self.count));
                //res = Some(Kminmer::new(&self.curr_sk[0..0], 0, 0,0));

                // new_from_hash is currently only good when KH is u64 and not Simd: otherwise slows things down at the
                // Dashmap level due to many kminmer hash collision
                res = Some(KminmerHash::new_from_hash(nthash, self.curr_pos[self.count], j + self.l - 1, self.count, rev));
                self.count += 1;
                break; 
            }
            else
            {
                self.kminmer_fhash ^= hash.rotate_left((self.k-1-(self.curr_sk.len()-1)) as u32);
                self.kminmer_rhash ^= hash.rotate_left((self.curr_sk.len()-1) as u32);
            }
        }
        res
    }
}


// a test implementation to hash just a single k-min-mer, hence it is non-rolling 
#[allow(dead_code)]
fn nthash1_minimizer_space(kminmer: &[H]) -> (KH, bool)
{
    let mut kminmer_fhash :KH= 0;
    let mut kminmer_rhash :KH= 0;
    let k = kminmer.len();
    for (i,hash) in kminmer.iter().enumerate()
    {
        kminmer_fhash ^= (*hash as KH).rotate_left((k-1-i) as u32);
        kminmer_rhash ^= (*hash as KH).rotate_left(i     as u32);
    }
    let hash = if kminmer_fhash < kminmer_rhash { kminmer_fhash } else { kminmer_rhash } as KH;
    let rev = kminmer_rhash < kminmer_fhash;
    (hash, rev)
}
