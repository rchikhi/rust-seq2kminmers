// kminmer.rs
// Contains the "Kminmer" struct, a k-mer in minimizer-space.

use std::hash::{Hash, Hasher};
use std::vec::Vec;
use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;
use fxhash::{hash, hash32, hash64};

use crate::{H,KH};

pub trait Kminmer {
    fn new(mers: &[H], start: usize, end: usize, offset: usize) -> Self;
    fn get_hash(&self) -> KH;
}

#[derive(Clone, Debug)]
pub struct KminmerVec {
    mers: Vec<H>, // Raw Vec of minimizer hashes
    pub start: usize, // Start location
    pub end: usize, // End location
    pub offset: usize, // Offset (index in the k-min-mer array)
    pub rev: bool, // Strand direction
}


impl Kminmer for KminmerVec {
    // Create a new Kminmer object.
    fn new(mers: &[H], start: usize, end: usize, offset: usize) -> Self {
        let mut obj = KminmerVec {
            mers: Vec::from(mers),
            start,
            end,
            offset,
            rev: false,
        };
        obj.normalize();
        obj     
    }
    
    // Hash the Vec of minimizer hashes to a u64 (this is used throughout the reference processing).
    fn get_hash(&self) -> KH {
        let mut hash = DefaultHasher::new();
        self.mers.hash(&mut hash);
        hash.finish() as KH
    }

}

impl KminmerVec {
    
    // Transform into canonical Kminmer for this object.
    pub fn normalize(&mut self) {
        let mut rev_mers = self.mers.clone();
        rev_mers.reverse();
        if rev_mers < self.mers {
            self.mers = rev_mers;
            self.rev = true;
        }
    }

    // Make sure a kminmer is normalized 
    pub fn is_normalized(&self) -> bool {
        let mut rev_mers = self.mers.clone();
        rev_mers.reverse();
        self.mers <= rev_mers
    }
    
    // Pretty-print a Kminmer.
    pub fn print(&self) -> String {
	// prints only the first 2 digits of each minimizer hash
	let mut s = String::new();
	for x in self.mers.iter() {
	    s = format!("{}{} ", s, x.to_string()[..2].to_string());
	}
	s
    }
    
    // Obtain a raw Vec of minimizer hashes.
    pub fn mers(&self) -> Vec<H> {
        self.mers.to_vec()
    }

    pub fn get_hash_usize(&self) -> usize {
        hash(&self.mers)
    }
    pub fn get_hash_u32(&self) -> u32 {
        hash32(&self.mers)
    }
    pub fn get_hash_u64(&self) -> u64 {
        hash64(&self.mers)
    }
}

// Various impls for Kminmer.
impl PartialEq for KminmerVec {
    fn eq(&self, other: &KminmerVec) -> bool {
        self.mers == other.mers
    }
}

impl Eq for KminmerVec {
}

impl Hash for KminmerVec {
    fn hash<G: Hasher>(&self, state: &mut G) {
        self.mers.hash(state);
    }
}

impl Default for KminmerVec {
    fn default() -> Self{KminmerVec{mers: vec![], start: 0, end: 0, offset: 0, rev: false}}
}

impl Ord for KminmerVec {
    fn cmp(&self, other: &Self) -> Ordering {
        self.mers.cmp(&other.mers)
    }
}

impl PartialOrd for KminmerVec{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug)]
pub struct KminmerHash {
    pub hash: KH, // hash from Vec of minimizer hashes
    pub start: usize, // Start location
    pub end: usize, // End location
    pub offset: usize, // Offset (index in the k-min-mer array)
    pub rev: bool, // Strand direction
}


impl Kminmer for  KminmerHash {
    // Create a new Kminmer object.
    fn new(mers: &[H], start: usize, end: usize, offset: usize) -> Self {
        let hash: KH;
        let rev;
        let mut rev_mers = mers.to_vec();
        rev_mers.reverse();
        let mers = mers.to_vec();
        if rev_mers < mers {
            hash = hash64(&rev_mers) as KH;
            rev = true;
        }
        else {
            hash = hash64(&mers) as KH; 
            rev = false;
        }
        KminmerHash {
            hash,
            start,
            end,
            offset,
            rev,
        }    
    }
    fn get_hash(&self) -> KH {
        self.hash as KH 
    }
}

impl KminmerHash {
    // Create a new Kminmer object.
    pub fn new_from_hash(hash :KH, start: usize, end: usize, offset: usize, rev: bool) -> Self {
        KminmerHash {
            hash,
            start,
            end,
            offset,
            rev,
        }    
    }
}

// Various impls for Kminmer.
impl PartialEq for KminmerHash {
    fn eq(&self, other: &Self) -> bool {
        self.hash == other.hash
    }
}

impl Eq for KminmerHash {
}

impl Default for KminmerHash {
    fn default() -> Self{KminmerHash{hash: 0, start: 0, end: 0, offset: 0, rev: false}}
}

impl Ord for KminmerHash {
    fn cmp(&self, other: &Self) -> Ordering {
        self.hash.cmp(&other.hash)
    }
}

impl PartialOrd for KminmerHash {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

