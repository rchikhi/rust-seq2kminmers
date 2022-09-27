// kminmer.rs
// Contains the "Kminmer" struct, a k-mer in minimizer-space.

use std::hash::{Hash, Hasher};
use std::vec::Vec;
use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;
use fxhash::{hash, hash32, hash64};

#[derive(Clone, Debug)]
pub struct Kminmer {
    mers: Vec<u64>, // Raw Vec of minimizer hashes
    pub start: usize, // Start location
    pub end: usize, // End location
    pub offset: usize, // Offset (index in the k-min-mer array)
    pub rev: bool, // Strand direction
}


impl Kminmer {
    // Create a new Kminmer object.
    pub fn new(mers: &[u64], start: usize, end: usize, offset: usize) -> Self {
        let mut obj = Kminmer {
            mers: Vec::from(mers),
            start,
            end,
            offset,
            rev: false,
        };
        obj.normalize();
        obj     
    }
    
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
    pub fn mers(&self) -> Vec<u64> {
        self.mers.to_vec()
    }

    // Hash the Vec of minimizer hashes to a u64 (this is used throughout the reference processing).
    pub fn get_hash(&self) -> u64 {
        assert!(self.is_normalized()); // TODO remove this if no issue
        let mut hash = DefaultHasher::new();
        self.mers.hash(&mut hash);
        hash.finish()
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
impl PartialEq for Kminmer {
    fn eq(&self, other: &Kminmer) -> bool {
        self.mers == other.mers
    }
}

impl Eq for Kminmer {
}

impl Hash for Kminmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.mers.hash(state);
    }
}

impl Default for Kminmer {
    fn default() -> Self{Kminmer{mers: vec![], start: 0, end: 0, offset: 0, rev: false}}
}

impl Ord for Kminmer {
    fn cmp(&self, other: &Self) -> Ordering {
        self.mers.cmp(&other.mers)
    }
}

impl PartialOrd for Kminmer{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug)]
pub struct KminmerHash {
    pub hash: u32, // hash from Vec of minimizer hashes
    pub start: usize, // Start location
    pub end: usize, // End location
    pub offset: usize, // Offset (index in the k-min-mer array)
    pub rev: bool, // Strand direction
}


impl KminmerHash {
    // Create a new Kminmer object.
    pub fn new(mers: &Vec<u64>, start: usize, end: usize, offset: usize) -> Self {
        let mut hash;
        let mut rev;
        let mut rev_mers = mers.clone();
        rev_mers.reverse();
        if &rev_mers < mers {
            hash = hash32(&rev_mers);
            rev = true;
        }
        else {
            hash = hash32(mers); 
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
    pub fn get_hash(&self) -> u32 {
        self.hash
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

impl PartialOrd for KminmerHash{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

