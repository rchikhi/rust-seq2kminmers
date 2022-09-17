use std::hash::{Hash, Hasher};
use std::vec::Vec;
use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;

#[derive(Clone, Debug)]
pub struct Kminmer {
    mers: Vec<u64>,
    pub start: usize,
    pub end: usize,
    pub offset: usize,
    pub rev: bool,
}

impl Kminmer {
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
    pub fn normalize(&mut self) {
        let mut rev_mers = self.mers.clone();
        rev_mers.reverse();
        if rev_mers < self.mers {
            self.mers = rev_mers.to_vec();
            self.rev = true;
        }
    }

    pub fn print(&self) -> String {
	// prints only the first 2 digits of each minimizer hash
	let mut s = String::new();
	for x in self.mers.iter() {
	    s = format!("{}{} ", s, x.to_string()[..2].to_string());
	}
	s
    }
    pub fn mers(&self) -> Vec<u64> {
        self.mers.to_vec()
    }
    pub fn get_hash(&self) -> u64 {
        let mut hash = DefaultHasher::new();
        self.mers.hash(&mut hash);
        hash.finish()
    }

}

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
