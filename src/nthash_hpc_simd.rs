
use crate::{encode_rle_simd, NtHashSIMDIterator};
use crate::{H}; // hash precision

//const TEST_CONST_K : usize = 32; 
// tested if giving a const instead of a variable k helps
// spoiler: just a bit (9% perf gain on a 40secs test, ie 3.x secs)
// TODO: revisit that

const BUFLEN :usize = 256;
const MAXIMUM_K_SIZE: usize = BUFLEN;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("K size {ksize} is out of range for the given sequence size {seq_size}")]
    KSizeOutOfRange { ksize: usize, seq_size: usize },
    #[error("K size {0} cannot exceed {MAXIMUM_K_SIZE}")]
    KSizeTooBig(usize),
}

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
pub struct NtHashHPCSIMDIterator<'a> {
    hpc_seq : String,
    hpc_pos :Vec<u32>,
    it : Option<NtHashSIMDIterator<'a>>
}

// this is just a wrapper over NtHashSIMDIterator except it passes the HPC string
// could be futher optimized if i basically combined HPC+NtHash-SIMD inside a single function
// would use a buffer instead of having the entire HPC string.
impl<'a> NtHashHPCSIMDIterator<'a> {
    pub fn new(seq: &'a [u8], l: usize, hash_bound: H) -> Result<NtHashHPCSIMDIterator<'a>> {
        let (hpc_seq, hpc_pos) = encode_rle_simd(seq);

        //let it = NtHashSIMDIterator::new(hpc_seq.as_bytes().as_ref(), l, hash_bound);

        let mut res = NtHashHPCSIMDIterator {
            hpc_seq: hpc_seq,
            hpc_pos,
            it: None
        };

        // OK I was trying to make a self referential structure.
        // And seems rust dislike them very much.
        // So look at what it made me do:
        let hpc_seq = res.hpc_seq.as_bytes().as_ptr() as *const u8;
        unsafe{
        res.it = Some(NtHashSIMDIterator::new(std::slice::from_raw_parts(hpc_seq,res.hpc_seq.len()), l, hash_bound));
        }
        Ok(res)
    }
}

impl<'a> Iterator for NtHashHPCSIMDIterator<'a> {
    type Item = (usize,H);


    fn next(&mut self) -> Option<(usize,H)> {
        let res = self.it.as_mut().unwrap_or_else(|| unsafe { std::hint::unreachable_unchecked() }).next();
        if let Some(res) = res {
            Some((self.hpc_pos[res.0] as usize,res.1))
        } else { 
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.it.as_ref().unwrap_or_else(|| unsafe { std::hint::unreachable_unchecked() }).size_hint() 
    }
}
