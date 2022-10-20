
/*unfinished code
 * */

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::alloc;

use crate::{H}; // hash precision


#[derive(Debug)]
#[allow(dead_code)]
pub struct NtHashHPCSIMDIterator {
    items: Vec<H>,
    pos: usize,
}


impl NtHashHPCSIMDIterator {
    /// Creates a new NtHashHPCSIMDIterator with internal state properly initialized.
    pub fn new(seq: &[u8], k: usize, hash_bound: H) -> NtHashHPCSIMDIterator{
        let _seq_len = seq.len();

        let items = nthash_hpc_simd_32(seq, k, hash_bound);

        NtHashHPCSIMDIterator {
            items,
            pos: 0,
        }
    }
}

impl Iterator for NtHashHPCSIMDIterator {
    type Item = (usize,H);


    fn next(&mut self) -> Option<(usize,H)> {
         //Some((prev_current_idx, hash))
         Some((0,0))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size_hint = self.items.len();
        (size_hint, Some(size_hint)) 
    }
}

// follows encode_rle_simd
// special case for ki=32, facilitates things
// 
// this is a WIP, completely unfinished (maybe for a long time)
//
// caution: only works for sequences lengths multiples of 32, otherwise truncates to 
// multiple of 32 below
#[allow(unused_assignments,unused_variables)]
pub fn nthash_hpc_simd_32(n: &[u8], k: usize, hash_bound: H) -> Vec<H>
{
    assert!(k == 32);
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5; // final vector size (len input string / (64/2) )
    let len = end_idx + if n.len() & 31 == 0 {0} else {1};

    unsafe {
        let hpc_layout = alloc::Layout::from_size_align_unchecked(len << 2, 8);
        let hpc_ptr = alloc::alloc(hpc_layout) as *mut H;

        // TODO need special treatment for first 256 bytes (32 nucleotides)

        // groups of 32 nucleotides at a time
        for i in 1..end_idx as isize {
            // loads 32 nucleotides into a simd vector
            let v = _mm256_loadu_si256(ptr.offset(i));

            // create a HPC mask (1 = keep nucleotide, 0 = should have been discarded by HPC)
            let v_shifted = _mm256_srli_epi16(v,32);
            let _v_cmp_mask = _mm256_cmpeq_epi16_mask(v,v_shifted);

            // TODO need a special treatment for the first byte given they may be
            // repeated from previous last 32 bytes

            //*hpc_ptr.offset(i) = ((a_cmp_mask as u32) << 16) | (b_cmp_mask as u32);
            
            // TODO draw the rest of the owl
        }

        if n.len() & 31 > 0 {
            // this needs to be added
            //*res_ptr.offset(end_idx as isize) = *n_to_bits_lut(&n[(end_idx << 5)..]).get_unchecked(0);
        }

        Vec::from_raw_parts(hpc_ptr, len, len)
    }

}
 
