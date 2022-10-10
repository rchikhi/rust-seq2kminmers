#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::alloc;

// original from rust-mdbg
pub fn encode_rle(s: &str) -> (String, Vec<usize>)
{
    let mut prev_char = '#';
    let mut hpc_seq = String::new();
    let mut pos_vec = Vec::<usize>::new();
    let mut prev_i = 0;
    for (i, c) in s.chars().enumerate() {
        if c == prev_char && "ACTGactgNn".contains(c) {continue;}
        if prev_char != '#' {
            hpc_seq.push(prev_char);
            pos_vec.push(prev_i);
            prev_i = i;
        }
        prev_char = c;
    }
    hpc_seq.push(prev_char);
    pos_vec.push(prev_i);
    (hpc_seq, pos_vec)
}

// plain version, just hpc the string
pub fn hpc(s: &str) -> String
{
    let mut prev_char = '#';
    let mut hpc_seq = String::with_capacity(s.len());
    for c in s.chars() {
        if c == prev_char {continue;}
        if prev_char != '#' {
            hpc_seq.push(prev_char);
        }
        prev_char = c;
    }
    hpc_seq.push(prev_char);
    hpc_seq
}

// first simd attempt following loosely the idea of
// https://hal.archives-ouvertes.fr/hal-02492824/document
// and based on the movemask nucleotide encoding of Daniel Liu
// https://github.com/Daniel-Liu-c0deb0t/cute-nucleotides/blob/master/src/n_to_bits.rs
//
// caution: only works for sequences lengths multiples of 32, otherwise truncates to 
// multiple of 32 below
pub fn encode_rle_simd(s: &str) -> (Vec<u64>,Vec<u32>)
{
    let n : &[u8] = s.as_bytes();
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5; // final vector size (len input string / (64/2) )
    let len = end_idx + if n.len() & 31 == 0 {0} else {1};
    unsafe {
    let nucl_mask = _mm256_set_epi8(
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7,
    );
        let layout = alloc::Layout::from_size_align_unchecked(len << 3, 8);
        let res_ptr = alloc::alloc(layout) as *mut u64;
        let hpc_layout = alloc::Layout::from_size_align_unchecked(len << 2, 8);
        let hpc_ptr = alloc::alloc(hpc_layout) as *mut u32;
        // groups of 32 nucleotides at a time
        for i in 0..end_idx as isize {
            // loads 32 nucleotides into a simd vector
            let v = _mm256_loadu_si256(ptr.offset(i));

            // permute because unpacks works on the low/high 64 bits in each lane
            // _mm256_permute4x64_epi64(a: __m256i, const IMM8: i32) -> __m256i:
            // Permutes 64-bit integers from a using control mask imm8.
            // Intel doc:
            /* DEFINE SELECT4(src, control) {
            CASE(control[1:0]) OF
            0:      tmp[63:0] := src[63:0]
            1:      tmp[63:0] := src[127:64]
            2:      tmp[63:0] := src[191:128]
            3:      tmp[63:0] := src[255:192]
            ESAC
            RETURN tmp[63:0]
            }
            dst[63:0] := SELECT4(a[255:0], imm8[1:0])
            dst[127:64] := SELECT4(a[255:0], imm8[3:2])
            dst[191:128] := SELECT4(a[255:0], imm8[5:4])
            dst[255:192] := SELECT4(a[255:0], imm8[7:6])
            dst[MAX:256] := 0
            */
            let v = _mm256_permute4x64_epi64(v, 0b11011000);

            // shift each group of two bits for each nucleotide to the end of each byte
            let lo = _mm256_slli_epi64(v, 6);
            let hi = _mm256_slli_epi64(v, 5);

            // interleave bytes then extract the bit at the end of each byte
            // _mm256_unpackhi_epi8(a: __m256i, b: __m256i) -> __m256i: Unpacks and interleave
            // 8-bit integers from the high half of each 128-bit lane in a and b.
            let a = _mm256_unpackhi_epi8(lo, hi);
            let b = _mm256_unpacklo_epi8(lo, hi);

            // at this point a contains a series of 32 bytes, to be considered by groups of two.
            // each pair of two bytes is a nucleotide, more specifically, only the msb of
            // each byte matters

            // zero extend after movemask
            // _mm256_movemask_epi8(a: __m256i) -> i32: Creates mask from the most significant bit
            // of each 8-bit element in a, return the result.
            let ascalar = (_mm256_movemask_epi8(a) as u32) as u64;
            let bscalar = (_mm256_movemask_epi8(b) as u32) as u64;
        
            let res :u64 = (ascalar << 32) | bscalar;
            *res_ptr.offset(i) = res;
           
            // in preparation for HPC mask, perform a AND on the a and b vector to only keep
            // the interesting bits
            let a = _mm256_and_si256(a, nucl_mask);
            let b = _mm256_and_si256(b, nucl_mask);

            // create a HPC mask (1 = keep nucleotide, 0 = should have been discarded by HPC)
            let a_shifted = _mm256_srli_epi16(a,16);
            let a_cmp_mask = _mm256_cmpeq_epi16_mask(a,a_shifted);
            let b_shifted = _mm256_srli_epi16(b,16);
            let b_cmp_mask = _mm256_cmpeq_epi16_mask(b,b_shifted);

            *hpc_ptr.offset(i) = ((a_cmp_mask as u32) << 16) | (b_cmp_mask as u32);
        }

        if n.len() & 31 > 0 {
            // this needs to be added
            //*res_ptr.offset(end_idx as isize) = *n_to_bits_lut(&n[(end_idx << 5)..]).get_unchecked(0);
        }

        (Vec::from_raw_parts(res_ptr, len, len), Vec::from_raw_parts(hpc_ptr, len, len))
    }

}
 
