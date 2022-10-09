#![allow(non_snake_case)]
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::alloc;

use crate::{H, FH}; // hash precision


#[derive(Debug)]
#[allow(dead_code)]
pub struct NtHashSIMDIterator {
    items: Vec<(usize,H)>,
    pos: usize,
}


impl NtHashSIMDIterator {
    /// Creates a new NtHashSIMDIterator with internal state properly initialized.
    pub fn new(seq: &[u8], k: usize, hash_bound: H) -> NtHashSIMDIterator{
        let _seq_len = seq.len();

        let items = nthash_simd_32bits(seq, k, hash_bound);

        NtHashSIMDIterator {
            items,
            pos: 0,
        }
    }
}

impl Iterator for NtHashSIMDIterator {
    type Item = (usize,H);

    fn next(&mut self) -> Option<(usize,H)> {
         let res;
         if self.pos < self.items.len()
         {
             res = Some(self.items[self.pos]);
             self.pos += 1;
         }
         else {
             res = None;
         }
         res
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size_hint = self.items.len();
        (size_hint, Some(size_hint)) 
    }
}



// convert kmers 8 -> 3 bit representation, N character is mapped to 4
//
unsafe fn _mm_CKX_epu8(_kmerSeq: __m128i) -> __m128i{
    // _mm_set1_epi8(a: i8) -> __m128i: Broadcasts 8-bit integer a to all elements.
    let _mask = _mm_set1_epi8(0x0f);

    // _mm_set_epi8([16 args]) -> __m128i: Sets packed 8-bit integers with the supplied values.
    let _table = _mm_set_epi8(
        4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 3, 1, 4, 0, 4);

    // _mm_shuffle_epi8(a: __m128i, b: __m128i) -> __m128i: Shuffles bytes from a according to the content of b.
    _mm_shuffle_epi8(
        _table,
        // _mm_and_si128(a: __m128i, b: __m128i) -> __m128i: Computes the bitwise AND of 128 bits (representing integer data) in a and b.
        _mm_and_si128(
            _kmerSeq,
            _mask))
}


unsafe fn _mm512_rori31_epu32<const IMM: u32, const THIRTYTWO_MINUS_IMM: u32>(_v: __m512i) -> __m512i {
    _mm512_or_epi32(
        _mm512_srli_epi32(
            _v,
            IMM),
        _mm512_srli_epi32(
            _mm512_slli_epi32(
                _v,
                THIRTYTWO_MINUS_IMM),
            1))
}



// rotate 31-right bits of "_v" to the right by _s position
// elements of _s must be less than 31
unsafe fn _mm512_rorv31_epu32(_v: __m512i, _s: __m512i) -> __m512i {
    let _32 = _mm512_set1_epi32(32);

    _mm512_or_epi32(
        _mm512_srlv_epi32(
            _v,
            _s),
        _mm512_srli_epi32(
            _mm512_sllv_epi32(_v,
                _mm512_sub_epi32(
                    _32,
                    _s)),
            1))
}


unsafe fn _mm512_LKX_epu32(kmerSeq: &[u8], offset :usize) -> __m512i {
    _mm512_cvtepu8_epi32(
        _mm_CKX_epu8(
            _mm_loadu_si128(
                kmerSeq.as_ptr().offset(offset as isize) as *const _)))
}

// 64-bit random seeds corresponding to bases and their complements
const SEEDA :u64 = 0x3c8bfbb395c60474;
const SEEDC :u64 = 0x3193c18562a02b4c;
const SEEDG :u64 = 0x20323ed082572324;
const SEEDT :u64 = 0x295549f54be24456;
//const SEEDN :u64 = 0x0000000000000000;

// load forward-strand kmers
unsafe fn _mm512_LKF_epu32(kmerSeq: &[u8], offset: usize) -> __m512i {
    let _seed = _mm512_set_epi32(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        (SEEDT >> 33) as i32,
        (SEEDG >> 33) as i32,
        (SEEDC >> 33) as i32,
        (SEEDA >> 33) as i32);

    let _kmer = _mm512_permutexvar_epi32(
        _mm512_LKX_epu32(
            kmerSeq, offset),
        _seed);

    _kmer
}

// load reverse-strand kmers
unsafe fn _mm512_LKR_epu32(kmerSeq: &[u8], offset: usize) -> __m512i {
    let _seed = _mm512_set_epi32(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        (SEEDA >> 33) as i32,
        (SEEDC >> 33) as i32,
        (SEEDG >> 33) as i32,
        (SEEDT >> 33) as i32);

    let _kmer = _mm512_permutexvar_epi32(
        _mm512_LKX_epu32(
            kmerSeq, offset),
        _seed);

    _kmer
}


// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
unsafe fn _mm512_NTF_epu32(kmerSeq: &[u8], k: usize) -> __m512i {
    let mut _hVal31 = _mm512_setzero_si512();

    for i in 0..k
    {
        _hVal31 = _mm512_rori31_epu32::<30,{32-30}>(_hVal31);

        let _kmer31 = _mm512_LKF_epu32(kmerSeq, i);

        _hVal31 = _mm512_xor_epi32(
            _hVal31,
            _kmer31);
    }

    _hVal31
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
unsafe fn _mm512_NTR_epu32(kmerSeq: &[u8], k: usize, _k: __m512i) -> __m512i {
    let _zero = _mm512_setzero_si512();

    let mut _hVal31 = _zero;

    for i in 0..k
    {
        let mut _kmer31 = _mm512_LKR_epu32(kmerSeq, i);

        _kmer31 = _mm512_rorv31_epu32(
            _kmer31,
            _k);

        _hVal31 = _mm512_xor_epi32(
            _hVal31,
            _kmer31);

        _hVal31 = _mm512_rori31_epu32::<1,{32-1}>(_hVal31);
    }

    _hVal31
}

// canonical ntHash
unsafe fn _mm512_NTC_epu32_initial(kmerSeq: &[u8], k: usize, _k: __m512i) -> (__m512i, __m512i, __m512i) {
    let _fhVal = _mm512_NTF_epu32(kmerSeq, k);
    let _rhVal = _mm512_NTR_epu32(kmerSeq, k, _k);

    let _hVal = _mm512_mask_blend_epi32(
        _mm512_cmpgt_epu32_mask(
            _fhVal,
            _rhVal),
        _fhVal,
        _rhVal);
    (_hVal, _fhVal, _rhVal)
}

// ----------------------- sliding function


// forward-strand ntHash for sliding k-mers
unsafe fn _mm512_NTF_epu32_sliding(_fhVal: __m512i, _k: __m512i, kmerSeq :&[u8], offsetOut :usize , offsetIn :usize) -> __m512i {

    // construct input kmers
    let mut _in31 = _mm512_LKF_epu32(kmerSeq, offsetIn);

    let mut _out31 = _mm512_LKF_epu32(kmerSeq, offsetOut);

    _out31 = _mm512_rorv31_epu32(
        _out31,
        _k);

    let mut _kmer31 = _mm512_xor_epi32(
        _in31,
        _out31);

    // scan-shift kmers
    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xfffe,
            _mm512_rori31_epu32::<30, {32 - 30}>(
                _kmer31)));

    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xfffc,
            _mm512_rori31_epu32::<29, {32 - 29}>(
                _kmer31)));

    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xfff0,
            _mm512_rori31_epu32::<27, {32 - 27}>(
                _kmer31)));

    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xff00,
            _mm512_rori31_epu32::<23, {32 - 23}>(
                _kmer31)));

    // var-shift the hash
    let mut _hVal31 = _mm512_permutexvar_epi32(
        _mm512_set1_epi32(15),
        _fhVal);

    let _shift31 = _mm512_set_epi32(
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30);

    _hVal31 = _mm512_rorv31_epu32(
        _hVal31,
        _shift31);

    // merge everything together
    _hVal31 = _mm512_xor_epi32(
        _hVal31,
        _kmer31);

    _hVal31
}




// reverse-complement ntHash for sliding k-mers
unsafe fn _mm512_NTR_epu32_sliding(_rhVal: __m512i, _k: __m512i, kmerSeq :&[u8], offsetOut :usize , offsetIn :usize) -> __m512i{
    // construct input kmers
    let mut _in31 = _mm512_LKR_epu32(kmerSeq, offsetIn);

    _in31 = _mm512_rorv31_epu32(
        _in31,
        _k);

    let _out31 = _mm512_LKR_epu32(kmerSeq, offsetOut);

    let mut _kmer31 = _mm512_xor_epi32(
        _in31,
        _out31);

    // scan-shift kmers
    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xfffe,
            _mm512_rori31_epu32::<1, {32 - 1}>(
                _kmer31)));

    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xfffc,
            _mm512_rori31_epu32::<2, {32 - 2}>(
                _kmer31)));

    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xfff0,
            _mm512_rori31_epu32::<4, {32 - 4}>(
                _kmer31)));

    _kmer31 = _mm512_xor_epi32(
        _kmer31,
        _mm512_maskz_expand_epi32(
            0xff00,
            _mm512_rori31_epu32::<8, {32 - 8}>(
                _kmer31)));

    // var-shift the hash
    let mut _hVal31 = _mm512_permutexvar_epi32(
        _mm512_set1_epi32(15),
        _rhVal);

    let _shift31 = _mm512_set_epi32(
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

    _hVal31 = _mm512_rorv31_epu32(
        _hVal31,
        _shift31);

    // merge everything together
    _hVal31 = _mm512_xor_epi32(
        _hVal31,
        _kmer31);

    _hVal31 = _mm512_rori31_epu32::<1, {32 - 1}>(_hVal31);

    _hVal31
}


// canonical ntHash for sliding k-mers
unsafe fn _mm512_NTC_epu32_sliding(kmerSeq: &[u8], offsetOut: usize, offsetIn: usize, _k: __m512i, _fhVal: __m512i, _rhVal: __m512i) -> (__m512i, __m512i, __m512i){
    let _fhVal = _mm512_NTF_epu32_sliding(_fhVal, _k, kmerSeq, offsetOut, offsetIn);
    let _rhVal = _mm512_NTR_epu32_sliding(_rhVal, _k, kmerSeq, offsetOut, offsetIn);

    let _hVal = _mm512_mask_blend_epi32(
        _mm512_cmpgt_epu32_mask(
            _fhVal,
            _rhVal),
        _fhVal,
        _rhVal);

    (_hVal, _fhVal, _rhVal)
}


// follows ntHash-AVX2, this is the AVX512 32bit hashes version
pub fn nthash_simd_32bits(n: &[u8], k: usize, hash_bound: H) -> Vec<(usize,H)>
{
    unsafe { 
    let debug = false;
    let length = n.len();
    // maybe need to pad n
    let hashes_layout = alloc::Layout::from_size_align_unchecked(length*32, 8);
    let hashes_ptr = alloc::alloc(hashes_layout) as *mut u32;
    let pos_layout = alloc::Layout::from_size_align_unchecked(length*32, 8);
    let pos_ptr = alloc::alloc(pos_layout) as *mut u32;
    let mut hashes_offset :u32 = 0;
    
        let _zero = _mm512_setzero_si512();
        let mut _pos = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
        let mut _16 = _mm512_set_epi32(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16);

        let _k = _mm512_set1_epi32((31 - (k % 31)) as i32);
        let density = (hash_bound as FH) /(H::max_value() as FH); 
        let hash_bound = ((density as f32) * (u32::max_value() as f32)) as u32;
        let _hashBound = _mm512_set1_epi32(hash_bound as i32);

        let (mut _hVal, mut _fhVal, mut _rhVal) = _mm512_NTC_epu32_initial(n, k, _k);

        let mut lo= _mm512_extracti64x4_epi64(_hVal, 0);
        let mut hval0 = _mm256_extract_epi32(lo, 0);
        if debug { println!("first hash AVX512x32 {:x}", hval0); }

        let sentinel = length - k + 1;
        let mut i = 16;
        while i < sentinel
        {
            (_hVal, _fhVal, _rhVal) = _mm512_NTC_epu32_sliding(n, i-1+k, i-1, _k, _fhVal, _rhVal);

            let mask = _mm512_cmp_epi32_mask(_hVal, _hashBound, 1 /*LT*/);
            _mm512_mask_compressstoreu_epi32(hashes_ptr.offset(hashes_offset as isize) as *mut u8, mask, _hVal);
            _mm512_mask_compressstoreu_epi32(pos_ptr.offset(hashes_offset as isize) as *mut u8, mask, _pos);
            let nb_ones = mask.count_ones();
            hashes_offset += nb_ones;

            _pos = _mm512_add_epi32(_pos, _16);
            i += 16;
        }

	if (length - k) % 16 < 8{
	    lo = _mm512_extracti64x4_epi64(_hVal, 0); }
	else {
	    lo = _mm512_extracti64x4_epi64(_hVal, 1); }
	if (length - k) % 8 == 0 {
	    hval0 = _mm256_extract_epi32(lo, 0); }
	else if (length - k) % 8 == 1 {
	    hval0 = _mm256_extract_epi32(lo, 1); }
	else if (length - k) % 8 == 2 {
	    hval0 = _mm256_extract_epi32(lo, 2); } 
	else if (length - k) % 8 == 3 { 
	    hval0 = _mm256_extract_epi32(lo, 3); }
	else if (length - k) % 8 == 4 {
	    hval0 = _mm256_extract_epi32(lo, 4); }
	else if (length - k) % 8 == 5 {
	    hval0 = _mm256_extract_epi32(lo, 5); }
	else if (length - k) % 8 == 6 { 
	    hval0 = _mm256_extract_epi32(lo, 6); }
	else if (length - k) % 8 == 7 {
	    hval0 = _mm256_extract_epi32(lo, 7); }
        if debug { println!("final hash AVX512x32 {:x}", hval0); }


    let vec_hashes : Vec<u32> = Vec::from_raw_parts(hashes_ptr as *mut u32, hashes_offset as usize, hashes_offset as usize);
    let vec_pos : Vec<u32> = Vec::from_raw_parts(pos_ptr as *mut u32, hashes_offset as usize, hashes_offset as usize);
    
    vec_pos.iter().map(|&p| p as usize).zip(vec_hashes.iter().map(|&p| p as u64)).collect()
    }
}
 
