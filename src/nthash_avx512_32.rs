// This one is faithful to ntHash 1

#![allow(non_snake_case)]
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::alloc;

use crate::{H,FH};


#[derive(Debug)]
#[allow(dead_code)]
pub struct NtHashSIMDIterator<'a> {
    string :&'a [u8],
    _positions :__m512i,
    _hVal :__m512i,
    _fhVal :__m512i,
    _rhVal :__m512i,
    hashes_ptr :*mut u32,
    pos_ptr :*mut u32,
    pos_in_hash :usize,
    nb_ones :usize,
    hash_bound :u32,
    length: usize,
    i: usize,
    k: usize,
}

impl<'a> NtHashSIMDIterator<'a> {
    /// Creates a new NtHashSIMDIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize, hash_bound: u32) -> NtHashSIMDIterator<'a>{
        assert!(k<=31);
        unsafe {
        let length = seq.len();
        // maybe need to pad seq so that it's a multiple of 16 nucleotides. let's see.

        let hashes_layout = alloc::Layout::from_size_align_unchecked(32*16, 8);
        let hashes_ptr = alloc::alloc(hashes_layout) as *mut u32;
        let pos_layout = alloc::Layout::from_size_align_unchecked(32*16, 8);
        let pos_ptr = alloc::alloc(pos_layout) as *mut u32;

        let mut _positions = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
        let mut _16 = _mm512_set_epi32(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16);

        let _ck = _mm512_set1_epi32((32 - (k % 32)) as i32);
        let density = (hash_bound as FH) /(H::MAX as FH); 
        let hash_bound = ((density as f32) * (u32::MAX as f32)) as u32;
        let _hashBound = _mm512_set1_epi32((hash_bound) as i32); // TODO need to figure out why _sometimes_ I need to divide hash_bound by 2 here to get values comparable to HashMode::Regular

        let (_hVal, _fhVal, _rhVal) = _mm512_NTC_epu32_initial(seq, k, _ck);

        //println!("after initial: hval {:#x?} fhval {:#x?} rhval {:#x?}", _hVal, _fhVal, _rhVal);

        let mask = _mm512_cmp_epu32_mask(_hVal, _hashBound, 1 /*LT*/);
        let nb_ones = mask.count_ones() as usize;
        _mm512_mask_compressstoreu_epi32(hashes_ptr as *mut u8, mask, _hVal);
        _mm512_mask_compressstoreu_epi32(pos_ptr as *mut u8, mask, _positions);
        _positions = _mm512_add_epi32(_positions, _16);

        NtHashSIMDIterator {
            string :seq,
            _positions,
            _hVal,
            _fhVal,
            _rhVal,
            hashes_ptr,
            pos_ptr,
            nb_ones,
            pos_in_hash: 0,
            hash_bound,
            length,
            i: 16,
            k,
        }
    }
    }
}


impl<'a> Iterator for NtHashSIMDIterator<'a> {
    type Item = (usize,u32);

    fn next(&mut self) -> Option<(usize,u32)> {
        unsafe {

        let res;
        let sentinel = self.length - self.k + 1;

        if self.pos_in_hash < self.nb_ones 
        {
            let pos =  std::slice::from_raw_parts(self.pos_ptr,    16*32)[self.pos_in_hash] as usize;
            let hash = std::slice::from_raw_parts(self.hashes_ptr, 16*32)[self.pos_in_hash] as u32;
            if pos >= sentinel {
                res = None
            } else { 
                res = Some((pos,hash));
            }
            self.pos_in_hash += 1;
        }
        else
        {
            let _16 = _mm512_set_epi32(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16);
            let _k  = _mm512_set1_epi32(self.k     as i32);
            let _km = _mm512_set1_epi32((self.k-1) as i32);
            let _hashBound = _mm512_set1_epi32((self.hash_bound) as i32); // TODO need to figure out why _sometimes I need to divide hash_bound by 2 here to get values comparable to HashMode::Regular
            let mut _hVal = self._hVal;
            let mut _fhVal = self._fhVal;
            let mut _rhVal = self._rhVal;
            let mut _positions = self._positions;
            let mut i = self.i;
            let k = self.k;
            let string = self.string;

            loop
            {
                if i >= sentinel
                {
                    res = None;
                    break;
                }

                (_hVal, _fhVal, _rhVal) = _mm512_NTC_epu32_sliding(string, i-1, i-1+k, _k, _km, _fhVal, _rhVal);
                //println!("at i {}: hval {:#x?} fhval {:#x?} rhval {:#x?}", i, _hVal, _fhVal, _rhVal);
                i += 16;
        

                let mut mask = _mm512_cmp_epu32_mask(_hVal, _hashBound, 1 /*LT*/);
                let mut nb_ones = mask.count_ones();
                if nb_ones > 0
                {
                    if i >= sentinel {
                        // avoid including info beyond the buffer
                        mask &= (((1 as u32) << (sentinel % 16) as u32) - 1) as u16; // ignore the rest of the buffer, the string end earlier
                        nb_ones = mask.count_ones();
                    }
                    _mm512_mask_compressstoreu_epi32(self.hashes_ptr  as *mut u8, mask, _hVal);
                    _mm512_mask_compressstoreu_epi32(self.pos_ptr  as *mut u8, mask, _positions);                
                }
                _positions = _mm512_add_epi32(_positions, _16);
                if nb_ones > 0
                {
                    res = Some((std::slice::from_raw_parts(self.pos_ptr,    16*32)[0] as usize, 
                                std::slice::from_raw_parts(self.hashes_ptr, 16*32)[0] as u32));
                    self.pos_in_hash = 1;
                    self.nb_ones = nb_ones as usize;
                    break;
                }
            }
            if res != None
            {
                self._hVal = _hVal;
                self._fhVal = _fhVal;
                self._rhVal = _rhVal;
                self._positions = _positions;
                self.i = i;
            }
        }

        res
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size_hint = (self.length as FH * (self.hash_bound as FH) /(u32::max_value() as FH)) as usize;
        (size_hint, Some(size_hint)) 
    }
}


// -------- follows ntHash-AVX2 (old code still in comments), this is a AVX512 32bit hashes version
//

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
const SHIFT :u64 = 0; // used to be 33, as in my vanilla nthash.hpp implementation in the nthash-AVX repo
// but then, it's not compatible with Luiz's rust implementation when hashes are casted as u64

// load forward-strand kmers
unsafe fn _mm512_LKF_epu32(kmerSeq: &[u8], offset: usize) -> __m512i {
    let _seed = _mm512_set_epi32(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        (SEEDT >> SHIFT) as i32,
        (SEEDG >> SHIFT) as i32,
        (SEEDC >> SHIFT) as i32,
        (SEEDA >> SHIFT) as i32);

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
        (SEEDA >> SHIFT) as i32,
        (SEEDC >> SHIFT) as i32,
        (SEEDG >> SHIFT) as i32,
        (SEEDT >> SHIFT) as i32);

    let _kmer = _mm512_permutexvar_epi32(
        _mm512_LKX_epu32(
            kmerSeq, offset),
        _seed);

    _kmer
}


// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
unsafe fn _mm512_NTF_epu32(kmerSeq: &[u8], k: usize) -> __m512i {
    let mut _hVal = _mm512_setzero_si512();

    for i in 0..k
    {
        //_hVal = _mm512_rori31_epu32::<30,{32-30}>(_hVal);
        _hVal = _mm512_rol_epi32(_hVal, 1);

        let _kmer = _mm512_LKF_epu32(kmerSeq, i);

        _hVal = _mm512_xor_epi32(
            _hVal,
            _kmer);
    }

    _hVal
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
unsafe fn _mm512_NTR_epu32(kmerSeq: &[u8], k: usize, _ck: __m512i) -> __m512i {
    let _zero = _mm512_setzero_si512();

    let mut _hVal = _zero;

    for i in 0..k
    {

        let mut _kmer = _mm512_LKR_epu32(kmerSeq, i);
        
        /*_kmer = _mm512_rorv31_epu32(
            _kmer,
            _ck);*/

        _kmer = _mm512_rorv_epi32(
            _kmer,
            _ck);

        _hVal = _mm512_xor_epi32(
            _hVal,
            _kmer);

        //_hVal = _mm512_rori31_epu32::<1,{32-1}>(_hVal);
        _hVal = _mm512_ror_epi32(_hVal, 1);
    }

    _hVal
}

// canonical ntHash
unsafe fn _mm512_NTC_epu32_initial(kmerSeq: &[u8], k: usize, _ck: __m512i) -> (__m512i, __m512i, __m512i) {
    let _fhVal = _mm512_NTF_epu32(kmerSeq, k);
    let _rhVal = _mm512_NTR_epu32(kmerSeq, k, _ck);

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

    /*_out31 = _mm512_rorv31_epu32(
        _out31,
        _ck);
    */
    _out31 = _mm512_rolv_epi32(
        _out31,
        _k);

    let mut _kmer = _mm512_xor_epi32(
        _in31,
        _out31);

    // scan-shift kmers
    // follows
    // https://upload.wikimedia.org/wikipedia/commons/thumb/e/ec/Hillis-Steele_Prefix_Sum.svg/450px-Hillis-Steele_Prefix_Sum.svg.png
    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xfffe,
//            _mm512_rori31_epu32::<30, {32 - 30}>(
            _mm512_rol_epi32(
                _kmer, 1)));

    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xfffc,
//            _mm512_rori31_epu32::<29, {32 - 29}>(
            _mm512_rol_epi32(
                _kmer, 2)));

    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xfff0,
//            _mm512_rori31_epu32::<27, {32 - 27}>(
            _mm512_rol_epi32(
                _kmer, 4)));

    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xff00,
//            _mm512_rori31_epu32::<23, {32 - 23}>(
            _mm512_rol_epi32(
                _kmer, 8)));


    // var-shift the hash
    let mut _hVal = _mm512_permutexvar_epi32(
        _mm512_set1_epi32(15),
        _fhVal);

    let _shift31 = _mm512_set_epi32(
   //     15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30);
        16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);

    /*
    _hVal = _mm512_rorv31_epu32(
        _hVal,
        _shift31);
    */
    
    _hVal = _mm512_rolv_epi32(
        _hVal,
        _shift31);

    // merge everything together
    _hVal = _mm512_xor_epi32(
        _hVal,
        _kmer);

    _hVal
}


// reverse-complement ntHash for sliding k-mers
unsafe fn _mm512_NTR_epu32_sliding(_rhVal: __m512i, _km: __m512i, kmerSeq :&[u8], offsetOut :usize , offsetIn :usize) -> __m512i{
    // construct input kmers
    let mut _in31 = _mm512_LKR_epu32(kmerSeq, offsetIn);

    //_in31 = _mm512_rorv31_epu32(
    _in31 = _mm512_rolv_epi32(
        _in31,
        _km);

    let mut _out31 = _mm512_LKR_epu32(kmerSeq, offsetOut);

    _out31 = _mm512_ror_epi32(
        _out31,
        1);

    let mut _kmer = _mm512_xor_epi32(
        _in31,
        _out31);

    // scan-shift kmers
    // follows
    // https://upload.wikimedia.org/wikipedia/commons/thumb/e/ec/Hillis-Steele_Prefix_Sum.svg/450px-Hillis-Steele_Prefix_Sum.svg.png
    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xfffe,
//            _mm512_rori31_epu32::<1, {32 - 1}>(
            _mm512_ror_epi32(
                _kmer, 1)));

    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xfffc,
//            _mm512_rori31_epu32::<2, {32 - 2}>(
            _mm512_ror_epi32(
                _kmer, 2)));

    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xfff0,
//            _mm512_rori31_epu32::<4, {32 - 4}>(
            _mm512_ror_epi32(
                _kmer, 4)));

    _kmer = _mm512_xor_epi32(
        _kmer,
        _mm512_maskz_expand_epi32(
            0xff00,
//            _mm512_rori31_epu32::<8, {32 - 8}>(
            _mm512_ror_epi32(
                _kmer, 8)));

    // var-shift the hash
    let mut _hVal = _mm512_permutexvar_epi32(
        _mm512_set1_epi32(15),
        _rhVal);
    

    let _shift31 = _mm512_set_epi32(
   //     16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,31);
        16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);

    //_hVal = _mm512_rorv31_epu32(
    _hVal = _mm512_rorv_epi32(
        _hVal,
        _shift31);

    // merge everything together
    _hVal = _mm512_xor_epi32(
        _hVal,
        _kmer);

    //_hVal = _mm512_rori31_epu32::<1, {32 - 1}>(_hVal);

    _hVal
}


// canonical ntHash for sliding k-mers
unsafe fn _mm512_NTC_epu32_sliding(kmerSeq: &[u8], offsetOut: usize, offsetIn: usize, _k: __m512i, _km: __m512i, _fhVal: __m512i, _rhVal: __m512i) -> (__m512i, __m512i, __m512i){
    let _fhVal = _mm512_NTF_epu32_sliding(_fhVal, _k,  kmerSeq, offsetOut, offsetIn);
    let _rhVal = _mm512_NTR_epu32_sliding(_rhVal, _km, kmerSeq, offsetOut, offsetIn);

    let _hVal = _mm512_mask_blend_epi32(
        _mm512_cmpgt_epu32_mask(
            _fhVal,
            _rhVal),
        _fhVal,
        _rhVal);

    (_hVal, _fhVal, _rhVal)
}

