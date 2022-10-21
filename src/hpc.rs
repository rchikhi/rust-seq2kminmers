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

// plain version, just hpc the string, no reported positions
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


// second simd attempt based on the first one
// this time a proper HPC string of chars is output
// 
// caution: 
// - unsure if the first nucleotide every 32th is correctly handled. probably not.
// - only works for sequences lengths multiples of 32, otherwise truncates to 
// multiple of 32 below
pub fn encode_rle_simd(s: &str) -> (String,Vec<u32>)
{
    let n : &[u8] = s.as_bytes();
    let ptr = n.as_ptr() as *const i8;
    let len = s.len();
    //let width = 32; // would be nice to be 32 but my machine doesn't support _mm512_mask_compressstoreu_epi8
    let width: usize = 16; // nb of nucleotides considered "in parallel"
    //let npos = width / 8;
    let end_idx = len / width;
    let mut hpc_len : usize= 0;
    unsafe {
        let mut _positions = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
        let mut _16 = _mm512_set_epi32(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16);
        //for i in 0..(npos-1) {
        //    _positions[i+1] = _mm512_add_epi32(_positions[i], _16);
       //}

        // len: worst-case size: hpc string has same length as original
        let res_layout = alloc::Layout::from_size_align_unchecked(len, 8);
        let res_ptr = alloc::alloc(res_layout) as *mut u8;
        let pos_layout = alloc::Layout::from_size_align_unchecked(len << 2, 8);
        let pos_ptr = alloc::alloc(pos_layout) as *mut u32;
        // The docs specify: "ptr needs to have been previously allocated via String/Vec<T> (at
        // least, it’s highly likely to be incorrect if it wasn’t)." And yet somehow.. allocation
        // with the lines below don't work:
        /*let mut res =     String::with_capacity(len);
        let mut pos = Vec::<u32>::with_capacity(len);
        let res_ptr = res.as_mut_ptr(); 
        let pos_ptr = pos.as_mut_ptr(); 
        */
        // groups of 'width' nucleotides at a time
        for i in 0..end_idx {
            // loads 'width' nucleotides into a simd vector
            // failed attempts, mainly because of incompat cpu on my cluster
            //let v = _mm512_loadu_epi8(ptr.offset(i));
            //let v = _mm512_cvtepu8_epi32(v128);
            let v128 = _mm_loadu_epi8(ptr.offset((i as isize) * (width as isize)));
            //println!("v128 {:#x?}",v128);

            // directly create a HPC mask 
            // (1 = keep nucleotide, 0 = should have been discarded by HPC)
            // this time skipping the 2bit conversion
            //let v_shifted = _mm512_bsrli_epi128(v,8); // the 0th and Xth nucl probably isnt well handled here
            let v128_shifted = _mm_slli_si128(v128,1);
            //let _v_cmp_mask = _mm512_cmpeq_epu8_mask(v,v_shifted);
            let mut mask = (! _mm_cmpeq_epi8_mask(v128,v128_shifted)) & 0xFFFE;
            // fix the mask with respect to previous nucleotide
            if i > 0 {
                mask |= (*ptr.offset((i*width) as isize) != *ptr.offset((i*width-1) as isize)) as u16;
            }
            else { 
                mask |= 1;
            }

            // store both the positions and the HPC sequence
            //_mm256_mask_compressstoreu_epi8(res_ptr.offset(hpc_len) as *mut u8, _v_cmp_mask, v256);
            //_mm512_mask_compressstoreu_epi8(res_ptr.offset(hpc_len) as *mut u8, _v_cmp_mask, v);
            //_mm_mask_compressstoreu_epi8(res_ptr.offset(hpc_len) as *mut u8, mask, v128);
            // would be nice to use these ^ but my machine doesnt support it
            let tmp_res_512 = _mm512_setzero_si512();
            let v128_512 = _mm512_cvtepi8_epi32(v128);
            let tmp_res_512 = _mm512_mask_compress_epi32(tmp_res_512, mask, v128_512);
            let tmp_res_128 = _mm512_cvtepi32_epi8(tmp_res_512);
            let mask128 = (1 << mask.count_ones()) - 1;
            _mm_mask_storeu_epi8            (res_ptr.offset(hpc_len as isize) as *mut i8, mask128, tmp_res_128);
            _mm512_mask_compressstoreu_epi32(pos_ptr.offset(hpc_len as isize) as *mut u8, mask,    _positions);
        
             _positions = _mm512_add_epi32(_positions, _16);
            hpc_len += mask.count_ones() as usize;
            //for j in 0..(npos-1) {
            //    _positions[j+1] = _mm512_add_epi32(_positions[j], _16);
            //}
        }

        if len & (width-1) > 0 {
            let i = end_idx;
            // do one last pass for the final bits, with a twist: 
            // read beyond the string buffer (really unsafe, huh).
            // made this code during halloween season so it's ok if spooky
            let v128 = _mm_loadu_epi8(ptr.offset((i as isize) * (width as isize)));
            let v128_shifted = _mm_slli_si128(v128,1);
            let mut mask = (! _mm_cmpeq_epi8_mask(v128,v128_shifted)) & 0xFFFE;
            mask |= (*ptr.offset((i*width) as isize) != *ptr.offset((i*width-1) as isize)) as u16;
            mask &= ((1 << (len & (width-1))) - 1) as u16; // ignore the rest of the buffer, the string end earlier
            let tmp_res_512 = _mm512_setzero_si512();
            let v128_512 = _mm512_cvtepi8_epi32(v128);
            let tmp_res_512 = _mm512_mask_compress_epi32(tmp_res_512, mask, v128_512);
            let tmp_res_128 = _mm512_cvtepi32_epi8(tmp_res_512);
            let mask128 = (1 << mask.count_ones()) - 1;
            _mm_mask_storeu_epi8            (res_ptr.offset(hpc_len as isize) as *mut i8, mask128, tmp_res_128);
            _mm512_mask_compressstoreu_epi32(pos_ptr.offset(hpc_len as isize) as *mut u8, mask,    _positions);
             _positions = _mm512_add_epi32(_positions, _16);
            hpc_len += mask.count_ones() as usize;
        }
       
        // corrupts memory if this is printed before returning the vec
        //println!("hpc str: {:?}",String::from_utf8(Vec::from_raw_parts(res_ptr.clone(), hpc_len, len)).unwrap().clone());
        //(String::from_utf16(&Vec::from_raw_parts(res_ptr, hpc_len, hpc_len)).unwrap(), Vec::from_raw_parts(pos_ptr, hpc_len, hpc_len))
        // also seems to corrupt memory (but only in the benchmark)
        //(String::from_utf8(Vec::from_raw_parts(res_ptr, hpc_len, len)).unwrap(), Vec::from_raw_parts(pos_ptr, hpc_len, len))
        (String::from_raw_parts(res_ptr, hpc_len, len), Vec::from_raw_parts(pos_ptr, hpc_len, len))
        //(res,pos)
    }
    //(String::new(),Vec::new())
}

