// A rust binding for ntHash through a C wrapper: https://github.com/rchikhi/ntHash-C

use std::ffi::CString;

use crate::{H};

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct nthashc {
        _unused: [u8; 0],
}
#[allow(non_camel_case_types)]
pub type nthashc_t = nthashc;
extern "C" {
    pub fn nthashc_create(
        s: *const ::std::os::raw::c_char,
        num_hashes: ::std::os::raw::c_uint,
        kmer_size: ::std::os::raw::c_uint,
    ) -> *mut nthashc_t;
}
extern "C" {
        pub fn nthashc_destroy(m: *mut nthashc_t);
}
extern "C" {
        pub fn nthashc_roll(m: *mut nthashc_t) -> ::std::os::raw::c_int;
}
extern "C" {
        pub fn nthashc_hashes(m: *mut nthashc_t, i: ::std::os::raw::c_uint) -> H; // unsure if i'm allowed to declare output as H as it's originally u64 (should work fine when h=u64 though)
}

pub fn nthash_c(seq : & [u8], kmer_size: u32)
{
    let seq_as_char = CString::new(seq).unwrap();
    let num_hashes = 1;
    unsafe {
        let n = nthashc_create(seq_as_char.as_ptr(), num_hashes, kmer_size);
        while nthashc_roll(n) != 0  {
            nthashc_hashes(n,0);
        }
        nthashc_destroy(n);
    }
}
