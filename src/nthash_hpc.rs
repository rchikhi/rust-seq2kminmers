// homopolymer compression version of ntHash1,
// which simulatenously HPC and record original offsets of sequences
//
//
// single-file version
// code adapted from ntHash1 and mainly from Luiz Irber's Rust crate port
// kept only the canonical version (no 'forward')

const BUFLEN :usize = 256;
// this implem is more limited
const MAXIMUM_K_SIZE: usize = BUFLEN;

use crate::{H, FH}; // hash precision

//const TEST_CONST_K : usize = 32; 
// tested if giving a const instead of a variable k helps
// spoiler: just a bit (9% perf gain on a 40secs test, ie 3.x secs)

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("K size {ksize} is out of range for the given sequence size {seq_size}")]
    KSizeOutOfRange { ksize: usize, seq_size: usize },
    #[error("K size {0} cannot exceed {MAXIMUM_K_SIZE}")]
    KSizeTooBig(usize),
}

pub type Result<T> = std::result::Result<T, Error>;

#[allow(overflowing_literals)]
const H_LOOKUP: [H; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x3c8b_fbb3_95c6_0474 as H;
    lookup[b'C' as usize] = 0x3193_c185_62a0_2b4c as H;
    lookup[b'G' as usize] = 0x2032_3ed0_8257_2324 as H;
    lookup[b'T' as usize] = 0x2955_49f5_4be2_4456 as H;
    lookup[b'N' as usize] = 0;
    lookup
};

#[allow(overflowing_literals)]
const RC_LOOKUP: [H; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x2955_49f5_4be2_4456 as H;
    lookup[b'C' as usize] = 0x2032_3ed0_8257_2324 as H;
    lookup[b'G' as usize] = 0x3193_c185_62a0_2b4c as H;
    lookup[b'T' as usize] = 0x3c8b_fbb3_95c6_0474 as H;
    lookup[b'N' as usize] = 0;
    lookup
};

#[inline(always)]
fn h(c: u8) -> H {
    unsafe{
    let val = H_LOOKUP.get_unchecked(c as usize);
    // this branch.. is unsatisfactory. but remarkably, having it doesn't seem to make a difference in benchs
  /*  if *val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", c as char)
    }*/
    *val
    }
}

#[inline(always)]
fn rc(nt: u8) -> H {
    unsafe{
    let val = RC_LOOKUP.get_unchecked(nt as usize);
/*    if *val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", nt as char)
    }*/
    *val
    }
}

/// An efficient iterator for returning rev-comp aware hashes in HPC space, keeping only those
/// below a hash bound. That's a very specialized application, useful for rust-mdbg.
///
/// We use a variable hash type H (default is u64 for classical ntHash, but u32 is also ok and
/// saves time)
///
/// Since it implements the `Iterator` trait it also
/// exposes many other useful methods. In this example we use `collect` to
/// generate all hashes and put them in a `Vec<H>`.
/// ```
///     use rust_seq2kminmers::NtHashHPCIterator;
///
///     # fn main() {
///     let seq = b"ACTGCACATGATGAGTAGATGATGATGATGATGATATGATGATAT";
///     let density = 0.1;
///     let hash_bound = ((density as FH) * (H::max_value() as FH)) as H;
///     let iter = NtHashHPCIterator::new(seq, 4, hash_bound).unwrap();
///     let hashes: Vec<(usize,H)> = iter.collect();
///     assert_eq!(hashes,
///                vec![(0, 1693589515812555183), (6, 876319423165292601), (13,771890730643629033), (16, 826464090118103095), (33, 1245321008145464903), (34,1193606442387228521)]);
///     # }
/// ```

#[derive(Debug)]
pub struct NtHashHPCIterator<'a> {
    seq: &'a [u8],
    k: usize, // TODO: I wonder what would happen if i set k has a const=32, in terms of performance
    fh: H,
    rh: H,
    current_idx_plus_k: usize,
    seq_len: usize,
    hash_bound: H,
    h_rc_buffer: [(H,H);  BUFLEN],
    idx_buffer: [usize;BUFLEN],
    buffer_pos: usize,
    first_element: bool,
}

impl<'a> NtHashHPCIterator<'a> {
    /// Creates a new NtHashHPCIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize, hash_bound: H) -> Result<NtHashHPCIterator<'a>> {
        let seq_len = seq.len();
        if k > seq_len {
            return Err(Error::KSizeOutOfRange {
                ksize: k,
                seq_size: seq.len(),
            });
        }
        if k > MAXIMUM_K_SIZE {
            return Err(Error::KSizeTooBig(k));
        }
        let mut fh = 0;
        let mut j = 0; // current position in non-HPC space
        let mut i = 0; // current position in HPC space
        let mut prev;
        let mut prev_j = 0;
        let mut v;

        assert!(k < 256); // who uses large minimizers anyway?
        let mut h_rc_buffer:   [(H,H);  BUFLEN]  = [(0,0); BUFLEN];
        let mut idx_buffer: [usize;BUFLEN]  = [0; BUFLEN];

        // pre-compute hash of first kmer in HPC space
        while i < k && j < seq_len
        {
            v = seq[j];
            let hv = h(v);
            h_rc_buffer[i] = (hv,0);
            idx_buffer[i] = j;
            fh ^= hv.rotate_left((k - i - 1) as u32);
            i += 1;
            // HPC on the fly
            prev = v;
            prev_j = j;
            while j < seq_len &&  seq[j] == prev { j += 1};
        }

        // if sequence is shorter than k in HPC space: hash will be incorrect and i<k , but that's ok, 
        // hash is only returned later in the iterator where proper length checks are performed
        assert!( (j >= seq_len && i < k) || (j <= seq_len && i==k) );

        // at this point just assume we're at position i==k in HPC space, j in non-HPC space,
        // read the sequence backwards to get hash of reverse complement
        i -= 1; 
        j = prev_j;
        let k_in_non_hpc_space = j;
        let mut rh = 0;
        loop
        {
            v = seq[j];
            let rcv = rc(v);
            let (hv, _osef) = h_rc_buffer[i];
            h_rc_buffer[i] = (hv,rcv);
            rh ^= rcv.rotate_left(i as u32);
            if i == 0 { break;}
            i -= 1;
            // HPC on the fly
            prev = v;
            while j > 0 && seq[j] == prev { j -= 1};
        }

        Ok(NtHashHPCIterator {
            seq,
            k,
            fh,
            rh,
            current_idx_plus_k: k_in_non_hpc_space,
            seq_len,
            hash_bound,
            h_rc_buffer,
            idx_buffer,
            buffer_pos : 0,
            first_element: true
        })
    }
}

impl<'a> Iterator for NtHashHPCIterator<'a> {
    type Item = (usize,H);


    fn next(&mut self) -> Option<(usize,H)> {
        unsafe {
            //let k = TEST_CONST_K;
            let k = self.k;
            let mut hash;
            let prev_current_idx;
            let mut prev :u8 = *self.seq.get_unchecked(self.current_idx_plus_k);
            let mut cur : u8 = prev;
            let mut h_seqk :H;
            let mut rc_seqk :H;

            // special case for very first element
            if std::intrinsics::unlikely(self.first_element)
            {
                self.first_element = false;
                loop
                {
                    self.current_idx_plus_k += 1;
                    cur = *self.seq.get_unchecked(self.current_idx_plus_k) ;
                    if self.current_idx_plus_k >= self.seq_len || cur != prev {
                        break
                    }
                }

                if self.current_idx_plus_k >= self.seq_len  {
                    return None;
                };

                // prepare everything for the next char
                (h_seqk, rc_seqk) = (h(cur),rc(cur));
                let pos = std::intrinsics::unchecked_rem(self.buffer_pos+k,BUFLEN);
                *self.h_rc_buffer.get_unchecked_mut(pos) = (h_seqk, rc_seqk);
                *self.idx_buffer. get_unchecked_mut(pos) = self.current_idx_plus_k;
                self.buffer_pos += 1;

                hash = H::min(self.rh, self.fh);
                if hash <= self.hash_bound {
                    prev_current_idx = *self.idx_buffer.get_unchecked(std::intrinsics::unchecked_rem(self.buffer_pos+BUFLEN-1,BUFLEN));
                    return Some((prev_current_idx, hash))
                }
            }
            else  {
                (h_seqk, rc_seqk)  = *self.h_rc_buffer .get_unchecked(std::intrinsics::unchecked_rem(self.buffer_pos+k-1,BUFLEN));
            }

            loop
            {
                let (h_seqi, rc_seqi)  = self.h_rc_buffer .get_unchecked(std::intrinsics::unchecked_rem(self.buffer_pos-1,BUFLEN));

                self.fh = self.fh.rotate_left(1) ^ h_seqi.rotate_left(k as u32) ^ h_seqk;

                self.rh = self.rh.rotate_right(1)
                    ^ rc_seqi.rotate_right(1)
                    ^ rc_seqk.rotate_left(k as u32 - 1);

                // move to next char
                prev = cur;
                loop
                {
                    self.current_idx_plus_k += 1;
                    if  std::intrinsics::unlikely(self.current_idx_plus_k >= self.seq_len) {
                        break;
                    }
                    cur = *self.seq.get_unchecked(self.current_idx_plus_k);
                    if cur != prev {
                        break
                    }
                }

                if self.current_idx_plus_k >= self.seq_len {
                    return None;
                };

                // prepare everything for the next char
                (h_seqk, rc_seqk) = (h(cur),rc(cur));
                let pos = std::intrinsics::unchecked_rem(self.buffer_pos+k,BUFLEN);
                *self.h_rc_buffer.get_unchecked_mut(pos) = (h_seqk, rc_seqk);
                *self.idx_buffer. get_unchecked_mut(pos) = self.current_idx_plus_k;
                self.buffer_pos += 1;

                hash = H::min(self.rh, self.fh);
                if hash <= self.hash_bound { break; }
            }

        prev_current_idx = *self.idx_buffer.get_unchecked(std::intrinsics::unchecked_rem(self.buffer_pos+BUFLEN-1,BUFLEN));
        Some((prev_current_idx, hash))
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size_hint = (self.seq_len - self.k + 1) * (((H::max_value() as FH) / (self.hash_bound as FH)) as H) as usize; // rough estimation
        (size_hint, Some(size_hint)) 
    }
}
