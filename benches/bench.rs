// based on https://raw.githubusercontent.com/luizirber/nthash/latest/benches/nthash.rs written by Luiz Irber
// still benchmarks the original crate, but..
// extended to benchmark this crate too,
// and also.. the C++ NtHash2 version!

#[macro_use]
extern crate criterion;

use criterion::{Bencher, Criterion, Fun};
use rand::distributions::{Distribution, Uniform};

use nthash::{nthash, NtHashIterator};
use rust_seq2kminmers::{KminmersIterator,Kminmer, NtHashHPCIterator};

fn nthash_bench(c: &mut Criterion) {
    let range = Uniform::from(0..4);
    let mut rng = rand::thread_rng();
    let seq = (0..10000)
        .map(|_| match range.sample(&mut rng) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        })
        .collect::<String>();

    let nthash_orig_it = Fun::new("nthash_orig_iterator", |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let iter = NtHashIterator::new(i.as_bytes(), 5).unwrap();
            //  iter.for_each(drop);
            let _res = iter.collect::<Vec<u64>>();
        })
    });

    let nthash_orig_simple = Fun::new("nthash_orig_simple", |b: &mut Bencher, i: &String| {
        b.iter(|| {
            nthash(i.as_bytes(), 5);
        })
    });

   let nthash_new_it_hpc = Fun::new("nthash_new_iterator_hpc", |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let density : f64 = 1.0;
            let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;
            let iter = NtHashHPCIterator::new(i.as_bytes(), 5, hash_bound).unwrap();
            let _res = iter.collect::<Vec<(usize, u64)>>();
        })
    });

    let kminmers = Fun::new("kminmers", |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let iter = KminmersIterator::new(i.as_bytes(), 10, 5, 0.01, false).unwrap();
            let _res = iter.collect::<Vec<Kminmer>>();
        })
    });

   let kminmers_hpc = Fun::new("kminmers_hpc", |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let iter = KminmersIterator::new(i.as_bytes(), 10, 5, 0.01, true).unwrap();
            let _res = iter.collect::<Vec<Kminmer>>();
        })
    });

    let functions = vec![nthash_orig_it, nthash_orig_simple, nthash_new_it_hpc, kminmers, kminmers_hpc];
    c.bench_functions("nthash", functions, seq);
}

criterion_group!(benches, nthash_bench);
criterion_main!(benches);
