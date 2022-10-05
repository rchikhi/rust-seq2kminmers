#![feature(bench_black_box)]
use flate2::read::GzDecoder;
use std::io::{BufRead, BufReader};
use lzzzz::lz4f::{BufReadDecompressor};
use std::time::Instant;
use rust_seq2kminmers::{KminmersIterator};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};

fn get_reader(filename_str: &str) -> Box<dyn BufRead + Send> {
    let mut filetype = "unzip";
    let file = match std::fs::File::open(filename_str) {
            Ok(file) => file,
            Err(error) => panic!("Error opening compressed file: {:?}.", error),
        };
    if filename_str.ends_with(".gz")  {filetype = "zip";}
    if filename_str.ends_with(".lz4") {filetype = "lz4";}
    let reader :Box<dyn BufRead + Send> = match filetype { 
        "zip" => Box::new(BufReader::new(GzDecoder::new(file))), 
        "lz4" => Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),
        _ =>     Box::new(BufReader::new(file)), 
    }; 
    reader
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    
    // A simple example given just a sequence in a string
    if args.len() < 2
    {
        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
        println!("Demonstrating how to construct k-min-mers (k=10, l=5, d=0.1) out of a test sequence: {}",std::str::from_utf8(seq).unwrap());
        let iter = KminmersIterator::new(seq, 10, 5, 0.1, true).unwrap();
        for kminmer in iter
        {
            println!("kminmer: {:?}",kminmer);
        }
    }


    // A complete example with FASTA parsing
    else
    {
        let k = 7;
        let l = 31;
        let d = 0.01;
        let queue_len = 200; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
                             // also: controls how many reads objects are buffered during fasta/fastq
                             // parsing
        let filename = &args[1];
        let nb_threads = args[2].parse::<i32>().unwrap();
        println!("Enumerating k-min-mers for the input file {} in parallel ({} threads)", filename, nb_threads);

        let mut is_fasta = false;
        if filename.contains(".fasta.") || filename.contains(".fa.") || filename.ends_with(".fa") || filename.ends_with(".fasta") { // not so robust
            is_fasta = true;
        }

        let process_seq_aux = |seq_str: &[u8], _seq_id: &str| -> Option<u64> {
            let iter = KminmersIterator::new(seq_str, k, l, d, true).unwrap();
            for kminmer in iter
            {
                //println!("seq_id {} kminmer: {:?}",seq_id,kminmer);
                
                // prevent from being optimized away
                std::hint::black_box(&kminmer);
            }
            return Some(1)
        };

        let process_seq_fasta = |record: seq_io::fasta::RefRecord, found: &mut Option<u64>| {
            let seq_str = record.seq().to_vec(); 
            let seq_id = record.id().unwrap().to_string();
            *found = process_seq_aux(&seq_str, &seq_id);

        };
        let process_seq_fastq = |record: seq_io::fastq::RefRecord, found: &mut Option<u64>| {
            let seq_str = record.seq(); 
            let seq_id = record.id().unwrap().to_string();
            *found = process_seq_aux(&seq_str, &seq_id);
        };
        let main_thread = |_found: &mut Option<u64>| { // runs in main thread
            None::<()>
        };

        let start = Instant::now();
        let buf = get_reader(&filename);
        if is_fasta {
            let reader = seq_io::fasta::Reader::new(buf);
            let _res = read_process_fasta_records(reader, nb_threads as u32, queue_len, process_seq_fasta, |_record, found| {main_thread(found)});
        }
        else {
            let reader = seq_io::fastq::Reader::new(buf);
            let _res = read_process_fastq_records(reader, nb_threads as u32, queue_len, process_seq_fastq, |_record, found| {main_thread(found)});
        }
        let duration = start.elapsed();
        println!("FASTA to kminmers in {:?}.", duration);

    }
}
