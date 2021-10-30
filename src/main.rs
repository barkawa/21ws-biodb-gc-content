extern crate clap;
extern crate bio;
use std::fs::File;

use bio::io::fasta;

#[derive(Debug, Default)]
struct BaseCount {
    a: u64,
    c: u64,
    t: u64,
    g: u64,
    other: u64,
}

fn count_bases(record: &fasta::Record) -> BaseCount {
    let mut count = BaseCount::default();

    for base in record.seq() {
        match base {
            b'A' => count.a += 1,
            b'C' => count.c += 1,
            b'T' => count.t += 1,
            b'G' => count.g += 1,
            _    => count.other += 1
        }
    }

    count
}

fn get_gc_content(record: &fasta::Record) -> f64 {
    let count = count_bases(record);
    
    let at = (count.a + count.t) as f64;
    let gc = (count.g + count.c) as f64;

    gc / (at + gc)
}

fn main() {
    let matches = clap::App::new("gc-content")
        .version("0.1")
        .about("Analyzes the GC-Content of a genome")
        .arg(clap::Arg::with_name("FILE")
            .help("Input file in FASTA format")
            .required(true)
            .index(1))
        .get_matches();
    
    let filename = matches.value_of("FILE").unwrap();
    
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(err) => {
            eprintln!("Couldn't open {}: {}", filename, err);
            return;
        }
    };
    
    let fasta_reader = fasta::Reader::new(file);

    match fasta_reader.records().collect::<Result<Vec<_>, _>>() {
        Ok(records) => {
            for record in records {
                println!("[{}] {}", record.id(), record.desc().unwrap_or_default());
                println!("  - GC Content: {:.2}%", get_gc_content(&record) * 100.);
            }
        },
        Err(err) => {
            eprintln!("Error parsing FASTA file: {}", err);
            return;
        }
    }
    
}