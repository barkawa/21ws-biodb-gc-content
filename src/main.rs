use std::fs::File;
use bio::io::fasta;
use plotters::prelude::*;


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

fn get_total_gc_content(record: &fasta::Record) -> f64 {
    let count = count_bases(record);
    
    let at = (count.a + count.t) as f64;
    let gc = (count.g + count.c) as f64;

    gc / (at + gc)
}


struct SlidingWindowAverage<'a> {
    data: &'a [u8],
    sum: u64,
    idx: usize,
    step: usize,
    size: usize,
}

impl SlidingWindowAverage<'_> {
    fn new(data: &[u8], size: usize, step: usize) -> SlidingWindowAverage {
        let mut first_window_sum = 0;

        for i in 0..size {
            first_window_sum += match data[i] {
                b'G' | b'C' => 1,
                _ => 0,
            };
        }

        SlidingWindowAverage { 
            data, 
            sum: first_window_sum, 
            idx: 0,
            step,
            size
        }
    }

}

impl Iterator for SlidingWindowAverage<'_> {
    type Item = f32;

    // its 2AM please dont judge 
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx + self.size + self.step < self.data.len() - 1 {
            // Subtract first n=step items
            for i in 0..self.step {
                self.sum -= match self.data[self.idx + i] {
                    b'G' | b'C' => 1,
                    _ => 0,
                };
            }

            // Add next n=step items
            for i in 0..self.step {
                self.sum += match self.data[self.idx + self.size + i] {
                    b'G' | b'C' => 1,
                    _ => 0,
                };
            }

            self.idx += self.step;

            // Calculate average
            Some(self.sum as f32 / self.size as f32)
        } else {
            None
        }
    }

}


fn plot(filename: &str, title: &str, record: &fasta::Record) {
    const SIZE: usize = 100000;
    const STEP: usize = 10000;

    let mut filename = filename.to_string();
    filename.push_str(".svg");

    let root_area = SVGBackend::new(
        &filename, (2000, 500)
        ).into_drawing_area();
    
    root_area.fill(&WHITE).unwrap();

    let mut context = ChartBuilder::on(&root_area)
        .margin(20)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(title, ("sans-serif", 20))
        .build_cartesian_2d(0..(record.seq().len()), 0f32..1f32)
        .unwrap();
    
    context
        .configure_mesh()
        .disable_mesh()
        .x_label_formatter(&|l| format!("{}M", l / 1000000))
        .x_labels(30)
        .draw()
        .unwrap();

    let window_avg_iter = SlidingWindowAverage::new(record.seq(), SIZE, STEP);
    context.draw_series(
        AreaSeries::new(
            (0..).step_by(STEP).zip(window_avg_iter),
            0.0,
            &BLACK.mix(0.1)
        ).border_style(&BLACK)
    ).unwrap();
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
                // println!("  - GC Content: {:.2}%", get_gc_content(&record) * 100.);
                plot(record.id(), record.desc().unwrap_or_default(), &record);
            }
        },
        Err(err) => {
            eprintln!("Error parsing FASTA file: {}", err);
            return;
        }
    }

    
}