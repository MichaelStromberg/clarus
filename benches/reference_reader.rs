use std::fs::File;
use std::io::BufReader;

use criterion::{Criterion, criterion_group, criterion_main};

use clarus::reference::reader::ReferenceReader;

const REFERENCE_PATH: &str = "data/reference/GRCh38.p14.dat";

fn bench_read_header(c: &mut Criterion) {
    c.bench_function("read_header (705 chromosomes)", |b| {
        b.iter(|| {
            let file = File::open(REFERENCE_PATH).unwrap();
            let mut reader = BufReader::new(file);
            let ref_reader = ReferenceReader::from_reader(&mut reader).unwrap();
            assert_eq!(ref_reader.chromosomes.len(), 705);
        });
    });
}

fn bench_load_chr1(c: &mut Criterion) {
    // Read header once, then benchmark loading chr1 repeatedly
    let file = File::open(REFERENCE_PATH).unwrap();
    let mut reader = BufReader::new(file);
    let ref_reader = ReferenceReader::from_reader(&mut reader).unwrap();
    let chr1_idx = ref_reader.get_index("chr1").unwrap();

    c.bench_function("load_chr1 (249 MB decompressed)", |b| {
        b.iter(|| {
            let file = File::open(REFERENCE_PATH).unwrap();
            let mut reader = BufReader::new(file);
            let _ = ReferenceReader::from_reader(&mut reader).unwrap();
            let data = ref_reader.load_chromosome(&mut reader, chr1_idx).unwrap();
            assert_eq!(data.sequence.len(), 248_956_422);
        });
    });
}

criterion_group!(benches, bench_read_header, bench_load_chr1);
criterion_main!(benches);
