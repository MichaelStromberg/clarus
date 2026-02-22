//! GFF3 file parser: converts GFF3 annotations into gene-transcript-exon hierarchies.

pub mod entry;
pub mod hierarchy;
pub mod parser;

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use flate2::read::GzDecoder;

use crate::error::Error;

use entry::Gff3Result;
use hierarchy::HierarchyBuilder;
use parser::ParsedLine;

/// Parse a gzip-compressed GFF3 file into a structured hierarchy.
pub fn parse_gff3_gz<R: Read>(
    reader: R,
    name_to_index: &HashMap<String, usize>,
) -> Result<Gff3Result, Error> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    parse_gff3(buf_reader, name_to_index)
}

/// Parse GFF3 from a buffered reader.
pub fn parse_gff3<R: BufRead>(
    reader: R,
    name_to_index: &HashMap<String, usize>,
) -> Result<Gff3Result, Error> {
    let mut builder = HierarchyBuilder::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line_num = line_num + 1;
        let line = line?;
        match parser::parse_line(&line, name_to_index)
            .map_err(|e| Error::Parse(format!("{e} (line {line_num}: {line})")))?
        {
            ParsedLine::Entry(entry) => builder.add(*entry)?,
            ParsedLine::Discarded | ParsedLine::Comment => continue,
            ParsedLine::EndOfSection => break,
        }
    }

    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn test_name_map() -> HashMap<String, usize> {
        let mut m = HashMap::new();
        m.insert("NC_000001.11".to_string(), 0);
        m
    }

    #[test]
    fn parse_worked_example() {
        let gff3 = "\
##gff-version 3
NC_000001.11\tBestRefSeq\tgene\t11874\t14409\t.\t+\t.\tID=gene-DDX11L1;Dbxref=GeneID:100287102,HGNC:HGNC:37102;gene=DDX11L1
NC_000001.11\tBestRefSeq\tmRNA\t11874\t14409\t.\t+\t.\tID=rna-NR_046018.2;Parent=gene-DDX11L1;Name=NR_046018.2
NC_000001.11\tBestRefSeq\texon\t11874\t12227\t.\t+\t.\tID=exon-NR_046018.2-1;Parent=rna-NR_046018.2;transcript_id=NR_046018.2
NC_000001.11\tBestRefSeq\texon\t12613\t12721\t.\t+\t.\tID=exon-NR_046018.2-2;Parent=rna-NR_046018.2;transcript_id=NR_046018.2
###";

        let reader = Cursor::new(gff3.as_bytes());
        let result = parse_gff3(std::io::BufReader::new(reader), &test_name_map()).unwrap();

        assert_eq!(result.genes.len(), 1);
        let gene = &result.genes[0];
        assert_eq!(gene.entry.start, 11874);
        assert_eq!(gene.entry.end, 14409);
        assert_eq!(gene.transcripts.len(), 1);

        let tx = &gene.transcripts[0];
        assert_eq!(tx.entry.attributes.name.as_deref(), Some("NR_046018.2"));
        assert_eq!(tx.exons.len(), 2);
        assert_eq!(tx.exons[0].start, 11874);
        assert_eq!(tx.exons[1].start, 12613);
    }
}
