//! GFF3 file parser: converts GFF3 annotations into gene-transcript-exon hierarchies.

pub mod ensembl_hierarchy;
pub mod ensembl_parser;
pub mod entry;
pub mod refseq_hierarchy;
pub mod refseq_parser;

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use flate2::read::GzDecoder;

use crate::biotype::BioTypeCategory;
use crate::error::Error;

use entry::{Gff3Entry, Gff3Result, ParsedLine};
use refseq_hierarchy::RefSeqHierarchyBuilder;

/// Parse a gzip-compressed RefSeq GFF3 file into a structured hierarchy.
pub fn parse_refseq_gff3_gz<R: Read>(
    reader: R,
    name_to_index: &HashMap<String, usize>,
) -> Result<Gff3Result, Error> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    parse_refseq_gff3(buf_reader, name_to_index)
}

/// Parse RefSeq GFF3 from a buffered reader.
pub fn parse_refseq_gff3<R: BufRead>(
    reader: R,
    name_to_index: &HashMap<String, usize>,
) -> Result<Gff3Result, Error> {
    let mut builder = RefSeqHierarchyBuilder::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line_num = line_num + 1;
        let line = line?;
        match refseq_parser::parse_line(&line, name_to_index)
            .map_err(|e| Error::Parse(format!("{e} (line {line_num}: {line})")))?
        {
            ParsedLine::Entry(entry) => builder.add(*entry)?,
            ParsedLine::Discarded | ParsedLine::Comment => {}
            ParsedLine::EndOfSection => break,
        }
    }

    builder.build()
}

/// Parse a gzip-compressed Ensembl GFF3 file into a structured hierarchy.
pub fn parse_ensembl_gff3_gz<R: Read>(
    reader: R,
    name_to_index: &HashMap<String, usize>,
) -> Result<Gff3Result, Error> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    parse_ensembl_gff3(buf_reader, name_to_index)
}

/// Parse Ensembl GFF3 from a buffered reader.
pub fn parse_ensembl_gff3<R: BufRead>(
    reader: R,
    name_to_index: &HashMap<String, usize>,
) -> Result<Gff3Result, Error> {
    let mut builder = ensembl_hierarchy::EnsemblHierarchyBuilder::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line_num = line_num + 1;
        let line = line?;
        match ensembl_parser::parse_line(&line, name_to_index)
            .map_err(|e| Error::Parse(format!("{e} (line {line_num}: {line})")))?
        {
            ParsedLine::Entry(entry) => builder.add(*entry)?,
            ParsedLine::Discarded | ParsedLine::Comment | ParsedLine::EndOfSection => {}
        }
    }

    builder.build()
}

/// Parse a gzip-compressed Ensembl regulatory GFF3 file.
///
/// Regulatory regions are flat intervals (no hierarchy). Returns entries
/// where the biotype category is `Regulatory`, sorted by chromosome/start/end.
pub fn parse_regulatory_gff3_gz<R: Read>(
    reader: R,
    name_to_index: &HashMap<String, usize>,
) -> Result<Vec<Gff3Entry>, Error> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    let mut regions = Vec::new();

    for (line_num, line) in buf_reader.lines().enumerate() {
        let line_num = line_num + 1;
        let line = line?;
        match ensembl_parser::parse_line(&line, name_to_index)
            .map_err(|e| Error::Parse(format!("{e} (line {line_num}: {line})")))?
        {
            ParsedLine::Entry(entry) => {
                if entry.biotype.category() == BioTypeCategory::Regulatory {
                    regions.push(*entry);
                }
            }
            ParsedLine::Discarded | ParsedLine::Comment | ParsedLine::EndOfSection => {}
        }
    }

    regions.sort_by(|a, b| {
        a.chromosome_index
            .cmp(&b.chromosome_index)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    Ok(regions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn refseq_name_map() -> HashMap<String, usize> {
        let mut m = HashMap::new();
        m.insert("NC_000001.11".to_string(), 0);
        m
    }

    #[test]
    fn parse_refseq_worked_example() {
        let gff3 = "\
##gff-version 3
NC_000001.11\tBestRefSeq\tgene\t11874\t14409\t.\t+\t.\tID=gene-DDX11L1;Dbxref=GeneID:100287102,HGNC:HGNC:37102;gene=DDX11L1
NC_000001.11\tBestRefSeq\tmRNA\t11874\t14409\t.\t+\t.\tID=rna-NR_046018.2;Parent=gene-DDX11L1;Name=NR_046018.2
NC_000001.11\tBestRefSeq\texon\t11874\t12227\t.\t+\t.\tID=exon-NR_046018.2-1;Parent=rna-NR_046018.2;transcript_id=NR_046018.2
NC_000001.11\tBestRefSeq\texon\t12613\t12721\t.\t+\t.\tID=exon-NR_046018.2-2;Parent=rna-NR_046018.2;transcript_id=NR_046018.2
###";

        let reader = Cursor::new(gff3.as_bytes());
        let result =
            parse_refseq_gff3(std::io::BufReader::new(reader), &refseq_name_map()).unwrap();

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

    fn ensembl_name_map() -> HashMap<String, usize> {
        let mut m = HashMap::new();
        m.insert("17".to_string(), 16);
        m.insert("1".to_string(), 0);
        m
    }

    #[test]
    fn parse_ensembl_worked_example() {
        let gff3 = "\
##gff-version 3
1\tensembl\tgene\t69091\t70008\t.\t+\t.\tID=gene:ENSG00000186092;Name=OR4F5;gene_id=ENSG00000186092;biotype=protein_coding;version=14
1\tensembl\tmRNA\t69091\t70008\t.\t+\t.\tID=transcript:ENST00000335137;Parent=gene:ENSG00000186092;Name=OR4F5-201;transcript_id=ENST00000335137;biotype=protein_coding;version=4
1\tensembl\texon\t69091\t70008\t.\t+\t.\tParent=transcript:ENST00000335137;exon_id=ENSE00002319515;version=1;ID=exon:ENSE00002319515
1\tensembl\tCDS\t69091\t69984\t.\t+\t0\tParent=transcript:ENST00000335137;protein_id=ENSP00000334393;ID=CDS:ENSP00000334393
###";

        let reader = Cursor::new(gff3.as_bytes());
        let result =
            parse_ensembl_gff3(std::io::BufReader::new(reader), &ensembl_name_map()).unwrap();

        assert_eq!(result.genes.len(), 1);
        let gene = &result.genes[0];
        assert_eq!(gene.entry.attributes.name.as_deref(), Some("OR4F5"));
        assert_eq!(
            gene.entry.attributes.gene_id.as_deref(),
            Some("ENSG00000186092")
        );
        assert_eq!(gene.transcripts.len(), 1);

        let tx = &gene.transcripts[0];
        assert_eq!(
            tx.entry.attributes.transcript_id.as_deref(),
            Some("ENST00000335137")
        );
        assert_eq!(tx.entry.attributes.version, Some(4));
        assert_eq!(tx.exons.len(), 1);
        assert_eq!(tx.cds_entries.len(), 1);
        assert_eq!(
            tx.cds_entries[0].attributes.protein_id.as_deref(),
            Some("ENSP00000334393")
        );
    }

    #[test]
    fn ensembl_artifact_transcript_skipped() {
        let gff3 = "\
1\tensembl\tgene\t100\t500\t.\t+\t.\tID=gene:ENSG00000000001;Name=TEST;gene_id=ENSG00000000001;biotype=protein_coding;version=1
1\tensembl\tmRNA\t100\t500\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000000001;transcript_id=ENST00000000001;biotype=artifact;version=1
###";
        let reader = Cursor::new(gff3.as_bytes());
        let result =
            parse_ensembl_gff3(std::io::BufReader::new(reader), &ensembl_name_map()).unwrap();
        assert_eq!(result.genes[0].transcripts.len(), 0);
    }

    #[test]
    fn ensembl_unconfirmed_transcript_skipped() {
        let gff3 = "\
1\tensembl\tgene\t100\t500\t.\t+\t.\tID=gene:ENSG00000000001;Name=TEST;gene_id=ENSG00000000001;biotype=protein_coding;version=1
1\tensembl\tmRNA\t100\t500\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000000001;transcript_id=ENST00000000001;biotype=unconfirmed_transcript;version=1
###";
        let reader = Cursor::new(gff3.as_bytes());
        let result =
            parse_ensembl_gff3(std::io::BufReader::new(reader), &ensembl_name_map()).unwrap();
        assert_eq!(result.genes[0].transcripts.len(), 0);
    }
}
