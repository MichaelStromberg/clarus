//! RefSeq cache creation pipeline.

use std::sync::Arc;

use anyhow::Result;

use crate::cache::writer::DataSourceVersion;
use crate::cli;
use crate::config::RefSeqCacheConfig;
use crate::evaluation;
use crate::genbank;
use crate::gff3;
use crate::sequence::{ProteinSequences, RnaSequences};
use crate::strand::Strand;
use crate::transcript::construction::{build_transcript_regions, detect_coding_region};
use crate::transcript::types::{IntermediateTranscript, Source};

use super::{
    PipelineContext, build_gene, inject_mt_cdna, local_path, open_file, parse_release_date,
    print_designation_summary, print_eval_stats, print_transcript_stats, write_cache_files,
};

/// Run the RefSeq transcript construction and cache writing pipeline.
#[allow(clippy::too_many_lines)]
pub fn run(ctx: &mut PipelineContext<'_>, config: &RefSeqCacheConfig) -> Result<()> {
    cli::section("RefSeq Parsing");

    let gff3_path = local_path(ctx.base_dir, "refseq", &config.gff3.url)?;
    let genbank_path = local_path(ctx.base_dir, "refseq", &config.genbank.url)?;
    let rna_path = local_path(ctx.base_dir, "refseq", &config.rna_fasta.url)?;
    let protein_path = local_path(ctx.base_dir, "refseq", &config.protein_fasta.url)?;

    // 1. Parse GFF3
    let gff3_result =
        gff3::parse_refseq_gff3_gz(open_file(&gff3_path)?, &ctx.reference.name_to_index)?;
    cli::kv(
        "GFF3",
        &format!(
            "{} genes, {} regulatory regions",
            cli::num(gff3_result.genes.len()),
            cli::num(gff3_result.regulatory_regions.len())
        ),
    );

    // 2. Parse GenBank
    let genbank_records = genbank::parse_genbank_gz(open_file(&genbank_path)?)?;
    cli::kv(
        "GenBank",
        &format!("{} CDS records", cli::num(genbank_records.len())),
    );

    // 3. Parse RNA FASTA
    let mut rna_sequences = RnaSequences::from_gz(open_file(&rna_path)?)?;
    cli::kv(
        "RNA FASTA",
        &format!("{} sequences", cli::num(rna_sequences.len())),
    );

    // 4. Parse Protein FASTA
    let protein_sequences = ProteinSequences::from_gz(open_file(&protein_path)?)?;
    cli::kv(
        "Protein FASTA",
        &format!("{} sequences", cli::num(protein_sequences.len())),
    );

    // 5. MT cDNA injection
    let mt_count = inject_mt_cdna(
        &mut rna_sequences,
        &gff3_result.genes,
        ctx.reference,
        ctx.ref_reader,
        |entry| {
            entry
                .attributes
                .name
                .as_deref()
                .unwrap_or(&entry.attributes.id)
                .to_string()
        },
    )?;
    if mt_count > 0 {
        cli::kv("MT cDNA injected", &cli::num(mt_count));
    }

    eprintln!();

    // ── Transcript Construction ──────────────────────────
    cli::section("RefSeq Transcript Construction");

    let mut transcripts: Vec<IntermediateTranscript> = Vec::new();
    let mut skipped_no_rna: u32 = 0;

    for gene_record in &gff3_result.genes {
        let gene_entry = &gene_record.entry;
        let on_reverse_strand = gene_entry.strand == Strand::Reverse;

        // Gene symbol resolution with HGNC overrides
        let gff3_symbol = gene_entry
            .attributes
            .gene_symbol
            .clone()
            .unwrap_or_default();
        let ncbi_gene_id = gene_entry.attributes.gene_id.clone();

        let (symbol, hgnc_id, ensembl_id) = if let Some(ref nid) = ncbi_gene_id {
            if let Some(hgnc_entry) = ctx.hgnc_db.get_by_ncbi_id(nid) {
                (
                    hgnc_entry.symbol.clone(),
                    Some(hgnc_entry.hgnc_id),
                    hgnc_entry.ensembl_gene_id.clone(),
                )
            } else {
                (gff3_symbol, gene_entry.attributes.hgnc_id, None)
            }
        } else {
            (gff3_symbol, gene_entry.attributes.hgnc_id, None)
        };

        let gene = build_gene(
            gene_entry.chromosome_index,
            symbol,
            ncbi_gene_id,
            ensembl_id,
            hgnc_id,
            on_reverse_strand,
        );

        for tx_record in &gene_record.transcripts {
            let tx_entry = &tx_record.entry;
            let tx_name = tx_entry
                .attributes
                .name
                .as_deref()
                .unwrap_or(&tx_entry.attributes.id);

            let Some(cdna_seq) = rna_sequences.remove(tx_name) else {
                skipped_no_rna += 1;
                continue;
            };

            let regions =
                build_transcript_regions(&tx_record.exons, &tx_record.matches, tx_entry.strand)?;

            let is_coding = tx_entry.biotype.is_coding() && !tx_record.cds_entries.is_empty();
            let coding_region = if is_coding {
                let protein_id = tx_record
                    .cds_entries
                    .first()
                    .and_then(|cds| cds.attributes.name.as_deref())
                    .or_else(|| protein_sequences.protein_id_for_transcript(tx_name))
                    .unwrap_or("");

                if protein_id.is_empty() {
                    None
                } else {
                    let protein_seq = protein_sequences
                        .get_by_protein_id(protein_id)
                        .map(<[u8]>::to_vec)
                        .unwrap_or_default();

                    let genbank_cds = genbank_records.get(tx_name);

                    match detect_coding_region(
                        &tx_record.cds_entries,
                        &regions,
                        tx_entry.strand,
                        protein_id,
                        protein_seq,
                        genbank_cds,
                    ) {
                        Ok(cr) => Some(cr),
                        Err(e) => {
                            cli::warning(&format!(
                                "coding region detection failed for {tx_name}: {e}"
                            ));
                            None
                        }
                    }
                }
            } else {
                None
            };

            transcripts.push(IntermediateTranscript {
                chromosome_index: tx_entry.chromosome_index,
                start: tx_entry.start,
                end: tx_entry.end,
                id: tx_name.to_string(),
                biotype: tx_entry.biotype,
                source: Source::RefSeq,
                strand: tx_entry.strand,
                gene: Arc::clone(&gene),
                transcript_regions: regions,
                coding_region,
                cdna_seq,
                designations: tx_entry.attributes.designations,
            });
        }
    }

    print_transcript_stats(&transcripts, skipped_no_rna);

    eprintln!();

    // ── Evaluation ───────────────────────────────────────
    cli::section("RefSeq Evaluation");

    let eval_stats =
        evaluation::evaluate_transcripts(&mut transcripts, &ctx.reference.chromosomes)?;
    print_eval_stats("RefSeq", &eval_stats)?;

    eprintln!();

    // ── Designation Summary ──────────────────────────────
    print_designation_summary("RefSeq", &mut transcripts);

    eprintln!();

    // ── Write Cache Files ────────────────────────────────
    let (year, month, day) = parse_release_date(&config.release_date)?;
    write_cache_files(
        &DataSourceVersion {
            name: &config.name,
            version: &config.version,
            release_year: year,
            release_month: month,
            release_day: day,
        },
        ctx.assembly,
        ctx.reference.reference_id,
        &ctx.reference.chromosomes,
        &transcripts,
        &gff3_result.regulatory_regions,
        ctx.out_dir,
    )?;

    Ok(())
}
