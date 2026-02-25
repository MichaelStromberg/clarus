//! Ensembl cache creation pipeline.

use std::sync::Arc;

use anyhow::Result;

use crate::cache::writer::DataSourceVersion;
use crate::cli;
use crate::config::EnsemblConfig;
use crate::evaluation;
use crate::gff3;
use crate::gff3::entry::Gff3Entry;
use crate::sequence::{ProteinSequences, RnaSequences};
use crate::strand::Strand;
use crate::transcript::construction::{build_transcript_regions, detect_coding_region};
use crate::transcript::types::{IntermediateTranscript, Source};

use super::{
    PipelineContext, build_gene, inject_mt_cdna, local_path, open_file, parse_release_date,
    print_designation_summary, print_eval_stats, print_transcript_stats, write_cache_files,
};

/// Run the Ensembl transcript construction and cache writing pipeline.
#[allow(clippy::too_many_lines)]
pub fn run(ctx: &mut PipelineContext<'_>, config: &EnsemblConfig) -> Result<()> {
    cli::section("Ensembl Parsing");

    let gff3_path = local_path(ctx.base_dir, "ensembl", &config.gff3.url)?;
    let reg_gff_path = local_path(ctx.base_dir, "ensembl", &config.regulatory_gff.url)?;
    let cdna_path = local_path(ctx.base_dir, "ensembl", &config.cdna_fasta.url)?;
    let ncrna_path = local_path(ctx.base_dir, "ensembl", &config.ncrna_fasta.url)?;
    let peptide_path = local_path(ctx.base_dir, "ensembl", &config.peptide_fasta.url)?;

    // 1. Parse main GFF3
    let gff3_result =
        gff3::parse_ensembl_gff3_gz(open_file(&gff3_path)?, &ctx.reference.name_to_index)?;
    cli::kv(
        "GFF3",
        &format!(
            "{} genes, {} regulatory regions",
            cli::num(gff3_result.genes.len()),
            cli::num(gff3_result.regulatory_regions.len())
        ),
    );

    // 2. Parse regulatory GFF3
    let regulatory_regions =
        gff3::parse_regulatory_gff3_gz(open_file(&reg_gff_path)?, &ctx.reference.name_to_index)?;
    cli::kv(
        "Regulatory GFF",
        &format!("{} regions", cli::num(regulatory_regions.len())),
    );

    // 3. Parse RNA FASTA (merge cdna + ncrna)
    let mut rna_sequences = RnaSequences::from_gz(open_file(&cdna_path)?)?;
    let cdna_count = rna_sequences.len();

    let ncrna_sequences = RnaSequences::from_gz(open_file(&ncrna_path)?)?;
    let ncrna_count = ncrna_sequences.len();
    rna_sequences.merge(ncrna_sequences)?;
    cli::kv(
        "RNA FASTA",
        &format!(
            "{} sequences (cdna: {}, ncrna: {})",
            cli::num(rna_sequences.len()),
            cli::num(cdna_count),
            cli::num(ncrna_count)
        ),
    );

    // 4. Parse Protein FASTA
    let protein_sequences = ProteinSequences::from_gz(open_file(&peptide_path)?)?;
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
            let stable_id = entry
                .attributes
                .transcript_id
                .as_deref()
                .unwrap_or_else(|| {
                    entry
                        .attributes
                        .id
                        .find(':')
                        .map_or(entry.attributes.id.as_str(), |pos| {
                            &entry.attributes.id[pos + 1..]
                        })
                });
            let version = entry.attributes.version.unwrap_or(1);
            format!("{stable_id}.{version}")
        },
    )?;
    if mt_count > 0 {
        cli::kv("MT cDNA injected", &cli::num(mt_count));
    }

    eprintln!();

    // ── Transcript Construction ──────────────────────────
    cli::section("Ensembl Transcript Construction");

    let mut transcripts: Vec<IntermediateTranscript> = Vec::new();
    let mut skipped_no_rna: u32 = 0;

    for gene_record in &gff3_result.genes {
        let gene_entry = &gene_record.entry;
        let on_reverse_strand = gene_entry.strand == Strand::Reverse;

        let ensembl_gene_id = gene_entry.attributes.gene_id.clone();

        // Gene symbol resolution: Ensembl uses Name attribute, falls back to gene_id
        let gff3_symbol = gene_entry
            .attributes
            .name
            .clone()
            .or_else(|| ensembl_gene_id.clone())
            .unwrap_or_default();

        // HGNC lookup by Ensembl gene ID
        let (ncbi_gene_id, hgnc_id) = if let Some(ref eid) = ensembl_gene_id {
            if let Some(hgnc_entry) = ctx.hgnc_db.get_by_ensembl_id(eid) {
                (hgnc_entry.entrez_id.clone(), Some(hgnc_entry.hgnc_id))
            } else {
                (None, None)
            }
        } else {
            (None, None)
        };

        let gene = build_gene(
            gene_entry.chromosome_index,
            gff3_symbol,
            ncbi_gene_id,
            ensembl_gene_id,
            hgnc_id,
            on_reverse_strand,
        );

        for tx_record in &gene_record.transcripts {
            let tx_entry = &tx_record.entry;

            let tx_stable_id = tx_entry
                .attributes
                .transcript_id
                .as_deref()
                .unwrap_or_else(|| {
                    tx_entry
                        .attributes
                        .id
                        .find(':')
                        .map_or(tx_entry.attributes.id.as_str(), |pos| {
                            &tx_entry.attributes.id[pos + 1..]
                        })
                });
            let version = tx_entry.attributes.version.unwrap_or(1);
            let versioned_id = format!("{tx_stable_id}.{version}");

            let Some(cdna_seq) = rna_sequences.remove(&versioned_id) else {
                skipped_no_rna += 1;
                continue;
            };

            let regions =
                build_transcript_regions(&tx_record.exons, &tx_record.matches, tx_entry.strand)?;

            let is_coding = tx_entry.biotype.is_coding() && !tx_record.cds_entries.is_empty();
            let coding_region = if is_coding {
                let protein_id = protein_sequences
                    .protein_id_for_transcript(&versioned_id)
                    .or_else(|| {
                        tx_record
                            .cds_entries
                            .first()
                            .and_then(|cds| cds.attributes.protein_id.as_deref())
                    })
                    .unwrap_or("");

                if protein_id.is_empty() {
                    None
                } else {
                    let protein_seq = protein_sequences
                        .get_by_protein_id(protein_id)
                        .map(<[u8]>::to_vec)
                        .unwrap_or_default();

                    match detect_coding_region(
                        &tx_record.cds_entries,
                        &regions,
                        tx_entry.strand,
                        protein_id,
                        protein_seq,
                        None,
                    ) {
                        Ok(cr) => Some(cr),
                        Err(e) => {
                            cli::warning(&format!(
                                "coding region detection failed for {versioned_id}: {e}"
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
                id: versioned_id,
                biotype: tx_entry.biotype,
                source: Source::Ensembl,
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
    cli::section("Ensembl Evaluation");

    let eval_stats =
        evaluation::evaluate_transcripts(&mut transcripts, &ctx.reference.chromosomes)?;
    print_eval_stats("Ensembl", &eval_stats)?;

    eprintln!();

    // ── Designation Summary ──────────────────────────────
    print_designation_summary("Ensembl", &mut transcripts);

    eprintln!();

    // ── Write Cache Files ────────────────────────────────
    // Combine regulatory regions from main GFF3 and regulatory GFF3
    let mut all_regulatory: Vec<Gff3Entry> = gff3_result.regulatory_regions;
    all_regulatory.extend(regulatory_regions);
    all_regulatory.sort_by(|a, b| {
        a.chromosome_index
            .cmp(&b.chromosome_index)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

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
        &all_regulatory,
        ctx.out_dir,
    )?;

    Ok(())
}
