//! BioType enum covering all 102 Sequence Ontology terms from GFF3 files.

use std::fmt;

use crate::error::Error;

/// Category that determines how a GFF3 feature is processed during hierarchy aggregation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioTypeCategory {
    Gene,
    Transcript,
    Regulatory,
    ChildFeature,
    AlignmentMatch,
    NotUseful,
    Unsupported,
}

/// All 102 recognized BioType values from GFF3 feature type column 3.
/// Encoded as a single byte for binary cache serialization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum BioType {
    AberrantProcessedTranscript = 0,
    AntisenseRna = 1,
    BiologicalRegion = 2,
    CaatSignal = 3,
    CageCluster = 4,
    CdnaMatch = 5,
    Cds = 6,
    Centromere = 7,
    Chromosome = 8,
    ChromosomeBreakpoint = 9,
    ConservedRegion = 10,
    CtcfBindingSite = 11,
    CGeneSegment = 12,
    DirectRepeat = 13,
    DispersedRepeat = 14,
    DnaseIHypersensitiveSite = 15,
    DGeneSegment = 16,
    DLoop = 17,
    Enhancer = 18,
    EnhancerBlockingElement = 19,
    EpigeneticallyModifiedRegion = 20,
    Exon = 21,
    FivePrimeUtr = 22,
    GcRichPromoterRegion = 23,
    Gene = 24,
    GuideRna = 25,
    ImprintingControlRegion = 26,
    Insulator = 27,
    InvertedRepeat = 28,
    JGeneSegment = 29,
    LincRna = 30,
    LincRnaGene = 31,
    LncRna = 32,
    LocusControlRegion = 33,
    Match = 34,
    MatrixAttachmentSite = 35,
    MeioticRecombinationRegion = 36,
    Microsatellite = 37,
    Minisatellite = 38,
    MiRna = 39,
    MiRnaGene = 40,
    MitoticRecombinationRegion = 41,
    MobileGeneticElement = 42,
    MRna = 43,
    MtGene = 44,
    NcRna = 45,
    NcRnaGene = 46,
    NcPrimaryTranscript = 47,
    NmdTranscriptVariant = 48,
    NonAllelicHomologousRecombinationRegion = 49,
    NucleotideCleavageSite = 50,
    NucleotideMotif = 51,
    OpenChromatinRegion = 52,
    OriginOfReplication = 53,
    PrimaryTranscript = 54,
    ProcessedPseudogene = 55,
    ProcessedTranscript = 56,
    Promoter = 57,
    PromoterFlankingRegion = 58,
    ProteinBindingSite = 59,
    Pseudogene = 60,
    PseudogenicTranscript = 61,
    RecombinationFeature = 62,
    Region = 63,
    RegulatoryRegion = 64,
    RepeatInstabilityRegion = 65,
    RepeatRegion = 66,
    ReplicationRegulatoryRegion = 67,
    ReplicationStartSite = 68,
    ResponseElement = 69,
    Rna = 70,
    RnaseMrpRna = 71,
    RnasePRna = 72,
    RRna = 73,
    RRnaGene = 74,
    Scaffold = 75,
    ScaRna = 76,
    ScRna = 77,
    SequenceAlteration = 78,
    SequenceAlterationArtifact = 79,
    SequenceComparison = 80,
    SequenceFeature = 81,
    SequenceSecondaryStructure = 82,
    Silencer = 83,
    SnoRna = 84,
    SnoRnaGene = 85,
    SnRna = 86,
    SnRnaGene = 87,
    Supercontig = 88,
    TandemRepeat = 89,
    TataBox = 90,
    TelomeraseRna = 91,
    TfBindingSite = 92,
    ThreePrimeUtr = 93,
    Transcript = 94,
    TranscriptionalCisRegulatoryRegion = 95,
    TRna = 96,
    UnconfirmedTranscript = 97,
    VaultRna = 98,
    VdGeneSegment = 99,
    VGeneSegment = 100,
    YRna = 101,
}

impl BioType {
    #[must_use]
    pub fn category(self) -> BioTypeCategory {
        match self {
            // Gene types
            Self::Gene
            | Self::Pseudogene
            | Self::LincRnaGene
            | Self::MiRnaGene
            | Self::MtGene
            | Self::NcRnaGene
            | Self::ProcessedPseudogene
            | Self::RRnaGene
            | Self::SnoRnaGene
            | Self::SnRnaGene => BioTypeCategory::Gene,

            // Transcript types
            Self::AberrantProcessedTranscript
            | Self::AntisenseRna
            | Self::CGeneSegment
            | Self::DGeneSegment
            | Self::GuideRna
            | Self::JGeneSegment
            | Self::LincRna
            | Self::LncRna
            | Self::MiRna
            | Self::MRna
            | Self::NcRna
            | Self::NcPrimaryTranscript
            | Self::NmdTranscriptVariant
            | Self::PrimaryTranscript
            | Self::ProcessedTranscript
            | Self::PseudogenicTranscript
            | Self::Rna
            | Self::RnaseMrpRna
            | Self::RnasePRna
            | Self::RRna
            | Self::ScaRna
            | Self::ScRna
            | Self::SnoRna
            | Self::SnRna
            | Self::TelomeraseRna
            | Self::Transcript
            | Self::TRna
            | Self::UnconfirmedTranscript
            | Self::VaultRna
            | Self::VdGeneSegment
            | Self::VGeneSegment
            | Self::YRna => BioTypeCategory::Transcript,

            // Regulatory types
            Self::CaatSignal
            | Self::CtcfBindingSite
            | Self::DnaseIHypersensitiveSite
            | Self::Enhancer
            | Self::EnhancerBlockingElement
            | Self::EpigeneticallyModifiedRegion
            | Self::GcRichPromoterRegion
            | Self::ImprintingControlRegion
            | Self::Insulator
            | Self::LocusControlRegion
            | Self::MatrixAttachmentSite
            | Self::OpenChromatinRegion
            | Self::Promoter
            | Self::PromoterFlankingRegion
            | Self::RegulatoryRegion
            | Self::ReplicationRegulatoryRegion
            | Self::ResponseElement
            | Self::Silencer
            | Self::TataBox
            | Self::TfBindingSite
            | Self::TranscriptionalCisRegulatoryRegion => BioTypeCategory::Regulatory,

            // Child feature types
            Self::Exon | Self::Cds => BioTypeCategory::ChildFeature,

            // Alignment match
            Self::CdnaMatch => BioTypeCategory::AlignmentMatch,

            // Not useful (silently discarded)
            Self::BiologicalRegion
            | Self::Centromere
            | Self::Chromosome
            | Self::FivePrimeUtr
            | Self::InvertedRepeat
            | Self::Match
            | Self::Region
            | Self::Scaffold
            | Self::SequenceAlterationArtifact
            | Self::Supercontig
            | Self::ThreePrimeUtr => BioTypeCategory::NotUseful,

            // Unsupported
            Self::CageCluster
            | Self::ChromosomeBreakpoint
            | Self::ConservedRegion
            | Self::DirectRepeat
            | Self::DispersedRepeat
            | Self::DLoop
            | Self::MeioticRecombinationRegion
            | Self::Microsatellite
            | Self::Minisatellite
            | Self::MitoticRecombinationRegion
            | Self::MobileGeneticElement
            | Self::NonAllelicHomologousRecombinationRegion
            | Self::NucleotideCleavageSite
            | Self::NucleotideMotif
            | Self::OriginOfReplication
            | Self::ProteinBindingSite
            | Self::RecombinationFeature
            | Self::RepeatInstabilityRegion
            | Self::RepeatRegion
            | Self::ReplicationStartSite
            | Self::SequenceAlteration
            | Self::SequenceComparison
            | Self::SequenceFeature
            | Self::SequenceSecondaryStructure
            | Self::TandemRepeat => BioTypeCategory::Unsupported,
        }
    }

    /// Returns true for coding transcript biotypes (mRNA + immunoglobulin gene segments).
    #[must_use]
    pub fn is_coding(self) -> bool {
        matches!(
            self,
            Self::MRna
                | Self::CGeneSegment
                | Self::DGeneSegment
                | Self::JGeneSegment
                | Self::VGeneSegment
                | Self::VdGeneSegment
        )
    }

    #[must_use]
    pub fn to_byte(self) -> u8 {
        self as u8
    }
}

/// All 102 BioType variants in discriminant order, enabling safe u8 → BioType conversion.
const ALL_BIOTYPES: [BioType; 102] = [
    BioType::AberrantProcessedTranscript,
    BioType::AntisenseRna,
    BioType::BiologicalRegion,
    BioType::CaatSignal,
    BioType::CageCluster,
    BioType::CdnaMatch,
    BioType::Cds,
    BioType::Centromere,
    BioType::Chromosome,
    BioType::ChromosomeBreakpoint,
    BioType::ConservedRegion,
    BioType::CtcfBindingSite,
    BioType::CGeneSegment,
    BioType::DirectRepeat,
    BioType::DispersedRepeat,
    BioType::DnaseIHypersensitiveSite,
    BioType::DGeneSegment,
    BioType::DLoop,
    BioType::Enhancer,
    BioType::EnhancerBlockingElement,
    BioType::EpigeneticallyModifiedRegion,
    BioType::Exon,
    BioType::FivePrimeUtr,
    BioType::GcRichPromoterRegion,
    BioType::Gene,
    BioType::GuideRna,
    BioType::ImprintingControlRegion,
    BioType::Insulator,
    BioType::InvertedRepeat,
    BioType::JGeneSegment,
    BioType::LincRna,
    BioType::LincRnaGene,
    BioType::LncRna,
    BioType::LocusControlRegion,
    BioType::Match,
    BioType::MatrixAttachmentSite,
    BioType::MeioticRecombinationRegion,
    BioType::Microsatellite,
    BioType::Minisatellite,
    BioType::MiRna,
    BioType::MiRnaGene,
    BioType::MitoticRecombinationRegion,
    BioType::MobileGeneticElement,
    BioType::MRna,
    BioType::MtGene,
    BioType::NcRna,
    BioType::NcRnaGene,
    BioType::NcPrimaryTranscript,
    BioType::NmdTranscriptVariant,
    BioType::NonAllelicHomologousRecombinationRegion,
    BioType::NucleotideCleavageSite,
    BioType::NucleotideMotif,
    BioType::OpenChromatinRegion,
    BioType::OriginOfReplication,
    BioType::PrimaryTranscript,
    BioType::ProcessedPseudogene,
    BioType::ProcessedTranscript,
    BioType::Promoter,
    BioType::PromoterFlankingRegion,
    BioType::ProteinBindingSite,
    BioType::Pseudogene,
    BioType::PseudogenicTranscript,
    BioType::RecombinationFeature,
    BioType::Region,
    BioType::RegulatoryRegion,
    BioType::RepeatInstabilityRegion,
    BioType::RepeatRegion,
    BioType::ReplicationRegulatoryRegion,
    BioType::ReplicationStartSite,
    BioType::ResponseElement,
    BioType::Rna,
    BioType::RnaseMrpRna,
    BioType::RnasePRna,
    BioType::RRna,
    BioType::RRnaGene,
    BioType::Scaffold,
    BioType::ScaRna,
    BioType::ScRna,
    BioType::SequenceAlteration,
    BioType::SequenceAlterationArtifact,
    BioType::SequenceComparison,
    BioType::SequenceFeature,
    BioType::SequenceSecondaryStructure,
    BioType::Silencer,
    BioType::SnoRna,
    BioType::SnoRnaGene,
    BioType::SnRna,
    BioType::SnRnaGene,
    BioType::Supercontig,
    BioType::TandemRepeat,
    BioType::TataBox,
    BioType::TelomeraseRna,
    BioType::TfBindingSite,
    BioType::ThreePrimeUtr,
    BioType::Transcript,
    BioType::TranscriptionalCisRegulatoryRegion,
    BioType::TRna,
    BioType::UnconfirmedTranscript,
    BioType::VaultRna,
    BioType::VdGeneSegment,
    BioType::VGeneSegment,
    BioType::YRna,
];

impl TryFrom<u8> for BioType {
    type Error = Error;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        ALL_BIOTYPES
            .get(value as usize)
            .copied()
            .ok_or_else(|| Error::Parse(format!("invalid biotype byte: {value}")))
    }
}

impl std::str::FromStr for BioType {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "aberrant_processed_transcript" => Ok(Self::AberrantProcessedTranscript),
            "antisense_RNA" => Ok(Self::AntisenseRna),
            "biological_region" => Ok(Self::BiologicalRegion),
            "CAAT_signal" => Ok(Self::CaatSignal),
            "CAGE_cluster" => Ok(Self::CageCluster),
            "cDNA_match" => Ok(Self::CdnaMatch),
            "CDS" => Ok(Self::Cds),
            "centromere" => Ok(Self::Centromere),
            "chromosome" => Ok(Self::Chromosome),
            "chromosome_breakpoint" => Ok(Self::ChromosomeBreakpoint),
            "conserved_region" => Ok(Self::ConservedRegion),
            "CTCF_binding_site" => Ok(Self::CtcfBindingSite),
            "C_gene_segment" => Ok(Self::CGeneSegment),
            "direct_repeat" => Ok(Self::DirectRepeat),
            "dispersed_repeat" => Ok(Self::DispersedRepeat),
            "DNaseI_hypersensitive_site" | "DNAseI_hypersensitive_site" => {
                Ok(Self::DnaseIHypersensitiveSite)
            }
            "D_gene_segment" => Ok(Self::DGeneSegment),
            "D_loop" => Ok(Self::DLoop),
            "enhancer" => Ok(Self::Enhancer),
            "enhancer_blocking_element" => Ok(Self::EnhancerBlockingElement),
            "epigenetically_modified_region" => Ok(Self::EpigeneticallyModifiedRegion),
            "exon" => Ok(Self::Exon),
            "five_prime_UTR" => Ok(Self::FivePrimeUtr),
            "GC_rich_promoter_region" => Ok(Self::GcRichPromoterRegion),
            "gene" => Ok(Self::Gene),
            "guide_RNA" => Ok(Self::GuideRna),
            "imprinting_control_region" => Ok(Self::ImprintingControlRegion),
            "insulator" => Ok(Self::Insulator),
            "inverted_repeat" => Ok(Self::InvertedRepeat),
            "J_gene_segment" => Ok(Self::JGeneSegment),
            "lincRNA" => Ok(Self::LincRna),
            "lincRNA_gene" => Ok(Self::LincRnaGene),
            "lnc_RNA" | "lncRNA" => Ok(Self::LncRna),
            "locus_control_region" => Ok(Self::LocusControlRegion),
            "match" => Ok(Self::Match),
            "matrix_attachment_site" => Ok(Self::MatrixAttachmentSite),
            "meiotic_recombination_region" => Ok(Self::MeioticRecombinationRegion),
            "microsatellite" => Ok(Self::Microsatellite),
            "minisatellite" => Ok(Self::Minisatellite),
            "miRNA" => Ok(Self::MiRna),
            "miRNA_gene" => Ok(Self::MiRnaGene),
            "mitotic_recombination_region" => Ok(Self::MitoticRecombinationRegion),
            "mobile_genetic_element" => Ok(Self::MobileGeneticElement),
            "mRNA" => Ok(Self::MRna),
            "mt_gene" => Ok(Self::MtGene),
            "ncRNA" => Ok(Self::NcRna),
            "ncRNA_gene" => Ok(Self::NcRnaGene),
            "nc_primary_transcript" => Ok(Self::NcPrimaryTranscript),
            "NMD_transcript_variant" => Ok(Self::NmdTranscriptVariant),
            "non_allelic_homologous_recombination_region" => {
                Ok(Self::NonAllelicHomologousRecombinationRegion)
            }
            "nucleotide_cleavage_site" => Ok(Self::NucleotideCleavageSite),
            "nucleotide_motif" => Ok(Self::NucleotideMotif),
            "open_chromatin_region" => Ok(Self::OpenChromatinRegion),
            "origin_of_replication" => Ok(Self::OriginOfReplication),
            "primary_transcript" => Ok(Self::PrimaryTranscript),
            "processed_pseudogene" => Ok(Self::ProcessedPseudogene),
            "processed_transcript" => Ok(Self::ProcessedTranscript),
            "promoter" => Ok(Self::Promoter),
            "promoter_flanking_region" => Ok(Self::PromoterFlankingRegion),
            "protein_binding_site" => Ok(Self::ProteinBindingSite),
            "pseudogene" => Ok(Self::Pseudogene),
            "pseudogenic_transcript" => Ok(Self::PseudogenicTranscript),
            "recombination_feature" => Ok(Self::RecombinationFeature),
            "region" => Ok(Self::Region),
            "regulatory_region" => Ok(Self::RegulatoryRegion),
            "repeat_instability_region" => Ok(Self::RepeatInstabilityRegion),
            "repeat_region" => Ok(Self::RepeatRegion),
            "replication_regulatory_region" => Ok(Self::ReplicationRegulatoryRegion),
            "replication_start_site" => Ok(Self::ReplicationStartSite),
            "response_element" => Ok(Self::ResponseElement),
            "RNA" => Ok(Self::Rna),
            "RNase_MRP_RNA" => Ok(Self::RnaseMrpRna),
            "RNase_P_RNA" => Ok(Self::RnasePRna),
            "rRNA" => Ok(Self::RRna),
            "rRNA_gene" => Ok(Self::RRnaGene),
            "scaffold" => Ok(Self::Scaffold),
            "scaRNA" => Ok(Self::ScaRna),
            "scRNA" => Ok(Self::ScRna),
            "sequence_alteration" => Ok(Self::SequenceAlteration),
            "sequence_alteration_artifact" => Ok(Self::SequenceAlterationArtifact),
            "sequence_comparison" => Ok(Self::SequenceComparison),
            "sequence_feature" => Ok(Self::SequenceFeature),
            "sequence_secondary_structure" => Ok(Self::SequenceSecondaryStructure),
            "silencer" => Ok(Self::Silencer),
            "snoRNA" => Ok(Self::SnoRna),
            "snoRNA_gene" => Ok(Self::SnoRnaGene),
            "snRNA" => Ok(Self::SnRna),
            "snRNA_gene" => Ok(Self::SnRnaGene),
            "supercontig" => Ok(Self::Supercontig),
            "tandem_repeat" => Ok(Self::TandemRepeat),
            "TATA_box" => Ok(Self::TataBox),
            "telomerase_RNA" => Ok(Self::TelomeraseRna),
            "TF_binding_site" => Ok(Self::TfBindingSite),
            "three_prime_UTR" => Ok(Self::ThreePrimeUtr),
            "transcript" => Ok(Self::Transcript),
            "transcriptional_cis_regulatory_region" => Ok(Self::TranscriptionalCisRegulatoryRegion),
            "tRNA" => Ok(Self::TRna),
            "unconfirmed_transcript" => Ok(Self::UnconfirmedTranscript),
            "vault_RNA" => Ok(Self::VaultRna),
            "VD_gene_segment" => Ok(Self::VdGeneSegment),
            "V_gene_segment" => Ok(Self::VGeneSegment),
            "Y_RNA" => Ok(Self::YRna),
            _ => Err(Error::Parse(format!("unrecognized biotype: '{s}'"))),
        }
    }
}

impl fmt::Display for BioType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Self::AberrantProcessedTranscript => "aberrant_processed_transcript",
            Self::AntisenseRna => "antisense_RNA",
            Self::BiologicalRegion => "biological_region",
            Self::CaatSignal => "CAAT_signal",
            Self::CageCluster => "CAGE_cluster",
            Self::CdnaMatch => "cDNA_match",
            Self::Cds => "CDS",
            Self::Centromere => "centromere",
            Self::Chromosome => "chromosome",
            Self::ChromosomeBreakpoint => "chromosome_breakpoint",
            Self::ConservedRegion => "conserved_region",
            Self::CtcfBindingSite => "CTCF_binding_site",
            Self::CGeneSegment => "C_gene_segment",
            Self::DirectRepeat => "direct_repeat",
            Self::DispersedRepeat => "dispersed_repeat",
            Self::DnaseIHypersensitiveSite => "DNaseI_hypersensitive_site",
            Self::DGeneSegment => "D_gene_segment",
            Self::DLoop => "D_loop",
            Self::Enhancer => "enhancer",
            Self::EnhancerBlockingElement => "enhancer_blocking_element",
            Self::EpigeneticallyModifiedRegion => "epigenetically_modified_region",
            Self::Exon => "exon",
            Self::FivePrimeUtr => "five_prime_UTR",
            Self::GcRichPromoterRegion => "GC_rich_promoter_region",
            Self::Gene => "gene",
            Self::GuideRna => "guide_RNA",
            Self::ImprintingControlRegion => "imprinting_control_region",
            Self::Insulator => "insulator",
            Self::InvertedRepeat => "inverted_repeat",
            Self::JGeneSegment => "J_gene_segment",
            Self::LincRna => "lincRNA",
            Self::LincRnaGene => "lincRNA_gene",
            Self::LncRna => "lncRNA",
            Self::LocusControlRegion => "locus_control_region",
            Self::Match => "match",
            Self::MatrixAttachmentSite => "matrix_attachment_site",
            Self::MeioticRecombinationRegion => "meiotic_recombination_region",
            Self::Microsatellite => "microsatellite",
            Self::Minisatellite => "minisatellite",
            Self::MiRna => "miRNA",
            Self::MiRnaGene => "miRNA_gene",
            Self::MitoticRecombinationRegion => "mitotic_recombination_region",
            Self::MobileGeneticElement => "mobile_genetic_element",
            Self::MRna => "mRNA",
            Self::MtGene => "mt_gene",
            Self::NcRna => "ncRNA",
            Self::NcRnaGene => "ncRNA_gene",
            Self::NcPrimaryTranscript => "nc_primary_transcript",
            Self::NmdTranscriptVariant => "NMD_transcript_variant",
            Self::NonAllelicHomologousRecombinationRegion => {
                "non_allelic_homologous_recombination_region"
            }
            Self::NucleotideCleavageSite => "nucleotide_cleavage_site",
            Self::NucleotideMotif => "nucleotide_motif",
            Self::OpenChromatinRegion => "open_chromatin_region",
            Self::OriginOfReplication => "origin_of_replication",
            Self::PrimaryTranscript => "primary_transcript",
            Self::ProcessedPseudogene => "processed_pseudogene",
            Self::ProcessedTranscript => "processed_transcript",
            Self::Promoter => "promoter",
            Self::PromoterFlankingRegion => "promoter_flanking_region",
            Self::ProteinBindingSite => "protein_binding_site",
            Self::Pseudogene => "pseudogene",
            Self::PseudogenicTranscript => "pseudogenic_transcript",
            Self::RecombinationFeature => "recombination_feature",
            Self::Region => "region",
            Self::RegulatoryRegion => "regulatory_region",
            Self::RepeatInstabilityRegion => "repeat_instability_region",
            Self::RepeatRegion => "repeat_region",
            Self::ReplicationRegulatoryRegion => "replication_regulatory_region",
            Self::ReplicationStartSite => "replication_start_site",
            Self::ResponseElement => "response_element",
            Self::Rna => "RNA",
            Self::RnaseMrpRna => "RNase_MRP_RNA",
            Self::RnasePRna => "RNase_P_RNA",
            Self::RRna => "rRNA",
            Self::RRnaGene => "rRNA_gene",
            Self::Scaffold => "scaffold",
            Self::ScaRna => "scaRNA",
            Self::ScRna => "scRNA",
            Self::SequenceAlteration => "sequence_alteration",
            Self::SequenceAlterationArtifact => "sequence_alteration_artifact",
            Self::SequenceComparison => "sequence_comparison",
            Self::SequenceFeature => "sequence_feature",
            Self::SequenceSecondaryStructure => "sequence_secondary_structure",
            Self::Silencer => "silencer",
            Self::SnoRna => "snoRNA",
            Self::SnoRnaGene => "snoRNA_gene",
            Self::SnRna => "snRNA",
            Self::SnRnaGene => "snRNA_gene",
            Self::Supercontig => "supercontig",
            Self::TandemRepeat => "tandem_repeat",
            Self::TataBox => "TATA_box",
            Self::TelomeraseRna => "telomerase_RNA",
            Self::TfBindingSite => "TF_binding_site",
            Self::ThreePrimeUtr => "three_prime_UTR",
            Self::Transcript => "transcript",
            Self::TranscriptionalCisRegulatoryRegion => "transcriptional_cis_regulatory_region",
            Self::TRna => "tRNA",
            Self::UnconfirmedTranscript => "unconfirmed_transcript",
            Self::VaultRna => "vault_RNA",
            Self::VdGeneSegment => "VD_gene_segment",
            Self::VGeneSegment => "V_gene_segment",
            Self::YRna => "Y_RNA",
        };
        write!(f, "{s}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn byte_round_trip_all() {
        for byte in 0..=101u8 {
            let bt = BioType::try_from(byte).unwrap();
            assert_eq!(bt.to_byte(), byte);
        }
    }

    #[test]
    fn invalid_byte() {
        assert!(BioType::try_from(102).is_err());
        assert!(BioType::try_from(255).is_err());
    }

    #[test]
    fn from_str_all_canonical() {
        let cases = [
            (
                "aberrant_processed_transcript",
                BioType::AberrantProcessedTranscript,
            ),
            ("antisense_RNA", BioType::AntisenseRna),
            ("biological_region", BioType::BiologicalRegion),
            ("CAAT_signal", BioType::CaatSignal),
            ("CAGE_cluster", BioType::CageCluster),
            ("cDNA_match", BioType::CdnaMatch),
            ("CDS", BioType::Cds),
            ("centromere", BioType::Centromere),
            ("chromosome", BioType::Chromosome),
            ("gene", BioType::Gene),
            ("mRNA", BioType::MRna),
            ("exon", BioType::Exon),
            ("lncRNA", BioType::LncRna),
            ("rRNA", BioType::RRna),
            ("tRNA", BioType::TRna),
            ("pseudogene", BioType::Pseudogene),
            ("enhancer", BioType::Enhancer),
            ("promoter", BioType::Promoter),
            ("region", BioType::Region),
            ("Y_RNA", BioType::YRna),
        ];
        for (s, expected) in cases {
            assert_eq!(s.parse::<BioType>().unwrap(), expected, "failed for '{s}'");
        }
    }

    #[test]
    fn from_str_aliases() {
        // lnc_RNA → lncRNA
        assert_eq!("lnc_RNA".parse::<BioType>().unwrap(), BioType::LncRna);
        // DNAseI_hypersensitive_site → DNaseI_hypersensitive_site
        assert_eq!(
            "DNAseI_hypersensitive_site".parse::<BioType>().unwrap(),
            BioType::DnaseIHypersensitiveSite
        );
        assert_eq!(
            "DNaseI_hypersensitive_site".parse::<BioType>().unwrap(),
            BioType::DnaseIHypersensitiveSite
        );
    }

    #[test]
    fn from_str_unknown() {
        assert!("unknown_type".parse::<BioType>().is_err());
    }

    #[test]
    fn category_gene() {
        assert_eq!(BioType::Gene.category(), BioTypeCategory::Gene);
        assert_eq!(BioType::Pseudogene.category(), BioTypeCategory::Gene);
        assert_eq!(BioType::MiRnaGene.category(), BioTypeCategory::Gene);
    }

    #[test]
    fn category_transcript() {
        assert_eq!(BioType::MRna.category(), BioTypeCategory::Transcript);
        assert_eq!(BioType::LncRna.category(), BioTypeCategory::Transcript);
        assert_eq!(
            BioType::CGeneSegment.category(),
            BioTypeCategory::Transcript
        );
    }

    #[test]
    fn category_regulatory() {
        assert_eq!(BioType::Enhancer.category(), BioTypeCategory::Regulatory);
        assert_eq!(BioType::Promoter.category(), BioTypeCategory::Regulatory);
    }

    #[test]
    fn category_child_feature() {
        assert_eq!(BioType::Exon.category(), BioTypeCategory::ChildFeature);
        assert_eq!(BioType::Cds.category(), BioTypeCategory::ChildFeature);
    }

    #[test]
    fn is_coding() {
        assert!(BioType::MRna.is_coding());
        assert!(BioType::CGeneSegment.is_coding());
        assert!(BioType::DGeneSegment.is_coding());
        assert!(BioType::JGeneSegment.is_coding());
        assert!(BioType::VGeneSegment.is_coding());
        assert!(BioType::VdGeneSegment.is_coding());
        assert!(!BioType::LncRna.is_coding());
        assert!(!BioType::Gene.is_coding());
    }

    #[test]
    fn display_round_trip() {
        for byte in 0..=101u8 {
            let bt = BioType::try_from(byte).unwrap();
            let s = bt.to_string();
            let parsed: BioType = s.parse().unwrap();
            assert_eq!(bt, parsed, "display round-trip failed for byte {byte}");
        }
    }

    #[test]
    fn not_useful_region() {
        assert_eq!(BioType::Region.category(), BioTypeCategory::NotUseful);
    }
}
