"""
04_wgs_extension_pipeline.py
=============================
Conceptual pipeline for extending PCR-based gene detection to
WGS-based AMR surveillance — the next phase of this research.

This script outlines the bioinformatics workflow for:
  1. Quality control of raw Illumina reads
  2. Assembly and annotation
  3. Resistance gene identification via CARD/ResFinder
  4. MLST and plasmid typing
  5. Whole genome phylogenetics

Status: PLANNED — requires access to Illumina sequencing data.
        Structured here as a reproducible template for when
        sequencing data becomes available.

Author: Okeke Chiamaka Juliet
        Wellcome Connecting Science — AMR Databases & Genotype Prediction, 2024
        Genomac Institute Inc. — Bioinformatics Course 201, 2024
"""

# ─── PHASE 2: WGS WORKFLOW TEMPLATE ──────────────────────────────────────────
#
# This pipeline is designed to run on raw FASTQ files from Illumina short-read
# sequencing of the ESBL-positive clinical isolates characterised in Phase 1.
#
# Required tools (to be installed in pipeline environment):
#   - FastQC / MultiQC      → Read quality assessment
#   - Trimmomatic           → Adapter trimming and quality filtering
#   - SPAdes                → De novo genome assembly
#   - Prokka                → Genome annotation
#   - CARD / RGI            → Resistance gene identification from assembly
#   - ResFinder             → Phenotype prediction from sequence
#   - mlst                  → Multi-locus sequence typing
#   - PlasmidFinder         → Plasmid replicon identification
#   - Snippy                → Core SNP phylogeny
#   - IQ-TREE               → Maximum likelihood tree construction

WGS_PIPELINE_STEPS = [
    {
        "phase": "01",
        "name": "Read QC",
        "tool": "FastQC + MultiQC",
        "command_template": "fastqc {sample}_R1.fastq.gz {sample}_R2.fastq.gz -o qc/",
        "purpose": "Assess raw read quality — per-base quality scores, adapter contamination, GC content",
        "amr_relevance": "Poor quality reads produce assembly errors that can create false resistance gene calls"
    },
    {
        "phase": "02",
        "name": "Adapter Trimming",
        "tool": "Trimmomatic",
        "command_template": "trimmomatic PE {sample}_R1.fastq.gz {sample}_R2.fastq.gz "
                            "{sample}_R1_trim.fastq.gz /dev/null "
                            "{sample}_R2_trim.fastq.gz /dev/null "
                            "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
        "purpose": "Remove Illumina adapters and low-quality bases",
        "amr_relevance": "Clean reads improve assembly contiguity — critical for detecting full-length resistance genes"
    },
    {
        "phase": "03",
        "name": "De Novo Assembly",
        "tool": "SPAdes",
        "command_template": "spades.py -1 {sample}_R1_trim.fastq.gz -2 {sample}_R2_trim.fastq.gz "
                            "-o assembly/{sample}/ --careful",
        "purpose": "Assemble trimmed reads into contigs",
        "amr_relevance": "Assembly quality determines whether resistance genes on plasmids are fully recovered"
    },
    {
        "phase": "04",
        "name": "Genome Annotation",
        "tool": "Prokka",
        "command_template": "prokka assembly/{sample}/contigs.fasta --outdir annotation/{sample}/ "
                            "--genus Escherichia --species coli --prefix {sample}",
        "purpose": "Annotate coding sequences, rRNAs, tRNAs",
        "amr_relevance": "Provides genomic context for resistance genes — identifies mobile genetic elements"
    },
    {
        "phase": "05",
        "name": "Resistance Gene Identification",
        "tool": "CARD / RGI",
        "command_template": "rgi main --input_sequence annotation/{sample}/{sample}.faa "
                            "--output_file card_results/{sample} --input_type protein "
                            "--alignment_tool BLAST",
        "purpose": "Identify all AMR genes using CARD Resistance Gene Identifier",
        "amr_relevance": "Moves beyond PCR gene family detection to exact variant identification "
                         "(e.g. blaCTX-M-15 vs blaCTX-M-1) with CARD ontology classification"
    },
    {
        "phase": "06",
        "name": "Phenotype Prediction",
        "tool": "ResFinder",
        "command_template": "python resfinder.py -i assembly/{sample}/contigs.fasta "
                            "-o resfinder/{sample}/ -s 'Escherichia coli' -l 0.6 -t 0.9",
        "purpose": "Predict resistance phenotype from sequence — computationally infer MIC from genotype",
        "amr_relevance": "Genotype-phenotype correlation: validates PCR/AST findings and predicts "
                         "resistance to antibiotics not tested in the clinical panel"
    },
    {
        "phase": "07",
        "name": "MLST Typing",
        "tool": "mlst (Torsten Seemann)",
        "command_template": "mlst assembly/{sample}/contigs.fasta",
        "purpose": "Multi-locus sequence typing — identify clonal lineages (STs)",
        "amr_relevance": "Determines whether ESBL spread is driven by a single epidemic clone "
                         "(e.g. E. coli ST131) or multiple independent acquisition events"
    },
    {
        "phase": "08",
        "name": "Plasmid Typing",
        "tool": "PlasmidFinder",
        "command_template": "python plasmidfinder.py -i assembly/{sample}/contigs.fasta "
                            "-o plasmidfinder/{sample}/ -p enterobacteriaceae",
        "purpose": "Identify plasmid replicon types",
        "amr_relevance": "Determines whether blaCTX-M, blaTEM, blaSHV are chromosomally integrated "
                         "or plasmid-borne — essential for understanding horizontal gene transfer dynamics"
    },
    {
        "phase": "09",
        "name": "Core Genome SNP Phylogeny",
        "tool": "Snippy + IQ-TREE",
        "command_template": "snippy-multi samples.txt --ref reference.fasta --cpus 8 > snippy_cmd.sh\n"
                            "iqtree -s core.full.aln -m GTR+G -bb 1000 -nt AUTO",
        "purpose": "Construct maximum likelihood phylogenetic tree from core genome SNPs",
        "amr_relevance": "Reconstructs transmission networks — identifies if hospital isolates are "
                         "clonally related (outbreak) or polyclonal (independent acquisition)"
    },
    {
        "phase": "10",
        "name": "WHO GLASS Submission",
        "tool": "NCBI Submission Portal + AMRFinderPlus",
        "command_template": "amrfinder -n assembly/{sample}/contigs.fasta -o amrfinder/{sample}.tsv "
                            "--organism Escherichia --plus",
        "purpose": "Prepare sequences for submission to WHO GLASS / NCBI Pathogen portal",
        "amr_relevance": "Adds West African clinical AMR data to globally underrepresented repositories — "
                         "directly addressing the surveillance gap this research is designed to close"
    }
]

if __name__ == "__main__":
    print("=" * 70)
    print("WGS-BASED AMR SURVEILLANCE PIPELINE — PHASE 2 TEMPLATE")
    print("Owerri Clinical ESBL Isolates — Extension of Nsofor et al. 2023")
    print("=" * 70)
    print()

    for step in WGS_PIPELINE_STEPS:
        print(f"[PHASE {step['phase']}] {step['name']} — {step['tool']}")
        print(f"  Purpose:       {step['purpose']}")
        print(f"  AMR relevance: {step['amr_relevance']}")
        print()

    print("─" * 70)
    print("TRANSITION FROM PHASE 1 → PHASE 2")
    print("─" * 70)
    print()
    print("Phase 1 (Complete — Published 2023):")
    print("  PCR detected blaCTX-M, blaTEM, blaSHV GENE FAMILIES")
    print("  → Answers: Which resistance genes are present?")
    print()
    print("Phase 2 (Planned — WGS):")
    print("  CARD/RGI will identify exact GENE VARIANTS (e.g. blaCTX-M-15)")
    print("  ResFinder will predict resistance phenotype from sequence")
    print("  MLST will identify clonal lineages (e.g. ST131, ST258)")
    print("  PlasmidFinder will characterise mobile genetic elements")
    print("  IQ-TREE will reconstruct hospital transmission networks")
    print("  → Answers: How are these genes spreading? Who transmitted to whom?")
    print("             Which variants? On which plasmids? In which lineages?")
    print()
    print("This is the full AMR bioinformatics pipeline that Phase 1 made necessary.")
    print("Training: Wellcome Connecting Science (2024) + Genomac Institute (2024)")
    print("Author:   Okeke Chiamaka Juliet | okekechiamakajuls@gmail.com")
