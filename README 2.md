# ESBL Gene Detection in Clinical Isolates from Southeast Nigeria
### Molecular Epidemiology of Antimicrobial Resistance — Owerri, Imo State

[![Research](https://img.shields.io/badge/Published-Reviews_%26_Research_in_Medical_Microbiology-blue)](https://journals.lww.com/rrmicrobiology)
[![Institute](https://img.shields.io/badge/Institution-FUTO-green)](https://futo.edu.ng)
[![Wellcome](https://img.shields.io/badge/Training-Wellcome_Connecting_Science-red)](https://www.wellcomeconnectingscience.org)
[![Genomac](https://img.shields.io/badge/Training-Genomac_Institute_Inc-orange)](https://genomac.org)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)
[![GitHub Pages](https://img.shields.io/badge/Dashboard-Live_View-brightgreen)](https://chiamaka-okeke.github.io/amr-nigeria-esbl)

---

## 🔬 Overview

This repository documents the molecular epidemiology analysis underpinning the peer-reviewed study:

> **Nsofor CA, Moses A, Onyeakazi CM, Okeke CJ, Ikegbunam MN. (2023).** Detection of *blaCTX-M*, *blaTEM*, and *blaSHV* genes in clinical isolates of *Escherichia coli* and *Klebsiella pneumoniae.* *Reviews and Research in Medical Microbiology*, 34(2):66–72. LWW/Wolters Kluwer.

Southeast Nigeria is one of the most under-surveilled regions in sub-Saharan Africa for antimicrobial resistance (AMR). This study investigated Extended-Spectrum Beta-Lactamase (ESBL)-producing *Escherichia coli* and *Klebsiella pneumoniae* from three major hospitals in Owerri, contributing the first PCR-confirmed ESBL gene characterisation data from this clinical region to the global AMR evidence base.

**This repository extends the published findings** with interactive visualisation, reproducible analysis scripts, and a structured framework for integrating this molecular data with whole-genome sequencing (WGS) pipelines — the natural next step in the AMR bioinformatics trajectory this research initiated.

---

## 📊 Key Findings at a Glance

| Metric | Value |
|--------|-------|
| Total clinical samples collected | 433 |
| Samples yielding bacterial growth | 249 (57.5%) |
| ESBL-positive isolates confirmed | 92 (36.9%) |
| *blaCTX-M* prevalence | 78.3% of ESBL+ isolates |
| *blaTEM* prevalence | 65.2% of ESBL+ isolates |
| *blaSHV* prevalence | 41.3% of ESBL+ isolates |
| Ceftazidime resistance rate | **90.8%** |
| Co-resistance (≥2 gene families) | 54.3% of ESBL+ isolates |

> ⚠️ **90.8% ceftazidime resistance** — among the highest documented in Nigerian clinical literature at time of publication. This has direct implications for empirical prescribing in tertiary hospitals in Imo State.

---

## 🏥 Study Sites

Clinical samples were collected from three hospitals in Owerri, Imo State, Nigeria:

- **Imo State Specialist Hospital** — tertiary public referral centre
- **Federal Medical Centre Owerri** — federal tertiary hospital
- **Private clinical facilities, Owerri** — supplementary community-level sampling

All isolates were consecutively collected over the study period. Sample collection followed standard clinical microbiology protocols and institutional ethics frameworks.

---

## 🧬 Methodology

### 1. Bacterial Culture & Identification
- Standard culture on MacConkey agar and blood agar
- Biochemical identification using conventional microbiological tests
- Species confirmed as *E. coli* and *K. pneumoniae*

### 2. Phenotypic ESBL Confirmation
- **Double-Disk Synergy Test (DDST)** — clavulanic acid-augmented disk pair (ceftazidime + amoxicillin-clavulanate)
- **Antimicrobial Susceptibility Testing** — Kirby-Bauer disk diffusion per CLSI breakpoints
- Panel of 12 antibiotics including 3rd-generation cephalosporins, fluoroquinolones, aminoglycosides, carbapenems

### 3. PCR-Based Molecular Detection
```
Target genes:    blaCTX-M  |  blaTEM  |  blaSHV
Method:          Conventional PCR
DNA extraction:  Boiling lysis method
Visualisation:   Agarose gel electrophoresis (1.5%)
Controls:        Positive (ATCC reference strains) + Negative (no-template)
```

### 4. Data Analysis
- Resistance gene prevalence calculated per isolate and per bacterial species
- Co-resistance patterns computed across all three gene families
- Fisher's exact test for species-specific gene distribution comparisons
- Minimum inhibitory concentration (MIC) proxy analysis via disk zone diameters

---

## 📁 Repository Structure

```
amr-nigeria-esbl/
│
├── README.md                          ← This file
├── index.html                         ← Interactive dashboard (GitHub Pages)
│
├── data/
│   ├── isolate_metadata.csv           ← Sample IDs, source hospital, species
│   ├── susceptibility_results.csv     ← AST zone diameters per antibiotic
│   ├── pcr_gene_results.csv           ← Binary gene detection results per isolate
│   └── esbl_summary_stats.csv         ← Aggregated prevalence statistics
│
├── scripts/
│   ├── 01_data_cleaning.py            ← Raw data QC and formatting
│   ├── 02_resistance_prevalence.py    ← Gene prevalence calculations
│   ├── 03_cooccurrence_analysis.py    ← Multi-gene co-resistance patterns
│   ├── 04_ast_analysis.py             ← Antibiotic susceptibility analysis
│   └── 05_visualisation.py            ← Figures and summary plots
│
├── results/
│   ├── gene_prevalence_table.csv      ← Final prevalence by species and gene
│   ├── cooccurrence_matrix.csv        ← Gene co-occurrence frequency matrix
│   └── ast_resistance_profile.csv     ← Resistance rates per antibiotic class
│
├── figures/
│   ├── fig1_gene_prevalence.png       ← Bar chart: gene frequency by species
│   ├── fig2_cooccurrence_heatmap.png  ← Heatmap: gene family co-resistance
│   ├── fig3_ast_profile.png           ← Antibiotic resistance spectrum
│   └── fig4_hospital_distribution.png ← Sample source breakdown
│
└── requirements.txt                   ← Python dependencies
```

---

## 🚀 Running the Analysis

### Prerequisites
```bash
pip install -r requirements.txt
```

### Execute the full pipeline
```bash
# Step 1: Clean and format raw data
python scripts/01_data_cleaning.py

# Step 2: Calculate resistance gene prevalence
python scripts/02_resistance_prevalence.py

# Step 3: Analyse co-occurrence patterns
python scripts/03_cooccurrence_analysis.py

# Step 4: Antibiotic susceptibility analysis
python scripts/04_ast_analysis.py

# Step 5: Generate all figures
python scripts/05_visualisation.py
```

All outputs are written to `/results/` and `/figures/`.

---

## 🔭 Next Steps: From PCR to WGS

This study answers: *"Which resistance genes are present?"*

The natural next question — *"How did these genes arrive, how are they spreading, and what evolutionary pressures are shaping them?"* — requires Whole Genome Sequencing (WGS)-based analysis. Planned extensions of this work include:

| Phase | Method | Question Addressed |
|-------|--------|-------------------|
| **Phase 1** (complete) | PCR gene detection | Presence/absence of blaCTX-M, blaTEM, blaSHV |
| **Phase 2** (planned) | Short-read WGS (Illumina) | Full resistome characterisation; plasmid typing |
| **Phase 3** (planned) | Phylogenetic reconstruction | Transmission networks; clonal lineage tracing |
| **Phase 4** (planned) | Population genomics | Gene dissemination dynamics across Nigerian hospitals |
| **Phase 5** (planned) | WHO GLASS contribution | West African sequence data submission to global repositories |

> This repository is the foundation for a WGS-based AMR surveillance pipeline for Owerri clinical settings — contributing to Nigeria's representation in global genomic databases, which are currently dominated by European and North American sequences.

---

## 🌍 Why This Matters

- Nigeria has **14,165 hospitals** and one of the highest burdens of drug-resistant infections in sub-Saharan Africa
- West Africa is **critically underrepresented** in global AMR genomic databases (WHO GLASS, NCBI, CARD)
- 90.8% ceftazidime resistance means **frontline antibiotics are failing** in Owerri hospitals with no computational surveillance infrastructure to track or predict spread
- This study is one of the **few PCR-confirmed ESBL gene characterisation datasets** from this region in the published literature

---

## 📚 References & Resources

- **Published paper:** Nsofor et al. (2023), *Reviews and Research in Medical Microbiology*, 34(2):66–72
- **AMR Databases:** [CARD](https://card.mcmaster.ca) | [ResFinder](https://cge.food.dtu.dk/services/ResFinder/) | [NDARO](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/) | [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/)
- **WHO GLASS:** [Global Antimicrobial Resistance and Use Surveillance System](https://www.who.int/initiatives/glass)
- **Training:** [Wellcome Connecting Science](https://www.wellcomeconnectingscience.org) | [Genomac Institute Inc.](https://genomac.org)

---

## 👩🏾‍🔬 Author

**Okeke Chiamaka Juliet**
Published AMR Researcher | Wellcome & Genomac Bioinformatics Trained | Biotechnologist

📧 okekechiamakajuls@gmail.com
🔬 Federal University of Technology Owerri (FUTO), B.Tech Biotechnology, 2022
📄 Wellcome Connecting Science — AMR Databases & Genotype Prediction, 2024
💻 Genomac Institute Inc. — Bioinformatics Course 201, 2024
🌍 Lagos, Nigeria

---

## 📄 License

This repository is licensed under the MIT License. The underlying research data belongs to the original study authors. Scripts and visualisation code are freely reusable with attribution.

---

*"PCR-based gene detection tells you which resistance genes are present in a hospital at a given moment. It cannot tell you how those genes arrived, where they are going, or what evolutionary pressures are shaping them. Those questions require whole genome sequencing data, phylogenetic analysis, and computational pipelines. This repository is the first step toward answering them."*
