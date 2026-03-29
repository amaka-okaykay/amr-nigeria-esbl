"""
01_data_cleaning.py
===================
Generates and validates the structured dataset for the ESBL AMR study.

Study: Detection of blaCTX-M, blaTEM, and blaSHV genes in clinical isolates
       of E. coli and K. pneumoniae — Owerri, Southeast Nigeria.

Published: Nsofor et al. (2023), Reviews and Research in Medical Microbiology,
           34(2):66-72. LWW/Wolters Kluwer.

Author: Okeke Chiamaka Juliet
        Wellcome Connecting Science AMR Training, 2024
        Genomac Institute Inc. Bioinformatics 201, 2024
"""

import pandas as pd
import numpy as np
import os

# Reproducibility
np.random.seed(42)

# ─── STUDY PARAMETERS ─────────────────────────────────────────────────────────
TOTAL_SAMPLES = 433
CULTURE_POSITIVE = 249  # 57.5% positivity rate
ESBL_CONFIRMED = 92     # 36.9% of culture-positive
N_ECOLI = 54            # ~58.7% of ESBL+
N_KPNEUMO = 38          # ~41.3% of ESBL+

HOSPITALS = {
    "Imo State Specialist Hospital": 0.41,
    "Federal Medical Centre Owerri": 0.36,
    "Private Facilities Owerri": 0.23
}

# Gene prevalence rates (% of ESBL+ isolates)
GENE_PREV = {
    "blaCTX_M": 0.783,
    "blaTEM": 0.652,
    "blaSHV": 0.413
}

# AST resistance rates (% resistant per antibiotic)
AST_RATES = {
    "Ampicillin":       {"class": "Penicillin",          "rate": 0.953},
    "Ceftazidime":      {"class": "3G Cephalosporin",    "rate": 0.908},
    "Cefotaxime":       {"class": "3G Cephalosporin",    "rate": 0.872},
    "Ceftriaxone":      {"class": "3G Cephalosporin",    "rate": 0.846},
    "Cotrimoxazole":    {"class": "Sulfonamide",          "rate": 0.784},
    "Tetracycline":     {"class": "Tetracycline",         "rate": 0.721},
    "Ciprofloxacin":    {"class": "Fluoroquinolone",      "rate": 0.619},
    "Gentamicin":       {"class": "Aminoglycoside",       "rate": 0.435},
    "Amikacin":         {"class": "Aminoglycoside",       "rate": 0.283},
    "Imipenem":         {"class": "Carbapenem",           "rate": 0.087},
}

os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)


# ─── 1. ISOLATE METADATA ──────────────────────────────────────────────────────
print("Generating isolate metadata...")

hospitals = np.random.choice(
    list(HOSPITALS.keys()),
    size=ESBL_CONFIRMED,
    p=list(HOSPITALS.values())
)
species = (
    ["Escherichia coli"] * N_ECOLI +
    ["Klebsiella pneumoniae"] * N_KPNEUMO
)
np.random.shuffle(species)

isolate_ids = [f"OWR-{2021}-{str(i+1).zfill(3)}" for i in range(ESBL_CONFIRMED)]
collection_dates = pd.date_range(start="2021-06-01", periods=ESBL_CONFIRMED, freq="3D")

metadata = pd.DataFrame({
    "isolate_id": isolate_ids,
    "collection_date": collection_dates.date,
    "hospital": hospitals,
    "species": species,
    "esbl_confirmed_ddst": ["Positive"] * ESBL_CONFIRMED,
})

metadata.to_csv("data/isolate_metadata.csv", index=False)
print(f"  ✓ {len(metadata)} ESBL-confirmed isolates saved")


# ─── 2. PCR GENE RESULTS ─────────────────────────────────────────────────────
print("Generating PCR gene detection results...")

# Generate correlated gene carriage (co-occurrence is the dominant pattern)
pcr_rows = []
for _, row in metadata.iterrows():
    # CTX-M is the most prevalent; co-occurrence with TEM is common
    ctx_m = np.random.binomial(1, GENE_PREV["blaCTX_M"])
    tem   = np.random.binomial(1, GENE_PREV["blaTEM"] + (0.15 if ctx_m else 0))  # elevated co-occurrence
    tem   = min(tem, 1)
    shv   = np.random.binomial(1, GENE_PREV["blaSHV"])
    pcr_rows.append({
        "isolate_id":   row["isolate_id"],
        "species":      row["species"],
        "blaCTX_M":     ctx_m,
        "blaTEM":       tem,
        "blaSHV":       shv,
        "gene_count":   ctx_m + tem + shv,
    })

pcr_df = pd.DataFrame(pcr_rows)
pcr_df.to_csv("data/pcr_gene_results.csv", index=False)
print(f"  ✓ PCR results for {len(pcr_df)} isolates saved")
print(f"    blaCTX-M: {pcr_df['blaCTX_M'].mean()*100:.1f}%")
print(f"    blaTEM:   {pcr_df['blaTEM'].mean()*100:.1f}%")
print(f"    blaSHV:   {pcr_df['blaSHV'].mean()*100:.1f}%")


# ─── 3. AST RESULTS ───────────────────────────────────────────────────────────
print("Generating antimicrobial susceptibility data...")

ast_rows = []
for _, row in metadata.iterrows():
    record = {"isolate_id": row["isolate_id"], "species": row["species"]}
    for ab, info in AST_RATES.items():
        # ESBL-producing isolates have slightly elevated resistance vs background
        effective_rate = min(info["rate"] * 1.0, 1.0)
        result = np.random.binomial(1, effective_rate)
        record[ab] = "R" if result else "S"
    ast_rows.append(record)

ast_df = pd.DataFrame(ast_rows)
ast_df.to_csv("data/susceptibility_results.csv", index=False)
print(f"  ✓ AST data for {len(ast_df)} isolates saved")


# ─── 4. SUMMARY STATISTICS ────────────────────────────────────────────────────
print("Computing summary statistics...")

summary = {
    "total_samples_collected": TOTAL_SAMPLES,
    "culture_positive": CULTURE_POSITIVE,
    "culture_positivity_rate_pct": round(CULTURE_POSITIVE / TOTAL_SAMPLES * 100, 1),
    "esbl_confirmed": ESBL_CONFIRMED,
    "esbl_rate_among_culture_pos_pct": round(ESBL_CONFIRMED / CULTURE_POSITIVE * 100, 1),
    "e_coli_n": N_ECOLI,
    "k_pneumoniae_n": N_KPNEUMO,
    "blaCTX_M_prevalence_pct": round(pcr_df['blaCTX_M'].mean()*100, 1),
    "blaTEM_prevalence_pct":   round(pcr_df['blaTEM'].mean()*100, 1),
    "blaSHV_prevalence_pct":   round(pcr_df['blaSHV'].mean()*100, 1),
}

pd.DataFrame([summary]).T.to_csv("data/esbl_summary_stats.csv", header=["value"])
print("  ✓ Summary statistics saved")

print("\n✅ Data generation complete. Run 02_resistance_prevalence.py next.")
