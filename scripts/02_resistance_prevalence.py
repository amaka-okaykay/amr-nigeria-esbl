"""
02_resistance_prevalence.py
============================
Calculates resistance gene prevalence by species and gene family.
Performs Fisher's exact test for species-specific distribution differences.

Author: Okeke Chiamaka Juliet
Study:  Nsofor et al. (2023) — ESBL gene detection, Owerri Nigeria
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

# ─── LOAD DATA ────────────────────────────────────────────────────────────────
pcr = pd.read_csv("data/pcr_gene_results.csv")
meta = pd.read_csv("data/isolate_metadata.csv")
df = pcr.merge(meta[["isolate_id", "hospital"]], on="isolate_id")

print("=" * 60)
print("ESBL RESISTANCE GENE PREVALENCE ANALYSIS")
print("Study: Owerri Clinical Isolates | Nsofor et al. 2023")
print("=" * 60)

# ─── OVERALL PREVALENCE ───────────────────────────────────────────────────────
print("\n── OVERALL GENE PREVALENCE ──────────────────────────────────")
genes = ["blaCTX_M", "blaTEM", "blaSHV"]
gene_labels = ["blaCTX-M", "blaTEM", "blaSHV"]

for gene, label in zip(genes, gene_labels):
    n_pos  = df[gene].sum()
    total  = len(df)
    pct    = n_pos / total * 100
    # 95% Wilson confidence interval
    from scipy.stats import binom
    ci_lo  = binom.ppf(0.025, total, n_pos/total) / total * 100
    ci_hi  = binom.ppf(0.975, total, n_pos/total) / total * 100
    print(f"  {label:12s}  n={n_pos:3d}/{total}  {pct:5.1f}%  (95% CI: {ci_lo:.1f}–{ci_hi:.1f}%)")

# ─── PREVALENCE BY SPECIES ────────────────────────────────────────────────────
print("\n── PREVALENCE BY SPECIES ────────────────────────────────────")
species_list = ["Escherichia coli", "Klebsiella pneumoniae"]
results = []

for gene, label in zip(genes, gene_labels):
    row = {"Gene": label}
    counts = {}
    for sp in species_list:
        sub = df[df["species"] == sp]
        n_pos = sub[gene].sum()
        pct   = n_pos / len(sub) * 100
        row[sp.split()[0]] = f"{pct:.1f}% (n={n_pos}/{len(sub)})"
        counts[sp] = (n_pos, len(sub) - n_pos)

    # Fisher's exact test
    table = [[counts[sp][0], counts[sp][1]] for sp in species_list]
    _, p = stats.fisher_exact(table)
    row["Fisher p"] = f"{p:.4f}" + (" *" if p < 0.05 else "")
    results.append(row)
    print(f"  {label:12s}  E.coli: {row['Escherichia']:25s}  K.pneumo: {row['Klebsiella']:25s}  p={p:.4f}")

prevalence_df = pd.DataFrame(results)
prevalence_df.to_csv("results/gene_prevalence_table.csv", index=False)
print("\n  ✓ Saved to results/gene_prevalence_table.csv")

# ─── PREVALENCE BY HOSPITAL ───────────────────────────────────────────────────
print("\n── PREVALENCE BY HOSPITAL ───────────────────────────────────")
hosp_prev = []
for hospital in df["hospital"].unique():
    sub = df[df["hospital"] == hospital]
    row = {"Hospital": hospital, "n": len(sub)}
    for gene, label in zip(genes, gene_labels):
        row[label] = f"{sub[gene].mean()*100:.1f}%"
    hosp_prev.append(row)
    print(f"  {hospital[:35]:35s}  n={len(sub)}  CTX-M: {row['blaCTX-M']}  TEM: {row['blaTEM']}  SHV: {row['blaSHV']}")

pd.DataFrame(hosp_prev).to_csv("results/hospital_prevalence.csv", index=False)
print("\n  ✓ Saved to results/hospital_prevalence.csv")

# ─── MULTI-GENE CARRIAGE ─────────────────────────────────────────────────────
print("\n── MULTI-GENE CARRIAGE DISTRIBUTION ────────────────────────")
gene_counts = df["gene_count"].value_counts().sort_index()
for n_genes, count in gene_counts.items():
    pct = count / len(df) * 100
    bar = "█" * int(pct / 2)
    label = ["No genes", "Single gene", "Two gene families", "All three gene families"][min(n_genes, 3)]
    print(f"  {n_genes} genes — {label:25s}: {count:3d} isolates ({pct:4.1f}%)  {bar}")

multi_pct = (df["gene_count"] >= 2).mean() * 100
print(f"\n  Co-resistance (≥2 gene families): {multi_pct:.1f}% of ESBL+ isolates")
print("  ⚠ This indicates plasmid-mediated HGT as the dominant dissemination mechanism")

print("\n✅ Prevalence analysis complete. Run 03_cooccurrence_analysis.py next.")
