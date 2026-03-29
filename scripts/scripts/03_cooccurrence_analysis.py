"""
03_cooccurrence_analysis.py
============================
Analyses multi-gene co-resistance patterns among ESBL-positive isolates.
Computes gene family co-occurrence frequencies and clinical implications.

Author: Okeke Chiamaka Juliet
Study:  Nsofor et al. (2023) — ESBL gene detection, Owerri Nigeria
"""

import pandas as pd
import numpy as np
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")

# ─── LOAD DATA ────────────────────────────────────────────────────────────────
pcr = pd.read_csv("data/pcr_gene_results.csv")
ast = pd.read_csv("data/susceptibility_results.csv")
genes = ["blaCTX_M", "blaTEM", "blaSHV"]
gene_labels = {"blaCTX_M": "blaCTX-M", "blaTEM": "blaTEM", "blaSHV": "blaSHV"}

print("=" * 60)
print("GENE CO-OCCURRENCE & AST ANALYSIS")
print("Study: Owerri Clinical Isolates | Nsofor et al. 2023")
print("=" * 60)

# ─── SINGLE-GENE PATTERNS ────────────────────────────────────────────────────
print("\n── GENE CARRIAGE PATTERNS ───────────────────────────────────")

patterns = []

# Single gene only
for g in genes:
    others = [x for x in genes if x != g]
    mask = (pcr[g] == 1) & (pcr[others[0]] == 0) & (pcr[others[1]] == 0)
    n = mask.sum()
    pct = n / len(pcr) * 100
    label = f"{gene_labels[g]} only"
    patterns.append({"pattern": label, "n": n, "pct": round(pct, 1), "type": "single"})
    print(f"  {label:30s}: {n:3d} ({pct:.1f}%)")

# Two-gene combinations
for g1, g2 in combinations(genes, 2):
    other = [x for x in genes if x not in [g1, g2]][0]
    mask = (pcr[g1] == 1) & (pcr[g2] == 1) & (pcr[other] == 0)
    n = mask.sum()
    pct = n / len(pcr) * 100
    label = f"{gene_labels[g1]} + {gene_labels[g2]}"
    patterns.append({"pattern": label, "n": n, "pct": round(pct, 1), "type": "double"})
    print(f"  {label:30s}: {n:3d} ({pct:.1f}%)")

# All three
mask_all = (pcr[genes[0]] == 1) & (pcr[genes[1]] == 1) & (pcr[genes[2]] == 1)
n = mask_all.sum()
pct = n / len(pcr) * 100
patterns.append({"pattern": "All three gene families", "n": n, "pct": round(pct, 1), "type": "triple"})
print(f"  {'All three gene families':30s}: {n:3d} ({pct:.1f}%)")

# No gene detected (edge case)
mask_none = (pcr[genes[0]] == 0) & (pcr[genes[1]] == 0) & (pcr[genes[2]] == 0)
n_none = mask_none.sum()
if n_none > 0:
    print(f"  {'No gene detected (PCR-neg)':30s}: {n_none:3d}")

patterns_df = pd.DataFrame(patterns)
patterns_df.to_csv("results/cooccurrence_matrix.csv", index=False)
print("\n  ✓ Saved to results/cooccurrence_matrix.csv")

# ─── CO-OCCURRENCE SUMMARY ────────────────────────────────────────────────────
print("\n── CLINICAL CO-RESISTANCE SUMMARY ───────────────────────────")
multi = patterns_df[patterns_df["type"].isin(["double", "triple"])]["n"].sum()
total = len(pcr)
print(f"  Isolates with ≥2 gene families: {multi}/{total} ({multi/total*100:.1f}%)")
print(f"  Most common multi-gene pattern: {patterns_df[patterns_df['type']=='double'].nlargest(1,'n')['pattern'].values[0]}")
print(f"  Triple-gene carriage: {patterns_df[patterns_df['type']=='triple']['n'].values[0]} isolates")
print("\n  ⚠ Dominant co-resistance pattern suggests:")
print("    → Conjugative plasmid transfer as primary mechanism")
print("    → Class 1 integron-mediated multi-gene cassette carriage likely")
print("    → WGS plasmid typing needed to confirm mobile genetic element context")

# ─── AST CO-RESISTANCE WITH GENE CARRIAGE ────────────────────────────────────
print("\n── AST RESISTANCE RATES IN ESBL+ ISOLATES ───────────────────")
antibiotics = [c for c in ast.columns if c not in ["isolate_id", "species"]]
ast_rates = []
for ab in antibiotics:
    n_r = (ast[ab] == "R").sum()
    pct = n_r / len(ast) * 100
    ast_rates.append({"antibiotic": ab, "n_resistant": n_r, "total": len(ast), "resistance_pct": round(pct, 1)})
    level = "CRITICAL" if pct > 80 else "HIGH" if pct > 60 else "MODERATE" if pct > 40 else "LOW"
    print(f"  {ab:20s}: {pct:5.1f}%  [{level}]")

ast_summary = pd.DataFrame(ast_rates).sort_values("resistance_pct", ascending=False)
ast_summary.to_csv("results/ast_resistance_profile.csv", index=False)
print("\n  ✓ Saved to results/ast_resistance_profile.csv")

# ─── CEFTAZIDIME HIGHLIGHT ────────────────────────────────────────────────────
print("\n── CEFTAZIDIME RESISTANCE — CLINICAL SIGNIFICANCE ──────────")
cef_r = ast_summary[ast_summary["antibiotic"] == "Ceftazidime"]["resistance_pct"].values[0]
print(f"  Ceftazidime resistance: {cef_r}%")
print(f"  This antibiotic is the primary empirical treatment for gram-negative")
print(f"  bacteraemia in Nigerian tertiary hospitals. A {cef_r}% resistance rate")
print(f"  means frontline therapy is likely ineffective in ~9 of 10 confirmed cases.")
print(f"\n  Imipenem (carbapenem) resistance: {ast_summary[ast_summary['antibiotic']=='Imipenem']['resistance_pct'].values[0]}%")
print(f"  → Carbapenems remain active but must be protected as last-resort agents.")

print("\n✅ Co-occurrence and AST analysis complete. Run 05_visualisation.py for figures.")
