"""
For each sig perturbo (TF -> DEG) pair, find all E2G enhancer links for that DEG
and sort by E2G score. Outputs:
  - detailed_{label}.tsv.gz : one row per (TF, DEG, enhancer) triple, sorted by score desc
  - deg_overlap_result_{label}.tsv : one-line summary appended to shared summary file

Usage: python e2g_deg_overlap_worker.py <label> <e2g_file_path>
"""

import os
import sys
import pandas as pd

label    = sys.argv[1]
e2g_path = sys.argv[2]

OUT_DIR    = "/nfs/turbo/umms-welchjd-code/code/wzqian/WG4_jamboree/GRN"
PAIRS_PATH = os.path.join(OUT_DIR, "perturbo_sig_pairs.tsv")
SUMMARY    = os.path.join(OUT_DIR, "e2g_overlap_summary.tsv")

print(f"[{label}] Loading sig pairs...")
sig = pd.read_csv(PAIRS_PATH, sep="\t")
sig_degs = set(sig["gene_id"].unique())
total_pairs = len(sig)
print(f"  Sig (TF, DEG) pairs: {total_pairs:,}  |  Unique DEGs: {len(sig_degs):,}")

print(f"[{label}] Loading E2G...")
e2g = pd.read_csv(e2g_path, sep="\t", comment="#")
print(f"  E2G rows: {len(e2g):,}  |  Unique genes: {e2g['GeneEnsemblID'].nunique():,}")

has_score = "Score" in e2g.columns

# ── join: sig pairs <- E2G on DEG gene id ────────────────────────────────────
e2g_cols = ["GeneEnsemblID", "ElementChr", "ElementStart", "ElementEnd",
            "ElementName", "ElementClass", "GeneSymbol"]
if has_score:
    e2g_cols += ["Score"]

e2g_sub = e2g[e2g_cols].copy()
e2g_sub = e2g_sub.rename(columns={
    "GeneEnsemblID": "gene_id",
    "ElementChr":    "enhancer_chr",
    "ElementStart":  "enhancer_start",
    "ElementEnd":    "enhancer_end",
    "ElementName":   "enhancer_name",
    "ElementClass":  "enhancer_class",
    "GeneSymbol":    "gene_symbol_e2g",
})
if has_score:
    e2g_sub = e2g_sub.rename(columns={"Score": "e2g_score"})

# merge: each sig pair gets all E2G enhancer rows for that DEG
detailed = sig.merge(e2g_sub, on="gene_id", how="inner")
if has_score:
    detailed = detailed.sort_values(
        ["tf_symbol", "gene_id", "e2g_score"], ascending=[True, True, False]
    )
else:
    detailed = detailed.sort_values(["tf_symbol", "gene_id"])

# ── save detailed file ────────────────────────────────────────────────────────
detail_path = os.path.join(OUT_DIR, f"detailed_{label}.tsv.gz")
detailed.to_csv(detail_path, sep="\t", index=False, compression="gzip")
print(f"[{label}] Detailed file: {detail_path}  ({len(detailed):,} rows)")

# ── summary stats ─────────────────────────────────────────────────────────────
degs_in_e2g   = set(detailed["gene_id"].unique())
pairs_covered = detailed[["tf_element_id", "gene_id"]].drop_duplicates()
pct_degs  = 100 * len(degs_in_e2g) / len(sig_degs)  if sig_degs  else 0
pct_pairs = 100 * len(pairs_covered) / total_pairs   if total_pairs else 0
med_score = detailed["e2g_score"].median() if has_score else float("nan")

print(f"[{label}] DEGs covered:       {len(degs_in_e2g):,} / {len(sig_degs):,} ({pct_degs:.2f}%)")
print(f"[{label}] Pairs covered:      {len(pairs_covered):,} / {total_pairs:,} ({pct_pairs:.2f}%)")
if has_score:
    print(f"[{label}] Median E2G score:  {med_score:.4f}")

# ── append one row to shared summary ─────────────────────────────────────────
row = pd.DataFrame([{
    "label":               label,
    "e2g_rows":            len(e2g),
    "e2g_unique_genes":    e2g["GeneEnsemblID"].nunique(),
    "sig_pairs":           total_pairs,
    "unique_degs":         len(sig_degs),
    "degs_in_e2g":         len(degs_in_e2g),
    "pct_degs_covered":    round(pct_degs, 4),
    "pairs_with_e2g":      len(pairs_covered),
    "pct_pairs_covered":   round(pct_pairs, 4),
    "median_e2g_score":    med_score,
}])
# write header only if file doesn't exist yet
write_header = not os.path.exists(SUMMARY)
row.to_csv(SUMMARY, sep="\t", index=False, mode="a", header=write_header)
print(f"[{label}] Appended to summary: {SUMMARY}")
