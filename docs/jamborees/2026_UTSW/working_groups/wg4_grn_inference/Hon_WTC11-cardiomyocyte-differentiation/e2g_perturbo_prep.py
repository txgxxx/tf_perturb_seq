"""
Step 1: filter calibrated trans results (empirical_pval_adj < 0.05) and save:
  - perturbo_sig_pairs.tsv  (TF, DEG pairs)
Input: calibrated trans results file (is_cis=False already filtered upstream)
"""

import os
import pandas as pd

TRANS = (
    "/nfs/turbo/umms-welchjd-code/code/wzqian/WG4_jamboree/cardiomyocyte_calibration/"
    "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq_2026_04_19_no_spacer_calibrated_trans_results.tsv"
)
OUT_DIR = "/nfs/turbo/umms-welchjd-code/code/wzqian/WG4_jamboree/GRN"
PVAL_THRESHOLD = 0.05

os.makedirs(OUT_DIR, exist_ok=True)

print("Loading calibrated trans results...")
trans = pd.read_csv(TRANS, sep="\t")
print(f"  Total rows: {len(trans):,}")

sig = trans[trans["empirical_pval_adj"] < PVAL_THRESHOLD].copy()
print(f"  Significant (empirical_pval_adj < {PVAL_THRESHOLD}): {len(sig):,}")
print(f"  Unique perturbed TFs (element_symbol): {sig['element_symbol'].nunique()}")
print(f"  Unique DEGs (tested_gene_id):          {sig['tested_gene_id'].nunique()}")

# save sig pairs: TF gene id, TF symbol, DEG gene id, DEG symbol, log2fc
pairs = sig[["element_id", "element_symbol", "tested_gene_id", "tested_gene_symbol",
             "log2fc", "log2fc_se", "empirical_pval_adj"]].copy()
pairs.columns = ["tf_element_id", "tf_symbol", "gene_id", "gene_symbol",
                 "log2fc", "log2fc_se", "empirical_pval_adj"]

pairs_path = os.path.join(OUT_DIR, "perturbo_sig_pairs.tsv")
pairs.to_csv(pairs_path, sep="\t", index=False)
print(f"  Saved: {pairs_path}")
