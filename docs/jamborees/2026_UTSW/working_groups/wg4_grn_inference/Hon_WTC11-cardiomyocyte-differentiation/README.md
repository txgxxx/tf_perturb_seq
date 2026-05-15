# Hon WTC11 Cardiomyocyte Differentiation GRN Scripts

This folder records the scripts used for the WG4 GRN analysis connecting
Perturbo trans-effects with E2G enhancer-gene links in the Hon WTC11
cardiomyocyte differentiation dataset.

## Scripts

| Script | Used for |
|--------|----------|
| `e2g_perturbo_prep.py` | Filters calibrated Perturbo trans results to significant TF-to-DEG pairs using `empirical_pval_adj < 0.05`, then writes `perturbo_sig_pairs.tsv` with TF IDs, DEG IDs, gene symbols, log2 fold change, standard error, and adjusted p-values. |
| `e2g_deg_overlap_worker.py` | Takes one E2G file and joins significant Perturbo DEGs to all E2G enhancer-gene links for those genes, producing a detailed `(TF, DEG, enhancer)` table and appending summary overlap statistics. |

## Workflow

1. Run `e2g_perturbo_prep.py` to create the significant Perturbo TF-to-DEG pair table.
2. Run `e2g_deg_overlap_worker.py` once per E2G input file to generate detailed enhancer-link overlaps and summary coverage metrics.
