#!/usr/bin/env python3
"""
add_pvals_to_IPA_sheet.py
-------------------------
• logFC exactly as in the previous script (freq‑normalised raw counts, minus sign).
• Two‑sided Fisher’s exact test per gene.
• Benjamini–Hochberg FDR correction.
• Output has Gene Name, logFC, pval, qval (FDR).
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ─────────────────────────────── CONFIG ──────────────────────────────── #
INFILE  = "hyPB_cordycepin_200_significant_184_genes.xlsx"
OUTFILE = "cordycepin_vs_DMSO_logFC_pval.xlsx"

LIBSIZE_CORDY = 51_053
LIBSIZE_DMSO  = 37_273
PSEUDO        = 0.5   # same half‑insertion pseudocount for logFC only
# ─────────────────────────────────────────────────────────────────────── #


def main() -> None:
    df = pd.read_excel(INFILE)

    # logFC (library-size-normalised, direction flipped)
    freq_c = (df["Insertions"] + PSEUDO) / (LIBSIZE_CORDY + PSEUDO)
    freq_d = (
        df["hyPB_Control_DMSO Insertions"] + PSEUDO
    ) / (LIBSIZE_DMSO + PSEUDO)
    df["logFC"] = -np.log2(freq_c / freq_d)  # minus makes "more insertions" negative

    # Fisher's exact test
    pvals = []
    for a, b in zip(df["Insertions"], df["hyPB_Control_DMSO Insertions"]):
        table = [[a, b], [LIBSIZE_CORDY - a, LIBSIZE_DMSO - b]]
        pvals.append(fisher_exact(table, alternative="two-sided")[1])

    df["pval"] = pvals

    # Benjamini–Hochberg FDR
    df["qval"] = multipletests(df["pval"], method="fdr_bh")[1]

    # save the four columns IPA can ingest (it ignores extras)
    df[["Gene Name", "logFC", "pval", "qval"]].to_excel(OUTFILE, index=False)
    print(f"✓ Wrote {OUTFILE}")


if __name__ == "__main__":
    main()
