"""
Apply trained ARCHE-G surrogate models to a new cohort

Input is a genes x samples file (first column gene symbols), tab- or comma-separated.
Expression should be log2; pass "log2" as the 4th argument for linear values
(e.g. FPKM/TPM) so they are log2(x + 1) transformed first.
Genes are z-scored within the cohort to match how the models were trained.

usage: python scripts/1-Signatures/GenomicSurrogates/apply_ARCHE-G_surrogates.py <cohort> <model_dir> <out_dir> [log2]
"""
import os
import sys
import joblib
import numpy as np
import pandas as pd

dir = "/Users/julianguyen/Documents/BCaATAC"
os.chdir(dir)
sys.path.insert(0, dir)
from utils.gene_surrogate.helpers import zscore_genes


cohort, model_dir, out_dir = sys.argv[1:4]
log2 = len(sys.argv) > 4 and sys.argv[4] == "log2"
os.makedirs(out_dir, exist_ok=True)

model_dir = f"data/results/data/1-Signatures/geneSurrogates/{model_dir}"

# load in cohort (separator detected from the file), log2 if linear, and z-score genes within cohort
rna = pd.read_csv(cohort, sep=None, engine="python", index_col=0).T
if log2:
    rna = np.log2(rna + 1)
rna = zscore_genes(rna)

preds = pd.DataFrame(index=rna.index)
for i in range(1, 7):
    bundle = joblib.load(f"{model_dir}/arche{i}_model.joblib")
    features = bundle["features"]

    # missing or unexpressed (zero variance) genes and missing values are set to the cohort mean (0)
    missing = [f for f in features if f not in rna.columns]
    if missing:
        print(f"ARCHE{i}: {len(missing)}/{len(features)} model genes missing or not expressed in cohort: {missing}")
    X = rna.reindex(columns=features).fillna(0)

    preds[f"ARCHE{i}"] = bundle["model"].predict(X)

preds.index.name = "sample"
name = os.path.splitext(os.path.basename(cohort))[0]
preds.to_csv(f"{out_dir}/{name}_ARCHE_predictions.csv")
print(f"Saved predictions for {preds.shape[0]} samples to {out_dir}/{name}_ARCHE_predictions.csv")
