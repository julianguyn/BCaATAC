"""
Apply trained ARCHE-G surrogate models (from ARCHE-G_surrogates.ipynb) to a new cohort.

Input is a genes x samples .cct file (first column gene symbols, log2 expression).
Genes are z-scored within the cohort to match how the models were trained.

usage: python scripts/1-Signatures/apply_ARCHE-G_surrogates.py <cohort.cct> <out_dir>
"""
import os
import sys
import joblib
import pandas as pd

dir = "/Users/julianguyen/Documents/BCaATAC"
os.chdir(dir)
sys.path.insert(0, dir)
from utils.gene_surrogate.helpers import zscore_genes

model_dir = "data/results/data/1-Signatures/geneSurrogates"

cct, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

# load in cohort and z-score genes within cohort
rna = pd.read_csv(cct, sep='\t', index_col=0).T
rna = zscore_genes(rna)

preds = pd.DataFrame(index=rna.index)
for i in range(1, 7):
    bundle = joblib.load(f"{model_dir}/arche{i}_model.joblib")
    features = bundle["features"]

    # missing genes and missing values are set to the cohort mean (0)
    missing = [f for f in features if f not in rna.columns]
    if missing:
        print(f"ARCHE{i}: {len(missing)}/{len(features)} model genes missing from cohort: {missing}")
    X = rna.reindex(columns=features).fillna(0)

    preds[f"ARCHE{i}"] = bundle["model"].predict(X)

preds.index.name = "sample"
name = os.path.basename(cct).replace(".cct", "")
preds.to_csv(f"{out_dir}/{name}_ARCHE_predictions.csv")
print(f"Saved predictions for {preds.shape[0]} samples to {out_dir}/{name}_ARCHE_predictions.csv")
