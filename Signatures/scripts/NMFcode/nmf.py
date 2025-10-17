import sys
import numpy as np
import pandas as pd
import gzip
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF, non_negative_factorization

sys.path.append('./utils')
import OONMFhelpers
import OONMF

rank=sys.argv[1]
print("working for rank",int(rank))

A= pd.read_table('../../data/BCa_binary.2.matrix', header=0,index_col=0).T

# rename columns and remove genomic locations
for col in range(len(A.columns)):
    label = f"{A.columns[col]}:{A.iloc[0,col]}|{A.iloc[1,col]}"
    A.columns.values[col] = label
A = A.iloc[2:]

print(A.shape)
Nc = int(rank)
seed = 20 # (not very important for NNDSVD)
a = OONMF.NMFobject(theNcomps=Nc)
a.performNMF(data=A, randomseed=seed, theinit='nndsvd')
#a.writeNMF_CSV(Basis_foutname= str(Nc)+'Rank_NNDSVD_Basis.tsv', Mixture_foutname=str(Nc)+'Rank_NNDSVD_Mixture.tsv')
a.writeNMF(Basis_foutname= str(Nc)+'Rank_NNDSVD_Basis.npy', Mixture_foutname=str(Nc)+'Rank_NNDSVD_Mixture.npy')
print("basis")
print(a.Basis.shape)
print("mixture")
print(a.Mixture.shape)

print("precision1")
[sample_PR, total_PR] = a.precision_recall_curve(A.to_numpy())
print(total_PR.head())
total_PR['precision'] = total_PR['TP'] / (total_PR['TP'] + total_PR['FP'])
total_PR['recall'] = total_PR['TP'] / (total_PR['TP'] + total_PR['FN'])
total_PR['accuracy'] = (total_PR['TP'] + total_PR['TN']) / (total_PR['TP'] + total_PR['FP'] + total_PR['FN'] + total_PR['TN'])
total_PR['F1'] = 2*total_PR['precision']*total_PR['recall'] / (total_PR['precision']+total_PR['recall'])
total_PR=pd.DataFrame(total_PR)
total_PR.to_csv(str(Nc)+'total_pr1.csv',sep="\t")
AUPRC = np.trapz([1] + list(total_PR['recall'].values) + [0], [0] + list(total_PR['precision'].values) +[1])
print(AUPRC)
#AUPRC.to_csv(str(Nc)+'auc1.csv',sep="\t")

print("presion2")
[sample_PR, total_PR] = a.quick_precision_recall_curve(A.to_numpy())
print(total_PR.head())
total_PR['precision'] = total_PR['TP'] / (total_PR['TP'] + total_PR['FP'])
total_PR['recall'] = total_PR['TP'] / (total_PR['TP'] + total_PR['FN'])
total_PR['accuracy'] = (total_PR['TP'] + total_PR['TN']) / (total_PR['TP'] + total_PR['FP'] + total_PR['FN'] + total_PR['TN'])
total_PR['F1'] = 2*total_PR['precision']*total_PR['recall'] / (total_PR['precision']+total_PR['recall'])
total_PR=pd.DataFrame(total_PR)
total_PR.to_csv(str(Nc)+'total_pr2.csv',sep="\t")
AUPRC = np.trapz([1] + list(total_PR['recall'].values) + [0], [0] + list(total_PR['precision'].values) +[1])
print(AUPRC)
#AUPRC.to_csv(str(Nc)+'auc2.csv',sep="\t")

print("presion3")
[sample_PR, total_PR] = a.precision_recall_curveDHS(A.to_numpy())
print(total_PR.head())
total_PR['precision'] = total_PR['TP'] / (total_PR['TP'] + total_PR['FP'])
total_PR['recall'] = total_PR['TP'] / (total_PR['TP'] + total_PR['FN'])
total_PR['accuracy'] = (total_PR['TP'] + total_PR['TN']) / (total_PR['TP'] + total_PR['FP'] + total_PR['FN'] + total_PR['TN'])
total_PR['F1'] = 2*total_PR['precision']*total_PR['recall'] / (total_PR['precision']+total_PR['recall'])
total_PR=pd.DataFrame(total_PR)
total_PR.to_csv(str(Nc)+'total_pr3.csv',sep="\t")
AUPRC = np.trapz([1] + list(total_PR['recall'].values) + [0], [0] + list(total_PR['precision'].values) +[1])
print(AUPRC)
#AUPRC.to_csv(str(Nc)+'auc3.csv',sep="\t")
