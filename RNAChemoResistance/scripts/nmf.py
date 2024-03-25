import sys
import numpy as np
import pandas as pd
import gzip
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF, non_negative_factorization
import scipy.sparse as sp

sys.path.append('./utils')
import OONMFhelpers
import OONMF

rank=sys.argv[1]
print("working for rank",int(rank))

#A= pd.read_table('../procdata/BCa_binary.2.matrix', header=0,index_col=0).T ## MODIFY: add path to matrix
A= pd.read_table('../procdata/adjcounts.tsv', header=0,index_col=0).T

# Remove the last 6 rows containing genomic information
#A = A.iloc[:-6]

# Rename columns
#for col in range(len(A.columns)):
#    label = f"{A.columns[col]}:{A.iloc[0,col]}|{A.iloc[1,col]}"
#    A.columns.values[col] = label

# Remove first two rows containing genomic locations
#A = A.iloc[2:]
#/cluster/projects/bhklab/references/GENCODE/human/GRCh38_v45/genome.fa


print(A.shape)
Nc = int(rank)
seed = 20 # (not very important for NNDSVD)
a = OONMF.NMFobject(theNcomps=Nc)
a.performNMF_betadivergence(data=A, randomseed=seed, theinit='nndsvd', filename_addon=rank)
#a.writeNMF_CSV(Basis_foutname= str(Nc)+'Rank_NNDSVD_Basis.tsv', Mixture_foutname=str(Nc)+'Rank_NNDSVD_Mixture.tsv')
a.writeNMF(foutname= "../results/jocelyn/"+str(Nc))
print("basis")
print(a.Basis.shape)
print("mixture")
print(a.Mixture.shape)