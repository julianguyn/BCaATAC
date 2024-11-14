import sys
import numpy as np
import pandas as pd
import gzip
import matplotlib.pyplot as plt
sys.path.append('./utils')
import OONMFhelpers
import OONMF
from OONMFhelpers import get_barsortorder1

rank=int(sys.argv[1])
#rank=20
print("working for rank",rank)

decomp = OONMF.NMFobject(rank)
decomp.matrix_input_name(str(rank)+'Rank_NNDSVD_Basis.npy',str(rank)+'Rank_NNDSVD_Mixture.npy')
decomp.read_matrix_input(compressed=True)
print(decomp.Mixture.shape)
print("matrix read")

sampnamePD = pd.read_table('sample.names.txt',header=0)
sampnamePD['full_name'] = sampnamePD.names
fullnames = sampnamePD.full_name.values
print(len(fullnames))

coo=pd.read_table('coordinates',header=0)
nms=list(coo['coordinates'])
#barsortorder =[ 8 19 16 13 14 10 18  9  2  1 17  5  4  3 15 12  7  0  6 11]
#n=decomp.Mixture[0:20,0:1000]

n=decomp.Mixture[0:rank,np.sum(decomp.Mixture, axis=0)>0]
print("here")
print(n.shape)
barsortorder1 = get_barsortorder1(n.T,nms)
print(barsortorder1)
print("sorting done")

decomp.normalize_matrices()
#nm=decomp.NormedMixture[0:20,0:1000]
#decomp.make_stacked_bar_plot(2117181, decomp.Mixture.T, 'coordinate_rank'+str(rank)+'.png',names = [],official_order=False, barsortorder=barsortorder)
nm=decomp.NormedMixture
#decomp.make_stacked_bar_plot(1870245, nm, 'coordinate_rank'+str(rank)+'.png',names = nms,official_order=False, barsortorder=barsortorder1)
decomp.make_stacked_bar_plot(1501209, nm, 'coordinate_rank'+str(rank)+'.png',names = nms,official_order=False, barsortorder=barsortorder1)

