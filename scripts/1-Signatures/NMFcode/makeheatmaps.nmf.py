import sys
import numpy as np
import pandas as pd
import gzip
import matplotlib.pyplot as plt
sys.path.append('scripts/1-Signatures/NMFcode/utils')
import OONMFhelpers
import OONMF

rank=int(sys.argv[1])
print("working for rank",rank)

decomp = OONMF.NMFobject(rank)
decomp.matrix_input_name('data/procdata/NMF/'+str(rank)+'Rank_NNDSVD_Basis.npy','data/procdata/NMF/'+str(rank)+'Rank_NNDSVD_Mixture.npy')
decomp.read_matrix_input(compressed=True)
decomp.Basis.shape
print("matrix read")

sampnamePD = pd.read_table('sample.names.txt',header=0)
sampnamePD['full_name'] = sampnamePD.names
fullnames = sampnamePD.full_name.values
len(fullnames)

from OONMFhelpers import get_barsortorder
bar_graph_sort_order = get_barsortorder(decomp.Basis)

decomp.make_standard_heatmap_plot(410, decomp.Basis, 'rank'+str(rank)+'.png', names=np.array(fullnames), barsortorder= bar_graph_sort_order)
