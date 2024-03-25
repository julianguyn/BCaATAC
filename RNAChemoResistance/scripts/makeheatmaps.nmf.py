import sys
import numpy as np
import pandas as pd
import gzip
import matplotlib.pyplot as plt
sys.path.append('./utils')
import OONMFhelpers
import OONMF

rank=int(sys.argv[1])
print("working for rank",rank)

decomp = OONMF.NMFobject(rank)
decomp.matrix_input_name('../results/jocelyn/'+str(rank)+'Rank_NNDSVD_Basis.npy','../results/jocelyn/'+str(rank)+'Rank_NNDSVD_Mixture.npy') # MODIFY: set to path of npy files
decomp.read_matrix_input(compressed=True)
decomp.Basis.shape
print("matrix read")

sampnamePD = pd.read_table('../metadata/jocelyn/samples.txt',header=0) # MODIFY: set to a sample file
sampnamePD['full_name'] = sampnamePD.samples
fullnames = sampnamePD.full_name.values
len(fullnames)

from OONMFhelpers import get_barsortorder
bar_graph_sort_order = get_barsortorder(decomp.Basis)

decomp.make_standard_heatmap_plot(410, decomp.Basis, '../results/jocelyn/heatmap_rank'+str(rank)+'.png', names=np.array(fullnames), barsortorder= bar_graph_sort_order)
