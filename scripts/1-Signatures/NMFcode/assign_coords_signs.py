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

comp=int(sys.argv[2])
print("working for component",comp)

decomp = OONMF.NMFobject(rank)
decomp.matrix_input_name('data/procdata/NMF/'+str(rank)+'Rank_NNDSVD_Basis.npy','data/procdata/NMF/'+str(rank)+'Rank_NNDSVD_Mixture.npy')
decomp.read_matrix_input(compressed=True)
print(decomp.Basis.shape)
print(decomp.Mixture.shape)
print("matrix read")

sampnamePD = pd.read_table('sample.names.txt',header=0)
sampnamePD['full_name'] = sampnamePD.names
fullnames = sampnamePD.full_name.values
#len(fullnames)

coo=pd.read_table('coordinates',header=0)
nms=list(coo['coordinates'])

green_cut = np.argmax(decomp.Mixture, axis=0) == comp
print("coordinate",len(green_cut[green_cut]))
comps=np.array(nms)[green_cut]
pd.DataFrame(comps).to_csv("coordinate_comp"+str(comp)+".txt")

blue_cut = np.argmax(decomp.Basis, axis=1) == comp
print("samples",len(blue_cut[blue_cut]))
comps=np.array(fullnames)[blue_cut]
pd.DataFrame(comps).to_csv("sample_comp"+str(comp)+".txt")

print("done")
