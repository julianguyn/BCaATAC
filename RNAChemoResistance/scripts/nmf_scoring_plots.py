import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
NC_ar = list(range(1,12)) # MODIFY: change to number of ranks assessed
NC_ar = np.array(NC_ar)

sys.path.append('./utils')
import OONMF 
import OONMFhelpers as OH


# extract reconstruction error values
Fr_ar = []
KL_ar = []
IS_ar = []
for comp in NC_ar:
    finname = '../results/jocelyn/'+str(comp)+'TotalPR.txt'
    #finname = 'signatures/'+str(comp)+'total_pr1.csv' ## MODIFY: change to directory of where the results are
    DF = pd.read_table(finname)
    Fr_ar.append(np.abs(DF.Fr.values[0]))
    KL_ar.append(DF.KL.values[0])
    IS_ar.append(DF.IS.values[0])


plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, Fr_ar, 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel('Frobenius', fontsize=30)
plt.savefig('../results/jplots/Frobenius.png', bbox_inches='tight')
plt.close()

print("Frobenius done")

plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, np.gradient(Fr_ar, NC_ar), 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel(r'$\frac{\partial Frobenius}{\partial k}$', fontsize=30)
plt.savefig('../results/jplots/partial_Frobenius.png', bbox_inches='tight')
plt.close()

print("partial_Frobenius done")

plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, KL_ar, 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel('Kullback-Leibler', fontsize=30)
plt.savefig('../results/jplots/Kullback-Leibler.png', bbox_inches='tight')
plt.close()

print("Kullback-Leibler done")

plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, np.gradient(KL_ar, NC_ar), 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel(r'$\frac{\partial Kullback-Leibler}{\partial k}$', fontsize=30)
plt.savefig('../results/jplots/partial_Kullback-Leibler.png', bbox_inches='tight')
plt.close()

print("partial_Kullback-Leibler done")

plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, IS_ar, 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel('Itakura-Saito', fontsize=30)
plt.savefig('../results/jplots/Itakura-Saito.png', bbox_inches='tight')
plt.close()

print("Itakura-Saito done")

plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, np.gradient(IS_ar, NC_ar), 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel(r'$\frac{\partial Itakura-Saito}{\partial k}$', fontsize=30)
plt.savefig('../results/jplots/partial_Itakura-Saito.png', bbox_inches='tight')
plt.close()

print("partial_Itakura-Saito done")