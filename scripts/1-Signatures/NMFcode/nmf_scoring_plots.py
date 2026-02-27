import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
NC_ar = list(range(2,14))
NC_ar = np.array(NC_ar)

sys.path.append('scripts/1-Signatures/NMFcode/utils')
import OONMF 
import OONMFhelpers as OH

plt.clf()
plt.figure(figsize=(10,10))
prAUC_ar = []
for comp in NC_ar:
    finname = str(comp)+'total_pr1.csv'
    DF = pd.read_table(finname)
    prec = DF.TP.values / (DF.TP.values + DF.FP.values)
    recall = DF.TP.values / (DF.TP.values + DF.FN.values)
    F1 = 2*prec*recall / (prec+recall)
    prAUC_ar.append(np.trapz([1] + list(recall) + [0], [0] + list(prec) +[1]))
    plt.plot(DF.threshold.values, F1, '-k')

    
plt.plot([0.35, 0.35], [0.2, 0.8], '--r')
plt.plot([0.4, 0.4], [0.2, 0.8], ':r')
plt.plot([0.3, 0.3], [0.2, 0.8], ':r')
OH.increase_axis_fontsize()
plt.xlabel('decision boundary', fontsize=30)
plt.ylabel('F1 score', fontsize=30)
plt.savefig('F1_score.png', bbox_inches='tight')
plt.close()

print("f1 done")

plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, prAUC_ar, 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel('AUPRC', fontsize=30)
plt.savefig('AUPRC.png', bbox_inches='tight')
plt.close()

print("auprc done")

plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, np.gradient(prAUC_ar, NC_ar), 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel(r'$\frac{\partial AUPRC}{\partial k}$', fontsize=30)
plt.savefig('partial_AUPRC.png', bbox_inches='tight')
plt.close()

print("partial_AUPRC done")

F1_ar733 = []
for comp in NC_ar:
    finname = str(comp)+'total_pr1.csv'
    DF = pd.read_table(finname)
    prec = DF.TP.values / (DF.TP.values + DF.FP.values)
    recall = DF.TP.values / (DF.TP.values + DF.FN.values)
    F1 = 2*prec*recall / (prec+recall)
    plt.plot(DF.threshold.values, F1, '--b')
    plt.savefig('df_threshold.png', bbox_inches='tight')
    F1_ar733.append(F1[6])

F1_ar733 = np.array(F1_ar733)
plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, F1_ar733, 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel('F1', fontsize=30)
plt.savefig('F1_values.png', bbox_inches='tight')
plt.close()

print("F1 values done")


plt.clf()
plt.figure(figsize=(10,10))
plt.plot(NC_ar, np.gradient(F1_ar733, NC_ar), 'ok')
OH.increase_axis_fontsize()
plt.xlabel('k', fontsize=30)
plt.ylabel(r'$\frac{\partial F1}{\partial k}$', fontsize=30)
plt.savefig('F1_partial.png', bbox_inches='tight')
plt.close()

print("F1 partial done")
