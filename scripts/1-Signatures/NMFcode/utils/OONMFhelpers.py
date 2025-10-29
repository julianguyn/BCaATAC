ClusterMode = False

import datetime
import time
from datetime import date
import numpy as np
import pandas as pd
today = str(date.today())

'''
functions:

mytime - returns the time as a string

get_today - returns date as string

increase_axis_fontsize - increases the label and tick size on the current plot

get_barsortorder - this is for ordering our stacked bar graphs. It basically sorts the vectors by component, one at a time. 

get_barsortorder1 - this is for ordering our stacked bar graphs.

get_barsortorder_OfficialOrder - probably does not deserve to be a function on its own, but whatever. It illustrates how to make an arbitrary order vs. natural order.

'''

if (ClusterMode):
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

def mytime():
    ts = time.time()
    sttime = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H:%M:%S')
    return sttime
    
def get_today():
	today = str(date.today())
	return today

def increase_axis_fontsize(axis_fontsize=30):
    ax = plt.gca()
    ticklabels = ax.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(axis_fontsize)
        label.set_family('serif')
    ticklabels = ax.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(axis_fontsize)
        label.set_family('serif')

def get_barsortorder(relevantMatrix):
    # assumes rows are the data, columns are NMF components
    WinningComponent = np.argmax(relevantMatrix, axis=1)
    barsortorder = np.array([])
    for i in range(relevantMatrix.shape[1]):
        desired_order = np.argsort(-relevantMatrix[:,i])
        relevant_cut = WinningComponent[desired_order]==i
        barsortorder = np.append(barsortorder, desired_order[relevant_cut])
    barsortorder = barsortorder.astype(int)
    return barsortorder


def get_barsortorder1(relevantMatrix,nms):
    # assumes rows are the data, columns are NMF components
    print("from sorter")
    print(relevantMatrix[1:5,1:5])
    WinningComponent = np.argmax(relevantMatrix, axis=1)
    #print("wincomp",WinningComponent,len(WinningComponent))
    barsortorder = np.array([])
    #print(relevantMatrix.shape,"shape")
    for i in range(relevantMatrix.shape[1]):
        desired_order = np.argsort(-relevantMatrix[:,i])
        #print("desired_order",desired_order)
        relevant_cut = WinningComponent[desired_order]==i
        #print("relevant_cut",relevant_cut)
        #print(len(relevant_cut),"shape")
        #print(i,"number")
        #m=pd.DataFrame(data=relevantMatrix[desired_order[relevant_cut]])
        #print("mshape",m.shape)
       #m.to_csv(str(i)+'component.mat', sep='\t', header=True, float_format='%.2f', index=False)
        barsortorder = np.append(barsortorder, desired_order[relevant_cut])
        #print("barsortorder",barsortorder)
        n = [x for _,x in sorted(zip(barsortorder,nms))]
        #print(len(n),"n")
        mat=relevantMatrix[desired_order[relevant_cut]].T
        print(mat.shape,"component matrix shape")
        ncol=mat.shape[1]
        m=pd.DataFrame(data=mat,columns=n[0:ncol])
        m.to_csv(str(i)+'component.mat', sep='\t', header=True, float_format='%.2f', index=False)
    barsortorder = barsortorder.astype(int)
    return barsortorder


def get_barsortorder_OfficialOrder(relevantMatrix):
    # much more ugly but gets the job done
    # assumes rows are the data, columns are NMF components
    WSO = np.array([7,5,15,9,12,14,3,8,13,2,4,6,16,11,10,1]).astype(int) - 1
    WinningComponent = np.argmax(relevantMatrix, axis=1)
    barsortorder = np.array([])
    for i in range(relevantMatrix.shape[1]):
        desired_order = np.argsort(-relevantMatrix[:,WinningComponent[WSO[i]]])
        relevant_cut = WinningComponent[desired_order]==WSO[i]
        barsortorder = np.append(barsortorder, desired_order[relevant_cut])
    barsortorder = barsortorder.astype(int)
    return barsortorder
	
