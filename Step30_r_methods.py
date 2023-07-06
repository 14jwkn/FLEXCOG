#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize sample fMRI timeseries, dFC(t), and LE(t) for the methods description.
Output:
fMRI.jpg Contains two sample fMRI timeseries.
dFC.jpg Contains sample dFC(t).
LE1D.jpg Contains sample LE(t).
LE2D.jpg Contains sample LE(t)*LE(t)^T.

Usage: 
    Step30_r_methods.py <subject> 
    
Arguments:

    <subject> Subject ID

"""

import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from docopt import docopt

#Convert vectorized into sorted matrix.
def visvec(dFCwin,netlabels):
    
    #Format.
    corrmat = pd.DataFrame(dFCwin)
    
    #Add a column for the labels.
    rowlabelled = pd.concat([pd.Series(netlabels),pd.DataFrame(corrmat)],axis=1)
    
    #Add a row for the labels.
    colnetlabels = [0] + netlabels
    rowlabelled.loc[-1] = colnetlabels
    rowlabelled.index = rowlabelled.index + 1
    collabelled = rowlabelled.sort_index()
    collabelled.columns = range(361)
    
    #Adds axes labels to enable reference.
    collabelled = collabelled.rename_axis('Index')
    collabelled = collabelled.rename_axis('Columns',axis='columns')
    
    #Sort the rows and columns.w
    rowsort = collabelled.sort_values(by=[0,'Index'],axis=0)
    colsort = rowsort.sort_values(by=[0,'Columns'],axis=1)
    
    #Reset indices. Save the matrix to list.
    reformatted = colsort.reset_index(drop=True)
    reformatted.columns = range(reformatted.shape[1])
    return reformatted

if __name__ == '__main__':
    __spec__ = None

    #Gets subject ID.
    arguments = docopt(__doc__)
    subject = arguments['<subject>']
    print('Doing:',subject)
    
    #Mean timeseries.
    infile = '../outputs/r_meants/'+subject+'/demean_rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv'
    meants = pd.read_csv(infile,header=None)
    xscale = np.arange(0,((40)*0.72),0.72)
    region1 = meants.T.iloc[0:40,0]
    region2 = meants.T.iloc[0:40,1]
    plt.plot(xscale,region1,alpha=1,color='blue')
    plt.plot(xscale,region2,alpha=1,color='red')
    plt.xlim([-1,30])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Time Point',labelpad=5)
    plt.ylabel('Region',labelpad=5)
    outfile = '../outputs/outcollect/fMRI.jpg'
    plt.savefig(outfile,dpi=720)

    #Reads in network labels and converts them to integers.
    with open('colenetlabels.txt') as f:
        netlabels = [label.rstrip() for label in f] 
    netlabels = list(map(int,netlabels))
    
    #Plot dFC.
    infile = '../outputs/outcollect/dFC.csv'
    dfc = pd.read_csv(infile,header=None)
    dfc = visvec(dfc,netlabels)

    #Plot with labels.
    allmin = -1
    allmax = 1
    g = sns.clustermap(dfc.iloc[1:,1:],row_cluster=False,col_cluster=False,
                       cmap='RdBu_r',center=0,vmin=allmin,vmax=allmax)
    g.cax.set_position([1,.2,.03,.45])
    plt.setp(g.ax_heatmap.get_xticklabels(),rotation=90) 
    g.ax_heatmap.tick_params(left=True,bottom=True,right=False,
                             labelleft=True,labelbottom=True,labelright=False,
                             axis='both',which='major',labelsize=20)
    g.fig.suptitle('dFC',fontsize=70,ha='center',y=0.93,x=0.576)
    ax = g.ax_heatmap
    ax.set_xlabel('Region',fontsize=36,labelpad=20)
    ax.set_ylabel('Region',fontsize=36,labelpad=20)
    ax.yaxis.set_label_position('left')
    ax.axes.set_yticks([1,90,180,270,360])
    ax.axes.set_yticklabels([1,90,180,270,360])
    ax.axes.set_xticks([1,90,180,270,360])
    ax.axes.set_xticklabels([1,90,180,270,360])
    
    #Save image.
    outfile = '../outputs/outcollect/dFC.jpg'
    plt.savefig(outfile,bbox_inches='tight',dpi=720)
    plt.close() 
    
    #Plot LE.
    infile = '../outputs/outcollect/LE.csv'
    le = pd.read_csv(infile,header=None)
    le = pd.concat((le,pd.DataFrame(netlabels)),axis=1)
    le.columns = ['FC','Label']
    le = pd.DataFrame(le.rename_axis('Region').sort_values(by=['Label','Region'],axis=0).loc[:,'FC'])
    
    #Plot with labels.
    plt.rcParams["figure.figsize"] = 1,5
    x = [1,90,180,270,360]
    y = le.values
    fig, ax = plt.subplots(ncols=1)
    ax.imshow(y,cmap='RdBu_r',aspect='auto')
    ax.axes.set_yticks([1,90,180,270,360])
    ax.axes.set_yticklabels([1,90,180,270,360])
    ax.set_xticks([])
    plt.tight_layout()
    fig.suptitle('Leading Eigenvector',fontsize=30,ha='center',y=1.1,x=0.576)
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    
    #Save image.
    outfile = '../outputs/outcollect/LE1D.jpg'
    plt.savefig(outfile,bbox_inches='tight',dpi=720)
    plt.close() 
    
    #Plot with labels.
    g = sns.clustermap(le.dot(le.T),row_cluster=False,col_cluster=False,
                       cmap='RdBu_r',center=0)
    g.cax.set_position([1,.2,.03,.45])
    plt.setp(g.ax_heatmap.get_xticklabels(),rotation=90) 
    g.ax_heatmap.tick_params(left=True,bottom=True,right=False,
                             labelleft=True,labelbottom=True,labelright=False,
                             axis='both',which='major',labelsize=20)
    g.fig.suptitle('LE',fontsize=70,ha='center',y=0.93,x=0.576)
    ax = g.ax_heatmap
    ax.set_xlabel('Region',fontsize=36,labelpad=20)
    ax.set_ylabel('Region',fontsize=36,labelpad=20)
    ax.yaxis.set_label_position('left')
    ax.axes.set_yticks([1,90,180,270,360])
    ax.axes.set_yticklabels([1,90,180,270,360])
    ax.axes.set_xticks([1,90,180,270,360])
    ax.axes.set_xticklabels([1,90,180,270,360])
    
    #Save image.
    outfile = '../outputs/outcollect/LE2D.jpg'
    plt.savefig(outfile,bbox_inches='tight',dpi=720)
    plt.close() 
