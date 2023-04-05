#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the subject file group and selected k, plot the 2D dFC or LE 2D states to
display.
Output:
'k'+statetype+'_unscaled.jpg' Plot containing auto-scaled dFC or LE 2D states.
'k'+statetype+'.jpg' Plot containing cross-k scaled dFC or LE 2D states.
'netleg_k'+statetype+'.jpg' Plot containing the network labels.
'colleg_k'+statetype+'.jpg' Plot containing the color bar for cross-k scaling.
'k'+statetype+'_dist.jpg' Histograms of edge strength values for each state.

Usage: 
    31_r_LE_group_state_vis.py <subgroup> <k> <statetype>
    
Arguments:
    
    <subfile> Subject file
    <k> K for k-clustering
    <statetype> dFC or LE

"""

import os, h5py
import numpy as np
import pandas as pd
import seaborn as sns
from docopt import docopt
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.patches import Patch
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Convert vectorized into sorted matrix.
def visvec(cwin,netlabels):
    
    #Format.
    corrmat = pd.DataFrame(cwin)
    
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
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    subfile = args['<subfile>']
    k = args['<k>']
    statetype = args['<statetype>']
    print('Doing:',subfile,k,statetype)
    
    #Set up I/O.
    if (subfile == 'r_full_submain.txt'):
        subgroup = 'full'  
    elif (subfile == 'r_half_submain.txt'):
        subgroup = 'half'    
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/state_images/')               
    os.makedirs(outpath,exist_ok=True)
        
    #Read in network labels and convert them to integers.
    with open('colenetlabels.txt') as f:
        netlabels = [label.rstrip() for label in f] 
    netlabels = list(map(int,netlabels))
    
    #Generate colorbar components from the network labels.
    groups = sorted(netlabels)
    groupstr = []
    namelist = ['Primary Visual (VIS1)','Secondary Visual (VIS2)',
                'Somatomotor (SMN)','Cingulo-Opercular (CON)',
                'Dorsal Attention (DAN)','Language (LAN)',
                'Frontoparietal (FPN)','Auditory (AUD)','Default Mode (DMN)',
                'Posterior Multimodal (PMM)','Ventral Multimodal (VMM)',
                'Orbitoaffective (ORA)']
    for net in groups:
        groupstr.append(namelist[int(net)-1])
    lut_colors = list(sns.color_palette(palette='tab10',n_colors=len(namelist)))
    lut_dict = OrderedDict(zip(namelist,sns.color_palette(palette='tab10',n_colors=len(namelist))))
    cmap = colors.LinearSegmentedColormap.from_list('Networks',lut_colors,N=len(namelist))
    class_labels = np.array(groups).reshape((1,360))
    class_labels = np.subtract(class_labels,1)
    
    #Read file with sorting.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/'+
              subgroup+'/allk/')       
    infile = (inpath+subgroup+'_sortclass.csv')
    reorder_mat = pd.read_csv(infile,header=0,index_col=0)
    reorder_mat.index = [str(x) for x in reorder_mat.index]
    reorder_mat.columns = [str(x) for x in reorder_mat.columns] 
    
    #Read in the states and limits.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
             subgroup+'/'+k+'/') 
    states = []
    limvals = []
    for i in range(int(k)):

        #Read in the LE or dFC states.
        infile = (inpath+'k_'+statetype+'.h5')
        instore = h5py.File(infile,'r')
        inkey = ('/'+subgroup+'_'+str(i+1))
        inmat = np.array(instore[inkey]).T
        states.append(inmat)
        inmat = np.copy(inmat)
        inmat[np.diag_indices_from(inmat)] = float('nan')
        limvals.append(np.nanmin(inmat))
        limvals.append(np.nanmax(inmat))
        instore.close()
        
    #Read the shuffler and shuffle if it the order has been picked.
    reorder_mat.loc[k,:].dropna().values
    reorder_k = reorder_mat.loc[k,:].dropna().values.astype(int).tolist()
    if reorder_k != []:
        states = [states[i-1] for i in reorder_k]
      
    #For each state.
    state_collect = []
    for i in range(int(k)):
        
        #Select the state, reformat, and append.
        print('Doing state:',i+1)
        cwin = states[i]
        reformatted = visvec(cwin,netlabels)
        state_collect.append(reformatted.iloc[1:,1:].values)
    
    #Plot each k for state, unscaled.
    pidx = 0
    fig, axs = plt.subplots(nrows=int(k),ncols=1)
    for ax, dat in zip(axs.ravel(),state_collect):

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        divnorm=colors.TwoSlopeNorm(vcenter=0)
        
        #Plot.
        im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
        x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
        y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
        
        #Remove axes values.
        ax.axis('off')
        x_ax.set_axis_off()
        y_ax.set_axis_off()    
        
        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = outpath+'k'+statetype+'_unscaled.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close() 
    
    #Find limits.
    cmin = min(limvals)   
    cmax = max(limvals)       
    
    #Generate labels.
    plabs = []
    for i in range(int(k)):
        plabs.append(str(i+1))
        
    #Plot each state, scaled.
    pidx = 0
    fig, axs = plt.subplots(nrows=int(k),ncols=1)
    for ax, dat in zip(axs.ravel(),state_collect):

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        divnorm=colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
        
        #Plot.
        im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
        x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
        y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
        
        #Remove axes values.
        ax.axis('off')
        x_ax.set_axis_off()
        y_ax.set_axis_off()    
        
        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = outpath+'k'+statetype+'.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close() 
    
    #Save legend.
    handles = [Patch(facecolor=lut_dict[name]) for name in lut_dict]
    plt.legend(handles,lut_dict,title='Networks',loc='center') 
    plt.axis('off')
    outfile = outpath+'netleg_k'+statetype+'.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()  
    
    #Save color bar.
    plt.imshow(state_collect[0],cmap='RdBu_r',vmin=cmin,vmax=cmax)    
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    outfile = outpath+'colleg_k'+statetype+'.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()  

    #Re-vectorize.
    statevec_collect = []
    for i in range(int(k)):
        
        #Select the state, reformat, and append.
        print('Doing state:',i+1)
        cwin = states[i]
        statevec_collect.append(cwin[np.triu_indices_from(cwin,k=1)])    
    
    #Generate labels.
    plabs = []
    for i in range(int(k)):
        plabs.append(str(i+1))
    
    #Set bin width.
    binwidth = 0.001
    
    #Set value limits.
    custom_xlim = (-0.055,0.20)
    custom_ylim = (0,5000)
    rc('xtick',labelsize=20) 
    rc('ytick',labelsize=20) 
   
    #Plot each k for dFC.
    plt.rcParams["figure.figsize"] = (5,30)
    pidx = 0
    fig, axs = plt.subplots(nrows=int(k),ncols=1)
    for ax, dat in zip(axs.ravel(),state_collect):
        
        #Plot.
        plt.setp(ax,xlim=custom_xlim,ylim=custom_ylim)
        ax.hist(dat,bins=np.arange(min(dat),max(dat)+binwidth,binwidth))
        
        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = outpath+'k'+statetype+'_dist.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close() 
    