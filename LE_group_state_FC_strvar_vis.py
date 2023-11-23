#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the selected k, calculate the correlation between total FC strength and variability 
across states, plot the FC strength and variability within each state as 2D matrices with limits,
and plot the FC strength and variability within each state as histograms.
Output:
dFC_str_unscaled.jpg Plot containing auto-scaled dFC strength 2D states.
dFC_str.jpg Plot containing cross-state scaled dFC strength 2D states.
colleg_dFC_str.jpg  Color bar for cross-state scaling for dFC strength 2D states.
dFC_var_unscaled.jpg Plot containing auto-scaled dFC variability 2D states.
dFC_var.jpg Plot containing cross-state scaled dFC variability 2D states.
colleg_dFC_var.jpg Color bar for cross-state scaling for dFC variability 2D states.
dFC_str_dist.jpg Histogram of dFC strength distribution across edges.
dFC_var_dist.jpg Histogram of dFC variability distribution across edges.

Usage: 
    LE_group_state_FC_strvar_vis.py <k>
    
Arguments:

    <k> k for k-clustering

"""

import os, h5py
import numpy as np
import pandas as pd
import seaborn as sns
from docopt import docopt
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import colors, rc, rcParams, rcParamsDefault
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import spearmanr

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
    
    #Sort the rows and columns.
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
    k = args['<k>']

    #Set up I/O.
    subgroup = 'full'      
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/group_strvar/')               
    os.makedirs(outpath,exist_ok=True)
        
    #Read in network labels and convert them to integers.
    with open('colenetlabels.txt') as f:
        netlabels = [label.rstrip() for label in f] 
    netlabels = list(map(int,netlabels))
    
    #Generate colorbar components from the network labels.
    groups = sorted(netlabels)
    groupstr = []
    namelist = ['Primary Visual (VIS1)','Secondary Visual (VIS2)','Somatomotor (SMN)',
                'Cingulo-Opercular (CON)','Dorsal Attention (DAN)','Language (LAN)',
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

    #Save legend.
    handles = [Patch(facecolor=lut_dict[name]) for name in lut_dict]
    plt.legend(handles,lut_dict,title='Networks',loc='center') 
    plt.axis('off')
    outfile = outpath+'netleg_dFC.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()  

    #Set inpath.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/group_strvar/') 
    
    #Do strength.
    infile = (inpath+'dFC_strvar.h5')
    instore = h5py.File(infile,'r')
    inkey = ('/dFCstr')
    totalstr = np.array(instore[inkey]).T
    instore.close()
    totalstr = np.mean(totalstr,axis=0)

    #Do variability.
    instore = h5py.File(infile,'r')
    inkey = ('/dFCvar')
    totalvar = np.array(instore[inkey]).T
    instore.close()
    totalvar = np.mean(totalvar,axis=0)

    #Correlate.
    totalcorr = spearmanr(np.abs(totalstr),totalvar)
    totalcorr = round(totalcorr.statistic,3)
    print('Total correlation:',totalcorr)

    #Read in the states and limits.
    states_str = []
    states_var = []
    limvals_str = []
    limvals_var = []
    for i in range(int(k)):

        #Do strength.
        infile = (inpath+'dFC_strvar.h5')
        instore = h5py.File(infile,'r')
        inkey = ('/dFCstr_k'+str(i+1)+'_2D')
        inmat = np.array(instore[inkey]).T
        instore.close()
        states_str.append(inmat)
        inmat = np.copy(inmat)
        inmat[np.diag_indices_from(inmat)] = float('nan')
        limvals_str.append(np.nanmin(inmat))
        limvals_str.append(np.nanmax(inmat))

        #Do variability.
        instore = h5py.File(infile,'r')
        inkey = ('/dFCvar_k'+str(i+1)+'_2D')
        inmat = np.array(instore[inkey]).T
        instore.close()
        states_var.append(inmat)
        inmat = np.copy(inmat)
        inmat[np.diag_indices_from(inmat)] = float('nan')
        limvals_var.append(np.nanmin(inmat))
        limvals_var.append(np.nanmax(inmat))
    strmin = min(limvals_str)
    strmax = max(limvals_str)
    varmin = min(limvals_var)
    varmax = max(limvals_var)

    #For each state.
    str_collect = []
    var_collect = []
    for i in range(int(k)):
        
        #Select the state, reformat, and append.
        print('Doing state:',i+1)
        cmat = states_str[i]
        reformatted = visvec(cmat,netlabels)
        str_collect.append(reformatted.iloc[1:,1:].values)
        cmat = states_var[i]
        reformatted = visvec(cmat,netlabels)
        var_collect.append(reformatted.iloc[1:,1:].values)

    #Generate labels.
    plabs = []
    for i in range(int(k)):
        plabs.append(str(i+1))
    
    #Set label size.
    labsize = 20

    #Plot each k for dFC strength, matrix unscaled.
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=int(k))
    for ax, dat in zip(axs.ravel(),str_collect):

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
    outfile = outpath+'dFC_str_unscaled.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close() 

    #Plot each k for dFC strength, matrix.
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=int(k))
    for ax, dat in zip(axs.ravel(),str_collect):

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        divnorm=colors.TwoSlopeNorm(vmin=strmin,vcenter=0,vmax=strmax)
        
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
    outfile = outpath+'dFC_str.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()  
    
    #Save color bar.
    divnorm=colors.TwoSlopeNorm(vmin=strmin,vcenter=0,vmax=strmax)
    c = plt.pcolormesh(str_collect[0],cmap='RdBu_r',norm=divnorm)
    cbar = plt.colorbar(c)
    cbar.ax.tick_params(labelsize=labsize)
    cbar.ax.set_yscale('linear')
    outfile = outpath+'colleg_dFC_str.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()  

    #Plot each k for dFC variability, matrix unscaled.
    colortype = 'twilight'
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=int(k))
    for ax, dat in zip(axs.ravel(),var_collect):

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        
        #Plot.
        im = ax.imshow(dat,cmap=colortype)
        x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
        y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
        
        #Remove axes values.
        ax.axis('off')
        x_ax.set_axis_off()
        y_ax.set_axis_off()    

        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = outpath+'dFC_var_unscaled.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()  
   
    #Plot each k for dFC variability, matrix.
    colortype = 'twilight'
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=int(k))
    for ax, dat in zip(axs.ravel(),var_collect):

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        
        #Plot.
        im = ax.imshow(dat,cmap=colortype,vmin=varmin,vmax=varmax)
        x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
        y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
        
        #Remove axes values.
        ax.axis('off')
        x_ax.set_axis_off()
        y_ax.set_axis_off()    

        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = outpath+'dFC_var.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()    
    
    #Save color bar.
    divnorm=colors.TwoSlopeNorm(vmin=varmin,vcenter=(varmin+((varmax-varmin)/2)),vmax=varmax)
    c = plt.pcolormesh(var_collect[0],cmap=colortype,norm=divnorm)
    cbar = plt.colorbar(c)
    cbar.ax.tick_params(labelsize=labsize)
    outfile = outpath+'colleg_dFC_var.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()  

    #Generate vectors.
    str_vec = []
    var_vec = []
    for i in range(int(k)):
        din = str_collect[i]
        cvec = din[np.triu_indices(360,1)]
        str_vec.append(cvec)
        din = var_collect[i]
        cvec = din[np.triu_indices(360,1)]
        var_vec.append(cvec)   
 
    #Set bin width and histogram color.
    binwidth = 0.001
    histcolor = 'silver'

    #Plot for each k for dFC strength, histogram.
    custom_xlim = (strmin,strmax)
    custom_ylim = (0,600)
    rc('xtick',labelsize=labsize) 
    rc('ytick',labelsize=labsize) 
    plt.rcParams['figure.figsize'] = (40,5)
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=int(k))
    for ax, dat in zip(axs.ravel(),str_vec):
        
        #Plot.
        plt.setp(ax,xlim=custom_xlim,ylim=custom_ylim)
        ax.hist(dat,bins=np.arange(min(dat),max(dat)+binwidth,binwidth),color=histcolor)
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        
        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = outpath+'dFC_str_dist.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()    

    #Plot each k for dFC variability, histogram.
    custom_xlim = (varmin,varmax+0.05)
    custom_ylim = (0,7500)
    rc('xtick',labelsize=labsize) 
    rc('ytick',labelsize=labsize) 
    plt.rcParams["figure.figsize"] = (40,5)
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=int(k))
    for ax, dat in zip(axs.ravel(),var_vec):
        
        #Plot.
        plt.setp(ax,xlim=custom_xlim,ylim=custom_ylim)
        ax.hist(dat,bins=np.arange(min(dat),max(dat)+binwidth,binwidth),color=histcolor)
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))

        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = outpath+'dFC_var_dist.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()    
    