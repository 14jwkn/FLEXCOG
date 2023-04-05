#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the subject file group and selected k, plot the static FC and occurrence-weighted
dFC(t) matrices.
Output:
sFCdFC.jpg Plots of the static FC and occurrence-weighted average dFC(t) matrices.
colleg_sFC.jpg Color bar for the static FC matrix.
colleg_dFC.jpg Color bar for the dFC matrix.

Usage: 
    r_LE_group_state_sFCdFC_vis.py <subfile> <k> 
    
Arguments:
    
    <subfile> Subject file
    <k> K for k-clustering

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
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    subfile = args['<subfile>']
    k = args['<k>']
    print('Doing:',subfile,k)
    
    #Set up I/O.
    if (subfile == 'r_full_submain.txt'):
        subgroup = 'full'  
    elif (subfile == 'r_half_submain.txt'):
        subgroup = 'half'         
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               '/'+subgroup+'/'+k+'/state_images/')               
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
    
    #Read in the mean sFC and dFC matrices.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/') 
    inkey = '/mean_2D'     
    infile = (inpath+'sub_mean_sFC.h5')
    instore = h5py.File(infile,'r')
    mean_sFC = np.array(instore[inkey]).T
    infile = (inpath+'sub_mean_dFC.h5')
    instore = h5py.File(infile,'r')
    mean_dFC = np.array(instore[inkey]).T
    
    #Set limits.
    inmat = np.copy(mean_sFC)
    inmat[np.diag_indices_from(inmat)] = float('nan')
    min_sFC = np.nanmin(inmat)
    max_sFC = np.nanmax(inmat)
    inmat = np.copy(mean_dFC)
    inmat[np.diag_indices_from(inmat)] = float('nan')
    min_dFC = np.nanmin(inmat)
    max_dFC = np.nanmax(inmat)
    
    #Rearrange.
    state_collect = []
    reformatted = visvec(mean_sFC,netlabels)
    state_collect.append(reformatted.iloc[1:,1:].values)
    reformatted = visvec(mean_dFC,netlabels)
    state_collect.append(reformatted.iloc[1:,1:].values)
      
    #Plot for sFC.
    pidx = 0
    fig, axs = plt.subplots(nrows=2,ncols=1)
    for ax, dat in zip(axs.ravel(),state_collect):
        
        #Set min and max.
        if pidx == 0:
            minval = min_sFC
            maxval = max_sFC
        elif pidx == 1:
            minval = min_dFC
            maxval = max_dFC

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        divnorm=colors.TwoSlopeNorm(vmin=minval,vcenter=0,vmax=maxval)
        
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
    outfile = outpath+'sFCdFC.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close() 
    
    #Save colorbar.
    divnorm = colors.TwoSlopeNorm(vmin=min_sFC,vcenter=0,vmax=max_sFC)
    c = plt.pcolormesh(state_collect[0],cmap='RdBu_r',norm=divnorm)
    cbar = plt.colorbar(c)
    cbar.ax.tick_params(labelsize=20) 
    outfile = outpath+'colleg_sFC.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close() 
    
    #Save colorbar.
    divnorm = colors.TwoSlopeNorm(vmin=min_dFC,vcenter=0,vmax=max_dFC)
    c = plt.pcolormesh(state_collect[1],cmap='RdBu_r',norm=divnorm)
    cbar = plt.colorbar(c)
    cbar.ax.tick_params(labelsize=20) 
    outfile = outpath+'colleg_dFC.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close() 
