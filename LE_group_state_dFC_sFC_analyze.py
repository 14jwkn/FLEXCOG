#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the selected k, generate the state-wise sFC matrices, average dFC states, full sFC matrix, and
occurrence-weighted average dFC matrix. Save and plot the matrices, and find the correlations between 
sFC and dFC matrices.
Output:
'S'+str(kidx+1)+'_sFC.csv' State-wise sFC Pearson correlation matrices for the timepoints classified to one dFC state.
dFC_sFC.jpg Plots of the state-wise sFC matrices and average dFC states.
dFC_sFC_corr.csv Correlations between state-wise sFC matrices and average dFC states.
'dFC_sFC_'+str(kidx+1)+'_1D.csv' Vectorized state-wise sFC matrices and average dFC states for normality tests.
full_sFC.csv Full sFC matrix across all time points.
full_dFC.csv Occurrence-weighted average dFC matrix.
full_dFC_sFC.jpg Plots of the full sFC and occurrence-weighted average dFC matrices.
full_dFC_sFC_corr.csv Correlation between full sFC and occurrence-weighted average dFC matrices.
dFC_sFC_full_1D.csv Vectorized full sFC and occurrence-weighted average dFC matrices for normality tests.

Usage: 
    r_LE_group_state_dFC_sFC_analyze.py <k>
    
Arguments:

    <k> K for k-clustering

"""

import os, h5py, sys
import numpy as np
import pandas as pd
import seaborn as sns
from docopt import docopt
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Convert into sorted matrix.
def visform(dFCwin,netlabels):
    
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

#Convert from vectorized to 2D.
def vec2mat(invec,vecsplit):
    
    #Initialize.
    invec = invec.values
    outmat = np.zeros((360,360))
    
    #For each cell.
    for vidx,[vlab1,vlab2] in enumerate(vecsplit):
        
        #Convert to Python indexing.
        pyth1 = int(vlab1) - 1
        pyth2 = int(vlab2) - 1
        
        #Put the value in the current cells.
        outmat[pyth1,pyth2] = invec[vidx,0]
        outmat[pyth2,pyth1] = invec[vidx,0]
    
    #Return matrix.
    return outmat

if __name__ == '__main__':
    __spec__ = None

    #Catches command-line arguments.
    args = docopt(__doc__)
    k = args['<k>']

    #Set output file name. 
    subgroup = 'full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/dFC_sFC/') 

    #Set parameters.
    nk = int(k)
    nts = 1200 - 2
    nroi = 360
    nconn = int((nroi*(nroi-1))/2)
    labsize = 20

    #Read in vectorized labels.
    veclabels = pd.read_csv('atlasMvec2.csv',header=None).values.tolist()
    vecsplit = []
    for [vecstr] in veclabels:
        vecsplit.append(vecstr.split('-'))

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

    #Read subjects.
    subfile = 'r_full_submain.txt'   
    with open(subfile) as f:
       subjects = [subject.rstrip() for subject in f]
    nsub = len(subjects)

    #Set runs.
    runs = ['REST1_LR','REST1_RL','REST2_LR','REST2_RL']
    nruns = len(runs)
    subts = nts*nruns

    #Read in universal clustering.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/')
    infile = (inpath+'uni_subcent.h5')
    inkey = ('/'+subgroup)
    store = pd.HDFStore(infile,'r')
    fullclust = pd.DataFrame(store.select(inkey)).values
    store.close()

    #Read in all the time points.
    fullts = np.zeros(((subts*nsub),nroi))
    for subidx in range(nsub):
        csub = subjects[subidx]
        print(csub)

        #Set starting and ending subject index.
        sub_startidx = (subidx)*subts
        sub_endidx = sub_startidx + subts

        #Read each run.
        for runidx in range(nruns):
            crun = runs[runidx]

            #Set starting and ending run index.
            run_startidx = sub_startidx + (runidx)*nts
            run_endidx = run_startidx + nts

            #Use run label to read file and remove first and last time points to be consistent with dFC.
            inpath = ('../outputs/r_meants/'+csub+'/')
            infile = (inpath+'demean_rfMRI_'+crun+'_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv')
            inmat = pd.read_csv(infile,header=None).T
            inmat = inmat.tail(-1)
            inmat = inmat.head(-1)

            #Append.
            fullts[run_startidx:run_endidx,:] = inmat

    #Read in time points for each cluster and produce sorted Pearson's correlation matrices.
    sFC_collect = []
    sFC_lim = []
    for kidx in range(nk):

        #Extract the time points.
        kts = fullts[fullclust[:,0] == (kidx+1),:]

        #Find the sFC matrix.
        k_sFC = np.corrcoef(kts.T)

        #Sort it.
        k_sFC = visform(k_sFC,netlabels)
        k_sFC = k_sFC.iloc[1:,1:].values
        sFC_collect.append(k_sFC)

        #Save it.
        outfile = (outpath+'S'+str(kidx+1)+'_sFC.csv')
        pd.DataFrame(k_sFC).to_csv(outfile,index=None,header=None)

        #Find limits.
        k_sFC = np.copy(k_sFC)
        k_sFC[np.diag_indices_from(k_sFC)] = float('nan')
        sFC_lim.append(np.nanmin(k_sFC))
        sFC_lim.append(np.nanmax(k_sFC))
    
    #Find max and min.
    sFC_max = np.max(sFC_lim)
    sFC_min = np.min(sFC_lim)    

    #Plot each k.
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=int(k))
    for ax, dat in zip(axs.ravel(),sFC_collect):

        #Make the diagonal zero.
        dat[np.diag_indices_from(dat)] = 0

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        b_ax = divider.append_axes('bottom',size='5%',pad=0.05)
        divnorm=colors.TwoSlopeNorm(vmin=sFC_min,vcenter=0,vmax=sFC_max)
        
        #Plot.
        im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
        x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
        y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
        cbar = fig.colorbar(im,cax=b_ax,orientation='horizontal')
        ticks = [sFC_min,0,sFC_max]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{tick:.3f}' for tick in ticks])
        cbar.ax.tick_params(labelsize=5,rotation=270,width=0.5)
        cbar.outline.set_linewidth(0.5)
        
        #Remove axes values.
        ax.axis('off')
        x_ax.set_axis_off()
        y_ax.set_axis_off()    
        
        #Increment.
        pidx = pidx + 1

    #Save.
    outfile = (outpath+'dFC_sFC.jpg')
    plt.savefig(outfile,dpi=720)
    plt.close() 

    #Read file with sorting for actual average dFC states.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/'+
              subgroup+'/allk/')       
    infile = (inpath+subgroup+'_sortclass.csv')
    reorder_mat = pd.read_csv(infile,header=0,index_col=0)
    reorder_mat.index = [str(x) for x in reorder_mat.index]
    reorder_mat.columns = [str(x) for x in reorder_mat.columns] 

    #Extract the actual average dFC states.
    dFC_collect = []
    for kidx in range(nk):

        #Read in the actual average dFC state and sort.
        inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                  subgroup+'/'+k+'/') 
        infile = (inpath+'k_dFC.h5')
        instore = h5py.File(infile,'r')
        inkey = ('/'+subgroup+'_'+str(kidx+1))
        inmat = np.array(instore[inkey]).T
        curr_dFC = visform(inmat,netlabels).iloc[1:,1:].values
        dFC_collect.append(curr_dFC)
    
    #Reorder actual average dFC states.
    reorder_mat.loc[k,:].dropna().values
    reorder_k = reorder_mat.loc[k,:].dropna().values.astype(int).tolist()
    if reorder_k != []:
        dFC_collect = [dFC_collect[i-1] for i in reorder_k]
    
    #Find the correlations between the generated matrices and actual average dFC states.
    dFC_sFC_corr = pd.DataFrame(np.zeros((nk,2)),index=[('S'+str(x+1)) for x in range(nk)],
                                                 columns=['Pearson','Spearman'])
    for kidx in range(nk):

        #Extract current sorted sFC and dFC.
        curr_sFC = sFC_collect[kidx]
        curr_dFC = dFC_collect[kidx]

        #Vectorize both matrices, correlate, and append.
        curr_sFC_1D = curr_sFC[np.triu_indices(nroi,k=1)]
        curr_dFC_1D = curr_dFC[np.triu_indices(nroi,k=1)]
        dFC_sFC_corr.loc[('S'+str(kidx+1)),'Pearson'] = pd.DataFrame((curr_sFC_1D,curr_dFC_1D)).T.corr().values[0,1]
        dFC_sFC_corr.loc[('S'+str(kidx+1)),'Spearman'] = pd.DataFrame((curr_sFC_1D,curr_dFC_1D)).T.corr(method='spearman').values[0,1]

        #Save the vectorize matrices for normality tests.
        outfile = (outpath+'dFC_sFC_'+str(kidx+1)+'_1D.csv')
        outmat = pd.DataFrame((curr_dFC_1D,curr_sFC_1D))
        outmat.to_csv(outfile)

    #Save it.
    outfile = (outpath+'dFC_sFC_corr.csv')
    dFC_sFC_corr.to_csv(outfile)

    #Calculate full sFC and reformat.
    full_collect = []
    full_sFC = np.corrcoef(fullts.T)
    full_sFC = visform(full_sFC,netlabels).iloc[1:,1:].values
    full_collect.append(full_sFC)

    #Save it.
    outfile = (outpath+'full_sFC.csv')
    pd.DataFrame(full_sFC).to_csv(outfile,index=None,header=None)

    #Collect the occurrence-weighted dFC sum.
    dFC_sum = np.zeros((nconn))
    final = []
    for subidx in range(nsub):
        csub = subjects[subidx]

        #Read.
        infile = (outpath+csub+'_dFCsum.h5')
        store = pd.HDFStore(infile,'r')
        outkey = ('sub_'+csub)
        cmat = store.select(outkey)
        store.close()
 
        #Add.
        dFC_sum += cmat.values[:,0]
    
    #Generate average.
    dFC_mean = dFC_sum/(subts*nsub)

    #Unvectorize and reformat.
    dFC_2D = vec2mat(pd.DataFrame(dFC_mean),vecsplit)
    dFC_2D = visform(dFC_2D,netlabels).iloc[1:,1:].values
    full_collect.append(dFC_2D)

    #Save it.
    outfile = (outpath+'full_dFC.csv')
    pd.DataFrame(dFC_2D).to_csv(outfile,index=None,header=None)

    #Add empty matrices to list to keep sizing consistent.
    for eidx in range(nk-2):
        full_collect.append(np.ones((nroi,nroi)))

    #Plot each matrix.
    pidx = 0
    fig, axs = plt.subplots(nrows=1,ncols=nk)
    for ax, dat in zip(axs.ravel(),full_collect):

        #Make the diagonal zero.
        dat[np.diag_indices_from(dat)] = 0

        #Create color bar axes.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        b_ax = divider.append_axes('bottom',size='5%',pad=0.05)
        divnorm=colors.TwoSlopeNorm(vcenter=0)
        
        #Plot.
        im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
        x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
        y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
        cbar = fig.colorbar(im,cax=b_ax,orientation='horizontal')
        ticks = [np.min(dat),0,np.max(dat)]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{tick:.3f}' for tick in ticks])
        cbar.ax.tick_params(labelsize=5,rotation=270,width=0.5)
        cbar.outline.set_linewidth(0.5)
        
        #Remove axes values.
        ax.axis('off')
        x_ax.set_axis_off()
        y_ax.set_axis_off()    
        
        #Increment.
        pidx = pidx + 1
    
    #Save.
    outfile = (outpath+'full_dFC_sFC.jpg')
    plt.savefig(outfile,dpi=720)
    plt.close() 

    #Re-vectorize and save correlations between full sFC and occurrence-weighted dFC average.
    sFC_1D = full_sFC[np.triu_indices(nroi,k=1)]
    dFC_1D = dFC_2D[np.triu_indices(nroi,k=1)]
    full_dFC_sFC_corr = pd.DataFrame((pd.DataFrame((sFC_1D,dFC_1D)).T.corr().values[0,1],
                                   pd.DataFrame((sFC_1D,dFC_1D)).T.corr(method='spearman').values[0,1]),
                                   index=['Pearson','Spearman']).T
    outfile = (outpath+'full_dFC_sFC_corr.csv')
    full_dFC_sFC_corr.to_csv(outfile)

    #Save the vectorized matrices for normality tests.
    outfile = (outpath+'dFC_sFC_full_1D.csv')
    outmat = pd.DataFrame((dFC_1D,sFC_1D))
    outmat.to_csv(outfile)
 