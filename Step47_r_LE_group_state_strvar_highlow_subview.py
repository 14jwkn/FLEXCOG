#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the subject file group and selected k, plot the dFC states for the top
and bottom 10% of average strength and variability scores.
Output:
clim+'_'+statetype+'_'+strvar+'.jpg' Matrix plots for the dFC states for high and low strength/variability.
'global_colleg_'+statetype+'_'+strvar+'.jpg' Color bar for the matrix plots.

Usage: 
    Step47_r_LE_group_state_strvar_highlow_subview.py <subgroup> <k> <strvar>
    
Arguments:

    <subgroup> Subject group
    <k> k for k-means
    <strvar> Strength (mean) or variability (std)
    
"""

import h5py, os
import numpy as np
import pandas as pd
from docopt import docopt
import seaborn as sns
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Sort matrix by labels.
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
    
    #Read in parameters.
    args = docopt(__doc__)
    subgroup = args['<subgroup>']
    k = args['<k>']
    strvar = args['<strvar>']
    print(subgroup,k,strvar)

    #Set parameters.
    nk = int(k)
    limval = 0.1
    statetype = 'dFC'
    
    #Set outpath.
    basepath =  ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/allreconfig_subsets/strvar/')
    
    #Read in network labels and convert them to integers.
    with open('colenetlabels.txt') as f:
        netlabels = [label.rstrip() for label in f] 
    netlabels = list(map(int,netlabels))
    
    #Read in vectorized labels.
    veclabels = pd.read_csv('atlasMvec2.csv',header=None).values.tolist()
    vecsplit = []
    for [vecstr] in veclabels:
        vecsplit.append(vecstr.split('-'))
    
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
    reorder_mat.loc[k,:].dropna().values
    reorder_k = reorder_mat.loc[k,:].dropna().values.astype(int).tolist()
    
    #Generate large matrices for each state.
    high_list = []
    low_list = []
    for kidx in range(nk):
        ck = str(kidx+1)
        print(ck)

        #Read in all subjects of the high and low groups.
        targtype = (statetype+'_'+strvar+'_'+ck)
        infile = (basepath+'high_'+targtype+'.txt')
        with open(infile) as f:
            highsubs = [subject.rstrip() for subject in f]
        infile = (basepath+'low_'+targtype+'.txt')
        with open(infile) as f:
            lowsubs = [subject.rstrip() for subject in f]
        
        #Limit.
        nsubs = round(limval*len(highsubs))
        highsubs = highsubs[0:nsubs]
        lowsubs = lowsubs[0:nsubs]
        
        #Go through high and low.
        for clim in ['high','low']:
            if clim == 'high':
                sublist = highsubs
            elif clim == 'low':
                sublist = lowsubs
            print(clim)
 
            #Produce matrix.
            highmat = np.zeros((nsubs,64620))
            lowmat = np.zeros((nsubs,64620))
                
            #Go through subjects.
            for sidx in range(nsubs):
                
                #Read in states.    
                sub = sublist[sidx]
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                          subgroup+'/'+k+'/'+sub+'/')
                inkey = ('/'+subgroup)
                infile = (inpath+'med_'+statetype+'.h5')
                instore = h5py.File(infile,'r')
                statesall = np.array(instore[inkey]).T
                instore.close()
                
                #Sort.
                statesall = np.array([statesall[(i-1),:] for i in reorder_k])

                #Append.
                if clim == 'high':
                    highmat[sidx,:] = statesall[kidx,:]
                elif clim == 'low':
                    lowmat[sidx,:] = statesall[kidx,:]
            
            #Append.
            if clim == 'high':
                high_list.append(highmat)
            elif clim == 'low':    
                low_list.append(lowmat)
    
    #Find centroids.
    highcent = np.zeros((nk,64620))
    lowcent = np.zeros((nk,64620))
    for kidx in range(nk):
        
        #Extract.
        highmat = high_list[kidx]
        lowmat = low_list[kidx]
        
        #Centroid.
        highcent[kidx,:] = np.median(highmat,axis=0)
        lowcent[kidx,:] = np.median(lowmat,axis=0)
    
    #Save centroids.
    outpath = basepath
    outfile = (outpath+'high_'+statetype+'_'+strvar+'_cent.csv')
    pd.DataFrame(highcent).to_csv(outfile,header=False,index=False)
    outfile = (outpath+'low_'+statetype+'_'+strvar+'_cent.csv')
    pd.DataFrame(lowcent).to_csv(outfile,header=False,index=False)
    print('Saved centroids.')
    
    #Go through matrices for limits.
    limvals = []
    for clim in ['high','low']:
        if clim == 'high':
            ccent = highcent
        elif clim == 'low':
            ccent = lowcent
        limvals.append(np.max(ccent))
        limvals.append(np.min(ccent))
               
    #Find limits.
    maxval = max(limvals)
    minval = min(limvals)
    
    #Define outpath.
    outpath = basepath+'plots/'
    os.makedirs(outpath,exist_ok=True)
    
    #Go through top and bottom to plot.
    for clim in ['high','low']:
        if clim == 'high':
            ccents = highcent
        elif clim == 'low':
            ccents = lowcent
            
        #Convert to 2D.
        collect2D = []
        for kidx in range(nk):
            cstate = pd.DataFrame(ccents[kidx,:])
            cstate = vec2mat(cstate,vecsplit)
            reformatted = visvec(cstate,netlabels)
            collect2D.append(reformatted.iloc[1:,1:].values)
            
        #Plot the states with global weights.
        pidx = 0
        fig, axs = plt.subplots(nrows=1,ncols=nk)
        for ax, dat in zip(axs.ravel(),collect2D):
            
            #Define min and max according if concerned.
            cmax = maxval
            cmin = minval
    
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
        outfile = outpath+clim+'_'+statetype+'_'+strvar+'.jpg'
        plt.savefig(outfile,dpi=720)
        plt.close()   
    
    #Save global color bars.
    divnorm = colors.TwoSlopeNorm(vmin=minval,vcenter=0,vmax=maxval)
    c = plt.pcolormesh(collect2D[0],cmap='RdBu_r',norm=divnorm)
    cbar = plt.colorbar(c)
    cbar.ax.tick_params(labelsize=20) 
    outfile = outpath+'global_colleg_'+statetype+'_'+strvar+'.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close() 
    print('Saved plots.')
        