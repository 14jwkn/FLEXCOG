#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the selected k and subject ID, find the average 
transition distance between LE(t) for given clustering transitions within each 
run, then average across runs.
Output:
statejump.h5 Contains all transition distances for each subject.

Usage: 
    LE_group_statejump.py <k> <subject>
    
Arguments:

    <k> k for k-clustering
    <subject> Subject ID

"""

import os, h5py, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from docopt import docopt
from scipy.spatial.distance import cityblock

if __name__ == '__main__':
    __spec__ = None

    #Catches command-line arguments.
    args = docopt(__doc__)
    k = args['<k>']
    subject = args['<subject>']
    
    #Set output file name. If it already exists, quit.
    subgroup = 'full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/'+subject+'/') 
    os.makedirs(outpath,exist_ok=True)
    outfile = (outpath+'statejump.h5')
    
    #Set parameters.
    nk = int(k)
    nts = 1200 - 2
    
    #Do each run separately.
    runs = ['REST1_LR','REST1_RL','REST2_LR','REST2_RL']
    nruns = len(runs)
    subts = nts*nruns
    
    #Read subjects.
    subfile = 'r_full_submain.txt'   
    with open(subfile) as f:
       subjects = [subject.rstrip() for subject in f]
      
    #Read in the subject clustering.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/')
    infile = (inpath+'uni_subcent.h5')
    inkey = ('/'+subgroup)
    store = pd.HDFStore(infile,'r')
    fullclust = pd.DataFrame(store.select(inkey)).values
    store.close()
    subidx = subjects.index(subject)
    startidx = (subidx)*subts
    endidx = (subidx+1)*subts
    fullclust = fullclust[startidx:endidx,0]
    
    #Set transition labels.
    translabs = []
    for kidx in range(nk):
        for kidx2 in range(nk):
            translabs.append(str(kidx+1)+'-'+str(kidx2+1))
    
    #Initialize collectors across runs.
    run_transmean = np.zeros((nruns,nk*nk))
        
    #For each run.
    for ridx in range(nruns):
        
        #Extract.
        crun = runs[ridx]
    
        #Read in the subject windows.
        infile = ('../outputs/r_LE_dFC/'+crun+'/'+subject+
                  '/LE_dFC.h5')     
        inkey = '/LE_LE'
        store = h5py.File(infile,'r')
        submat = np.array(store[inkey]).T
        store.close()
        ndim = np.shape(submat)[1]
        
        #Extract clustering for the current run.
        startidx = (ridx)*nts
        endidx = (ridx+1)*nts
        subclust = fullclust[startidx:endidx]
    
        #Initialize.
        alldist = np.zeros((nts,1))
        transcount = pd.DataFrame(np.zeros((1,nk*nk)))
        transcount.columns = translabs
        
        #Find the distance between consecutive time points.
        for cidx in range(nts - 1):
            
            #Extract next time point.
            nidx = cidx + 1
            
            #Extract the FC at current and next time point.
            cFC = submat[cidx,:]
            nFC = submat[nidx,:]
            
            #Find the distance between and save.
            cdist = cityblock(nFC,cFC)
            alldist[nidx,0] = cdist
            
            #Extract the cluster at the current and next time point.
            cclust = subclust[cidx]
            nclust = subclust[nidx]
            ctrans = str(cclust)+'-'+str(nclust)
            
            #Add to the count of the specific transition.
            transcount.loc[0,ctrans] = transcount.loc[0,ctrans] + 1
        
        #Initialize transition distance collection based on maximum count.
        maxcount = int(np.max(np.max(transcount)))
        transall = np.zeros((maxcount,nk*nk))
        transall[:] = np.nan
        transall = pd.DataFrame(transall)
        transall.columns = translabs
        
        #Find the current ID to add in maxcount.
        ccounter = pd.DataFrame(np.zeros((1,nk*nk)))
        ccounter.columns = translabs
        
        #Find the distance between consecutive time points.
        for cidx in range(nts - 1):
            
            #Extract next time point.
            nidx = cidx + 1
            
            #Extract the FC at current and next time point.
            cFC = submat[cidx,:]
            nFC = submat[nidx,:]
            
            #Find the distance between and save.
            cdist = cityblock(nFC,cFC)
            
            #Extract the cluster at the current and next time point.
            cclust = subclust[cidx]
            nclust = subclust[nidx]
            ctrans = str(cclust)+'-'+str(nclust)
            
            #Extract current ID.
            cID = ccounter.loc[0,ctrans]
             
            #Add to the distance of the specific transition.
            transall.loc[cID,ctrans] = cdist
            
            #Update the current ID.
            ccounter.loc[0,ctrans] = cID + 1
        
        #Find the mean and append.
        transmean = transall.mean()
        run_transmean[ridx,:] = transmean
 
    #Average matrices across runs.
    avg_transmean = np.nanmean(run_transmean,0)
    
    #Package.
    avg_transmean = pd.DataFrame(avg_transmean).T
    avg_transmean.columns = translabs
    
    #Save.
    outkey = ('/transmean')
    store = pd.HDFStore(outfile)
    store.put(outkey,avg_transmean,format='table')
    store.close()
    print('Saved.')     
