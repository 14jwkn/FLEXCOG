#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the subject file group and selected k, do PCA on the LE(t) across subjects
and retrieve the PCs.
Output:
LEfull_PCreduce.h5 Contains the PCs for the LE(t) across subjects and labels.
LEfull_PCreduce_explained.txt Contains the variance explained for the PCs.

Usage: 
    39_r_LE_group_state_dimtaker_LEfull.py <subgroup> <k> 
    
Arguments:

    <subgroup> Subject group
    <k> k for k-means
    
"""

import h5py, os
import numpy as np
import pandas as pd
from docopt import docopt
from sklearn.decomposition import PCA

if __name__ == '__main__':
    __spec__ = None
    
    #Read in parameters.
    args = docopt(__doc__)
    subgroup = args['<subgroup>']
    k = args['<k>']
    print(subgroup,k)
    
    #Set outpath.
    basepath =  ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/lowdim_LEfull/')
    os.makedirs(basepath,exist_ok=True)
    
    #Set parameters.
    nwin = 4792
    nreg = 360
    nconn = int((nreg*(nreg-1))/2)
    nk = int(k)
    
    #Read in subjects.
    if subgroup == 'full':
        subfile = 'r_full_submain.txt'  
    elif subgroup == 'half':
        subfile = 'r_half_submain.txt'   
    with open(subfile) as f:
       subjects = [subject.rstrip() for subject in f]
    nsubs = len(subjects) 

    #Set dimensions.
    ndim = 3
    ndat = nwin
    nrows = ndat*nsubs

    #Set up an empty container. 
    dFCmat = np.empty((nk+nrows,nreg),dtype='float32')
     
    #Get best iteration.
    bestfile = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/'+
              'group/best_iter/'+subgroup+'/best_iter.csv')
    bestiter = pd.read_csv(bestfile,header=0,index_col=0)
    bestlabel = k
    iteration = str(bestiter.loc[subgroup,bestlabel])
    
    #Read in group centroids.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/')
    infile = (inpath+'subcent_'+iteration+'.h5')
    inkey = ('/'+subgroup)
    store = h5py.File(infile,'r')
    group_states = np.array(store[inkey]) 
    group_states = group_states.T
    store.close()
    
    #Place in container.
    dFCmat[0:nk,:] = group_states
    
    #For each subject.
    for subidx in range(nsubs):
    
        #Extract.
        csub = subjects[subidx]
        print(csub)
    
        #Set data indices.
        datstart = (ndat*(subidx)) + nk
        datend = (ndat*(subidx+1)) + nk
    
        #Read in data and append.
        inkey = '/LE_LE'
        Lfile1 = ('../outputs/r_LE_dFC/REST1_LR/'+csub+'/LE_dFC.h5')     
        Rfile1 = ('../outputs/r_LE_dFC/REST1_RL/'+csub+'/LE_dFC.h5') 
        Lfile2 = ('../outputs/r_LE_dFC/REST2_LR/'+csub+'/LE_dFC.h5') 
        Rfile2 = ('../outputs/r_LE_dFC/REST2_RL/'+csub+'/LE_dFC.h5') 
        instore = h5py.File(Lfile1,'r')
        dFCmatL1 = np.array(instore[inkey]).T
        instore.close()   
        instore = h5py.File(Rfile1,'r')
        dFCmatR1 = np.array(instore[inkey]).T
        instore.close()  
        instore = h5py.File(Lfile2,'r')
        dFCmatL2 = np.array(instore[inkey]).T
        instore.close()   
        instore = h5py.File(Rfile2,'r')
        dFCmatR2 = np.array(instore[inkey]).T
        instore.close() 
        dFCmat[datstart:datend,:] = np.vstack((dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2))
        del dFCmatL1, dFCmatR1, dFCmatL2, dFCmatR2

    #Do dimension reduction.
    print('Reducing dimensions.')
    pcaobj = PCA(n_components=ndim)
    pcaobj.fit(dFCmat)
    dFClow = pd.DataFrame(pcaobj.transform(dFCmat))

    #Generate labels.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/')
    infile = (inpath+'uni_subcent.h5')
    inkey = ('/'+subgroup)
    store = pd.HDFStore(infile,'r')
    fullclust = pd.DataFrame(store.select(inkey)).values
    store.close()
    clulab = np.zeros((nk+nrows,2))
    for subidx in range(nsubs):
    
        #Extract.
        csub = int(subjects[subidx])
        print(csub)
    
        #Set data indices.
        datstart = (ndat*(subidx)) + nk
        datend = (ndat*(subidx+1)) + nk
        
        #Put in subject labels.
        clulab[datstart:datend,0] = [csub]*ndat
        
        #Put in cluster labels.
        cluidx = subjects.index(str(csub))
        clustart = (nwin*(cluidx))
        cluend = nwin*((cluidx+1))
        clulab[datstart:datend,1] = fullclust[clustart:cluend,0]
    
    #Format.
    clulab = pd.DataFrame(clulab)
    
    #Save the matrices.
    outfile = (basepath+'LEfull_PCreduce.h5')
    outstore = pd.HDFStore(outfile)
    outkey = ('/dFClow')
    outstore.put(outkey,dFClow,format='table')
    outstore.close()
    outstore = pd.HDFStore(outfile)
    outkey = ('/dFClow_lab')
    outstore.put(outkey,clulab,format='table')
    outstore.close()
    outfile = (basepath+'LEfull_PCreduce_explained.txt')
    pd.DataFrame(pcaobj.explained_variance_ratio_).to_csv(outfile,index=False,header=False)
    print('Saved.')
    