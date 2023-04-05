#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the subject file group, selected k, and subject ID, find the average 
idiosyncrasy of each LE(t) from group LE(t) across runs.
Output:
statesim.h5 Contains all idiosyncrasies for each subject.

Usage: 
    21_r_LE_group_statesim.py <subgroup> <k> <subject>
    
Arguments:
    
    <subgroup> Subject group
    <k> K for k-clustering
    <subject> Subject ID

"""

import os, h5py, sys
import numpy as np
import pandas as pd
from docopt import docopt
from scipy.spatial.distance import cityblock

if __name__ == '__main__':
    __spec__ = None

    #Catches command-line arguments.
    args = docopt(__doc__)
    subgroup = args['<subgroup>']
    k = args['<k>']
    subject = args['<subject>']
    print(subgroup,k,subject)
    
    #Set output file name. If it already exists, quit.
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/'+subject+'/') 
    os.makedirs(outpath,exist_ok=True)
    outfile = (outpath+'statesim.h5')
    
    #Read subjects.
    if (subgroup == 'full'):
        subfile = 'r_full_submain.txt' 
    elif (subgroup == 'half'):
        subfile = 'r_half_submain.txt'      
    with open(subfile) as f:
       subjects = [subject.rstrip() for subject in f]

    #Read file with sorting.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/'+
              subgroup+'/allk/')       
    infile = (inpath+subgroup+'_sortclass.csv')
    reorder_mat = pd.read_csv(infile,header=0,index_col=0)
    reorder_mat.index = [str(x) for x in reorder_mat.index]
    reorder_mat.columns = [str(x) for x in reorder_mat.columns] 

    #Get best iteration.
    bestfile = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/'+
                'group/best_iter/'+subgroup+'/best_iter.csv')
    bestiter = pd.read_csv(bestfile,header=0,index_col=0)
    bestlabel = k
    iteration = str(bestiter.loc[subgroup,bestlabel])
        
    #Read in the states.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/')
    infile = (inpath+'subcent_'+iteration+'.h5')
    inkey = ('/'+subgroup)
    store = h5py.File(infile,'r')
    group_states = np.array(store[inkey]) 
    group_states = group_states.T
        
    #Create shuffler with Python indexing and shuffle order.
    reorder_k = reorder_mat.loc[k,:].dropna().values.astype(int).tolist()
    reorder_k = [(x-1) for x in reorder_k]
    vecmat = group_states[reorder_k,]
    
    #Read in the subject windows.
    Lfile1 = ('../outputs/r_LE_dFC/REST1_LR/'+subject+
                '/LE_dFC.h5')     
    Rfile1 = ('../outputs/r_LE_dFC/REST1_RL/'+subject+
                '/LE_dFC.h5') 
    Lfile2 = ('../outputs/r_LE_dFC/REST2_LR/'+subject+
                '/LE_dFC.h5')
    Rfile2 = ('../outputs/r_LE_dFC/REST2_RL/'+subject+
                '/LE_dFC.h5')
    inkey = ('/LE_LE')
    store = h5py.File(Lfile1,'r')
    dFCmatL1 = np.array(store[inkey]).T
    store.close()
    store = h5py.File(Rfile1,'r')
    dFCmatR1 = np.array(store[inkey]).T
    store.close()
    store = h5py.File(Lfile2,'r')
    dFCmatL2 = np.array(store[inkey]).T
    store.close()
    store = h5py.File(Rfile2,'r')
    dFCmatR2 = np.array(store[inkey]).T
    store.close()
    submat = np.vstack((dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2))
    nts = np.shape(submat)[0]
    
    #Read in the subject clusterings.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/'+k+'/')
    infile = (inpath+'uni_subcent.h5')
    inkey = ('/'+subgroup)
    store = pd.HDFStore(infile,'r')
    subclust = pd.DataFrame(store.select(inkey)).values
    store.close()
    subidx = subjects.index(subject)
    startidx = (subidx)*nts
    endidx = (subidx+1)*nts
    subclust = subclust[startidx:endidx,0]
    
    #Find the distance of each clustered window from their group states.
    windist = np.zeros((nts,1))
    for widx in range(nts):
        
        #Extract the current state.
        cstate = subclust[widx]
        
        #Convert to Python indexing and extract the group state for the current state.
        pythstate = (cstate-1)
        cgroup = vecmat[pythstate,:]
        
        #Extract the current window.
        cwin = submat[widx,:]
        
        #Find the distance between the group state and the window.
        cdist = cityblock(cgroup,cwin)
        windist[widx,0] = cdist
    
    #Find the mean of these distances for each state.
    meandist = np.zeros((1,int(k)))
    for kidx in range(int(k)):
        
        #Convert to nonzero indexing.
        klab = kidx + 1
        
        #Find all values which belong to the state, save the average and std.
        cbelong = windist[np.where(subclust == klab),0]
        meandist[0,kidx] = np.mean(cbelong)
    
    #Package.
    outdist = meandist
    outdist = pd.DataFrame(outdist)
    outdist.index = ['Mean']
    
    #Save for the subject.
    outkey = ('/sub_'+subject)
    store = pd.HDFStore(outfile)
    store.put(outkey,outdist,format='table')
    store.close()
    print('Saved.')
