#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the selected k, gather the transition distances across all subjects across 
runs and within run.
Output:
transmean.csv Contains transition distances for all subjects, across runs.
crun+'_transmean.csv' Contains transition distances for all subjects, within run.

Usage: 
    LE_group_statejump_gather.py <k> 
    
Arguments:

    <k> k for k-clustering

"""

import os, h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None

    #Catches command-line arguments.
    args = docopt(__doc__)
    k = args['<k>']
    
    #Set parameters.
    nk = int(k)
    
    #Read subjects.
    subgroup = 'full'
    subfile = 'r_full_submain.txt'    
    with open(subfile) as f:
       subjects = [subject.rstrip() for subject in f]
    nsub = len(subjects)
    
    #Generate a matrix for transition distances.
    transmean = np.zeros((nsub,nk*nk))
    for sidx in range(nsub):
        
        #Extract.
        csub = subjects[sidx]
        
        #Read in.
        inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                  subgroup+'/'+k+'/'+csub+'/') 
        infile = (inpath+'statejump.h5')
        inkey = ('/transmean')
        store = pd.HDFStore(infile,'r')
        transmean[sidx,:] = store.select(inkey)
        store.close()
    
    #Set transition labels.
    translabs = []
    for kidx in range(nk):
        for kidx2 in range(nk):
            translabs.append(str(kidx+1)+'-'+str(kidx2+1))

    #Format.
    transmean = pd.DataFrame(transmean)
    transmean.index = subjects
    transmean.columns = translabs
  
    #Save the matrices.
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/statejump/')
    os.makedirs(outpath,exist_ok=True)
    outfile = (outpath+'transmean.csv')
    transmean.to_csv(outfile)

    #For each run.
    runs = ['REST1_LR','REST1_RL','REST2_LR','REST2_RL']
    nrun = len(runs)
    for ridx in range(nrun):
        crun = runs[ridx]

        #Generate a matrix for transition distances.
        transmean = np.zeros((nsub,nk*nk))
        for sidx in range(nsub):
            
            #Extract.
            csub = subjects[sidx]
     
            #Read in the transition distances for this run.
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                      subgroup+'/'+k+'/'+csub+'/') 
            infile = (inpath+'statejump.h5')
            inkey = ('/run_transmean')
            store = pd.HDFStore(infile,'r')
            transmean[sidx,:] = store.select(inkey).loc[:,crun]
            store.close()
        
        #Package.
        transmean = pd.DataFrame(transmean,index=subjects,columns=translabs)

        #Save.
        outfile = (outpath+crun+'_transmean.csv')
        transmean.to_csv(outfile)
    print('Saved.')
