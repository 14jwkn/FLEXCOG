#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the selected k, gather the idiosyncrasies across all subjects across runs
and within run.
Output:
meansim.csv Contains idiosyncrasies for all subjects, across runs.
crun+'_meansim.csv' Contains idiosyncrasies for all subjects, within run.

Usage: 
    LE_group_statesim_gather.py <k> 
    
Arguments:

    <k> K for k-clustering

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
    
    #Generate a matrix for idiosyncrasy.
    meandist = np.zeros((nsub,nk))
    for sidx in range(nsub):
        
        #Extract.
        csub = subjects[sidx]
        
        #Read in the matrix.
        inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                  subgroup+'/'+k+'/'+csub+'/') 
        infile = (inpath+'statesim.h5')
        inkey = ('/sub_'+csub)
        store = pd.HDFStore(infile,'r')
        outdist = store.select(inkey)
        store.close()
        
        #Save.
        meandist[sidx,:] = outdist.loc['Mean',:]
    
    #Format.
    meandist = pd.DataFrame(meandist)
    meandist.index = subjects
    meandist.columns = [('State_'+str(i+1)) for i in range(nk)]
    
    #Save the matrix.
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/statesim/')
    os.makedirs(outpath,exist_ok=True)
    outfile = (outpath+'meansim.csv')
    meandist.to_csv(outfile)

    #For each run.
    runs = ['REST1_LR','REST1_RL','REST2_LR','REST2_RL']
    nrun = len(runs)
    for ridx in range(nrun):
        crun = runs[ridx]

        #Generate a matrix for idiosyncrasy.
        meandist = np.zeros((nsub,nk))
        for sidx in range(nsub):
            
            #Extract.
            csub = subjects[sidx]
            
            #Read in the matrix for idiosyncrasy for the run.
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                     subgroup+'/'+k+'/'+csub+'/') 
            infile = (inpath+'statesim.h5')
            inkey = ('/runwise_'+csub)
            store = pd.HDFStore(infile,'r')
            meandist[sidx,:] = store.select(inkey).loc[crun,:]
            store.close()
        
        #Format.
        meandist = pd.DataFrame(meandist,index=subjects,columns=[('State_'+str(i+1)) for i in range(nk)])
   
        #Save the matrices.
        outfile = (outpath+crun+'_meansim.csv')
        meandist.to_csv(outfile)
    print('Saved.')
    