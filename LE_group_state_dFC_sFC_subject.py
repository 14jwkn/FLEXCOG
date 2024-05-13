#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the selected k and subject ID, calculate the occurrence-weighted sum of dFC matrices for the subject.
Output:
subject+'_dFCsum.h5' Contains the occurrence-weighted sum of dFC matrices.

Usage: 
    LE_group_state_dFC_sFC_subject.py <k> <subject>
    
Arguments:

    <k> K for k-clustering
    <subject> Subject ID

"""

import os, h5py, sys
import numpy as np
import pandas as pd
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None

    #Catches command-line arguments.
    args = docopt(__doc__)
    k = args['<k>']
    subject = args['<subject>']

    #Set output file name. 
    subgroup = 'full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/dFC_sFC/') 
    os.makedirs(outpath,exist_ok=True)

    #Set parameters.
    nk = int(k)
    nts = 1200 - 2
    nroi = 360
    nconn = int((nroi*(nroi-1))/2)

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

    #Calculate occurrence.
    ndFC = len(fullclust)
    occur = np.zeros((nk))
    for kidx in range(nk):
        occur[kidx] = np.sum(fullclust==(kidx+1))/len(fullclust)

    #Calculate occurrence-weighted mean of dFC states.
    dFC_mean = np.zeros((nconn))
    subidx = subjects.index(subject)

    #Set starting and ending subject index.
    sub_startidx = (subidx)*subts
    sub_endidx = sub_startidx + subts

    #Read each run.
    for runidx in range(nruns):
        crun = runs[runidx]

        #Set starting and ending run index.
        run_startidx = sub_startidx + (runidx)*nts
        run_endidx = run_startidx + nts
        run_allidx = [x for x in range(run_startidx,run_endidx)]

        #Use run label to read file.
        inkey = '/LE_dFC'
        infile = ('../outputs/r_LE_dFC/'+crun+'/'+subjects[subidx]+'/LE_dFC.h5')
        store = h5py.File(infile,'r')
        submat = np.array(store[inkey]).T
        store.close()

        #Add each window to the running mean.
        for winidx in range(nts):
            full_winidx = run_allidx[winidx]

            #Extract cluster and corresponding occurrence.
            cclust = fullclust[full_winidx,0]
            coccur = occur[cclust-1]

            #Add weighted dFC window.
            dFC_mean += submat[winidx]*coccur
    
    #Save the subject's sum.
    outfile = (outpath+subject+'_dFCsum.h5')
    store = pd.HDFStore(outfile)
    outkey = ('sub_'+subject)
    store.put(outkey,pd.DataFrame(dFC_mean),format='table')
    store.close()
