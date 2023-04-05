#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the subject file group and selected k, gather the transition distances across
all subjects.
Output:
transmean.csv Contains transition distances for all subjects.

Usage: 
    r_LE_group_statejump_gather.py <subgroup> <k> 
    
Arguments:

    <subgroup> Subject group
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
    subgroup = args['<subgroup>']
    k = args['<k>']
    
    #Set parameters.
    nk = int(k)
    
    #Read subjects.
    if (subgroup == 'full'):
        subfile = 'r_full_submain.txt'  
    elif (subgroup == 'half'):
        subfile = 'r_half_submain.txt'   
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
    print('Saved.')
