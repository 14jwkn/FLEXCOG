#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the subject file group and selected k, gather the idiosyncrasies across
all subjects.
Output:
meansim.csv Contains idiosyncrasies for all subjects.

Usage: 
    22_r_LE_group_statesim_gather.py <subgroup> <k> 
    
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
    
    #Generate a matrix for mean and SD.
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
    print('Saved.')
    