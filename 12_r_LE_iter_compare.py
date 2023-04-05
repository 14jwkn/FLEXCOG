#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the current subject group for the subject file, find the best iteration of
k-medians clustering for each k. Values with lower within-cluster sum of 
distance (WCS) scores are better.
Output:
best_iter.csv For all ks for this subject group, contains the best iteration.  

Usage: 
    r_LE_iter_compare.py <subfile>
    
Arguments:
    
    <subfile> File name of subject list

"""

import os, h5py
import numpy as np
import pandas as pd
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None
    
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    subfile = args['<subfile>']
    # subfile = 'r_full_submain.txt'
    print('Doing:',subfile)
    
    #Set output folder path. If output folder doesn't exist, creates it.
    if subfile == 'r_half_submain.txt':
        subgroup = 'half'
    elif subfile == 'r_full_submain.txt':
        subgroup = 'full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/best_iter/'+subgroup+'/')
    os.makedirs(outpath,exist_ok=True)
   
    #Set output file names. 
    iterfile = outpath + 'best_iter.csv'
    
    #Go through each iteration to make labels.
    row_labs = []
    for k in ['2','3','4','5','6']:    
        row_labs.append(k)
                                
    #Create empty matrix to save best iterations for the group.
    best_mat = pd.DataFrame(np.zeros((1,len(row_labs))))
    best_mat.columns = row_labs
    best_mat.index = [subgroup]
    
    #Go through each iteration and discover the best ones.               
    for k in ['2','3','4','5','6']:  
        niters = 500
        curr_lab = k
        bestwcs = float('inf')
        bestitx = 0
        for itx in range(1,niters+1):
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/'+
                      'ver_MATLAB/'+subgroup+'/'+k+'/')  
            infile = (inpath+'substat_'+str(itx)+'.h5')
            inkey = ('/'+subgroup)
            if not os.path.exists(infile):
                print(curr_lab,str(itx),'does not exist.')
                continue
            statstore = h5py.File(infile,'r')
            statmat = np.array(statstore[inkey]).T
            statstore.close()
            wcs = statmat[0,:].sum()
            if wcs < bestwcs:
                bestwcs = wcs
                bestitx = itx
        best_mat.loc[:,curr_lab] = bestitx

    #Save the table of best iterations.
    best_mat.to_csv(iterfile,header=True,index=True)
