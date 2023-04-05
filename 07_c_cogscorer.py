#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract 10 cognitive test scores of interest for each subject, then find the
first PC from PCA to use as g.
Output:
cogscores.csv Contains g calculated as PC1.
gPCAvalidity.csv Contains the variance explained from PCA.

Usage: 
    07_c_cogscorer.py <subgroup> 
    
Arguments:

    <subgroup> Subject group

"""


import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from docopt import docopt
    
if __name__ == '__main__':
    __spec__ = None

     #Catches command-line arguments.
    args = docopt(__doc__)
    subgroup = args['<subgroup>']
    
    #Read in subjects.
    if (subgroup == 'full'):
        subfile = 'r_full_submain.txt'  
    elif (subgroup == 'half'):
        subfile = 'r_half_submain.txt'

    #Read in subjects.
    with open(subfile) as f:
        subjects = [subject.rstrip() for subject in f]
    
    #Get cognitive scores.
    openpath = '../inputs/data/hcp/HCP1200_Data_Dictionary.csv'
    cogscores = pd.read_csv(openpath,header=0,usecols=[
                            'Subject',
                            'CardSort_Unadj',
                            'Flanker_Unadj',
                            'ListSort_Unadj',
                            'PicSeq_Unadj',
                            'PicVocab_Unadj',
                            'ProcSpeed_Unadj',
                            'ReadEng_Unadj',
                            'PMAT24_A_CR',
                            'IWRD_TOT',
                            'VSPLOT_TC'
                            ])
    
    #Converts subjects of interest gotten as strings into integers in order to
    #index dataframe which contains integers.
    numsubint = list(map(int,subjects))
    cogint = cogscores[cogscores.Subject.isin(numsubint)][:]
    
    #Set the names of the cognitive tests used.
    cogtests = ['CardSort_Unadj',
                'Flanker_Unadj',
                'ListSort_Unadj',
                'PicSeq_Unadj',
                'PicVocab_Unadj',
                'ProcSpeed_Unadj',
                'ReadEng_Unadj',
                'PMAT24_A_CR',
                'IWRD_TOT',
                'VSPLOT_TC']
    
    #Slice out corresponding columns and convert it to numpy for processing.
    pgtable = cogint.loc[:,cogtests]
    ngtable = pgtable.to_numpy()
    
    #Z-score standardization of scores.
    mgtable = ngtable - ngtable.mean(axis=0,keepdims=True)
    zgtable = mgtable / ngtable.std(axis=0,ddof=1,keepdims=True)
    
    #Conduct PCA to set the first principal component score as g-factor score.
    pca = PCA()
    gscore = pca.fit_transform(zgtable)[:,0]
    
    #Verifies PCA validity, finds factor loadings and explained variance of 
    #each component. Sets the names of the components.
    components = ['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']
    
    #Computes factor loadings and creates data frame from it.
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    loading_mat = pd.DataFrame(loadings,columns=components,index=cogtests)
    
    #Adds row for explained variance of the component.
    loading_mat.loc['Variance Explained'] = pca.explained_variance_ratio_
    
    #Add the score to the list and isolate the scores of interest. 
    #Given that all the loadings for each cognitive test are negative for PC1
    #and a higher score means better performance for each test, we flip the 
    #sign on the PC1 loadings for better interpretability - a higher score 
    #means better performance.
    outcogtable = cogint.assign(gPCA=-gscore)
    outcogtable = outcogtable.loc[:,['Subject','gPCA']]
    
    #Format.
    outcogtable.Subject = outcogtable.Subject.astype(str)
    
    #Set output folder path. If output folder doesn't exist, creates it.
    outpath = ('../outputs/c_cognition/'+subgroup+'/')
    os.makedirs(outpath,exist_ok=True)
    
    #Save.
    outcogtable.to_csv(outpath +'cogscores.csv',index=False)
    loading_mat.to_csv(outpath+'gPCAvalidity.csv',header=True,index=True) 
        