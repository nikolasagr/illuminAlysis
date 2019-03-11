#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:30:07 2019

@author: nicolasagrotis
"""

# Quality Control Function architecture:
            
# 1) remove unreliable functions 'remove_unreliable_samples'
            #- missing
            #-outliers
            #- ! Remove corresponding samples in other list entries

#2) infer sex 'infer_sex':
            # - get the F&M column
            # - maybe try and also disply the results
            # - ask for threshold
            
#3) call SNPS 'call_snps' :
            # - identify the snps using the thresholds
            # - plot? ASK tim
            
#4) snps distance matrix 'snp_distance':
            # - plot the heatmap
            # - display the matrix
            
#5) identify duplicates 'identify_replicates'
            # - print duplicates/replicates

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyreadr
from datetime import datetime
from math import sqrt
import timeit


from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix,classification_report
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import confusion_matrix,classification_report



#-----------------------------------------------------------------------------------------#
# Data Loading

control_beads = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/hm450_controls.Rds')
control_beads = control_beads[None]

# In functions its dnam
cpgs = pd.read_csv("/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/cpgs.csv")

controls = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_controls.Rds')
controls = controls[None]

#snps = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/eira_snps.Rds')
#snps = snps[None]

snps = pd.read_csv("/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/snps.csv")

samples = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_samples.Rds')
samples = samples[None]

# Load individual characteristics data

covars = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/Covariates.Rds')
covars = covars[None]


#-----------------------------------------------------------------------------------------#

# samples.shape
# Out[669]: (560, 28)

def remove_unreliable_samples(samples,threshold):
        
    # Apply the threshold for the missing data per row
    samples=samples.loc[samples['missing']<threshold]
    
    # Apply the outlier detection based on the values of  bc1.grn,bc1.red,bc2
    tech_vars=samples[['bc1.grn','bc1.red','bc2']]
    #samples.set_index('sample.id',inplace=True)
    
    
    mean=[]
    sd=[]
    m_1=[]
    m_2=[]
    
    # set the outlier threshold to be between mean + 2sd and mean -2sd
    for m in range(0,tech_vars.shape[1]):
        
        mean_col= tech_vars.iloc[:,m].mean()
        mean.append(mean_col)
    
        sd_col= tech_vars.iloc[:,m].std()
        sd.append(sd_col)
        
        m_1.append(mean[m] - (2 * sd[m]))
        m_2.append(mean[m] + (2 * sd[m]))    
        
          
    for i in range(0,len(mean)):
    
        tech_vars = tech_vars[tech_vars.iloc[:,i]>m_1[i]]
        tech_vars = tech_vars[tech_vars.iloc[:,i]<m_2[i]]
        
        # Print the indices common to the output of the missing data and the ouput of the tech_vars
        common = samples.index.intersection(tech_vars.index)

        samples=samples.loc[common]
    
    return samples

#-----------------------------------------------------------------------------------------#    
    










