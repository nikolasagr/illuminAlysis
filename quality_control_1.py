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
            #- Remove corresponding samples in other list entries

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
from math import sqrt
import timeit

#-----------------------------------------------------------------------------------------#
# Data Loading

control_beads = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/hm450_controls.Rds')
control_beads = control_beads[None]

# In functions its dnam
cpgs = pd.read_csv("/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/cpgs.csv",index_col='Unnamed: 0')

controls = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_controls.Rds')
controls = controls[None]

#snps = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/eira_snps.Rds')
#snps = snps[None]

snps = pd.read_csv("/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/snps.csv",index_col='Unnamed: 0')

samples = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_samples.Rds')
samples = samples[None]

# Load individual characteristics data

covars = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/Covariates.Rds')
covars = covars[None]


#-----------------------------------------------------------------------------------------#
# 1) remove unreliable functions 'remove_unreliable_samples'
            #- missing
            #-outliers
            
# samples.shape
# Out[669]: (560, 28)
# SMV (573,27)
# Simple(560,27)
# PCA (501)

def remove_unreliable_samples(samples,threshold):
        
    # Apply the threshold for the missing data per row
    samples=samples.loc[samples['missing']<threshold]
    cpgs=cpgs.loc[samples.index.intersection(cpgs.index)]
    
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

#s_out=samples_outliers.index
#s_simp=samples_simple.index

#s_simp.intersection(s_out)

#-----------------------------------------------------------------------------------------#    

def k_mean_sex_infer(samples):

    from sklearn import datasets
    from sklearn.cluster import KMeans
    import sklearn.metrics as sm
    
    x=samples[['median.chrX','missing.chrY']]
    
    # K Means Cluster
    model = KMeans(n_clusters=2)
    model.fit(x)
    
    # View the results
    # Set the size of the plot
    plt.figure(figsize=(14,7))
     
    # Create a colormap
    colormap = np.array(['red', 'lime'])
    
    # Plot the Original Classifications
    plt.subplot(1, 2, 1)
    plt.scatter(x['median.chrX'], x['missing.chrY'], s=40)
    plt.title('Real Classification')
     
    # Plot the Models Classifications
    plt.subplot(1, 2, 2)
    plt.scatter(x['median.chrX'], x['missing.chrY'], c=colormap[model.labels_], s=40)
    plt.title('K Mean Classification')

    
    samples['sex_Kmeans']=model.labels_
    
    samples.loc[(samples['sex_Kmeans']==1),'sex_Kmeans']='F'
    samples.loc[(samples['sex_Kmeans']==0),'sex_Kmeans']='M'



#-----------------------------------------------------------------------------------------#    
#2) infer sex 'infer_sex':
            # - get the F&M column
            # - ask for threshold
plt.scatter(x=samples['median.chrX'],y=samples['missing.chrY'])

def infer_sex(samples,threshold_chrX=0.37,threshold_chrY=0.39):
    
    #x=float(median_X)
    #y=float(missing_Y)
    
    # Subset the median.chrX and missing.chrY from the samples dataset
    
    # From bibliography set hard boundaries to descriminate between males and females
    # hard boundaries: if median.chrX' < 0.37 and missing Y chromosome is smaller than 0.39 then M
                      #if median.chrX' > 0.37 and missing Y chromosome is bigger than 0.39 then F
                      # Otherwise set the value to NaN
                      
    samples.loc[(samples['median.chrX'] < threshold_chrX) & (samples['missing.chrY'] < threshold_chrY), 'sex'] = 'M'
    samples.loc[(samples['median.chrX'] > threshold_chrX) & (samples['missing.chrY'] > threshold_chrY), 'sex'] = 'F'
    samples.loc[(samples['median.chrX'] < threshold_chrX) & (samples['missing.chrY'] > threshold_chrY), 'sex'] = np.nan
    samples.loc[(samples['median.chrX'] > threshold_chrX) & (samples['missing.chrY'] < threshold_chrY), 'sex'] = np.nan    
    
    #Count the number of males and females
    num_males=samples.loc[samples.sex == 'M', 'sex'].count()
    num_females=samples.loc[samples.sex == 'F', 'sex'].count()
    print("Number of Males:",num_males)
    print("Number of Females:",num_females)
    
    samples.set_index("sample.id",inplace=True)
    
    return samples

#-----------------------------------------------------------------------------------------#    

# snps distribution plot
# The allele at a SNP locus can be inferred from SNP intensities measured on the BeadChip.
            
def snps_distribution(snps):
    
    for i in range(0,snps.shape[0]):
        a=snps.iloc[i]
        b=a.index
        snp_vals=a.values


        snp_vals=pd.DataFrame(snp_vals)
        snp_vals['snps_name']=b

        snp_vals.drop(snp_vals.index[0],inplace=True)
        snp_vals.columns=['val','snps_name']
        plt.scatter(snp_vals['snps_name'],snp_vals['val'])
        #plt.show()
        
#-----------------------------------------------------------------------------------------#    

#3) call SNPS 'call_snps' :
            # - identify the snps using the thresholds
            # - plot? ASK tim

def identify_replicates(snps,threshold,samples):
    
    snps.set_index('Unnamed: 0',inplace=True)
    snps.index.name='sample_id'
    
    snps[snps<=0.2]=0
    snps[snps>=0.8]=2
    snps[(snps>0.2 ) & (snps<0.8)]=1
    
    dist_matrix = np.empty((snps.shape[0], snps.shape[0]))
    dist_matrix[:,:] = np.nan

    for i in range(0,snps.shape[0]):
        for j in range(i+1,snps.shape[0]):
            dist_matrix[j, i] = abs(snps.iloc[i,:]-snps.iloc[j,:]).sum()
            dist_m=pd.DataFrame(dist_matrix)
            dist_m.index=snps.index
            dist_m.columns=snps.index
            
    ax = sns.heatmap(dist_m, annot=True)
    ax.set_title('Distance Matrix for SNPSs')
    
    dist_m=dist_m < threshold
    
    rows=[]
    columns=[]
    
    for i in range(0,dist_m.shape[0]):
        for j in range(0,dist_m.shape[0]):
            
            if dist_m.iloc[i,j] == True:
                
                rows.append(dist_m.index[i])
                columns.append(dist_m.index[j])
                
    
    sex_ident=samples['sex']
    
    for n in range (0,len(rows)):
            
        if sex_ident[rows[n]] == sex_ident[columns[n]]:
                
            print('Replicate detected:',rows[n],columns[n]) 


#-----------------------------------------------------------------------------------------#

