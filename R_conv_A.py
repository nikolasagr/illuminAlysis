#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 17:43:55 2019

@author: nicolasagrotis
"""
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


#samplesheet = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illumiAlysis/illumiData/Sample_sheet.Rds')
#samplesheet = samplesheet[None]


#-----------------------------------------------------------------------------------------#
# Remove samples with >10% NAs

def missing_less_10(samples):

    samples=samples.set_index("sample.id")

    # Before 689, after 673
    samples.iloc[:,:].isnull().sum()
    samples= samples.loc[samples.isnull().mean(axis=1) < 0.1,]
    
    return samples
    
#-----------------------------------------------------------------------------------------#

# Remove samples with >10% NAs
# Correct
# Produced the same result as R
def missing_thresh(samples):
    
    samples=samples.loc[samples['missing']<0.1]
    
    return samples


#-----------------------------------------------------------------------------------------#
#OneClassSVM is an algorithm that specializes in learning the expected distributions in a dataset. OneClassSVM is especially useful as a novelty detector method if you can first provide data cleaned from outliers; otherwise, it’s effective as a detector of multivariate outliers. In order to have OneClassSVM work properly, you have two key parameters to fix:

# fit the model
def SMV_outliers(samples):
    
    from sklearn import preprocessing
    from sklearn import svm
    
    tech_vars=samples[['bc1.grn','bc1.red','bc2']]
    
    # Standardisation prior to outlier detection
    
    tech_vars_stand = preprocessing.scale(tech_vars)
    tech_vars_stand = pd.DataFrame(tech_vars_stand)
    tech_vars_stand.index=tech_vars.index
    tech_vars_stand.columns=tech_vars.columns
    
    # OneClassSVM oulier detector
    # nu, #nu, which can be calculated by the following formula:
    # nu_estimate = 0.95 * f + 0.05, where f is the percentage of expected outliers 
    # (a number from 1 to 0). If your purpose is novelty detection, f will be 0.

    # gamma, telling the algorithm whether to follow or approximate the 
    # dataset distributions. For novelty detection, it is better to have a value of 0 or superior 
    # (follow the distribution); for outlier detection values, smaller than 0 values are preferred 
    # (approximate the distribution).

    clf = svm.OneClassSVM(nu=0.05, kernel="rbf", gamma=0.01)
    clf.fit(tech_vars_stand)

    y_pred_test = clf.predict(tech_vars_stand)
    
    # Number of outliers stored in n_error_test
    n_error_test = y_pred_test[y_pred_test == -1].size
    print("Number of outliers removed:",n_error_test)

    y_pred_test_df=pd.DataFrame(y_pred_test)

    y_pred_test_df.index=tech_vars_stand.index

    # Subset the samples that are = 1 which are the non outliers
    
    y_pred_test_df=y_pred_test_df[y_pred_test_df==1]

    
    y_pred_test_df= y_pred_test_df.dropna()

    # Apply the subsetted indices to the samples dataframe
    common = samples.index.intersection(y_pred_test_df.index)

    samples=samples.loc[common]
    
    return samples

#-----------------------------------------------------------------------------------------#

# Outlier Method using SD
# Outliers are identified based by the criterion:
# Either they are higher than Mean + 2SD
# Or they are lower than Mean - 2SD
    
def simple_outlier(samples):

    tech_vars=samples[['bc1.grn','bc1.red','bc2']]
    samples.set_index('sample.id',inplace=True)
    
    
    mean=[]
    sd=[]
    m_1=[]
    m_2=[]
    
    for m in range(0,tech_vars.shape[1]):
        
        mean_col= tech_vars.iloc[:,m].mean()
        mean.append(mean_col)
    
        sd_col= tech_vars.iloc[:,s].std()
        sd.append(sd_col)
        
        m_1.append(mean[m] - (2 * sd[m]))
        m_2.append(mean[m] + (2 * sd[m]))    
        
          
    for i in range(0,len(mean)):
    
        tech_vars = tech_vars[tech_vars.iloc[:,i]>m_1[i]]
        tech_vars = tech_vars[tech_vars.iloc[:,i]<m_2[i]]
        
    return tech_vars

#571 Samples left after outlier removal
    

#-----------------------------------------------------------------------------------------#
# Infer the sex from the median X chromosome value and the missing Y chromosome values

def infer_sex(samples):
    
    # Subset the median.chrX and missing.chrY from the samples dataset
    
    sex_info=samples[['median.chrX','missing.chrY']]
    
    # From bibliography set hard boundaries to descriminate between males and females
    # hard boundaries: median.chrX' < 0.37
    conditions = [(sex_info['median.chrX'] < 0.37) & (sex_info['missing.chrY'] < 0.39),
              (sex_info['median.chrX'] > 0.37) & (sex_info['missing.chrY'] > 0.39)]
    choices = ['M', 'F']
    sex_info['sex'] = np.select(conditions, choices)
    
    # Set the correct indexing 
    samples=samples.set_index("sample.id")
    sex_info.index=samples.index
    
    # Count the number of males and females
    num_males=sex_info.loc[sex_info.sex == 'M', 'sex'].count()
    num_females=sex_info.loc[sex_info.sex == 'F', 'sex'].count()
    print("Number of Males:",num_males)
    print("Number of Females:",num_females)
    
    return sex_info

#-----------------------------------------------------------------------------------------#
# Infer the sex from the median X chromosome value and the missing Y chromosome values

def infer_sex_2(samples,median_X,missing_Y):
    
    x=float(median_X)
    y=float(missing_Y)
    
    # Subset the median.chrX and missing.chrY from the samples dataset
    
    # From bibliography set hard boundaries to descriminate between males and females
    # hard boundaries: median.chrX' < 0.37
    samples.loc[(samples['median.chrX'] < 0.37) & (samples['missing.chrY'] < 0.39), 'sex'] = 'M'
    samples.loc[(samples['median.chrX'] > 0.37) & (samples['missing.chrY'] > 0.39), 'sex'] = 'F'
    
    
    # Set the correct indexing
    #samples_A.index=samples.index
    
    #Count the number of males and females
    num_males=samples.loc[samples.sex == 'M', 'sex'].count()
    num_females=samples.loc[samples.sex == 'F', 'sex'].count()
    print("Number of Males:",num_males)
    print("Number of Females:",num_females)
    
    
    return samples

samples.set_index("sample.id",inplace=True)
    
#-----------------------------------------------------------------------------------------#
# Confirm that the sex inferal done earlier is correct by checking the covars
    
def sex_comp(covars,samples):
    
    f_num_covars=covars.loc[covars['gender']== 'f','gender'].count()
    m_num_covars=covars.loc[covars['gender']== 'm','gender'].count()
    
    sex_info=samples[['median.chrX','missing.chrY']]
    
    # From bibliography set hard boundaries to descriminate between males and females
    # hard boundaries: median.chrX' < 0.37
    sex_info.loc[(sex_info['median.chrX'] < 0.37) & (sex_info['missing.chrY'] < 0.39), 'sex'] = 'M'
    sex_info.loc[(sex_info['median.chrX'] > 0.37) & (sex_info['missing.chrY'] > 0.39), 'sex'] = 'F'
    samples=samples.set_index("sample.id")
    
    # Set the correct indexing
    sex_info.index=samples.index
    
    #Count the number of males and females
    num_males_samp=sex_info.loc[sex_info.sex == 'M', 'sex'].count()
    num_females_samp=sex_info.loc[sex_info.sex == 'F', 'sex'].count()
    
    if num_males_samp == m_num_covars:
        print("Male numbers of the infered function are the same")
    else:
        print("Male numbers of the infered function are not the same")
    
    if num_females_samp == f_num_covars:
        print("Female numbers of the infered function are the same")
    else:print("Male numbers of the infered function are not the same")
    
    return
        
#-----------------------------------------------------------------------------------------#   
# Identify replicates
    
snps.shape
snps.iloc[0]
snps.iloc[0].isnull().sum()
snps.isnull().sum().sum()
snps.dropna(inplace=True)

snps.shape[0]
#-----------------------------------------------------------------------------------------#   
# snps_distribution used to vicualise the distributions of the different snps

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
        plt.show()


#-----------------------------------------------------------------------------------------#   
# Categorise the snps into bands depending on their values

snps.set_index('Unnamed: 0',inplace=True)
snps.index.name='sample_id'    
    
snps_decomposition = {}
for n in range(0, snps.shape[0]):
    
    snps_decomposition['sample_%s' % n] = []
    
    
    for i in snps.iloc[n,:]:
        
        if i <= 0.2:
            
            i = 0
            
            #lista.append(i)
        
        if i>= 0.8:
            
            i = 2
            
            #lista.append(i)
            
        elif i>0.2 or i<0.8: 
            
            i = 1
            
            #lista.append(i)
            
        snps_decomposition['sample_' + str(n)].append(i)
    
snps_Dec=pd.DataFrame.from_dict(snps_decomposition)
snps_Dec.index=snps.columns
snps_Dec.columns=snps.index

#-----------------------------------------------------------------------------------------#   

def snps_DeComp(snps):
    
    snps.set_index('Unnamed: 0',inplace=True)
    snps.index.name='sample_id'
    
    snps[snps<=0.2]=0
    snps[snps>=0.8]=2
    snps[(snps>0.2 ) & (snps<0.8)]=1
    
    return snps

#-----------------------------------------------------------------------------------------#   
# if a sample has the same values in its snps
# && the same sex
# Then the samples are duplicates
# create a new Dataframe that includes all the duplicates

from scipy.spatial import distance_matrix

dist_snps=pd.DataFrame(distance_matrix(snps.values, snps.values), index=snps.index, columns=snps.index)

#-----------------------------------------------------------------------------------------#   

# Distance matrix to see similarities in SNPS
def similarity_snps_m(snps_call):

    dist_matrix = np.empty((snps_call.shape[0], snps_call.shape[0]))
    dist_matrix[:,:] = np.nan
    
    for i in range(0,snps_call.shape[0]):
        for j in range(i+1,snps_call.shape[0]):
            dist_matrix[j, i] = abs(snps_call.iloc[i,:]-snps_call.iloc[j,:]).sum()
            dist_m=pd.DataFrame(dist_matrix)
            dist_m.index=snps_call.index
            dist_m.columns=snps_call.index
            
            #ax = sns.heatmap(dist_m, annot=True, fmt="d")
            
    return dist_m

# Heatmap visualisation
    
ax = sns.heatmap(dist_m, annot=True)
ax.set_title('Distance Matrix for SNPSs')

#-----------------------------------------------------------------------------------------#   

# Replicates

def replicates_pullout(dist_m, threshold,sex_info):
    
    dist_m=dist_m < threshold
    
    rows=[]
    columns=[]
    
    for i in range(0,dist_m.shape[0]):
        for j in range(0,dist_m.shape[0]):
            
            if dist_m.iloc[i,j] == True:
                
                rows.append(dist_m.index[i])
                columns.append(dist_m.index[j])
                
    
    sex_ident=sex_info['sex']
    
    for n in range (0,len(rows)):
            
        if sex_ident[rows[n]] == sex_ident[columns[n]]:
                
            print('Replicate detected:',rows[n],columns[n])
            
#-----------------------------------------------------------------------------------------#   
 
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
            

            
        