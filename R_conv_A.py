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

# snps_distribution used to vicualise the distributions of the different snps

def snps_distribution(snps):
    
    for i in range(0,snps.shape[0]):
        a=snps.iloc[i]
        b=a.index
        c=a.values


        c=pd.DataFrame(c)
        c['snps_name']=b

        c.drop(c.index[0],inplace=True)
        c.columns=['val','snps_name']
        plt.scatter(c['snps_name'],c['val'])
        plt.show()












    
