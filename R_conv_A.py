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


import matplotlib.pyplot as plt
import seaborn as sns  # plot

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression

from sklearn.linear_model import LinearRegression
from datetime import datetime

from sklearn.metrics import mean_squared_error
from math import sqrt
from sklearn.metrics import confusion_matrix,classification_report
from sklearn.model_selection import train_test_split


#-----------------------------------------------------------------------------------------#
import pyreadr



control_beads = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/hm450_controls.Rds')
control_beads = control_beads[None]

# Load DNA Data
# too big to run
dnam = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_dnam.Rds')
dnam = dnam[None]

controls = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_controls.Rds')
controls = controls[None]

snps = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_snps.Rds')
snps = snps[None]

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
    samples = samples.loc[samples.isnull().mean(axis=1) < 0.1]
    
    return samples
    
#-----------------------------------------------------------------------------------------#

# Extract BC controls and perform multivariate outlier identification

from sklearn import preprocessing

# Subset the 'bc1.grn','bc1.red','bc2'
tech_vars=samples[['bc1.grn','bc1.red','bc2']]

# Standardise the tech_vars
tech_vars_stand = preprocessing.scale(tech_vars)
tech_vars_stand = pd.DataFrame(tech_vars_stand)
tech_vars_stand.index=tech_vars.index
tech_vars_stand.columns=tech_vars.columns

from sklearn.decomposition import PCA as sklearnPCA
import plotly.plotly as py

sklearn_pca = sklearnPCA(n_components=2)
tech_vars_pca = sklearn_pca.fit_transform(tech_vars_stand)


tech_vars_pca_DF=pd.DataFrame(tech_vars_pca)

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