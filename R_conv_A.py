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

# Remove samples with >10% NAs
# Function 


#
df.columns[df.isnull().mean() < 0.8]


#-----------------------------------------------------------------------------------------#

# set the columns and the index of the datasets

samples=samples.set_index("sample.id")

# Before 689, after 673
samples.iloc[:,:].isnull().sum()
samples = samples.loc[samples.isnull().mean(axis=1) < 0.1]

#-----------------------------------------------------------------------------------------#

# Extract BC controls and perform multivariate outlier identification

from sklearn import preprocessing

tech_vars=samples[['bc1.grn','bc1.red','bc2']]












