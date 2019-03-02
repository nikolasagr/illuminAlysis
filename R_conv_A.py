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

control_beads = pyreadr.read_r('/Users/nicolasagrotis/Desktop/CompEpi/PR/hm450_controls.Rds')
control_beads = control_beads[None]

# Load DNA Data

dnam = pyreadr.read_r('/Users/nicolasagrotis/Desktop/CompEpi/PR/eira_dnam.Rds')
dnam = dnam[None]

controls = pyreadr.read_r('/Users/nicolasagrotis/Desktop/CompEpi/PR/eira_controls.Rds')
controls = controls[None]

snps = pyreadr.read_r('/Users/nicolasagrotis/Desktop/CompEpi/PR/eira_snps.Rds')
snps = snps[None]

samples = pyreadr.read_r('/Users/nicolasagrotis/Desktop/CompEpi/PR/eira_samples.Rds')
samples = samples[None]

# Load individual characteristics data

covars = pyreadr.read_r('/Users/nicolasagrotis/Desktop/CompEpi/PR/Covariates.Rds')
covars = covars[None]

samplesheet = pyreadr.read_r('/Users/nicolasagrotis/Desktop/CompEpi/PR/Sample_sheet.Rds')
samplesheet = samplesheet[None]

# Remove samples with >10% NAs
# Function 


# Remove 
df.columns[df.isnull().mean() < 0.8]
