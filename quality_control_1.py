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