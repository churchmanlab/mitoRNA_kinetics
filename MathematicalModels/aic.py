#!/usr/bin/env python
# coding: utf-8

# ## aic.py
# Author: Robert Ietswaart, Gyan Prakash
# date last updated: 20221122
# License: BSD2.  
# Python v3.7.4 
# 
# For Erik McShane's project: 1 and 2 state RNA decay models: Akaike Information Criterium.
# 
# Starting point references: theory from McShane et al, Cell 2016. Sin et al, PLoS ONE, 2016.

import numpy as np

def _aic(n,k,rss):
    """first equation on page e7 of McShane et al
       k : number of parameters in model"""
    aic = 2.0*k + n*np.log(rss/n) + 2.0*k*(k+1)/(n-k-1)
    return aic
 
def calc_aic(n, rss1, rss2,k1,k2):
    """Calculate Akaike Information Criterium (AIC)
       from residuals for 1 and 2 state models
       n : number of data points
       rss1 : residual sum of squares of model 1
       rss2 : as rss1 for model 2
       k1 : number of free parameters of model 1
       k2 : as k1 for model 2
       """
    aic1 = _aic(n, k1, rss1) 
    aic2 = _aic(n, k2, rss2)
    aic_min = min(aic1, aic2)
    sum_of_probs = np.exp((aic_min - aic1) / 2) + np.exp((aic_min - aic2) / 2)
    prob_aic1 = np.exp((aic_min - aic1) / 2) / sum_of_probs
    prob_aic2 = np.exp((aic_min - aic2) / 2) / sum_of_probs
    return [prob_aic1, prob_aic2]

def _aic_nocr(n,k,rss):
    """Akaike Information Criterium without low number correction
       k : number of parameters in model"""
    aic = 2.0*k + n*np.log(rss/n)
    return aic

def calc_aic_nocr(n, rss1, rss2, k1, k2):
    aic1 = _aic_nocr(n, k1, rss1) 
    aic2 = _aic_nocr(n, k2, rss2)
    aic_min = min(aic1, aic2)
    sum_of_probs = np.exp((aic_min - aic1) / 2) + np.exp((aic_min - aic2) / 2)
    prob_aic1 = np.exp((aic_min - aic1) / 2) / sum_of_probs
    prob_aic2 = np.exp((aic_min - aic2) / 2) / sum_of_probs
    
    return [prob_aic1, prob_aic2]