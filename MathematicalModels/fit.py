#!/usr/bin/env python
# coding: utf-8

# ## fit.py
# Author: Robert Ietswaart
# date last updated: 20221122
# License: BSD2.  
# Python v3.7.4 
# 
# For Erik's mitoRNA project

import numpy as np
import new_total_ratio as ntr


def res_2models(k, times, model1, model2, lam_obs1, lam_obs2):
    """Calculate Residual
       for input models
       lam_obs : observed new to total ratio lambda
       k_fit : parameter to fit in model
       k_para: other parameters in model
       times : list with observed time points >0
       model1 : function describing mathematical model1 (mito total RNA)
       model2 : function describing mathematical model2 (mito ribo IP RNA)
       """
    lam_predict1 = model1(k, times)
    lam_predict2 = model2(k, times)
    res1 = lam_predict1 - lam_obs1
    res2 = lam_predict2 - lam_obs2
    return res1.append(res2, ignore_index=True)