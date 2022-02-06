#!/usr/bin/env python
# coding: utf-8

# ## new_total_ratio.py
# Author: Robert Ietswaart
# Date: 20220206
# License: BSD2.  
# Python v3.7.4 
# 
# MitoRibo Subcellular Timelapse Seq: fit time scales to his subcellular timelapse seq using observed new to total ratio (Lambda)
# of the mitoribo and mitototal subcellular RNA fractions


import numpy as np
from numba import jit, float64

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_total_one_step(kD, t):
    """ Fit total RNA degradation rate
        through a one step process
        t: time
        kD: (total) degradation rate
    """
    Lambda = 1 - np.exp(-kD * t)
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_total_one_step<0: ', Lambda, kD)
        Lambda = 0
    return Lambda

#MitoRibo model
@jit(float64(float64, float64, float64), nopython=True, cache=True)
def lam_ribo(kL, kD, t):
    """ Fit of ribosome entry rate using ribosomal compartment 
        t: time point (float)
        kL: mitoribosome entry rate
        kD: mitochondrial RNA Degradation rate
    """
    eps = 1e-16
    if np.absolute(kL) <= eps: #limit kL = 0 case
        Lambda = 1 - np.exp(-kD * t) * (1 + kD * t)
    else:
        Lambda = 1 - (kL + kD) / kL * np.exp(-kD * t) + kD / kL * np.exp(-(kD + kL) * t)
    if Lambda < 0:#numerical under/overflow
        Lambda = 0
    return Lambda
