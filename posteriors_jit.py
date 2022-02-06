#!/usr/bin/env python
# coding: utf-8

# ## posteriors_jit.py
# Author: Robert Ietswaart
# Date: 20220206
# License: BSD2.  
# Python v3.7.4 
# 
# MitoRibo Subcellular Timelapse Seq: Bayesian framework to fit time scales (Mean, MAP, 95%CIs) 
# using observed new to total ratio from GRAND-SLAM of the mitototal and mitoribo subcellular RNA fractions


import numpy as np
from scipy.special import beta
from numba import jit, float64


@jit(float64(float64, float64, float64), nopython=True) #NumbaWarning: Cannot cache compile: when you include cache=True
def beta_pdf(x, a, b):
    """defined on domain of x: [0,1]
    function structure as in scipy.stats.beta: 
    https://github.com/scipy/scipy/blob/v1.7.0/scipy/stats/_boost/include/func_defs.hpp
    https://www.boost.org/doc/libs/1_61_0/boost/math/distributions/beta.hpp
    
    Gaussian approximations in case of numerical underflow for beta function with large a ,b:
    https://stats.stackexchange.com/questions/2066/how-can-i-numerically-approximate-values-for-a-beta-distribution-with-large-al
    """
    if (x >= 0) and (x <= 1): 
        if (x == 1) and (b < 1):
            return np.inf
        elif (x == 0) and (a < 1):
            return np.inf
        _beta = beta(a, b)
        if _beta > 0:
            return x**(a - 1) * (1 - x)**(b - 1) / _beta
        else:
            #normal distribution approximation because of numerical underflow
            mu = a * (a + b)**(-1)
            sig = np.sqrt(a * b * ((a + b)**(2) * (1 + a + b))**(-1))
            return np.exp(-0.5 * ((x - mu) / sig)**(2)) * (sig * np.sqrt(2 * np.pi))**(-1) 
    else:
        return 0    

@jit(float64(float64, float64), nopython=True, cache=True)
def dlam_1var_da(a, t):
    """ Derivative of Lambda with respect to (only) var a at t
        t: time point (float)  
        a: rate a in domain [0,inf)
        der: derivative, function on a domain
    """
    der = t * np.exp(- a * t)
    return der

@jit(float64(float64, float64), nopython=True, cache=True)
def dlam_total_dkDtot(kDtot, t):
    """ Derivative of Lambda with respect to kDtot at t
        t: time point (float)
        kDtot: Cellular residence (Total Degradation) rate
        der: derivative (np.array), function on k domain
    """
    der = dlam_1var_da(kDtot, t)
    return der

@jit(float64(float64, float64, float64), nopython=True, cache=True)
def dlam_ribo_dkL(kL, kD, t):
    """ Derivative of Lambda with respect to kL at t
        t: time point (float)
        kL: ribosome entry rate
        kD: mitochondrial RNA degradation rate
        der: derivative (np.array), function on k domain
    """
    eps = 1e-16
    if np.absolute(kL) <= eps:
        der = 0 #limit approximation
    else:
        der = kD * kL**(-2) * np.exp(-kD * t) * (1 - (1 + t * kL) * np.exp(-kL * t))
    return der

@jit(float64(float64, float64, float64), nopython=True, cache=True)
def det_jac_mitoribo_from_tot(kL, kDtot, t):
    det = np.absolute(dlam_ribo_dkL(kL, kDtot, t) * dlam_total_dkDtot(kDtot, t))
    return det

@jit(float64(float64, float64), nopython=True, cache=True)
def det_jac_mitototal(kDtot, t):
    return dlam_total_dkDtot(kDtot, t)