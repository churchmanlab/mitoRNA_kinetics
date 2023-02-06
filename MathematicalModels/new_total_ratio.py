#!/usr/bin/env python
# coding: utf-8

# ## new_total_ratio.py
# Author: Robert Ietswaart
# date last updated: 20221122
# License: BSD2.  
# Python v3.7.4 
# 
# For Erik McShane's project: 1 and 2 state RNA decay models: new to total ratio Lambda.
# 
# Starting point references: theory from McShane et al, Cell 2016. Sin et al, PLoS ONE, 2016.

import numpy as np
k_bound_lo = 1e-4 #unit: min^-1: 1 per 7 days
k_bound_hi = 1e4 #unit: min^-1: 1 per 6 ms, if too restrictive: increase to 1e6

def lam1state_total(k,t):
    """ Model without 4sU dynamics
        t: time
        kD: mitochondrial RNA degradation rate of state A and B
        Derived in Sin et al, PLoS ONE 2016
    """
    kD = k['degradation']
    Lambda = 1 - np.exp(-kD * t)
    return Lambda

def lam_ribo_1deg(k, t):
    """ Fit of ribosome entry rate using ribosomal compartment 
        t: time point (float)
        kAB: mitoribosome entry rate
        kD: mitochondrial RNA Degradation rate: set equal in state A and B
    """
    kD = k['degradation']
    kAB = k['transitionAB']
    eps = k_bound_lo #1e-16
    if np.absolute(kAB) <= eps: #limit kAB = 0 case
        Lambda = 1 - np.exp(-kD * t) * (1 + kD * t)
    else:
        Lambda = 1 + kD / kAB * np.exp(-(kD + kAB) * t) - (kAB + kD) / kAB * np.exp(-kD * t) 
    return Lambda

def lam2state_total(k,t):
    """ Model without 4sU dynamics
        t: time
        kA: RNA degradation rate from state A
        kAB: transition rate from state A to state B
        kB: RNA degradation rate from state B
        Derived in Sin et al, PLoS ONE 2016: eq 19
    """
    kA = k['degradationA']
    kAB = k['transitionAB']
    kB = k['degradationB']
    Ap = kB * (kA - kB)
    Bp = kAB * (kA + kAB)
    Lambda = Ap / (Ap + Bp) * (1 - np.exp(-(kA + kAB) * t)) + \
             Bp/(Ap + Bp) * (1 - np.exp(-kB * t))
    return Lambda

def lam2state_ribo(k,t):
    """ Model without 4sU dynamics
        t: time
        kA: RNA degradation rate from state A
        kAB: transition rate from state A to state B
        kB: RNA degradation rate from state B
        Notebook N8p98
    """
    kA = k['degradationA']
    kAB = k['transitionAB']
    kB = k['degradationB'] 
    Lambda = 1 - kB / (kB - (kA + kAB)) * np.exp(-(kA + kAB) * t) + \
             (kA + kAB) / (kB - (kA + kAB)) * np.exp(-kB * t)
    return Lambda
