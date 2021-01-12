#!/usr/bin/env python
# coding: utf-8

# ## new_total_ratio.py
# Author: Robert Ietswaart, 20200402
# License: BSD2.  
# Python v3.7.4 
# 
# For Erik McShane's project: 1 and 2 state RNA decay models: new to total ratio Lambda.
# 
# Starting point references: theory from McShane et al, Cell 2016. Sin et al, PLoS ONE, 2016.

import numpy as np

def lam1state(k,t):
    """ Model includes 4sU dynamics: derivations by Robert Ietswaart, N6p63-80
        t: time
        kout: 4sU/U exchange out rate
        kA: RNA degradation rate from state A
    """
    kout = k[0]
    kA = k[1]
    Lambda = -kA/(kout-kA)*(1 - np.exp(-kout * t)) + \
             kout/(kout-kA)*(1 - np.exp(-kA * t))
    return Lambda

def lam2state(k,t):
    """ Model includes 4sU dynamics: derivations by Robert Ietswaart, N6p63-80
        t: time
        kout: 4sU/U exchange out rate
        kA: RNA degradation rate from state A
        kAB: transition rate from state A to state B
        kB: RNA degradation rate from state B
    """
    kout = k[0]
    kA = k[1]
    kAB = k[2]
    kB = k[3]
    Phi = kB*(kA+kAB)*(kB+kAB-kout)/((kB+kAB)*(kA+kAB-kout)*(kB-kout))
    Psi = kB*(kA-kB)*kout/((kB+kAB)*(kout-(kA+kAB))*(kA+kAB-kB))
    Omega = kAB*kout*(kA+kAB)/((kB+kAB)*(kA+kAB-kB)*(kout-kB)) #1-Phi-Psi
    Lambda = Phi*(1 - np.exp(-kout * t)) + \
             Psi*(1 - np.exp(-(kA+kAB) * t)) + \
             Omega*(1 - np.exp(-kB * t))
    return Lambda

def lam1state_toy_model(k,t):
    """ Model without 4sU dynamics
        t: time
        kA: RNA degradation rate from state A
        Derived in Sin et al, PLoS ONE 2016
    """
    kA = k[1]
    Lambda = 1 - np.exp(-kA * t)
    return Lambda

def lam2state_toy_model(k,t):
    """ Model without 4sU dynamics
        t: time
        kA: RNA degradation rate from state A
        kAB: transition rate from state A to state B
        kB: RNA degradation rate from state B
        Derived in Sin et al, PLoS ONE 2016: eq 19
    """
    kA = k[1]
    kAB = k[2]
    kB = k[3]
    Ap = kB*(kA-kB)
    Bp = kAB*(kA+kAB)
    Lambda = Ap/(Ap+Bp)*(1 - np.exp(-(kA+kAB) * t)) + \
             Bp/(Ap+Bp)*(1 - np.exp(-kB * t))
    return Lambda

def TCconv(k,t):
    """ 4sU dynamics
        t: time
        kin: 4sU in rate
        kout: 4sU out rate derivation N6p184: TC = int_{0}^{t} L(s)ds
    """
    kin = k[0]
    kout = k[1]
    TC = kin / kout *( 1 + (np.exp(-kout * t) - 1) / (kout * t) )
    #old: wrong derivation again N6p183: this is L(t) = 4sU(t)/U(t)
#     TC = kin * (1 - np.exp(-kout * t)) / kout 
    return TC
