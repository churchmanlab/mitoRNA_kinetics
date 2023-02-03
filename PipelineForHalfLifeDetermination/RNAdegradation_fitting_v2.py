#!/usr/bin/env python
# coding: utf-8
# USE
# python /Users/Mary/Desktop/Data/TimelapseSeq/Scripts/RNAdegradation_fitting.py /Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/ TL3_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew

# v2 has modified parameter values for the 2-state model, and uses the leastsquare method instead of ampgo for the minimize function

# In[1]:


import csv 
import ast
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import binom
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
from math import factorial as fac
from lmfit import minimize, Parameters, fit_report
import seaborn as sns
from mpl_toolkits.axes_grid1 import Divider, Size
import sys
import math
import argparse # to read in list from command line


CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--path",  # name on the CLI - drop the `--` for positional/required parameters
  type=str,
  default='',  # default if nothing is provided
)
CLI.add_argument(
  "--file",  # name on the CLI - drop the `--` for positional/required parameters
  type=str,
  default='',  # default if nothing is provided
)
CLI.add_argument(
  "--doubletime",  # name on the CLI - drop the `--` for positional/required parameters
  type=int,
  default=10000,  # default if nothing is provided
)
CLI.add_argument(
  "--tps",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=int,
  default=[],  # default if nothing is provided
)
CLI.add_argument(
  "--ng",  # name on the CLI - drop the `--` for positional/required parameters
  type=str,
  default='',  # default if nothing is provided
)
CLI.add_argument(
  "--toplot",  # name on the CLI - drop the `--` for positional/required parameters
  type=str,
  default='',  # default if nothing is provided
)

CLI.add_argument(
  "--genome",  # name on the CLI - drop the `--` for positional/required parameters
  type=str,
  default='',  # default if nothing is provided
)


# parse the command line
args = CLI.parse_args()

# In[17]:

path = args.path
filename = args.file
doubletime = args.doubletime # in minutes
tps_to_use = args.tps
normg = args.ng
toplot = args.toplot
genome = args.genome

print(path)
print(filename)
print(doubletime)
print(tps_to_use)
# loads the csv and removes the NaN values 
TC_MAP_in = pd.read_csv(path + "FracNew/" + filename + ".txt" ,sep='\t',header=(0)).dropna()
time_all = pd.Series(tps_to_use)
# time_all = pd.Series([0,15,30,60,120,240])


# Correction for cell doubling 
CDT = int(doubletime) ### Cell doubling time in minutes [constant]
appended_data=[]
col=[]
data=[]

for Gene_Name in TC_MAP_in.Symbol:
    data.append(Gene_Name)
    for TP in time_all:
        a=TC_MAP_in[TC_MAP_in.Symbol == Gene_Name].index[0]
        MFN=TC_MAP_in.loc[a].at[str(TP)] # measured fraction new at time point
        CFN = (MFN - (1-(1/(math.e**(np.log(2)/CDT))**TP)))/(1-(1-(1/(math.e**(np.log(2)/CDT))**TP))) 
        
#Corrected fraction new at timepoint
        data.append(CFN)
    appended_data.append(data)
    data = []
TC_MAP = pd.DataFrame(appended_data, columns=list(TC_MAP_in.columns))
TC_MAP.to_csv(path + "FracNew/" + filename + "_corr_" + str(doubletime) + "min.txt", sep="\t", index=False)




# In[18]:


#defining the RNA degradation functions

def lam1state_toy_model(params,t):
    """ Model without 4sU dynamics
        t: time
        kA: RNA degradation rate from state A
        Derived in Sin et al, PLoS ONE 2016
    """
    kA = params['degradationA']
    Lambda = 1 - np.exp(-kA * t)
    return Lambda

def lam2state_toy_model(params,t):
    """ Model without 4sU dynamics
        t: time
        kA: RNA degradation rate from state A
        kAB: transition rate from state A to state B
        kB: RNA degradation rate from state B
        Derived in Sin et al, PLoS ONE 2016: eq 19
    """
    kA = params['degradationA']
    kAB = params['transitionAB']
    kB = params['degradationB']
    Ap = kB*(kA-kB)
    Bp = kAB*(kA+kAB)
    Lambda = Ap/(Ap+Bp)*(1 - np.exp(-(kA+kAB) * t)) + Bp/(Ap+Bp)*(1 - np.exp(-kB * t))
    return Lambda


# In[19]:


def _aic1(n,k,rss):
    """first equation on page e7 of McShane et al
       k : number of parameters in model"""
    aic = 2.0*k + n*np.log(rss/n)+ 2.0*k*(k+1)/(n-k-1)
    return aic

#+ 2.0*k*(k+1)/(n-k-1) (As I was reading AIC under the Framework of Least Squares Estimation, it suggests to use this 
#correction parameter when n/K < 40 , so Iw asnt sure.Also suggests different AIC models for OLS optimization )

def _aic_nocr(n,k,rss):
    """first equation on page e7 of McShane et al
       k : number of parameters in model"""
    aic = 2.0*k + n*np.log(rss/n)
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
    aic1 = _aic1(n, k1, rss1) 
    aic2 = _aic1(n, k2, rss2)
    aic_min = min(aic1, aic2)
    sum_of_probs = np.exp((aic_min - aic1) / 2) + np.exp((aic_min - aic2) / 2)
    prob_aic1 = np.exp((aic_min - aic1) / 2) / sum_of_probs
    prob_aic2 = np.exp((aic_min - aic2) / 2) / sum_of_probs
    
    return [prob_aic1, prob_aic2]

def calc_aic_nocr(n, rss1, rss2,k1,k2):
    aic1 = _aic_nocr(n, k1, rss1) 
    aic2 = _aic_nocr(n, k2, rss2)
    aic_min = min(aic1, aic2)
    sum_of_probs = np.exp((aic_min - aic1) / 2) + np.exp((aic_min - aic2) / 2)
    prob_aic1 = np.exp((aic_min - aic1) / 2) / sum_of_probs
    prob_aic2 = np.exp((aic_min - aic2) / 2) / sum_of_probs
    
    return [prob_aic1, prob_aic2]


# In[20]:


def res_fitting(k,t, fit_model, data):
    predict_data =fit_model(k, t)
    return (predict_data - data)
   


# In[21]:


#defining the intial values for model 1 and model 2 in the Parameters. Out3 is the probability values from the AIC function for 
# Model1 vs Model 2. 
#Then calling dict in which we are storing all the three outcome on which we will perform the iteration. 

def final_outcome(input,gene):
    params = Parameters()
    params.add('degradationA', value=0.001 , min=0, max=np.inf)
    
    out1 = minimize(res_fitting, params, method='leastsquare', args=(time_all,lam1state_toy_model,input),nan_policy='omit')
    
    res1=out1.residual
    RSS_1 = sum(res1*res1)
    
    if len(time_all) > 5:    
        report1=fit_report(out1)
    
    
    params = Parameters()
    params.add('degradationA', value=0.001, min=0, max=np.inf)
    params.add('transitionAB', value=0.0005, min=0, max=np.inf)
    params.add('degradationB', value=0.00001, min=0, max=np.inf)
    
    
    out2 = minimize(res_fitting, params, method='leastsquare', args=(time_all,lam2state_toy_model,input),nan_policy='omit' )
    
    res2=out2.residual
    RSS_2 = sum(res2*res2)
    print(RSS_2)
    
    if len(time_all) > 5:
        report2=fit_report(out2)
    
    out2Arr=[out2.params['degradationA'].value,
    out2.params['transitionAB'].value,
    out2.params['degradationB'].value]
    
    out4= calc_aic_nocr(out2.ndata,out1.chisqr,out2.chisqr,out1.nvarys,out2.nvarys)
    state = ("one-state" if out4[0]>=0.50 else "two-state")

    # Only do correction when there are at least 5 time points (plus 0m)
    if len(time_all) > 5:
        out3 =calc_aic(out2.ndata,out1.chisqr,out2.chisqr,out1.nvarys,out2.nvarys)
    
    
        state = ("one-state" if out3[0]>=0.50 else "two-state")
    
        allOutcomes = {
          "Gene":gene,
          "Fit1 (Ka)": out1.params['degradationA'].value,
          "Fit2 (Ka, Kab, Kb)": out2Arr,
#           "AIC_correction (model1 vs model2)":out3, 
          "AIC_correction model1":out3[0],
          "AIC_correction model2":out3[1],
          "AIC model1":out4[0],
          "AIC model2":out4[1],
          "res1":res1,
          "res2":res2,
          "repor1":report1,
          "repor2":report2,
          "State":state,
          "RSS_1":RSS_1,
          "RSS_2":RSS_2
        }
    elif len(time_all) < 6:
        allOutcomes = {
          "Gene":gene,
          "Fit1 (Ka)": out1.params['degradationA'].value,
          "Fit2 (Ka, Kab, Kb)": out2Arr,
          "AIC model1":out4[0],
          "AIC model2":out4[1],
          "res1":res1,
          "res2":res2,
          "State":state,
          "RSS_1":RSS_1,
          "RSS_2":RSS_2
        }

    
    return allOutcomes


# In[22]:


# here we are defining a fucntion, which loops over different time points to provide MAP value for a gene. Then we pass 
#the output of the data to the final_outcome function to return a dataframe with all the genes. 

def ex_data(gene):

    col=[]
    data=[]
    Gene_Name=gene

    for t in time_all:
        a=TC_MAP[TC_MAP.Symbol == Gene_Name].index[0]
        col=TC_MAP[str(t)]
        data.append(col[a])

    input=pd.Series(data)
    return final_outcome(input,gene)

#ex_data('mt-co1')


# In[24]:


appended_data=[]
for symbol in (TC_MAP.Symbol):
    data = ex_data(symbol)
    appended_data.append(data)
    
# final_df = pd.DataFrame(appended_data, columns=["Gene", "Fit1 (Ka)","Fit2 (Ka, Kab, Kb)","AIC_correction (model1 vs model2)","State","AIC (model1 vs model2)"])
final_df = pd.DataFrame(appended_data, columns=["Gene", "Fit1 (Ka)","Fit2 (Ka, Kab, Kb)","AIC_correction model1", "AIC_correction model2","State","AIC model1", "AIC model2", "RSS_1", "RSS_2"])
final_df


# In[25]:


def data_plt(gene):

    col=[]
    data=[]
    Gene_Name=gene

    for t in time_all:
        a=TC_MAP[TC_MAP.Symbol == Gene_Name].index[0]
        col=TC_MAP[str(t)]
        data.append(col[a])

    return pd.Series(data)


# In[26]:


def plot_model (gene_name):
    
    x = time_all
    y = data_plt(gene_name)
    y1= final_outcome(data_plt(gene_name),gene_name)['res1'] + data_plt(gene_name)
#     y2= final_outcome(data_plt(gene_name),gene_name)['res2'] + data_plt(gene_name)

    fig = plt.figure(figsize=(8, 6))
    plt.plot(x, y, 'bo', label= gene_name)
    plt.plot(x, y1, 'r--', label='best fit_Model1')
#     plt.plot(x, y2, 'g-', label='best fit_Model2')    
    plt.legend(loc='best')
    plt.xlabel('Time (mins)')
    plt.ylabel('fraction 4sU labeled') 

#     return plt.show()
    return fig

def plot_models_1state (gene_name1, gene_name2, gene_name3, c1, c2, c3):
    
    x = time_all
    y = data_plt(gene_name1)
    y1= final_outcome(data_plt(gene_name1),gene_name1)['res1'] + data_plt(gene_name1)
    z= data_plt(gene_name2)
    z1= final_outcome(data_plt(gene_name2),gene_name2)['res1'] + data_plt(gene_name2)
    a= data_plt(gene_name3)
    a1= final_outcome(data_plt(gene_name3),gene_name3)['res1'] + data_plt(gene_name3)
    fig = plt.figure(figsize=(8, 6))
    plt.plot(x, y, 'o', label= gene_name1, color=c1)
    plt.plot(x, y1, '-', label='best fit_Model1', color=c1)
    plt.plot(x, z, 'o', label= gene_name2, color=c2)
    plt.plot(x, z1, '-', label='best fit_Model1', color=c2)
    plt.plot(x, a, 'o', label= gene_name3, color=c3)
    plt.plot(x, a1, '-', label='best fit_Model1', color=c3)
    plt.legend(loc='best')
    plt.xlabel('Time (mins)')
    plt.ylabel('fraction 4sU labeled') 
    return fig

def plot_models_2state (gene_name1, gene_name2, c1, c2):
    
    x = time_all
    y = data_plt(gene_name1)
    y1= final_outcome(data_plt(gene_name1),gene_name1)['res1'] + data_plt(gene_name1)
    y2= final_outcome(data_plt(gene_name1),gene_name1)['res2'] + data_plt(gene_name1)
    z= data_plt(gene_name2)
    z1= final_outcome(data_plt(gene_name2),gene_name2)['res1'] + data_plt(gene_name2)
    z2= final_outcome(data_plt(gene_name2),gene_name2)['res2'] + data_plt(gene_name2)

    fig = plt.figure(figsize=(8, 6))
    plt.plot(x, y, 'o', label= gene_name1, color=c1)
    plt.plot(x, y1, '--', label='best fit_Model1', color=c1)
    plt.plot(x, y2, '-', label='best fit_Model2', color=c1)
    plt.plot(x, z, 'o', label= gene_name2, color=c2)
    plt.plot(x, z1, '--', label='best fit_Model1', color=c2)
    plt.plot(x, z2, '-', label='best fit_Model2', color=c2)
    plt.legend(loc='best')
    plt.xlabel('Time (mins)')
    plt.ylabel('fraction 4sU labeled') 
    return fig

def plot_models_other (gene_name1, gene_name2, gene_name3, c1, c2, c3):
    
    x = time_all
    y = data_plt(gene_name1)
    y1= final_outcome(data_plt(gene_name1),gene_name1)['res1'] + data_plt(gene_name1)
    y2= final_outcome(data_plt(gene_name1),gene_name1)['res2'] + data_plt(gene_name1)
    z= data_plt(gene_name2)
    z1= final_outcome(data_plt(gene_name2),gene_name2)['res1'] + data_plt(gene_name2)
    z2= final_outcome(data_plt(gene_name2),gene_name2)['res2'] + data_plt(gene_name2)
    w= data_plt(gene_name3)
    w1= final_outcome(data_plt(gene_name3),gene_name3)['res1'] + data_plt(gene_name3)
    w2= final_outcome(data_plt(gene_name3),gene_name3)['res2'] + data_plt(gene_name3)

    fig = plt.figure(figsize=(8, 6))
    plt.plot(x, y, 'o', label= gene_name1, color=c1)
    plt.plot(x, y1, '--', label='best fit_Model1', color=c1)
    plt.plot(x, y2, '-', label='best fit_Model2', color=c1)
    plt.plot(x, z, 'o', label= gene_name2, color=c2)
    plt.plot(x, z1, '--', label='best fit_Model1', color=c2)
    plt.plot(x, z2, '-', label='best fit_Model2', color=c2)
    plt.plot(x, w, 'o', label= gene_name3, color=c3)
    plt.plot(x, w1, '--', label='best fit_Model1', color=c3)
    plt.plot(x, w2, '-', label='best fit_Model2', color=c3)
    plt.legend(loc='best')
    plt.xlabel('Time (mins)')
    plt.ylabel('fraction 4sU labeled') 
    return fig



# In[27]:

if toplot == '1state':
    g1 = 'MT-CO1'
    c1 = 'lightcoral'
    g2 = 'MT-ND1'
    c2 = 'lightskyblue'
    g3 = 'MT-RNR2'
    c3 = 'gold'
    if genome == 'mouse':
        g1 = 'mt-Co1'
        g2 = 'mt-Nd1'
        g3 = 'mt-Rnr2'
    fig = plot_models_1state(g1, g2, g3, c1, c2, c3)
    fig.savefig(path + "HalfLife/" + filename + "_halflives_corr_" + str(doubletime) + "_modelfit1state_v2.pdf")
        

if toplot == '2state':
    g1 = 'MT-ND6'
    c1 = 'darkgoldenrod'
    g2 = 'MT-antiND1'
    c2 = 'fuchsia'
    if genome == 'mouse':
        g1 = 'mt-Nd6'
        g2 = 'mt-antiNd1'
    fig = plot_models_2state(g1, g2, c1, c2)
    fig.savefig(path + "HalfLife/" + filename + "_halflives_corr_" + str(doubletime) + "_modelfit2state_v2.pdf")


if toplot == '2state' or toplot == '1state':
    g1 = 'MT-ND6'
    c1 = 'blue'
    g2 = 'MT-ND3'
    c2 = 'cornflowerblue'
    g3 = 'MT-CO3'
    c3 = 'lightcoral'
    if genome == 'mouse':
        g1 = 'mt-Nd6'
        g2 = 'mt-Nd3'
        g3 = 'mt-Co3'
    fig = plot_models_other(g1, g2, g3, c1, c2, c3)
    fig.savefig(path + "HalfLife/" + filename + "_halflives_corr_" + str(doubletime) + "_modelfitother_v2.pdf")

# In[38]:


# merge_plot=pd.merge(final_df, on="Gene")[["Gene","Fit1 (Ka)"]]
merge_plot=final_df[["Gene","Fit1 (Ka)"]]

merge_plot["Half Life"]= 0.693/merge_plot["Fit1 (Ka)"]

# Gene = ["MT-ATP8-6", "MT-CO1", "MT-CO2","MT-CO3", "MT-ND1",
#         "MT-ND2", "MT-ND3", "MT-ND4L-4", "MT-ND5" ,"MT-CYB"]#, "mt-rnr2"]#,"mt-rnr1"]

# merge_plot=merge_plot.reindex(merge_plot[merge_plot.Gene == gene].index[0] for gene in Gene).reset_index(drop=True)


# In[39]:


merge_plot["Normalized Half-Life"] = merge_plot["Half Life"]/merge_plot["Half Life"][merge_plot[merge_plot.Gene == normg].index[0]]
# merge_plot["Normalized Abundance"] = merge_plot["Abundance"]/merge_plot["Abundance"][merge_plot[merge_plot.Gene == "MT-CO1"].index[0]]
# merge_plot ["Expected abundance1"] = merge_plot["Normalized Abundance"]/merge_plot["Normalized Half-Life"]


# In[113]:


merge_plot.to_csv(path + "HalfLife/" + filename + "_halflives_corr_" + str(doubletime) + "min_v2.txt", sep="\t", index=False)

# Write probabilities for each model fit

prob_df = final_df.loc[0:,["Gene", "AIC_correction model1", "AIC_correction model2", "AIC model1", "AIC model2", "RSS_1", "RSS_2"]]
prob_df.to_csv(path + "HalfLife/" + filename + "_halflives_corr_" + str(doubletime) + "min_modelprobs_v2.txt", sep="\t", index=False, quoting=None)