import csv
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import binom
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from math import factorial as fac
from lmfit import minimize, Parameters, fit_report, Minimizer

import matplotlib.pyplot as plot
import matplotlib.patches as mpatches
import sys

######### USE
### Uncomment this section, update rates and paste entire block into terminal. Directory structure: run from directory called 'FracNew' (with logs folder), with matrices in the directory 'MMfrequency' which is on the same level as FracNew


# source ~/Python_Venv/bin/activate
# unset PYTHONPATH # to keep it from using the old version of numpy


##########

Exp = sys.argv[1]  #'TL8'
bgLibName = sys.argv[2]
LibName = sys.argv[3]  #'60m_tot'
MapMethod = sys.argv[4]
rateRegion = sys.argv[5]
useConv = sys.argv[6]
dist = sys.argv[7]
RatesOnly = sys.argv[8]
stringentFilter = sys.argv[9]

############# First section for getting rates #####################
#################### from background matrices #####################

# Step one is to use the 0 min / negative control sample matrix. Then you run that through the first three steps. 1) a simple one-rate model. 2) a model using two background rates. 3) the AIC which decides if we should continue with the one-rate or two-rate estimates
# 0 time point, to get background rates

if dist == 'orig':
    Tmatrix = '_TcountANDTCconv_'
elif dist == 'mod':
    Tmatrix = '_TcountANDTCconv_modTcountDist_'

# TC matrix
dfTC = pd.read_csv(
    "../MMfrequency/"
    + Exp
    + "_"
    + MapMethod
    + "_frag_"
    + rateRegion
    + Tmatrix
    + Exp
    + "_"
    + bgLibName
    + "_otherMMfilt_"
    + stringentFilter
    + ".txt",
    sep="\t",
    header=(0),
)
dfTC = dfTC.rename(columns={'Unnamed: 0':'Tcount'})


if useConv != 'TC':
    # CT matrix
    dfCT = pd.read_csv(
        "../MMfrequency/"
        + Exp
        + "_"
        + MapMethod
        + "_frag_"
        + rateRegion
        + "_CcountANDCTconv_"
        + Exp
        + "_"
        + bgLibName
        + "_otherMMfilt_"
        + stringentFilter
        + ".txt",
        sep="\t",
        header=(0),
    )
    dfCT = dfCT.rename(columns={'Unnamed: 0':'Ccount'})




# print(dfTC)
# print(dfCT)

# For TC
n=len(dfTC['Tcount'])
k=len(dfTC.columns)
# Make numpy array from the data frame
dataTC=np.zeros((n,k-1))

for n in range(len(dfTC['Tcount'])):
    for k in range(1,len(dfTC.columns)):
        dataTC[n,k-1]= dfTC.iloc[n,k]
print('TC background n: ', n)
print('TC background k: ', k)

if useConv != 'TC':
    # For CT
    n=len(dfCT['Ccount'])
    k=len(dfCT.columns)
    # Make numpy array from the data frame
    dataCT=np.zeros((n,k-1))

    for n in range(len(dfCT['Ccount'])):
        for k in range(1,len(dfCT.columns)):
            dataCT[n,k-1]= dfCT.iloc[n,k]

    print('CT background n: ', n)
    print('CT background k: ', k)

# Function for one-state model for TC
def input_conv(params, df, data, count): # ex. params, dfTC, dataTC, 'Tcount'
    pe1 = params['Error_1']

    n=len(df[count])
    k=len(df.columns)
    model_res=np.zeros((n,k-1))

    def gs_func(x,y):
        return binom.pmf(y,x,pe1) # (range from 0 to n, n number of trials, probability of success)
    
    for x in range(1,len(df[count])):
        for y in range(len(df.columns)-1):
            if x>=y:
                model_res[x-1,y] = gs_func(x,y) * sum(data[x-1,])
                   
    return model_res

# params = Parameters()
# Set initial value?
# params.add('Error_1', value=0.5, min=0, max=1)

# Function for two-state model
def Two_conv_background(params, df, data, count): #params, dfTC, dataTC, 'Tcount'
    pe1 = params['Error_1']
    pe2 = params['Error_2']
    pi = params['frac_err']
    
    n=len(df[count])
    k=len(df.columns)
    model_Two=np.zeros((n,k-1))
    
    def gs_func(x,y):
        return pi*(binom.pmf(y,x,pe1)) + (1-pi)*(binom.pmf(y,x,pe2))
    
    for x in range(1,len(df[count])):
        for y in range(len(df.columns)-1):
            if x>=y:
                model_Two[x-1,y] = gs_func(x,y) * sum(data[x-1,])

    return model_Two


def residual_conv(params, df, data, count):
    return np.concatenate((input_conv(params, df, data, count) - (data)), axis=None)

# Run one-state model with brute method
params = Parameters()
params.add('Error_1', min=0, max=0.05, brute_step=.0001)

fitterTC = Minimizer(residual_conv, params, fcn_args=(dfTC, dataTC, 'Tcount'))
if useConv != 'TC':
    fitterCT = Minimizer(residual_conv, params, fcn_args=(dfCT, dataCT, 'Ccount'))

result_bruteTC = fitterTC.minimize(method='brute', keep=25)
if useConv != 'TC':
    result_bruteCT = fitterCT.minimize(method='brute', keep=25)

#display(result_bruteTC)
#display(result_bruteCT)
err1startTC = result_bruteTC.brute_x0
if useConv != 'TC':
    err1startCT = result_bruteCT.brute_x0

print('TC start value:' + str(err1startTC)) # value
print('TC chi-square:' + str(result_bruteTC.brute_fval)) # chi-square
if useConv != 'TC':
    print('CT start value:' + str(err1startCT)) # value
    print('CT chi-square:' + str(result_bruteCT.brute_fval)) # chi-square

# Run one-state model
paramsTC = Parameters()
paramsTC.add('Error_1', value=err1startTC, min=0, max=1)
if useConv != 'TC':
    paramsCT = Parameters()
    paramsCT.add('Error_1', value=err1startCT, min=0, max=1)

out1stateTC = minimize(residual_conv, paramsTC, args=(dfTC, dataTC, 'Tcount'), method='leastsquare')
err1_1state_TC = out1stateTC.params['Error_1'].value
print('TC 1state Error_1 value:' + str(err1_1state_TC)) # value

if useConv != 'TC':
    out1stateCT = minimize(residual_conv, paramsCT, args=(dfCT, dataCT, 'Ccount'), method='leastsquare')
    err1_1state_CT = out1stateCT.params['Error_1'].value
    print('CT 1state Error_1 value:' + str(err1_1state_CT)) # value


def residual_conv2(params, df, data, count): #params, dfTC, dataTC, 'Tcount'):
    return np.concatenate((Two_conv_background(params, df, data, count) - (data)), axis=None)

# Run two-state model
params = Parameters()
params.add('Error_1', value=0.001, min=0, max=1)
params.add('Error_2', value= 0.0001, min=0, max=1)
params.add('frac_err', value=0.001, min=0, max=1)

def residual_conv2(params, df, data, count): #params, dfTC, dataTC, 'Tcount'):
    return np.concatenate((Two_conv_background(params, df, data, count) - (data)), axis=None)

out2stateTC = minimize(residual_conv2, params, args=(dfTC, dataTC, 'Tcount'), method='leastsquare')
if useConv != 'TC':
    out2stateCT = minimize(residual_conv2, params, args=(dfCT, dataCT, 'Ccount'), method='leastsquare')


# Store background error rates
err1_2state_TC = out2stateTC.params['Error_1'].value
err2_2state_TC = out2stateTC.params['Error_2'].value
frac_2state_TC = out2stateTC.params['frac_err'].value
print('TC 2state Error values:' + str(err1_2state_TC) + str(err2_2state_TC) + str(frac_2state_TC)) # value

if useConv != 'TC':
    err1_2state_CT = out2stateCT.params['Error_1'].value
    err2_2state_CT = out2stateCT.params['Error_2'].value
    frac_2state_CT = out2stateCT.params['frac_err'].value

    print('CT 2state Error values:' + str(err1_2state_CT) + str(err2_2state_CT) + str(frac_2state_CT)) # value

def _aic1(n,k,rss):
    """first equation on page e7 of McShane et al
       k : number of parameters in model"""
    aic = 2.0*k + n*np.log(rss/n)
    return aic

#+ 2.0*k*(k+1)/(n-k-1) (As I was reading AIC under the Framework of Least Squares Estimation, it suggests to use this 
#correction parameter when n/K < 40 , so Iw asnt sure.Also suggests different AIC models for OLS optimization )


def calc_aic(n, rss1, rss2,k1,k2, conv):
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
    
    print((aic_min - aic1), (aic_min - aic2))
    print('Prob1state Prob2state:')
    print(prob_aic1, prob_aic2, '\n')
    print(conv)

    if conv == 'TC' and prob_aic1 > prob_aic2:
        err1TC = err1_1state_TC
        err2TC = err1_1state_TC
        fracerr2TC = 0.5
    elif conv == 'TC' and prob_aic2 > prob_aic1:
        err1TC = err1_2state_TC
        err2TC = err2_2state_TC
        fracerr2TC = frac_2state_TC
    elif conv == 'CT' and prob_aic1 > prob_aic2:
        err1CT = err1_1state_CT
        err2CT = err1_1state_CT
        fracerr2CT = 0.5
    elif conv == 'CT' and prob_aic2 > prob_aic1:
        err1CT = err1_2state_CT
        err2CT = err2_2state_CT
        fracerr2CT = frac_2state_CT
    else:
        raise NotImplementedError()
        
    if conv == 'TC':
        return err1TC, err2TC, fracerr2TC
    elif conv == 'CT':
        return err1CT, err2CT, fracerr2CT


print('TC info')
err1TC, err2TC, fracerr2TC = calc_aic(out1stateTC.ndata, out1stateTC.chisqr,out2stateTC.chisqr,out1stateTC.nvarys,out2stateTC.nvarys, 'TC')
if useConv != 'TC':
    print('CT info')
    err1CT, err2CT, fracerr2CT = calc_aic(out1stateCT.ndata, out1stateCT.chisqr,out2stateCT.chisqr,out1stateCT.nvarys,out2stateCT.nvarys, 'CT')

print('TC:', err1TC, err2TC, fracerr2TC)
if useConv != 'TC':
    print('CT:', err1CT, err2CT, fracerr2CT)



df_RawTC = pd.read_csv(
    "../MMfrequency/"
    + Exp
    + "_"
    + MapMethod
    + "_frag_"
    + rateRegion
    + Tmatrix
    + Exp
    + "_"
    + LibName
    + "_otherMMfilt_"
    + stringentFilter
    + ".txt",
    sep="\t",
    header=(0),
)
df_RawTC = df_RawTC.rename(columns={'Unnamed: 0':'Tcount'})

if useConv != 'TC':
    df_RawCT = pd.read_csv(
        "../MMfrequency/"
        + Exp
        + "_"
        + MapMethod
        + "_frag_"
        + rateRegion
        + "_CcountANDCTconv_"
        + Exp
        + "_"
        + LibName
        + "_otherMMfilt_"
        + stringentFilter
        + ".txt",
        sep="\t",
        header=(0),
    )
    df_RawCT = df_RawCT.rename(columns={'Unnamed: 0':'Ccount'})




n=len(df_RawTC['Tcount'])
k=len(df_RawTC.columns)
RawTC=np.zeros((n,k-1))

for x in range(0,len(df_RawTC['Tcount'])):
      for y in range(0,len(df_RawTC.columns)):
        RawTC[x,y-1]= df_RawTC.iloc[x,y]
print('TC matrix length:')
print(len(RawTC))

if useConv != 'TC':
    n=len(df_RawCT['Ccount'])
    k=len(df_RawCT.columns)
    RawCT=np.zeros((n,k-1))

    for x in range(0,len(df_RawCT['Ccount'])):
          for y in range(0,len(df_RawCT.columns)):
            RawCT[x,y-1]= df_RawCT.iloc[x,y]
    print('CT matrix length:')
    print(len(RawCT))



# Function for getting TC conversion rate on sample of interest (input background rates from above)
def Two_conv(params, conv, err1, err2, fracerr, df_Raw, count, Raw):
    pc = params[conv]
    pe1 = err1 # Fill this in based on Error_1 and Error_2 above	
    pe2 = err2 # These three input values all refer to the background mismatch rate and you 
    pi_err = fracerr	# get them from the first part of the script using the 0 min matrix. If model 1 is correct, use error1 for pe1 and pe2 and pi_err (frac) doesn't matter
    pi = params['frac']    

    n=len(df_Raw[count])
    k=len(df_Raw.columns)
    model_Two=np.zeros((n,k-1))
    
    def gs_func(x,y):
        return pi*(binom.pmf(y,x,pc)) + (1-pi)*(pi_err*binom.pmf(y,x,pe1)+(1-pi_err)*binom.pmf(y,x,pe2))
    
    for x in range(1,len(df_Raw[count])):
        for y in range(len(df_Raw.columns)-1):
            if x>=y:
                model_Two[x-1,y] = gs_func(x,y) * sum(Raw[x-1,])

    return model_Two


paramsTC = Parameters()
paramsTC.add('TC', value=0.01, min=0, max=1)
paramsTC.add('frac', value=0.3, min=0, max=1)
if useConv != 'TC':
    paramsCT = Parameters()
    paramsCT.add('CT', value=0.01, min=0, max=1)
    paramsCT.add('frac', value=0.3, min=0, max=1)


def Two_residual_conv(params, conv, err1, err2, fracerr, df_Raw, count, Raw):
    return np.concatenate((Two_conv(params, conv, err1, err2, fracerr, df_Raw, count, Raw) - (Raw)), axis=None)


out2stateTC = minimize(Two_residual_conv, paramsTC, args=('TC', err1TC, err2TC, fracerr2TC, df_RawTC, 'Tcount', RawTC), method='leastsquare')
if useConv != 'TC':
    out2stateCT = minimize(Two_residual_conv, paramsCT, args=('CT', err1CT, err2CT, fracerr2CT, df_RawCT, 'Ccount', RawCT), method='leastsquare')


# Store TC and frac new 
TCrate = out2stateTC.params['TC'].value
TCfracnew = out2stateTC.params['frac'].value

if useConv != 'TC':
    CTrate = out2stateCT.params['CT'].value
    CTfracnew = out2stateCT.params['frac'].value

print('TCrate:', TCrate, 'fracNew:', TCfracnew)
if useConv != 'TC':
    print('CTrate:', CTrate, 'fracNew:', CTfracnew)

RateFile = open(Exp + "_" + LibName + "_" + MapMethod + "_ratesFrom_" + rateRegion + ".txt", "w")

if useConv == 'TC':
    RateFile.write("Rate" + "\t" + Exp + "_" + LibName + "\n" + "TC" + "\t" + str(TCrate) + "\n" + "TCfracnew" + "\t" + str(TCfracnew) )

if useConv == 'both':
    RateFile.write("Rate" + "\t" + Exp + "_" + LibName + "\n" + "TC" + "\t" + str(TCrate) + "\n" + "TC" + "\t" + str(CTrate) + "\n" + "CTfracnew" + "\t" + str(CTfracnew))

###################################################################


###################################################################
################## Now run through on each sample #####################
################## using background rates from above ##################
###################################################################
if RatesOnly != 'TRUE':
    # Create a file to write to
    if useConv == 'both':
        OutputFile = open(Exp + "_" + LibName + "_" + MapMethod + "_ratesFrom_" + rateRegion + "_FracNew_TCandCT_"+ dist+"dist.txt", "w")
    elif useConv == 'TC':
        OutputFile = open(Exp + "_" + LibName + "_" + MapMethod + "_ratesFrom_" + rateRegion + "_FracNew_TConly_"+ dist+"dist.txt", "w")
    elif useConv == 'CT':
        OutputFile = open(Exp + "_" + LibName + "_" + MapMethod + "_ratesFrom_" + rateRegion + "_FracNew_CTonly_"+ dist+"dist.txt", "w")


    regions = [

        "ND1",
        "ND2",
        "CO1",
        "CO2",
        "ATP8_6",
        "CO3",
        "ND3",
        "ND4L_4",
        "ND5",
        "CYB",
        "ND6",
    ]


    OutputFile.write("region\t" + Exp + "_" + LibName + "\n")


    for region in regions:
        print(region)
        # TC matrix
        df_RawTC = pd.read_csv(
            "../MMfrequency/"
            + Exp
            + "_"
            + MapMethod
            + "_frag_"
            + region
            + Tmatrix
            + Exp
            + "_"
            + LibName
            + "_otherMMfilt_"
            + stringentFilter
            + ".txt",
            sep="\t",
            header=(0),
        )
        df_RawTC = df_RawTC.rename(columns={"Unnamed: 0": "Tcount"})
        # CT matrix
        if useConv != 'TC':
            df_RawCT = pd.read_csv(
                "../MMfrequency/"
                + Exp
                + "_"
                + MapMethod
                + "_frag_"
                + region
                + "_CcountANDCTconv_"
                + Exp
                + "_"
                + LibName
                + "_otherMMfilt_"
                + stringentFilter
                + ".txt",
                sep="\t",
                header=(0),
            )
            df_RawCT = df_RawCT.rename(columns={"Unnamed: 0": "Ccount"})

        # Make numpy arrays from the data frames
        n = len(df_RawTC["Tcount"])
        k = len(df_RawTC.columns)
        RawTC = np.zeros((n, k - 1))

        for x in range(0, len(df_RawTC["Tcount"])):
            for y in range(0, len(df_RawTC.columns)):
                RawTC[x, y - 1] = df_RawTC.iloc[x, y]

        if useConv != 'TC':
            n = len(df_RawCT["Ccount"])
            k = len(df_RawCT.columns)
            RawCT = np.zeros((n, k - 1))

            for x in range(0, len(df_RawCT["Ccount"])):
                for y in range(0, len(df_RawCT.columns)):
                    RawCT[x, y - 1] = df_RawCT.iloc[x, y]

        # Function for getting fraction new from input TC conversion rate. Run Binomial_Conversion_Rate_Estimation-v1 interactively (in jupyter notebook)first to get inputs including TC rate.
        def Two_conv(params, conv, err1, err2, fracerr, df_Raw, count, Raw):
            pc = params[conv]
            pe1 = err1  # Fill this in based on Error_1 and Error_2 above
            pe2 = err2  # These three input values all refer to the background mismatch rate and you
            pi_err = fracerr  # get them from the first part of the script using the 0 min matrix. If model 1 is correct, use error1 for pe1 and pe2 and pi_err (frac) doesn't matter
            pi = params["frac"]

            n = len(df_Raw[count])
            k = len(df_Raw.columns)
            model_Two = np.zeros((n, k - 1))

            def gs_func(x, y):
                return pi * (binom.pmf(y, x, pc)) + (1 - pi) * (
                    pi_err * binom.pmf(y, x, pe1) + (1 - pi_err) * binom.pmf(y, x, pe2)
                )

            for x in range(1, len(df_Raw[count])):
                for y in range(len(df_Raw.columns) - 1):
                    if x >= y:
                        model_Two[x - 1, y] = gs_func(x, y) * sum(
                            Raw[
                                x - 1,
                            ]
                        )

            return model_Two

        # 	paramsTC = Parameters()
        # 	paramsTC.add('TC', value=TCrate, min=TCrate, max=TCrate+0.00000001)
        # 	paramsTC.add('frac', value=0.4, min=0, max=1)
        # 	paramsCT = Parameters()
        # 	paramsCT.add('CT', value=CTrate, min=CTrate, max=CTrate+0.00000001)
        # 	paramsCT.add('frac', value=0.4, min=0, max=1)


        # Function to get difference between model and raw data
        if useConv == 'both':
            def Two_residual_conv(
                params,
                err1TC,
                err2TC,
                fracerrTC,
                df_RawTC,
                RawTC,
                err1CT,
                err2CT,
                fracerrCT,
                df_RawCT,
                RawCT,
            ):

                res1 = Two_conv(
                    params, "TC", err1TC, err2TC, fracerrTC, df_RawTC, "Tcount", RawTC
                )
                res2 = Two_conv(
                    params, "CT", err1CT, err2CT, fracerrCT, df_RawCT, "Ccount", RawCT
                )

                return np.concatenate((res1 - RawTC, res2 - RawCT), axis=None)
        else:
            def Two_residual_conv(params, conv, err1, err2, fracerr, df_Raw, count, Raw):
                return np.concatenate((Two_conv(params, conv, err1, err2, fracerr, df_Raw, count, Raw) - (Raw)), axis=None)


        # Run function to get fraction new
    
        if useConv == 'both':
            params = Parameters()
            params.add("TC", value=TCrate, min=TCrate, max=TCrate + 0.00000001)
            params.add("CT", value=CTrate, min=CTrate, max=CTrate + 0.00000001)
            params.add("frac", value=0.4, min=0, max=1)
            out2 = minimize(
                Two_residual_conv,
                params,
                args=(
                    err1TC,
                    err2TC,
                    fracerr2TC,
                    df_RawTC,
                    RawTC,
                    err1CT,
                    err2CT,
                    fracerr2CT,
                    df_RawCT,
                    RawCT,
                ),
                method="leastsquare",
            )
        elif useConv == 'TC':
            params = Parameters()
            params.add("TC", value=TCrate, min=TCrate, max=TCrate + 0.00000001)
            params.add("frac", value=0.4, min=0, max=1)
            out2 = minimize(
                Two_residual_conv,
                params,
                args=(
                    'TC',
                    err1TC,
                    err2TC,
                    fracerr2TC,
                    df_RawTC,
                    'Tcount',
                    RawTC,
                ),
                method="leastsquare",
            )
        elif useConv == 'CT':
            params = Parameters()
            params.add("CT", value=CTrate, min=CTrate, max=CTrate + 0.00000001)
            params.add("frac", value=0.4, min=0, max=1)
            out2 = minimize(
                Two_residual_conv,
                params,
                args=(
                    'CT',
                    err1CT,
                    err2CT,
                    fracerr2CT,
                    df_RawCT,
                    'Ccount',
                    RawCT,
                ),
                method="leastsquare",
            )

        # Get fraction new out of lmfit minimizer object
        print(fit_report(out2))
        fnew = out2.params.valuesdict()["frac"]
        # Write to file
        OutputFile.write(region + "\t" + str(fnew) + "\n")
    
    
####### Combine files after running for all libs
# Combine samples
# lib=("0m_tot_A" "30m_tot_A" "60m_tot_A" "0m_IP_A" "30m_IP_A" "60m_IP_A" "0m_tot_B" "30m_tot_B" "60m_tot_B" "0m_IP_B" "30m_IP_B" "60m_IP_B")
# Exp='NHC1'
# MapMethod='MTnorRNA_t5MTMMinformed6_withDups'
# rateRegion='Bin1_2_3_4' # Bin1_2_3_4 MTnorRNA
# 
# 
# use="TConly" # TConly CTonly TCandCT
# paste <(awk '{print $1"\t"$2}' ${Exp}_${lib[0]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[1]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[2]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[3]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[4]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[5]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[6]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[7]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[8]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[9]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[10]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) <(awk '{print $2}' ${Exp}_${lib[11]}_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt) > ${Exp}_all_${MapMethod}_ratesFrom_${rateRegion}_FracNew_${use}.txt

