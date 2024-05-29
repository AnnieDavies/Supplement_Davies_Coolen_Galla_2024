'''
Code to perform data analysis using the delayed kernel models for the simulation study in Davies, Coolen and Galla (2024)
The same code can be used for scenario 1 and 2 (simply change the data files that are read in)
Written by A. Davies (2023).

The code is run from R using the reticulate package so that data generated in R can be analysed
The delayed kernel models are fitted to the R data using maximum likelihood.
Minimisation of the negative log likelihood is performed using Powell's method via the function scipy.optimize()

Then cross validation is used to estimate prediction error for one prediction window and five base times 

The models fitted here specify a fixed assocation parameter for individuals with s=0 (n=1)
For the decaying parameter model, the relevant changes to the code are commented out and labelled 's=0 Decaying Association Model'

Description of the delayed kernel models can be found in the manuscriupt 'Delayed kernels for longitudinal survival analysis and dynamic predictions' Davies, Coolen and Galla (2024). 
Detailed derivations of the functions evaluated by this code can be found in the corresponding Supplementary Material.
'''

import numpy as np
import os
import math
from numpy import exp, arange, sqrt, pi, log
import scipy.optimize as optimize
import time
import pandas as pd
import pyarrow.feather as feather
import multiprocessing as mp
from multiprocessing import Pool

#############################
#Maximum likelihood functions
#############################

'''
To perform maximum likelihood inference of delayed kernel models A and B we minimise the negative log likelihood 
The functions funcA() and funcB() output the negative log likelihood for models A and B respectively calling to functions mu_sum1(), j_sum() and mu_sum2() (versions A and B)
The equation that correpsonds to this function is Eq (22) in Davies, Coolen and Galla (2024).
We make use of Eqs (S14) and (S16) in the correpsonding supplementary material which give the results of the integration in Eq (22) when using the LOCF interpolation procedure for Models A and B.
'''

## MODEL A =======================================================================================

def musum1A(V, T, Z, Z_fix, X, i, j, p, q, n):
    mu_sum = 0.0
    #work out sum over n[] up to j-1 for use in Z and T indicators
    sum_nj = 0
    for k in range(0,j):
        sum_nj = sum_nj + n[k]
    
    s_j = T[sum_nj + n[j] - 1] #final observation time 
    min_XS = min(s_j, X[i])
    
    #sum over p longitudinal covariates
    for mu in range(0, p):
        tau = p + mu
        if n[j]==1: 
            a = 0
            mu_sum = mu_sum + V[mu] * Z[sum_nj*p + a * p + mu]                  #s=0 Fixed Association Model
            #mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]*exp(-X[i]/V[tau]) 	#s=0 Decaying Association Model
        else: 
            for a in range(0, n[j]-1): 
                if V[tau] == 0: 
                    mu_sum = mu_sum
                elif X[i] <= T[sum_nj + a]: 
                    mu_sum = mu_sum
                else:
                    min_XT = min(X[i], T[sum_nj + a + 1])
                    numer = exp(-(min_XS-min_XT)/V[tau])-exp(-(min_XS-T[sum_nj+a])/V[tau])
                    denom = 1.0 - exp(-min_XS / V[tau])
                    mu_sum = mu_sum + (V[mu] * Z[sum_nj*p + a * p + mu] * numer) / denom
    #sum over q fixed covariates
    for nu in range(0, q):
        mu_sum = mu_sum + V[2 * p + nu] * Z_fix[j*q + nu]
    return mu_sum

def j_sumA(V, T, Z, Z_fix, X, i, N_train, p, q, n):
    numer = 0.0
    for j in range(0, N_train):
        #work out the sum over n[] up to j-1 for use in Z and T indicators
        sum_nj = 0
        for k in range(0,j):
            sum_nj = sum_nj + n[k]
            
        s_j = T[sum_nj + n[j] - 1] #the final observation time = final switch time (when a=n[j]-1)

        if X[j] - X[i] < 0.0 or X[i] < 0.0: 
            numer = numer        
        else:
            e = exp(musum1A(V, T, Z, Z_fix, X, i, j, p, q, n))
            numer = numer + e
    return numer
    
def musum2A(V, T, Z, Z_fix, X, i, p, q, n):
    mu_sum = 0.0
    #work out sum over n[] up to i-1 for use in Z and T indicators
    sum_ni = 0
    for k in range(0,i):
        sum_ni = sum_ni + n[k]
    
    s_i = T[sum_ni + n[i] - 1] #final observation time 
    min_XS = min(s_i, X[i])
    
    #sum over p longitudinal covariates
    for mu in range(0, p):
        tau = p + mu
        if n[i]==1:
            a = 0
            mu_sum = mu_sum + V[mu] * Z[sum_ni*p + a * p + mu]
        else: 
            for a in range(0, n[i]-1): 
                if V[tau] == 0: 
                    mu_sum = mu_sum
                elif X[i] <= T[sum_ni + a]:
                    mu_sum = mu_sum
                else:
                    min_XT = min(X[i], T[sum_ni + a + 1])
                    numer = exp(-(min_XS - min_XT) / V[tau]) - exp(-(min_XS - T[sum_ni + a]) / V[tau])
                    denom = 1.0 - exp(-min_XS / V[tau])
                    mu_sum = mu_sum + (V[mu] * Z[sum_ni*p + a * p + mu] * numer) / denom
    #sum over q fixed covariates
    for nu in range(0, q):
        mu_sum = mu_sum + V[2 * p + nu] * Z_fix[i*q + nu]
    return mu_sum

def funcA(V, T, Z, Z_fix, X, N_train, p, q, n, Censor):
    func_val = 0.0
    if V[1] < 0: #penalise for negative tau
        func_val = func_val + 1000000000000
    else:
        for i in range(0, N_train):
            if Censor[i] == 0: 
                func_val = func_val
            else:
                func_val = func_val + log(j_sumA(V, T, Z, Z_fix, X, i, N_train, p, q, n)) - musum2A(V, T, Z, Z_fix, X, i, p, q, n)
    return func_val

## MODEL B =======================================================================================

def mu_sum1B(V, T, Z, Z_fix, X, i, j, p, q, n):
    mu_sum = 0.0

    #work out the sum over n[] up to j-1 for use in Z and T indicators
    sum_nj = 0
    for k in range(0,j):
        sum_nj = sum_nj + n[k]

    s_j = T[sum_nj + n[j] - 1] #the final observation time
    min_XS = min(s_j, X[i])

    #sum over p longitudinal covariates
    for mu in range(0, p):
        tau = p + mu
        if n[j]==1: 
            a = 0
            mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]                      #s=0 Fixed Association Model
            #mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]*exp(-X[i]/V[tau])   #s=0 Decaying Association Model
        else: 
            #sum over n[j]-1 longitudinal observations
            for a in range(0, n[j]-1): 
                if V[tau] == 0: 
                    mu_sum = mu_sum
                elif X[i] <= T[sum_nj + a]: 
                    mu_sum = mu_sum
                else:
                    min_XT = min(X[i], T[sum_nj + a + 1])
                    #first term
                    mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]*(min_XT - T[sum_nj + a])/min_XS #CHANGE
                    #second term
                    #term 2a:
                    diffa = X[i] - T[sum_nj + a] 
                    diffa1 = X[i] - min_XT
                    term2a = exp(-diffa1/V[tau]) - exp(-diffa/V[tau])
                    #term 2b:
                    diffs = X[i] - min_XS
                    term2b = (min_XT-T[sum_nj + a])*(exp(-diffs/V[tau])-exp(-X[i]/V[tau]))/min_XS #CHANGE
                    #result:
                    mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]*(term2a - term2b)
    #sum over q fixed covariates
    for nu in range(0, q):
        mu_sum = mu_sum + V[2*p+nu]*Z_fix[j*q+nu]
    return mu_sum

def j_sumB(V, T, Z, Z_fix, X, i, N_train, p, q, n):
    j_sum = 0.0
    #sum over j (individuals in the training data)
    for j in range(0, N_train):
        #work out the sum over n[] up to j-1 for use in Z and T indicators
        sum_nj = 0
        for k in range(0,j):
            sum_nj = sum_nj + n[k]
            
        s_j = T[sum_nj + n[j] - 1]  #the final observation time 

        if X[j] - X[i] < 0.0 or X[i] < 0.0:
            j_sum = j_sum
        else:
            e = exp(mu_sum1B(V, T, Z, Z_fix, X, i, j, p, q, n))
            j_sum = j_sum + e
    return j_sum


def mu_sum2B(V, T, Z, Z_fix, X, i, p, q, n):
    mu_sum = 0.0
    #work out the sum over n[] up to i-1 for use in Z and T indicators
    sum_ni = 0
    for k in range(0,i):
        sum_ni = sum_ni + n[k]

    s_i = T[sum_ni + n[i] - 1] #the final observation time = final switch time (when a=n[i]) #CHANGE
    min_XS = min(s_i,X[i]) #should always be s_i

    #sum over p longitudinal covariates
    for mu in range(0, p):
        tau = p + mu
        if n[i]==1: 
            a = 0
            mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu]                      #s=0 Fixed Association Model
            #mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu]*exp(-X[i]/V[tau])   #s=0 Decaying Association Model
        else: 
            #sum over n[i] longitudinal observations
            for a in range(0, n[i]-1):
                if V[tau] == 0: 
                    mu_sum = mu_sum
                elif X[i] <= T[sum_ni + a]: #this shouldn't ever be the case 
                    print("oh no mu sum2")
                    mu_sum = mu_sum
                else:
                    min_XT = min(X[i], T[sum_ni + a + 1])   #should always be T 
                    #first term
                    mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu]*(min_XT - T[sum_ni + a])/min_XS 
                    #second term
                    #term 2a:
                    diffa = X[i] - T[sum_ni + a] 
                    diffa1 = X[i] - min_XT
                    term2a = exp(-diffa1/V[tau]) - exp(-diffa/V[tau])
                    #term 2b:
                    diffs = X[i] - min_XS
                    term2b = (min_XT-T[sum_ni + a])*(exp(-diffs/V[tau])-exp(-X[i]/V[tau]))/min_XS 
                    #Result:
                    mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu]*(term2a - term2b)
    #sum over q fixed covariates
    for nu in range(0, q):
        mu_sum = mu_sum + V[2*p+nu]*Z_fix[i*q+nu]
    return mu_sum

def funcB(V, T, Z, Z_fix, X, N_train, p, q, n, Censor):
    function = 0.0
    #penalise for negative tau (NB: V[1] represents tau)
    if V[1]<=0: 
        function = function + 10e15
    else:
        #sum over i individuals in data
        for i in range(0, N_train):
            if Censor[i] == 0: 
                function = function 
            else:
                function = function + log(j_sumB(V, T, Z, Z_fix, X, i, N_train, p, q, n)) - mu_sum2B(V, T, Z, Z_fix, X, i, p, q, n)
    return function


#######################
# Pred error functions
#######################

'''
Functions to calculate the prediction error for the delayed kernel models 
Function PE_SQ (A and B) returns the prediction error for given base time, prediction time, and training and test data vectors.
PE_SQ (A and B) calls to the function S_t (A and B) which calculates survival probability.
Function S_t (A and B) calls to function EXP (A and B).
Within function PE_SQ, 'new data' vectors are created from the test data vectors restricting observations to be <= t (base time).
The Eq for prediction error is given in Eq (26) of Davies, Coolen and Galla (2024). 
The Eq for survival probability (of the DK model) for a given base time and prediction time is given in Eq (27). 
For the integrals we again we make use of Eqs (S14) and (S16) in the correpsonding Supplementary Material.
'''

## MODEL A =======================================================================================

def EXP_A(ti, kappa, tau, gamma, Z_t, T_t, Zf_t, n_t, j, p, q):
    summ = 0.0
    #work out the sum over n[] up to j-1 for use in Z and T indicators
    sum_nj = 0
    for k in range(0,j):
        sum_nj = sum_nj + n_t[k]
    sum_nj = int(sum_nj)
    s_j = T_t[sum_nj + int(n_t[j]) - 1] #the final obs time
    min_XS = min(s_j, ti)
    
    #sum over p longitudinal covariates
    for mu in range(0,p):
        #sum over n[i]-1 longitudinal observations        
        #print("in for mu")
        if int(n_t[j])==1: 
            a = 0
            summ = summ + kappa[mu] * Z_t[sum_nj*p + a * p + mu]                    #s=0 Fixed Association Model
            #summ = summ + kappa[mu]*Z_t[sum_nj*p + a*p + mu]*exp(-ti/tau[mu])      #s=0 Decaying Association Model
        else: 
            for a in range(0,int(n_t[j])-1): 
                if tau[mu] == 0.0: 
                    summ = summ
                elif ti <= T_t[sum_nj + a]:
                    summ = summ
                else:
                    min_XT = min(ti, T_t[sum_nj + a + 1])
                    numer = exp(-(min_XS - min_XT) / tau[mu]) - exp(-(min_XS - T_t[sum_nj + a]) / tau[mu])#CHANGE
                    denom = 1.0 - exp(-min_XS / tau[mu])
                    summ = summ + (kappa[mu] * Z_t[sum_nj*p + a * p + mu] * numer) / denom
    #sum over q fixed covariates
    for nu in range(0,q):
        summ = summ + gamma[nu] * Zf_t[j*q + nu]
    return exp(summ)

def S_t_A(kappa, tau, gamma, predTime, baseTime, X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train, t_obs, p, q):
    #j refers to the indiviudal in the test data
    #i and k label individuals in the training data
    
    #sum over individuals (event times) in the training data 
    sum_i = 0.0
    for i in range(0,N_train):
        if Censor[i] == 0 or X[i] < baseTime or X[i] > predTime:
            sum_i = sum_i
        else:
            #base hazard - sum over individuals in the training data:
            denom = 0.0
            for k in range(0,N_train):
                sum_nk = 0
                for l in range(0,k):
                    sum_nk = sum_nk + n[l]

                s_k = t_obs[sum_nk + n[k] - 1] #the final observation time 
                if X[i] < 0.0 or X[i] > X[k]:
                    denom = denom
                else:
                    denom = denom + EXP_A(X[i], kappa, tau, gamma, Z, T, Z_fix, n, k, p, q)
            #for individual j (test) we only have data up to baseTime hence 'new' vectors (from test data)
            numer = EXP_A(X[i], kappa, tau, gamma, Z_new, T_new, Z_fix_test, n_new, j, p, q)
            sum_i = sum_i + numer / denom
    S = exp(-sum_i)
    return S

def PE_SQ_A(predTime, baseTime, kappa, tau, gamma, X, X_test, Censor, Censor_test, Z, Z_test, Z_fix, Z_fix_test, t_obs, t_obs_test, n, n_test, N_train, N_test, T, p, q):
    #create Z_new, n_new, t_obs_new and T_new from baseTime and TEST vectors
    #keep only observations <= baseTime
    
    #No. of observations
    n_new = np.zeros(N_test)
    sum_new = 0
    for j in range(0,N_test):
        sum_nj = 0
        for i in range(0,j):
            sum_nj = sum_nj + n_test[i]

        count_n = 0
        for a in range(0,n_test[j]):
            if t_obs_test[sum_nj + a] <= baseTime:
                count_n = count_n + 1
            else:
                count_n = count_n
        n_new[j] = count_n #new n defined 
        sum_new = sum_new + n_new[j]

    sum_new = int(sum_new)
    Z_new = [0.0]*p*sum_new    #longitudinal covariates 
    t_obs_new = [0.0]*sum_new    #observation times 
    
    for j in range(0,N_test):
        new_nj = 0
        sum_nj = 0
        for i in range(0,j):
            new_nj = new_nj + n_new[i]
            sum_nj = sum_nj + n_test[i]
        new_nj = int(new_nj)
        sum_nj = int(sum_nj)

        #t_obs new************************************
        for a in range(0,int(n_new[j])):
            t_obs_new[new_nj + a] = t_obs_test[sum_nj + a]

        #Z_new****************************************
        for mu in range(0,p):
            for a in range(0,int(n_new[j])):
                Z_new[p*new_nj + p * a + mu] = Z_test[p*sum_nj + p * a + mu]
    T_new = t_obs_new 
    
    r = 0.0 #no. of subjects in test data at risk at baseTime
    summ = 0.0
    #sum over individuals in test data still at risk at baseTime (Xj>=baseTime)
    for j in range(0,N_test):
        if X_test[j] < baseTime:
            summ = summ
        else:
            r = r + 1.0
            #Prob of j's survival to predTime given covariates up to baseTime (new vectors) and survival to baseTime
            SprobTau = S_t_A(kappa, tau, gamma, predTime, baseTime, X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train, t_obs, p, q)
            #term 1: still alive at predTime
            if X_test[j] >= predTime:
                summ = summ + pow((1.0 - SprobTau), 2.0)
            #term 2: experienced an event by predTime (and after baseTime)
            elif Censor_test[j] == 1:
                summ = summ + pow(SprobTau, 2.0)
            #term 3: censored by predTime (and after baseTime)
            else:
                #Prob of j's survival to predTime given covariates up to baseTime (new vectors) and survival to X_test[j]
                SprobTj = S_t_A(kappa, tau, gamma, predTime, X_test[j], X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train, t_obs, p, q)
                summ = summ + pow((1.0 - SprobTau), 2.0)*SprobTj + pow(SprobTau, 2.0)*(1.0 - SprobTj)
    M_Y_sq=0.0
    if r>0:
        M_Y_sq = summ / r
    return M_Y_sq

## MODEL B =======================================================================================

def EXP_B(ti, kappa, tau, gamma, Z_t, T_t, Zf_t, n_t, j, p, q):
    summ = 0.0

    #work out the sum over n[] up to j-1 for use in Z and T indicators
    sum_nj = 0
    for k in range(0,j):
        sum_nj = sum_nj + n_t[k]
    sum_nj = int(sum_nj)
    s_j = T_t[sum_nj + int(n_t[j]) - 1] #the final observation time 
    min_XS = min(s_j,ti)

    #sum over p longitudinal covariates
    for mu in range(0,p):
        if int(n_t[j])==1: 
            a = 0
            summ = summ + kappa[mu]*Z_t[sum_nj*p + a*p + mu]                    #s=0 Fixed Association Model
            #summ = summ + kappa[mu]*Z_t[sum_nj*p + a*p + mu]*exp(-ti/tau[mu])  #s=0 Decaying Association Model
        else: 
            #sum over n[j] longitudinal observations
            for a in range(0,int(n_t[j])-1): 
                if ti <= T_t[sum_nj + a]: 
                    summ = summ
                else:
                    min_XT = min(ti, T_t[sum_nj + a + 1]) 
                    #first term
                    summ = summ + kappa[mu]*Z_t[sum_nj*p + a*p + mu]*(min_XT - T_t[sum_nj + a])/min_XS 
                    #second term
                    #term 2a:
                    diffa = ti - T_t[sum_nj + a] 
                    diffa1 = ti - min_XT
                    term2a = exp(-diffa1/tau[mu]) - exp(-diffa/tau[mu])
                    #term 2b:
                    diffs = ti - min_XS; 
                    term2b = (min_XT-T_t[sum_nj + a])*(exp(-diffs/tau[mu])-exp(-ti/tau[mu]))/min_XS 

                    summ = summ + kappa[mu]*Z_t[sum_nj*p + a*p + mu]*(term2a - term2b)
    #sum over q fixed covariates
    for nu in range(0,q):
        summ = summ + gamma[nu]*Zf_t[j*q+nu]
    e = exp(summ) 
    return e


def S_t_B(kappa, tau, gamma, predTime, baseTime, X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train, t_obs, p, q):
    #j refers to the indiviudal in the test data
    #i and k label individuals in the training data

    #sum over individuals (event times) in the training data
    sum_i = 0.0
    for i in range(0,N_train):
        if Censor[i] == 0 or X[i] < baseTime or X[i] > predTime:
            sum_i = sum_i
        else:
            #base hazard - sum over individuals in the training data:
            denom = 0.0
            for k in range(0,N_train):
                sum_nk = 0
                for l in range(0,k):
                    sum_nk = sum_nk + n[l]
                s_k = t_obs[sum_nk + n[k] - 1] #the final observation time = final switch time (when a=n[k]-1)
                if X[i] < 0.0 or X[i] > X[k]:
                    denom = denom
                else:
                    denom = denom + EXP_B(X[i], kappa, tau, gamma, Z, T, Z_fix, n, k, p, q)
            #for individual j (test) we only have data up to baseTime hence 'new' vectors (from test data)
            numer = EXP_B(X[i], kappa, tau, gamma, Z_new, T_new, Z_fix_test, n_new, j, p, q)
            sum_i = sum_i + numer/denom
    S = exp(-sum_i)
    return S


def PE_SQ_B(predTime, baseTime, kappa, tau, gamma, X, X_test, Censor, Censor_test, Z, Z_test, Z_fix, Z_fix_test, t_obs, t_obs_test, n, n_test, N_train, N_test, T, p, q):
    #create Z_new, n_new, t_obs_new and T_new from baseTime and TEST vectors
    #keep only observations <= baseTime
    
    #No. of observations
    n_new = np.zeros(N_test)
    sum_new = 0
    for j in range(0,N_test): 
        sum_nj = 0
        for i in range(0,j):
            sum_nj = sum_nj + n_test[i]

        count_n = 0
        for a in range(0,n_test[j]):
            if t_obs_test[sum_nj + a] <= baseTime:
                count_n = count_n + 1
            else:
                count_n = count_n
        n_new[j] = count_n  #new n defined 
        sum_new = sum_new + n_new[j]

    sum_new = int(sum_new)
    Z_new = [0.0]*p*sum_new             #longitudinal covariates (create vector of zeroes of size p*sum_new)
    t_obs_new = [0.0]*sum_new           #observation times (create vector of zeroes of size sum_new)

    for j in range(0,N_test): 
        new_nj = 0
        sum_nj = 0
        for i in range(0,j):
            new_nj = new_nj + n_new[i]
            sum_nj = sum_nj + n_test[i]
        new_nj = int(new_nj)
        sum_nj = int(sum_nj)
        
        #t_obs new************************************
        for a in range(0,int(n_new[j])):
            t_obs_new[new_nj+a] = t_obs_test[sum_nj+a]

        #Z_new****************************************
        for mu in range(0,p):
            for a in range(0,int(n_new[j])):
                Z_new[p*new_nj + p*a + mu] = Z_test[p*sum_nj + p*a + mu]

        #T_new****************************************
        #define vector of observation times for individual j (length n_new[j])
        t_o = [0]*int(n_new[j])
        for i in range(0,int(n_new[j])):
            t_o[i] = t_obs_test[sum_nj+i]

    T_new = t_obs_new 

    r = 0.0     #no. of subjects in test data at risk at baseTime
    summ = 0.0
    #sum over individuals in test data still at risk at baseTime (Xj>=baseTime)
    for j in range(0,N_test):
        if X_test[j] < baseTime:
            summ = summ 
        else:
            r = r + 1.0
            #Prob of j's survival to predTime given covariates up to baseTime (new vectors) and survival to baseTime
            SprobTau = S_t_B(kappa, tau, gamma, predTime, baseTime, X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train, t_obs, p, q)
            #term 1: still alive at predTime
            if X_test[j] >= predTime:
                summ = summ + pow((1.0-SprobTau), 2.0)
            #term 2: experienced an event by predTime (and after baseTime)
            elif Censor_test[j]==1:
                summ = summ + pow(SprobTau, 2.0)
            #term 3: censored by predTime (and after baseTime)
            else:
                #Prob of j's survival to predTime given covariates up to baseTime (new vectors) and survival to X_test[j]
                SprobTj = S_t_B(kappa, tau, gamma, predTime, X_test[j], X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train, t_obs, p, q)
                summ = summ + pow((1.0-SprobTau),2.0)*SprobTj + pow(SprobTau,2.0)*(1.0-SprobTj)
    M_Y_sq=0.0
    if r>0:
        M_Y_sq = summ / r
    
    return M_Y_sq




######################
# Functions to fit models and perform PE
#####################

def fit_DK(data, data_id, model):
    p = 1
    q = 1

    ## define data
    # Number of individuals
    id_counts = data['id'].value_counts()
    n_obs = [id_counts.get(i, 0) for i in data_id['id']]
    n_obs = [int(ind) for ind in n_obs]
    N = len(n_obs)

    # Fixed variables - use data.id
    X = data_id['Time'].tolist()
    Censor = data_id['event'].tolist()
    Z_fix = data_id['group'].tolist()
    Z_fix = [float(ind) for ind in Z_fix]
    
    # Longitudinal variables - use data
    t_obs = data['time'].tolist()
    Z = data['y'].tolist()

    #MINMISE
    init1 = np.random.uniform(-100,100)
    init2 = np.random.uniform(0,100)
    init3 = np.random.uniform(-100,100)
    initial_guess = [init1, init2, init3]

    
    st = time.time()
    if model=="A":
        result = optimize.minimize(funcA, initial_guess, args = (t_obs, Z, Z_fix, X, N, p, q, n_obs, Censor), method = 'Powell', options = {'xtol':3e-8, 'ftol':3e-8, 'return_all':True, 'disp':True}) 
    elif model=="B":
        result = optimize.minimize(funcB, initial_guess, args = (t_obs, Z, Z_fix, X, N, p, q, n_obs, Censor), method = 'Powell', options = {'xtol':3e-8, 'ftol':3e-8, 'return_all':True, 'disp':True}) 
    else:
        print("Model must be A or B")
    
    kappa = []
    for i in range(0,p):
        kappa.append(result.x[i])

    tau = []
    for i in range(0,p):
        tau.append(result.x[p + i])

    gamma = []
    for i in range(0,q):
        gamma.append(result.x[2 * p + i])
    
    return [kappa[0], tau[0], gamma[0], elap]#return results

def run_PE_DK(test_data, train_data, test_data_id, train_data_id, DKA_params, DKB_params):
    
    # define data variables-------------------------------------
    p = 1
    q = 1

    ## TEST
    id_counts_test = test_data['id'].value_counts()
    n_test = [id_counts_test.get(i, 0) for i in test_data_id['id']]
    n_test = [int(ind) for ind in n_test]
    N_test = len(n_test)
    X_test = test_data_id['Time'].tolist()
    Censor_test = test_data_id['event'].tolist()
    Z_fix_test = test_data_id['group'].tolist()
    Z_fix_test = [float(ind) for ind in Z_fix_test]
    t_obs_test = test_data['time'].tolist()
    Z_test = test_data['y'].tolist()
    
    
    ## TRAIN    
    id_counts_train = train_data['id'].value_counts()
    n_train = [id_counts_train.get(i, 0) for i in train_data_id['id']]
    n_train = [int(ind) for ind in n_train]
    N_train = len(n_train)
    X_train = train_data_id['Time'].tolist()
    Censor_train = train_data_id['event'].tolist()
    Z_fix_train = train_data_id['group'].tolist()
    Z_fix_train = [float(ind) for ind in Z_fix_train]
    t_obs_train = train_data['time'].tolist()
    Z_train = train_data['y'].tolist()
    T_train = t_obs_train
    
    # define parameters (as lists)--------------------------------------------
    kappa_A = [DKA_params[0]]
    tau_A = [DKA_params[1]]
    gamma_A = [DKA_params[2]]
    
    kappa_B = [DKB_params[0]]
    tau_B = [DKB_params[1]]
    gamma_B = [DKB_params[2]]
    
    
    #Calculate prediction error for different base times-----------------------
    PEA = []
    PEB = []
    
    win = 2.0
    for l in range(0,5):
        baseTime = 1.5+(l*2) 
        predTime = baseTime + win
        predErr_A = PE_SQ_A(predTime, baseTime, kappa_A, tau_A, gamma_A, X_train, X_test, Censor_train, Censor_test, Z_train, Z_test, Z_fix_train, Z_fix_test, t_obs_train, t_obs_test, n_train, n_test, N_train, N_test, T_train, p, q)
        predErr_B = PE_SQ_B(predTime, baseTime, kappa_B, tau_B, gamma_B, X_train, X_test, Censor_train, Censor_test, Z_train, Z_test, Z_fix_train, Z_fix_test, t_obs_train, t_obs_test, n_train, n_test, N_train, N_test, T_train, p, q)
        PEA.append(predErr_A)
        PEB.append(predErr_B)
    
    return [PEA, PEB]

######################
# Function to perform main tasks
#####################

def main_func(dat_split):
    dat = dat_split[0] #first element of list is the df
    splits = dat_split[1] #second element is the list of split ids
    
    #remove duplicated 'id'
    dat_id = dat[~dat.duplicated('id')]
    
    ## fit models to full data
    DKA = fit_DK(dat, dat_id, "A")
    DKB = fit_DK(dat, dat_id, "B")
    
    ##create dataframe with results
    cols_res = ["Method", "alpha", "gamma", "tau", "time"]
    rowA = ["DKA", DKA[0], DKA[2], DKA[1], DKA[3]]
    rowB = ["DKB", DKB[0], DKB[2], DKB[1], DKB[3]]
    data_res = [rowA, rowB]
    df_res = pd.DataFrame(data_res, columns=cols_res)
    
    # Create df for group results
    cols_gres = ["Group","kappa", "tau", "gamma", "PE1", "PE2", "PE3", "PE4", "PE5"]
    df_DKA = pd.DataFrame(columns=cols_gres)
    df_DKB = pd.DataFrame(columns=cols_gres)
    
    #Perform cross validation
    for g in range(0,10):
        #define test and training data for each group
        test_id = splits[g]
        test_data = dat[dat['id'].isin(test_id)]
        test_data_id = test_data[~test_data.duplicated('id')]
        train_data = dat[~dat['id'].isin(test_id)]
        train_data_id = train_data[~train_data.duplicated('id')]
        
        #fit models to training data
        DKA_train = fit_DK(train_data, train_data_id, "A")
        DKB_train = fit_DK(train_data, train_data_id, "B")
        
        #perform PE
        PE_list = run_PE_DK(test_data, train_data, test_data_id, train_data_id, DKA_train, DKB_train)
        PEA = PE_list[0]
        PEB = PE_list[1]
        
        rowAg = [g, DKA_train[0], DKA_train[1], DKA_train[2], PEA[0], PEA[1], PEA[2], PEA[3], PEA[4]]
        rowBg = [g, DKB_train[0], DKB_train[1], DKB_train[2], PEB[0], PEB[1], PEB[2], PEB[3], PEB[4]]
        df_DKA.loc[g] = rowAg
        df_DKB.loc[g] = rowBg
        
    # return data frames
    return [df_res, df_DKA, df_DKB]


def process_dat_split(dat_split):
    return main_func(dat_split)


###################
# Run main_func parallel
###################


# Now, let's run the multiprocessing code
if __name__ == '__main__':
    print("In python main")
    
    # Access the list of DataFrames from R using reticulate
    df_list_r = r.dat_split_list_r ## dat_split_list is the list of dataframes & splits in R with simulated data: [[data1, splits1], [data2, splits2],...]
    
    # Extract DataFrames into a list of pandas DataFrames
    df_list_python = [pd.DataFrame(element[0]) for element in df_list_r]

    # Combine DataFrames and associated lists into one list
    dat_split_list_py = [[df, element[1]] for df, element in zip(df_list_python, df_list_r)]
    print("Created python dat_split_list")
    
    # Run the multiprocessing code
    print("Begin multi-processing")
    with Pool() as pool:
        DK_results = pool.map(process_dat_split, dat_split_list_py)
    print(type(DK_results))
