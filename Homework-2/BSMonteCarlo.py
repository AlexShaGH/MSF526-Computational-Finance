# -*- coding: utf-8 -*-

# BSMonteCarlo.py - The function to price vanilla european-exercise options
# MSF 526
# Illinois Institute of Technology
# Homework 2
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import numpy as np
from numpy import exp,sqrt,maximum,mean,std
from numpy.random import randn
from scipy.stats import sem

# these are standard tenors for US Treasury Yield Curve as of 10/1/2020
RATE_CURVE_TENORS = [1/12,2/12,3/12,6/12,1,2,3,5,7,10,20,30]

def InterpolateRateCurve(curve,T,tenors=None):
    """ Interpolates the value of Yield Rate based on the provided rate curve
   
    Parameters
    ----------
    curve : array of floats
         represents rate curve. The length of curve and tenors must be the same
    T : float
        Time to interpolate rate for expressed in years
    tenors : array of floats, optional
        Tenors for Yield Rate curve . The default is None.

    Returns
    -------
    rate : float
        interpolated value of rate
    """
    if tenors is None:
        tenors = RATE_CURVE_TENORS #assume it's constant
    if len(tenors) != len(curve):
        return None
    return np.interp(T,tenors,curve)


def BSMonteCarlo(S0, K, T, sigma, checkpoints, rateCurve, samples=None):
    """ Estimates value of vanilla european-exercise options usinf Monte Carlo
    method in Black-Scholes-Merton model
    
    ### Call option for the purpose of programming assignment ###

    Parameters
    ----------
    S0 : float
        current price of underlying asset
    K : float
        option strike price
    T : float
        time to expiration expressed in years (1 month = 1/12 -> T = 0.08(3))
    sigma : float
        volatility of underlying expressed as a fraction
        (i.e. sigma = 0.5 stands for volatility of 50%)
    checkpoints : ordered list
        is an ordered list of integer sample counts at which to return
        the running mean, standard deviation, and estimated error
    rateCurve : numpy array
        is an InterestRateCurve stored as a numpy array
    samples : numpy array of floats, optional
        is a numpy array of uniform random samples to use. 
        The default is None

    Returns
    -------
    dict { 'TV': , # The final value ( i.e. mean at checkpoints[-1] )
            'Means': , # The running mean at each checkpoint
            'StdDevs': , # The running standard deviation at each checkpoint
            'StdErrs': , # The running standard error at each checkpoint
            }
    """
    
    #TODO: check parameters values
    
    M = checkpoints[-1]
    
    # check for samples and generate them if needed
    if samples is None:
        samples = randn(M,1)
    elif len(samples) < M:
        raise ValueError('Not enough samples: {0}'.format(len(samples)))
    
    # find the value of "r" for given T using rate curve
    r = InterpolateRateCurve(rateCurve,T)
    
    running_means = []
    running_stds = []
    running_st_errs = []
    vals = []
    
    # need these totals for efficient calculation of running standar deviation
    # see Professor Dixon's evening lecture from 09/29/2020 and 
    # https://stackoverflow.com/questions/1174984/how-to-efficiently-calculate-a-running-standard-deviation
    running_total_x = 0
    running_total_x2 = 0
    
    start = 0
    
    for i in range(len(checkpoints)):
        end = checkpoints[i]
        price_samples = S0*exp((r-0.5*sigma**2)*T + sigma*sqrt(T)*samples[start:end])
        running_total_x = running_total_x + sum(price_samples)
        running_total_x2 = running_total_x2 + sum(price_samples**2)
        running_means.append(running_total_x/end)
        running_stds.append(sqrt(running_total_x2/end - running_means[-1]**2))
        running_st_errs.append(running_stds[-1]/sqrt(end))
        vals.extend( exp(-r*T) * maximum(0,price_samples-K))
        start = end
        
    #vals = exp(-r*T) * maximum(0,price_samples-K)
    V0=mean(vals)
    #D=std(vals)
    #error_est=D/sqrt(M)
    #V0 = exp(-r*T) * maximum(0,running_means[-1]-K)
    
    
    # return all the calculations
    return { 'TV': V0, # The final value ( i.e. mean at checkpoints[-1] )
            'Means':  running_means,# The running mean at each checkpoint
            'StdDevs':  running_stds,# The running standard deviation at each checkpoint
            'StdErrs':  running_st_errs# The running standard error at each checkpoint
            }
    

K=110.0;
S0=100.0
sigma=0.4
T=2.5
M=1000000
samples=randn(M-1,1)
checkpoints = [100,500,1000,10000,100000,M]

rate_curve = np.array([0.08,0.08,0.10,0.11,0.12,0.13,0.16,0.28,0.47,0.69,1.23,1.46])

print(BSMonteCarlo(S0, K, T, sigma, checkpoints, rate_curve))
print(InterpolateRateCurve(rate_curve,T))

# matlab: [Call, Put] = blsprice(100,110,0.145,2.5,0.4)
# Call = 35.4805
# Put = 12.0333

# bsformula results for Call option:
# Price = 35.4805
# Delta = 0.7700
# Vega = 48.0165


