# -*- coding: utf-8 -*-

# BSImplVol.py
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import math

from BS import bsformula

from Bisect import bisect
from Newton import newton

def bsimpvol(callput, S0, K, r, T, price, q=0.,
             priceTolerance=0.01, method='bisect' , reportCalls=False):
    """
    Parameters
    ----------
    callput : integer
        option type: 1 = call, -1 = put
    S0 : float
        current price of underlying asset
    K : float
        option strike price
    r : float
        risk-free interest rate expressed as a fraction
        (i.e. r = 0.1 stands for 10%)
    T : float
        time to expiration expressed in years (1 month = 1/12 -> T = 0.08(3))
    sigma : float
        volatility expressed as a fraction
        (i.e. sigma = 0.5 stands for volatility of 50%)
    q : float, optional, the default is 0
        Yield, continuous return rate on the underlying expressed as a fraction
        (i.e. q = 0.05 stands for 5%)    
    price : float
        observed price of the option
    priceTolerance : float, optional
        The default is 0.01.
    method : method's name, optional
        indicates which numerical method to use for finding a root
        of the function. The default is 'bisect'.
    reportCalls : boolean, optional
        if set to true, then the function must return 2-tuple consisting of
        the volatility found and the number of times it made a call
        to the function bsformula. The default is False

    Returns
    -------
    impvol : float
        volatility found or NaN if it can not be found
    numberofcalls : integer
        the number of calls to bsformula if roportCalls is set to True
    """
    
    # check input parameters for NaN
    if math.isnan(callput):
        raise ValueError("callput argument can not be NaN")
    if math.isnan(K):
        raise ValueError("Strike price can not be NaN")
    if math.isnan(S0):
        raise  ValueError("Underlying price can not be NaN")
    if math.isnan(r):
        raise  ValueError("Interest rate can not be NaN")
    if math.isnan(T):
        raise  ValueError("Time to expiration can not be NaN")
    if math.isnan(price):
        raise  ValueError("Volatility can not be NaN")
    if math.isnan(q):
        raise ValueError("Yield rate can not be NaN")
    
    # check values of input parameters
    if abs(callput) != 1:
        raise ValueError("Unable to parse option type {0}".format(callput))
    if S0 <= 0:
        raise ValueError("Underlying price can not be zero or negative")
    if K <= 0:
        raise ValueError("Strike price can not be zero or negative")
    if T <= 0:
        raise ValueError("Expiration time can not be zero or negative")
    if price <= 0:
        raise ValueError("Price can not be zero or negative")    
    
    
    xvals=[]
    fdiffs=[]
    impvol = math.nan
    numberofcalls = 0
    start_sigma = 1 # we need starting point
    
    # we need to find sigma that would correspond to minimal difference
    # between value of the option calculated with BS formula
    # and real price of the option on the market
    f = lambda sigma: bsformula(callput,S0,K,r,T,sigma,q)[0] - price
    df = lambda sigma: bsformula(callput,S0,K,r,T,sigma,q)[2]    
    
    if method == 'bisect':
        xvals, fdiffs = bisect(0,f,start_sigma,tols=[0.001,priceTolerance])
    elif method == 'newton':
        xvals, fdiffs = newton(0,f,df,start_sigma,tols=[0.001,priceTolerance])
    else:
        raise ValueError("unknown method")

    impvol = xvals[-1]
    numberofcalls = len(xvals)

    if reportCalls:
        return impvol, numberofcalls
    else:
        return impvol