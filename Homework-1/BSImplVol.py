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
    xvals=[]
    fdiffs=[]
    impvol = math.nan
    numberofcalls = 0
    start_sigma = 0.35
    
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