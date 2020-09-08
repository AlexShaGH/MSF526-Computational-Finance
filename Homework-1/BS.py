# -*- coding: utf-8 -*-

# BS.py - Black - Scholes value, delta, vega calculations
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import math
from math import log, sqrt, exp
from scipy.stats import norm

def bsformula(callput, S0, K, r, T, sigma, q=0.):
    """ Calculates option value, delta and vega using extended Black Scholes 
    formula for the assets with continuously compounded dividend yield
    
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

    Returns
    -------
    optionValue : float
        value of the option for a underlying asset
    delta : float
        Option's delta - rate of change of the option value with respect to
        changes in the underlying asset's price
    vega : float
        Option's vega  - rate of change of the option value with respect to
        the change in volatility of the underlying asset
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
    if math.isnan(sigma):
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
    if sigma <= 0:
        raise ValueError("Volatility can not be zero or negative")

    d1 = (log(S0/K) + (r - q + 0.5*sigma**2)*T)/(sigma*sqrt(T))
    d2 = d1 - sigma*sqrt(T)

    optionValue = callput*S0*exp(-q*T)*norm.cdf(callput*d1) \
                  - callput*K*exp(-r*T)*norm.cdf(callput*d2)

    delta = callput*exp(-q*T)*norm.cdf(callput*d1)
    
    vega = S0*exp(-q*T)*sqrt(T)*norm.pdf(d1)
    
    return optionValue, delta, vega