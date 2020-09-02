# -*- coding: utf-8 -*-

# BS.py - Black - Scholes value, delta, vega calculations
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

from math import log, sqrt, exp
from scipy.stats import norm

def bsformula(callput, S0, K, r, T, sigma, q=0.):
    """ Calculates option value, delta and vega using extended Black Scholes 
    formula for the assets with continuously compounded dividend yield
    https://www.macroption.com/black-scholes-formula/

    Parameters
    ----------
    callput : integer
        option type: 1 = call, -1 = put
    S0 : float
        current price of underlying asset
    K : float
        option strike price
    r : float
        risk-free interest rate
    T : float
        time to expiration expressed in years (1 month = 1/12 = 0.08(3))
    sigma : float
        volatility
    q : float, optional
        continuous return rate on the underlying. The default is 0

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
    # TODO: check imput parameters
    
    call = 1
    put = -1
    
    d1 = (log(S0/K) + (r - q + (sigma**2)/2)*T)/(sigma*sqrt(T))
    d2 = d1 - sigma*sqrt(T)
    if callput == call:
        optionValue = S0*exp(-q*T)*norm.cdf(d1) - K*exp(-r*T)*norm.cdf(d2)
        delta = exp(-q*T)*norm.cdf(d1)
        vega = 0
    elif callput == put:
        optionValue =  K*exp(-r*T)*norm.cdf(-d2) - S0*exp(-q*T)*norm.cdf(-d1)
        delta = exp(-q*T)*(norm.cdf(d1) - 1)
        vega = 0        
    else:
        raise ValueError("Unable to parse option type {0}".format(callput))
    return (optionValue, delta, vega)


print(bsformula(1,50,49,0.02,0.25,0.33,0.02))
