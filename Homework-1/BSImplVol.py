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

from Bisect import newton, bisect

def bsimpvol( callput, S0, K, r, T, price, q=0.,
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
    price : TYPE
        DESCRIPTION.
    priceTolerance : TYPE, optional
        DESCRIPTION. The default is 0.01.
    method : TYPE, optional
        DESCRIPTION. The default is 'bisect'.
    reportCalls : TYPE, optional
        DESCRIPTION. The default is False.        

    Returns
    -------
    None.

    """
    xvals=[]
    fdiffs=[]
    
    
    
    if method == 'bisect':
        xvals, fdiffs = bisect(target,y,start,tols=tols,maxiter=maxiter)
        
    else if method == 'newton':
        xvals, fdiffs = newton(target,y,dy,start,tols=tols,maxiter=maxiter)
        
    else:
        raise ValueError("unknown method")
        
#unfinished
