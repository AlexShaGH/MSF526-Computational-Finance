# -*- coding: utf-8 -*-

# BinomialAmerican.py - binomial tree function for pricing American
# options in the Black-Scholes framework
# MSF 526
# Illinois Institute of Technology
# Homework 3
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import numpy as np
from numpy import exp,maximum,sqrt
import matplotlib.pyplot as plt

def binomialAmerican(S0, Ks, r, T, sigma, q, callputs, M):
    """ Calculates prices of American options in the Black-Scholes framework
    using binomial tree method
    
    Parameters
    ----------
    S0 : float
        the stock prices at time t0 (spot)
    Ks : float or array of floats
        strike price(s) - might be array or scalar, must be positive
    r : float
        risk-free interest rate expressed as a fraction
        (i.e. r = 0.1 stands for 10%)
    T : float
        the expiration date of the option
    sigma : float
        the constant volatility  
    q : float
        Yield, continuous return rate on the underlying expressed as a fraction
        (i.e. q = 0.05 stands for 5%)
    callputs : int or array of ints
        1=call;-1=put, otherwise will raise ValueError - array or scalar
    M : int
        the number of timesteps

    Raises
    ------
    ValueError
        when passed parameters are invalid

    Returns
    -------
    prices : array of floats
        array of calculated option prices

    """
    
    # check parameters
    if np.isnan(S0):
        raise  ValueError("Underlying price can not be NaN")
    if np.isnan(sigma):
        raise  ValueError("Volatility can not be NaN")
    if np.isnan(M):
        raise  ValueError("M can not be NaN")        
    if np.isnan(T):
        raise  ValueError("Time to expiration can not be NaN")

    # test if Ks and callputs are scalars and if so - convert them into vector
    if np.isscalar(Ks):
        Ks = [Ks]
    if np.isscalar(callputs):
        callputs = [callputs]        
        
    # check the values of input parameters
    if S0 <= 0:
        raise ValueError("Underlying price can not be zero or negative")
    if sigma <= 0:
        raise ValueError("Volatility can not be zero or negative")
    if T <= 0:
        raise ValueError("Expiration time can not be zero or negative")
    if M < 1:
        raise ValueError("M has to be 1 or more (at least 2 timesteps)")
    if np.shape(Ks) != np.shape(callputs):
        raise ValueError("Ks and callputs are of incompatible shapes")
    if not all(abs(x) == 1 for (x) in np.array(callputs)):
        raise ValueError("Invalid value in callputs, must be 1 or -1")
    if not all(y > 0 for (y) in np.array(Ks)):
        raise ValueError("Strike prices can not be zero or negative")
        
    # setup parameters
    dt = T/M #timestep
    u = exp(sigma*sqrt(dt))
    d = 1/u
    qu = (exp((r-q)*dt)-d)/(u-d)
    qd = 1-qu
    
    df = exp(-(r-q)*dt) # discount factor accounting for dividents
    
    N = M + 1 # number of leafs in recombining tree with M steps
    
    # create initial tree to hold stock prices, S0 is a root node
    stock_prices = [np.array([S0])] # [[S0],[Su,Sd],[Suu,Sud,Sdd],...]
    
    chain_lenght = len(callputs)
    
    # calculate stock prices over each step
    for i in range(M):
        prev_branches = stock_prices[-1]
        st = np.concatenate((prev_branches*u,[prev_branches[-1]*d]))
        stock_prices.append(st)
    #print(stock_prices)
    
    # traverse tree to compute spot price of the options
    
    payoffs = np.array([[0]*chain_lenght]*N,dtype='float32')
    #print(payoffs)
    print(stock_prices[-1])
    print(Ks)
    
    # payofs at terminal nodes
    for i in range(chain_lenght):
        if callputs[i] == 1:
            payoffs[:,i] =np.maximum(0,stock_prices[M]-Ks[i])
            #print(np.maximum(0,stock_prices[M]-Ks[i]))
        else:
            payoffs[:,i] = np.maximum(0,Ks[i]-stock_prices[M])
            #print(np.maximum(0,Ks[i]-stock_prices[M]))

    print(payoffs)
                    
    
    #
    
    #
    
    
        
    prices = []

    return prices



if __name__ == '__main__':
    # Problem 1 test:
    #For example, your code should return a 3-element array of prices at t0 if I call
    S0 = 40
    Ks = [30,35,40,45]
    callputs = [-1,1,-1,-1]    
    r = 0.05
    T = 1.2
    sigma =  0.15
    q = 0.01

    M = 10
    
    option_prices = binomialAmerican(S0, Ks, r, T, sigma, q, callputs, M)
    print(option_prices)
