# -*- coding: utf-8 -*-

# FiniteDifferenceEuropean.py - Finite difference scheme (explicit) for pricing
# European options in the Black-Scholes framework
# MSF 526
# Illinois Institute of Technology
# Homework 3
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu
# Code fragments from MSF526 lecture notes and notebooks used as well as 
# materails attached to the "Mastering Python for Finance, Second Edition" by 
# James Ma Weiming

__author__ = "oshashkov"

import numpy as np
from numpy import exp,maximum,linspace, interp

def fdEuropean(callput, S0, K, r, T, sigma, q, M, N, S_max):
    """ Calculates prices of European options in the Black-Scholes framework
    using explicit finite difference schema

    Parameters
    ----------
    callput : int
        1=call;-1=put, otherwise will raise ValueError
    S0 : float
        Underlying price
    K : float
        Strike price of the option
    r : float
        risk-free interest rate expressed as a fraction
        (i.e. r = 0.1 stands for 10%)
    T : float
        the expiration date of the option
        defines the pitch of the grid for the time axis
    sigma : float
        the constant volatility  
    q : float
        Yield, continuous return rate on the underlying expressed as a fraction
        (i.e. q = 0.05 stands for 5%)        
    M : int
        defines the pitch of the grid for the price axis
    N : int
        defines the pitch of the grid for the time axis
    S_max : float
        unrelaisticialy large price that security is unlikely to reach
        defines the pitch of the grid for the price axis

    Raises
    ------
    ValueError
        when passed parameters are invalid

    Returns
    -------
    float
        calculated option value
    """
    
    # check parameters
    if not abs(callput) == 1:
        raise ValueError("Invalid value of callput, must be 1 or -1")    
    if np.isnan(S0):
        raise  ValueError("Underlying price can not be NaN")
    if np.isnan(K):
        raise  ValueError("Strike price can not be NaN")        
    if np.isnan(sigma):
        raise  ValueError("Volatility can not be NaN")
    if np.isnan(M):
        raise  ValueError("M can not be NaN")
    if np.isnan(N):
        raise  ValueError("N can not be NaN")        
    if np.isnan(T):
        raise  ValueError("Time to expiration can not be NaN")

    # check the values of input parameters
    if S0 <= 0:
        raise ValueError("Underlying price can not be zero or negative")
    if S_max <= 0:
        raise ValueError("S_max can not be zero or negative")
    if K <= 0:
        raise ValueError("Strike price can not be zero or negative")        
    if sigma <= 0:
        raise ValueError("Volatility can not be zero or negative")
    if T <= 0:
        raise ValueError("Expiration time can not be zero or negative")
    if M < 1:
        raise ValueError("M has to be 1 or more")
    if N < 1:
        raise ValueError("N has to be 1 or more")        

    dt = T / float(N) # N + 1 time values

    # create datastructures
    i_values = np.arange(M)
    j_values = np.arange(N)
    grid = np.zeros((M + 1, N + 1))
    boundary_conds = linspace(0, S_max, M + 1)
    
    # setup boundary conditions
    if callput == 1:
        grid[:, -1] = maximum(0,boundary_conds - K)
        grid[-1, :-1] = (S_max-K)*exp(-(r-q)*dt*(N-j_values))
    else:
        grid[:, -1] = maximum(0,K-boundary_conds)
        grid[0, :-1] = (K-S_max)*exp(-(r-q)*dt*(N-j_values))
        
    # setup coefficients
    a = 0.5*dt*i_values*((sigma**2)*i_values - r + q)
    b = 1 - dt*((sigma**2)*(i_values**2) + r)
    c = 0.5*dt*i_values*((sigma**2)*i_values + r - q)
    #print('a = \n{0}'.format(a))
    #print('b = \n{0}'.format(b))
    #print('c = \n{0}'.format(c))
        
    # traverse grid
    for j in reversed(j_values):
        for i in range(M)[2:]:
            grid[i,j]=a[i]*grid[i-1,j+1]+b[i]*grid[i,j+1]+c[i]*grid[i+1,j+1]

    #print('grid = \n{0}'.format(grid))
    
    return interp(S0, boundary_conds, grid[:, 0])            


if __name__ == '__main__':
    # Problem 2 test:
    S0 = 40
    K = 30
    r = 0.05
    T = 1
    sigma =  0.15
    q = 0.01
    M = 100
    N = 1000
    S_max = 100
    
    print('Call option:')
    callput = 1
    print('Parameters:\nS0 = {0}, K = {1}, r = {2}, T = {3},\n\
Sigma = {4}, q = {5}, M = {6}, N = {7}, S_max = {8}\n'.format(S0,K,r,T,sigma,q,M,N, S_max))
    option = fdEuropean(callput, S0, K, r, T, sigma, q, M, N, S_max)
    print('Calculated option value = {0}\n'.format(option))
    
    M = 100
    N = 100    
    callput = 1    
    print('Call option:')    
    print('Parameters:\nS0 = {0}, K = {1}, r = {2}, T = {3},\n\
Sigma = {4}, q = {5}, M = {6}, N = {7}, S_max = {8}\n'.format(S0,K,r,T,sigma,q,M,N, S_max))
    option = fdEuropean(callput, S0, K, r, T, sigma, q, M, N, S_max)
    print('Calculated option value = {0}\n'.format(option))    

    S0 = 45
    K = 40
    r = 0.05
    T = 1
    sigma =  0.4
    q = 0.01
    M = 100
    N = 1000
    S_max = 100

    print('Put option:')
    callput = -1
    print('Parameters:\nS0 = {0}, K = {1}, r = {2}, T = {3},\n\
Sigma = {4}, q = {5}, M = {6}, N = {7}, S_max = {8}\n'.format(S0,K,r,T,sigma,q,M,N, S_max))
    option = fdEuropean(callput, S0, K, r, T, sigma, q, M, N, S_max)
    print('Calculated option value = {0}\n'.format(option))     
    
    M = 30
    N = 50     
    print('Put option:')
    callput = -1
    print('Parameters:\nS0 = {0}, K = {1}, r = {2}, T = {3},\n\
Sigma = {4}, q = {5}, M = {6}, N = {7}, S_max = {8}\n'.format(S0,K,r,T,sigma,q,M,N, S_max))
    option = fdEuropean(callput, S0, K, r, T, sigma, q, M, N, S_max)
    print('Calculated option value = {0}\n'.format(option))      
    