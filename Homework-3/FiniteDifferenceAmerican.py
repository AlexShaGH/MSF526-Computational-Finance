# -*- coding: utf-8 -*-

# FiniteDifferenceAmerican.py - Finite dierence scheme for pricing American 
# options in the Black-Scholes framework.
# MSF 526
# Illinois Institute of Technology
# Homework 3
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import numpy as np
import sys
from numpy import exp,maximum,linspace, interp

def fdAmerican(callput, S0, K, r, T, sigma, q, M, N, S_max):
    """ Calculates prices of American options in the Black-Scholes framework
    using Crank-Nicolson formulation of the finite difference scheme

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
    #dS = S_max / float(M) # not used

    # create datastructures
    omega = 1.2 # ask professor
    tol = 0.001 # ask professor
    i_values = np.arange(M+1)
    j_values = np.arange(N+1)    
    
    boundary_conds = linspace(0, S_max, M + 1)    
    
    # setup boundary conditions
    if callput == 1:
        payoffs = maximum(0, boundary_conds[1:M]-K)
    else:
        payoffs = maximum(0, K-boundary_conds[1:M])
        
    past_values = payoffs
    boundary_values = K * exp(-r*dt*(N-j_values))# r or r-q?
    
    # setup coefficients and matrixes
    alpha = 0.25*dt*((sigma**2)*(i_values**2) - (r-q)*i_values)
    beta = -dt*0.5*((sigma**2)*(i_values**2) + r)# r or r-q?
    gamma = 0.25*dt*((sigma**2)*(i_values**2) + (r-q)*i_values)
   
    # M1 = -np.diag(alpha[2:M], -1) + \
    #               np.diag(1-beta[1:M]) - \
    #               np.diag(gamma[1:M-1], 1)
    
    M2 = np.diag(alpha[2:M], -1) + \
                  np.diag(1+beta[1:M]) + \
                  np.diag(gamma[1:M-1], 1)    
    
    # traverse grid
    aux = np.zeros(M-1)
    new_values = np.zeros(M-1)

    for j in reversed(range(N)):
        aux[0] = alpha[1]*(boundary_values[j] + boundary_values[j+1])
        rhs = np.dot(M2, past_values) + aux
        old_values = np.copy(past_values)
        error = sys.float_info.max

        while tol < error:
            payoff = old_values[0] + omega/(1-beta[1]) * \
                (rhs[0] - (1-beta[1])*old_values[0] + gamma[1]*old_values[1])
            new_values[0] = max(payoffs[0], payoff)

            for k in range(M-2)[1:]:
                payoff = old_values[k] + omega/(1-beta[k+1]) * (rhs[k] + \
                    alpha[k+1]*new_values[k-1] - (1-beta[k+1])*old_values[k] + \
                    gamma[k+1]*old_values[k+1])
                
                new_values[k] = max(payoffs[k], payoff)
                
            payoff = old_values[-1] + omega/(1-beta[-2]) * (rhs[-1] + \
                alpha[-2]*new_values[-2] - (1-beta[-2])*old_values[-1])
                    
            new_values[-1] = max(payoffs[-1], payoff) 
                
            error = np.linalg.norm(new_values-old_values)
            old_values = np.copy(new_values)

            past_values = np.copy(new_values)

    values = np.concatenate(([boundary_values[0]], new_values, [0]))
    
    # linear interpolation on final values array
    return interp(S0, boundary_conds, values)

if __name__ == '__main__':
    # Problem 3 test:
    callput = 1
    S0 = 40
    K = 30
    r = 0.05
    T = 1
    sigma =  0.15
    q = 0.01
    M = 100
    N = 100
    S_max = 200

    option = fdAmerican(callput, S0, K, r, T, sigma, q, M, N, S_max)
    
    print(option)            
