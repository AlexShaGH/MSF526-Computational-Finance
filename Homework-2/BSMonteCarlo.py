# -*- coding: utf-8 -*-

# BSMonteCarlo.py - The function to price vanilla european-exercise options
# MSF 526
# Illinois Institute of Technology
# Homework 2
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import math
import numpy as np
from numpy import exp,sqrt,maximum,mean,std,log,cumsum
from numpy.random import randn,rand

def BSMonteCarlo(S0, K, T, sigma, checkpoints, rateCurve, samples=None):
    """ Calculates value of vanilla european-exercise options

    Parameters
    ----------
    S0 : float
        current price of underlying asset
    K : float
        option strike price
    T : float
        time to expiration expressed in years (1 month = 1/12 -> T = 0.08(3))
    sigma : float
        volatility expressed as a fraction
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
    
    '''    
    #
    # Monte Carlo valuation of European call option
    # in Black-Scholes-Merton model
    # bsm_mcs_euro.py
    #
    # Python for Finance, 2nd ed.
    # (c) Dr. Yves J. Hilpisch
    #
    
    
    # Parameter Values
    S0 = 100.  # initial index level
    K = 105.  # strike price
    T = 1.0  # time-to-maturity
    r = 0.05  # riskless short rate
    sigma = 0.2  # volatility
    
    I = 100000  # number of simulations
    
    # Valuation Algorithm
    z = np.random.standard_normal(I)  # pseudo-random numbers
    # index values at maturity
    ST = S0 * np.exp((r - 0.5 * sigma ** 2) * T + sigma * math.sqrt(T) * z)
    hT = np.maximum(ST - K, 0)  # payoff at maturity
    C0 = math.exp(-r * T) * np.mean(hT)  # Monte Carlo estimator
    
    # Result Output
    print('Value of the European call option %5.3f.' % C0)  '''
    
    #TODO: check parameters values
    
    #TODO: check for samples and generate them if needed
    
    #TODO: find out the value of "r" using rate curve
    

    
    
    #TODO: make all the calculations
    running means = []
    running_stds = []
    running_total_x = []
    running_total_x2 = []
    
    
    #TODO: return all the calculations
    '''return { 'TV': , # The final value ( i.e. mean at checkpoints[-1] )
            'Means': , # The running mean at each checkpoint
            'StdDevs': , # The running standard deviation at each checkpoint
            'StdErrs': , # The running standard error at each checkpoint
            }'''

    return { 'TV': , # The final value ( i.e. mean at checkpoints[-1] )
            'Means': , # The running mean at each checkpoint
            'StdDevs': , # The running standard deviation at each checkpoint
            'StdErrs': , # The running standard error at each checkpoint
            }
    

K=110.0;
S0=100.0
sigma=0.4
T=2.5
M=1000000
samples=randn(M,1)
checkpoints = [100,200,300,400, M]

rateCurve = np.array([0.08,0.08,0.10,0.11,0.12,0.13,0.16,0.28,0.47,0.69,1.23,1.46])

