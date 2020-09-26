# -*- coding: utf-8 -*-

# BSMonteCarlo.py - The function to price vanilla european-exercise options
# MSF 526
# Illinois Institute of Technology
# Homework 2
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"


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
    return { 'TV': , # The final value ( i.e. mean at checkpoints[-1] )
            'Means': , # The running mean at each checkpoint
            'StdDevs': , # The running standard deviation at each checkpoint
            'StdErrs': , # The running standard error at each checkpoint
            }