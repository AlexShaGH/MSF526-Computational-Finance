# -*- coding: utf-8 -*-

from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randn, rand
from BSMonteCarlo import BSMonteCarlo
from BSMonteCarlo import InterpolateRateCurve
from MCStockPrices import MCStockPrices

rate_curve = np.array(
    [0.08,0.08,0.10,0.11,0.12,0.13,0.16,0.28,0.47,0.69,1.23,1.46])

# Problem 1:
# S0=100.0
# K=110.0
# T=2.5
# sigma=0.4    
  
# checkpoints = 10**r_[2:9]
# M = checkpoints[-1]
# samples = rand(M,1)
# samples = randn(M,1)
# samples = None

# h=hist(samples,100)
# xlabel('z')
# ylabel('count')
# title('normal samples')

# matlab: [Call, Put] = blsprice(100,110,0.145,2.5,0.4)
# Call = 35.4805
# Put = 12.0333

# bsformula results for Call option:
# Price = 35.4805
# Delta = 0.7700
# Vega = 48.0165
# print(BSMonteCarlo(S0, K, T, sigma, checkpoints, rate_curve, samples))
# print(InterpolateRateCurve(rate_curve,T))

# Problem 2:
S0=100.0
sigma=0.4

T=2.5
 
time_steps = 100    
fixing_times = np.array((T/time_steps)*r_[1:time_steps+1])
N = len(fixing_times)

M = number_of_paths =20
samples = randn(N,M)
# h=hist(samples,100)
# xlabel('z')
# ylabel('count')
# title('normal samples')
# print(samples)

integrator = 'standard'

sim_stock_prices = MCStockPrices(S0, sigma, rate_curve, fixing_times, samples, integrator)
print(sim_stock_prices)
print(shape(sim_stock_prices))

