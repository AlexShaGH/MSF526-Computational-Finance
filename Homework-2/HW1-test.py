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
sigma=0.2

T=2.5
 
time_steps = 50
fixing_times = np.array((T/time_steps)*r_[1:time_steps+1])
N = len(fixing_times)
np.random.seed(20000)
M = number_of_paths = 2
samples = randn(N,M)
# h=hist(samples,100)
# xlabel('z')
# ylabel('count')
# title('normal samples')
# print(samples)



print('Fixing Times:\n{0}\n'.format(fixing_times))
rates = InterpolateRateCurve(rate_curve,fixing_times)
print('Rates for steps:\n{0}\n'.format(rates))

print('Samples:\n{0}\n'.format(samples))
#for index,sample in enumerate(samples):
#    plt.plot(fixing_times,sample,label=index, lw=1.5)
plt.plot(samples) 
plt.grid(True)    
plt.xlabel('time')
plt.ylabel('samples')
plt.legend()
plt.show()

integrator = 'standard'
sim_stock_prices = MCStockPrices(S0, sigma, rate_curve, fixing_times, samples, integrator)
print('Stock prices ({0}):\n{1}\n'.format(integrator,sim_stock_prices))
plt.plot(sim_stock_prices)    
plt.grid(True)    
plt.xlabel('time')
plt.ylabel('stock price')
plt.legend()
plt.show()

integrator = 'euler'
sim_stock_prices = MCStockPrices(S0, sigma, rate_curve, fixing_times, samples, integrator)
print('Stock prices ({0}):\n{1}\n'.format(integrator,sim_stock_prices))
plt.plot(sim_stock_prices[:, :10])     
plt.grid(True)    
plt.xlabel('time')
plt.ylabel('stock price')
plt.legend()
plt.show()

integrator = 'milstein'
sim_stock_prices = MCStockPrices(S0, sigma, rate_curve, fixing_times, samples, integrator)
print('Stock prices ({0}):\n{1}\n'.format(integrator,sim_stock_prices))
plt.plot(sim_stock_prices[:, :10])     
plt.grid(True)    
plt.xlabel('time')
plt.ylabel('stock price')
plt.legend()
plt.show()