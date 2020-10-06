# -*- coding: utf-8 -*-
#
# Monte Carlo valuation of European call options with NumPy (log version)
# mcs_full_vector_numpy.py
#
#import math
import numpy
from numpy import exp,sqrt,maximum,mean,std,cumsum, prod
from numpy.random import randn, rand
from time import time
import matplotlib.pyplot as plt

numpy.random.seed(20000)
t0=time()
S0 = 100.; T = 2.5; r = 0.08; sigma = 0.2
M = 50; dt = T / M; I = 2

samples = numpy.random.standard_normal((M, I))
samples = randn(M,2)

S=S0*exp(cumsum((r-0.5*sigma**2)*dt+sigma*sqrt(dt)*samples , axis=0))

print(samples)

plt.plot(S[:, :10])
plt.grid(True)
plt.xlabel('time step')
plt.ylabel('index level')
plt.show()

samples = randn(7, 10)
t = numpy.linspace(0.001,30,numpy.shape(samples)[1])
print(numpy.shape(samples))