# -*- coding: utf-8 -*-

# test_BS.py - Unit tests for bsformula, Newton and Bisect methods
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import math
from BS import bsformula
from Newton import newton
from Bisect import bisect
from BSImplVol import bsimpvol

print('***** Test case 1: Checking calculated results against matlab\n')

S0 = 100
K = 89
r = 0.05
T = 0.5
sigma = .5
q = 0.025

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}\n'.format(S0,K,r,T,sigma,q))

#mathlab results for the same initial values
callValue_ml = 19.8433
putValue_ml = 7.8881
callDelta_ml = 0.6972
putDelta_ml = -0.2903
vega_ml = 24.0568

print('Expected results:\nCall Price = {0}, Put Price = {1}\n\
Call Delta = {2}, Put Delta = {3}\n\
vega = {4}\n'.format(callValue_ml,putValue_ml,callDelta_ml,putDelta_ml,\
vega_ml))

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}\n'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}\n'.format(error))

print('***** Test case 2: Checking calculated results against matlab\n\
dividents = 0')

S0 = 90
K = 89
r = 0.03
T = 0.7
sigma = .3
q = 0

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}\n'.format(S0,K,r,T,sigma,q))

#mathlab results for the same initial values
callValue_ml = 10.3431
putValue_ml = 7.4936
callDelta_ml = 0.6001
putDelta_ml = -0.3999
vega_ml = 29.0889

print('Expected results:\nCall Price = {0}, Put Price = {1}\n\
Call Delta = {2}, Put Delta = {3}\n\
vega = {4}\n'.format(callValue_ml,putValue_ml,callDelta_ml,putDelta_ml,\
vega_ml))

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}\n'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}\n'.format(error))
    
print('***** Test case 3: Checking calculated results against matlab\n\
risk free rate = 0')

S0 = 70
K = 75
r = 0
T = 0.3
sigma = .5
q = 0

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}\n'.format(S0,K,r,T,sigma,q))

#mathlab results for the same initial values
callValue_ml = 5.6439
putValue_ml = 10.6439
callDelta_ml = 0.4542
putDelta_ml = -0.5458
vega_ml = 15.1949

print('Expected results:\nCall Price = {0}, Put Price = {1}\n\
Call Delta = {2}, Put Delta = {3}\n\
vega = {4}\n'.format(callValue_ml,putValue_ml,callDelta_ml,putDelta_ml,\
vega_ml))

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}\n'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}\n'.format(error))
    
print('***** Test case 4: Handling incorrect input parameters')
S0 = math.nan
K = 89
r = 0.05
T = 0.5
sigma = .5
q = 0.025

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q))
print('Expected results: Gracefully handled exception for both call and put')

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('call: bsformula returned ValueError: {0}'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('put: bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = math.nan
r = 0.05
T = 0.5
sigma = .5
q = 0.025

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q))
print('Expected results: Gracefully handled exception for both call and put')

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('call: bsformula returned ValueError: {0}'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('put: bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = 89
r = math.nan
T = 0.5
sigma = .5
q = 0.025

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q))
print('Expected results: Gracefully handled exception for both call and put')

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('call: bsformula returned ValueError: {0}'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('put: bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = 89
r = 0.05
T = math.nan
sigma = .5
q = 0.025

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q))
print('Expected results: Gracefully handled exception for both call and put')

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('call: bsformula returned ValueError: {0}'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('put: bsformula returned ValueError: {0}'.format(error))     
print('')

S0 = 100
K = 89
r = 0.05
T = 0.5
sigma = math.nan
q = 0.025

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q))
print('Expected results: Gracefully handled exception for both call and put')

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('call: bsformula returned ValueError: {0}'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('put: bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = 89
r = 0.05
T = 0.5
sigma = 0.5
q = math.nan

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q))
print('Expected results: Gracefully handled exception for both call and put')

callput = 1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('call: bsformula returned ValueError: {0}'.format(error))
callput = -1    
try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Put option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('put: bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = 89
r = 0.05
T = 0.5
sigma = 0.5
q = 0
callput = 0

print('Initial values:\n\
callput={6} S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,
sigma,q,callput))
print('Expected results: Gracefully handled exception')

try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}'.format(error))
print('')

S0 = -1
K = 89
r = 0.05
T = 0.5
sigma = 0.5
q = 0
callput = 1

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q))
print('Expected results: Gracefully handled exception')

try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = -1
r = 0.05
T = 0.5
sigma = 0.5
q = 0
callput = 1

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q,))
print('Expected results: Gracefully handled exception')

try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = 89
r = 0.05
T = -0.5
sigma = 0.5
q = 0
callput = 1

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q,))
print('Expected results: Gracefully handled exception')

try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}'.format(error))
print('')

S0 = 100
K = 89
r = 0.05
T = 0.5
sigma = -0.1
q = 0
callput = 1

print('Initial values:\n\
S0={0}, K={1}, r={2}, T={3}, sigma={4}, q={5}'.format(S0,K,r,T,sigma,q,))
print('Expected results: Gracefully handled exception')

try:
    optionValue, delta, vega = bsformula(callput,S0,K,r,T,sigma,q)
    print('bsformula results for Call option:\nPrice = {0:.4f}\n\
Delta = {1:.4f}\nVega = {2:.4f}\n'.format(optionValue,delta,vega))
except ValueError as error:    
    print('bsformula returned ValueError: {0}'.format(error))
    
    
    
target = 0
y = lambda x: x**3 + 2*x**2 - 5
dy = lambda x: 3*x**2 + 4*x
start = 5
tols = [0.00001,0.010]
maxiter = 1000
bounds = [0,3]

print('***** Test case 5: Root finding using')
print("********** Newton's Method *************")
print('Initial values:\n\
y=x**3+2*x**2-5\n\
dy=3*x**2+4*x\n\
target={0}, start={1}\n\
tolerances={2}, bounds={3}\n\
maxiter={4}\n'.format(target,start,tols,bounds,maxiter))
xvals, fdiffs = newton(target,y,dy,start,tols=tols,maxiter=maxiter)
print(xvals, fdiffs)
print('root = {0}'.format(xvals[-1]))
print('error = {0}'.format(fdiffs[-1]))
print('n = {0}\n'.format(len(xvals)))

print('***** Test case 6: Root finding using')
print("********** Bisect Method *************")
print('Initial values:\n\
y=x**3+2*x**2-5\n\
target={0}, start={1}\n\
tolerances={2}, bounds={3}\n\
maxiter={4}\n'.format(target,start,tols,bounds,maxiter))
xvals, fdiffs = bisect(target,y,start,tols=tols,maxiter=maxiter)
print(xvals, fdiffs)
print('root = {0}'.format(xvals[-1]))
print('error = {0}'.format(fdiffs[-1]))
print('n = {0}\n'.format(len(xvals)))    

callput = 1
S0 = 110
K = 100
r = 0.05
T = 0.5
q = 0.025
OptionPrice = 14

print('***** Test case 7: Finding Implied volatility')
print('Initial values:\n\
callput = {0}, S0 = {1}, K = {2}\n\
r = {3}, T = {4}, q = {5}\n\
OptionPrice = {6}\n'.format(callput,S0,K,r,T,q,OptionPrice))

print("********** Newton's Method *************")
try:
    impVol,callNum =  bsimpvol(callput,S0,K,r,T,OptionPrice,q,method='newton',reportCalls=True)
    print('Implied volatility calculation using Newton method:')
    print('implied volatility = {0},\nnumber of calls = {1}\n'.format(impVol,callNum))
except ValueError as error:
    print('bsimpvol returned ValueError: {0}\n'.format(error))

print("********** Bisect Method *************")
try:
    impVol,callNum =  bsimpvol(callput,S0,K,r,T,OptionPrice,q,method='bisect',reportCalls=True)
    print('Implied volatility calculation using Bisect method:')
    print('implied volatility = {0},\nnumber of calls = {1}\n'.format(impVol,callNum))
except ValueError as error:
    print('bsimpvol returned ValueError: {0}\n'.format(error))

callput = -1
print('Initial values:\n\
callput = {0}, S0 = {1}, K = {2}\n\
r = {3}, T = {4}, q = {5}\n\
OptionPrice = {6}\n'.format(callput,S0,K,r,T,q,OptionPrice))

print("********** Newton's Method *************")
try:
    impVol,callNum =  bsimpvol(callput,S0,K,r,T,OptionPrice,q,method='newton',reportCalls=True)
    print('Implied volatility calculation using Newton method:')
    print('implied volatility = {0},\nnumber of calls = {1}\n'.format(impVol,callNum))
except ValueError as error:
    print('bsimpvol returned ValueError: {0}\n'.format(error))

print("********** Bisect Method *************")
try:
    impVol,callNum =  bsimpvol(callput,S0,K,r,T,OptionPrice,q,method='bisect',reportCalls=True)
    print('Implied volatility calculation using Bisect method:')
    print('implied volatility = {0},\nnumber of calls = {1}\n'.format(impVol,callNum))
except ValueError as error:
    print('bsimpvol returned ValueError: {0}\n'.format(error))

