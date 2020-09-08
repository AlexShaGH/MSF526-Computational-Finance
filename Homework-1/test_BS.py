# -*- coding: utf-8 -*-

# test_BS.py - Unit tests for bsformula
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"
import math
from BS import bsformula

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