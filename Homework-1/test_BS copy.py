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
K = 110
r = 0.14500000000000002
T = 2.5
sigma = 0.4
q = 0


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