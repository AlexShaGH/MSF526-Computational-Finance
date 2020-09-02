# -*- coding: utf-8 -*-

# test_BS.py - Unit tests for bsformula
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

from BS import bsformula

callput = 1
S0 = 100
K = 89
r = 0.05
T = 0.5
sigma = .5
q = 0.05

try:
    print(bsformula(callput,S0,K,r,T,sigma,q))
except ValueError as error:    
    print('bsformula returned ValueError: {0}'.format(error))

