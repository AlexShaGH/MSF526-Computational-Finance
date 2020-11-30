# -*- coding: utf-8 -*-

# Short and long rates modeling
# MSF 526
# Illinois Institute of Technology
# Homework 4, Question 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu


__author__ = "oshashkov"



def BrennanSchwartzPaths(r0, theta0, a1,a2, b1, b2, b3, sigma1,sigma2, rho12, fixingTimes, samples):

that uses the usual Euler "step forward" method of forming Brennan Schwartz model short and long rate paths, where
'r0' is the time t value of the short rate
'theta0' is the time t value of the long rate
'sigma1' is the constant volatility of the short rate
'sigma2' is the constant volatility of the long rate
'a1' is the constant drift of the short rate
'a2' is the mean reversion rate of the short rate
'b1,b2,b3' are constant drift coecients of the long rate
'rho12' is the correlation between the wiener processes for the long and short rates
'fixingTimes' is the size N list of times relevant to the path
'samples' is a two-dimensional MxN array of random uniform (U[0; 1)) samples.
This function should return a dictionary with the following two keys:
'paths r' an M x N array consisting of M short rate paths.
'paths theta' an M x N array consisting of M long rate paths.



if __name__ == '__main__':