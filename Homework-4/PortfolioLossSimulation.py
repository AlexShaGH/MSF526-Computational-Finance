# -*- coding: utf-8 -*-

# Portfolio Loss Simulation
# MSF 526
# Illinois Institute of Technology
# Homework 4, Question 2
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu


__author__ = "oshashkov"



Evaluate the 99% quantile of the daily portfolio loss distribution for
an equally weighted portfolio with n correlated risky assets using MC
simulation. The correlation of asset returns is given by Rho ij = 0.2; i != j.
Each asset daily return is N(mu; sigma**2) where mu = -0.1 and sigma = 0.1. You
should implement the following function

def portfolioLossSimulation(mu, Sigma, w, alpha, samples):

where
'mu' is the n-vector of daily asset returns means
'Sigma' is the covariance matrix of daily asset returns
'samples' is a two-dimensional M x n array of random uniform
(U[0, 1)) samples.

The function should return the quantile estimate. Compare your MC
estimate with the exact quantile and plot the convergence of the MC
simulation to the exact quantile as the number of simulations M is
increased.



if __name__ == '__main__':