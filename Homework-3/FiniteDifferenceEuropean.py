# -*- coding: utf-8 -*-

# FiniteDifferenceEuropean.py - Finite dierence scheme for pricing European 
# options in the Black-Scholes framework.
# MSF 526
# Illinois Institute of Technology
# Homework 3
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

def fdEuropean(callput, S0, K, r, T, sigma, q, M, N, S_max):

    dS = Smax / float(M) # M+1 spatial values
    dt = T / float(N) #N+1 time values
    i_values = np.arange(M)
    j_values = np.arange(N)
    

    a = 0.5*dt*((sigma**2) *(i_values**2) -
                              r*i_values)
    b = 1 - dt*((sigma**2) * (i_values**2) + r)
    c = 0.5*dt*((sigma**2) *(i_values**2) +
                              r*i_values)
    
    # create datastructures
    grid = np.zeros(shape=(M+1, N+1)) #M+1 by N+1 (full) grid 
    boundary_conds = np.linspace(0, Smax, M+1) # M+1 payoff values (time boundary condition: last columns)
    
    # call payoffs at final value at spatial B.C.s
    if callput:
            grid[:, -1] = np.maximum(
                boundary_conds - K, 0)
            grid[-1, :-1] = np.maximum(Smax - self.K, 0) * \
                                 np.exp(-r *
                                        self.dt *
                                        (self.N-self.j_values))
        else:
            grid[:, -1] = \
                np.maximum(K-boundary_conds, 0)
            grid[0, :-1] = self.K * \
                               np.exp(-self.r *
                                      self.dt *
                                      (self.N-self.j_values))
    
    for j in reversed(self.j_values):
            for i in range(self.M)[1:]:
              self.grid[i,j] = self.a[i]*self.grid[i-1,j+1] +\
                                 self.b[i]*self.grid[i,j+1] + \
                                 self.c[i]*self.grid[i+1,j+1]
    
    
    return interp(S0, boundary_conds, grid[:, 0])


def fdEuropean(callput, S0, K, r, T, sigma, q, M, N, S_max):
    dS = Smax / float(M) # M+1 spatial values
    dt = T / float(N) #N+1 time values
    
    i_values = np.arange(M)
    j_values = np.arange(N)
    

    a = 0.5*dt*((sigma**2) *(i_values**2) -
                              r*i_values)
    b = 1 - dt*((sigma**2) * (i_values**2) + r)
    c = 0.5*dt*((sigma**2) *(i_values**2) +
                              r*i_values)    
    
	matval = np.zeros((M+1,N+1))
	vetS = np.linspace(0,S_max,M+1)
	
    veti = np.array(range(M+1))
    vetj = np.array(range(N+1))
    
	if callput == 1:
		matval[:,N] = np.maximum(vetS-K,0)
		matval[0,:] = 0
		matval[M,:] = S_max*np.exp(-q*dt*(N-vetj)) - K*np.exp(-r*dt*(N-vetj))
	else:
		matval[:,N] = np.maximum(K-vetS,0)
		matval[0,:] = K*np.exp(-r*dt*(N-vetj))
		matval[M,:] = 0
	a = 0.5*dt*(sigma**2*veti - (r-q))*veti
	b = 1 - dt*(sigma**2*veti**2 + (r-q))
	c = 0.5*dt*(sigma**2*veti + (r-q))*veti
	for j in range(N-1,-1,-1):
		for i in range(1,M):
			matval[i,j] = a[i]*matval[i-1,j+1] + b[i]*matval[i,j+1]+ \
			c[i]*matval[i+1,j+1]
	return np.exp(-q*T)*np.interp(S0, vetS, matval[:,0])


if __name__ == '__main__':
    # Problem 1 test:
    pass