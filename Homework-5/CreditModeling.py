# -*- coding: utf-8 -*-

# CreditModeling.py 
# MSF 526
# Illinois Institute of Technology
# Homework 5
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu
# Code fragments from MSF526 lecture notes and notebooks are used

__author__ = "oshashkov"

import numpy as np
from numpy import exp,sqrt,log,pi,sin,cos,dot,zeros
import sys
import inspyred
from inspyred import ec
from random import Random
from scipy.optimize import fmin_tnc
import scipy.stats as st
norminv = st.distributions.norm.ppf
norm = st.distributions.norm.cdf
import numpy.linalg
import matplotlib.pyplot as plt

class Option:
    def __init__(self, strike, maturity, bid, ask, underlying, typ, vol = 0):
        self.bid        = bid 
        self.ask        = ask 
        self.maturity   = maturity
        self.underlying = underlying
        self.strike     = strike
        self.vol        = vol
        self.typ        = typ
        
    def getMidPrice(self):
        return 0.5*(self.bid + self.ask)        
    
    def dump(self):
        return str(self.bid) +"," + str(self.ask) + "," +  str(self.strike) +"," + str(self.maturity) + "," + str(self.underlying) + "," + str(self.typ) + "," + str(self.vol)

def readTDAmeritradeData(fileName, max_chain_size = 1024):
    global W
    chain = []
    try:
      data = np.genfromtxt(str(fileName),delimiter=',',dtype=None,skip_header=1)	
    except IOError as e:
      print ("I/O error({0}): {1}".format(e.errno, e.strerror))
      raise
    except ValueError:
      print ("Could not convert data to an integer.")
      raise
    except:
      print ("Unexpected error:", sys.exc_info()[0])
      raise
    
    W = []    
 
    for i in range(len(data)):
      #if (i % 100 ==0):
        #print (i)
      tau = (data[i][1])/365.0 # there are approx 260 trading days in a year
      if data[i][11] != 0: # remove zero bids
        o = Option(float(data[i][16]),tau,float(data[i][11]),float(data[i][13]), float(data[i][0]),'C',float(data[i][7]))
        chain.append(o)
        if ((o.ask -o.bid) != 0.0):
          W.append(1.0/(o.ask-o.bid))
        else:
          W.append(1e-3)
      
      if data[i][17] != 0: # remove zero bids  
        o = Option(float(data[i][16]),tau,float(data[i][17]),float(data[i][19]), float(data[i][0]),'P',float(data[i][26]))
        chain.append(o)
        if ((o.ask -o.bid) != 0.0):
          W.append(1.0/(o.ask-o.bid))
        else:
          W.append(1e-3)
    
    W = np.array(W)
    chain = np.array(chain)
    indices = np.argsort(W) 
    
    W = W[indices[:max_chain_size]]/sum(W[indices[:max_chain_size]])
    return chain[indices[:max_chain_size]] 

def BatesCOS(S,K,T,r,sigma,lmbda,meanV,v0,rho,a_prime,b_prime,lmbda_prime,otype,N=256):
    c1 = r*T+(1-exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T
    c2 = 1.0/(8.0*lmbda**3)*(sigma*T*lmbda*exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+
        lmbda*rho*sigma*(1-exp(-lmbda*T))*(16.0*meanV-8.0*v0)+
        2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma**2+4.0*lmbda**2)+
        sigma**2*((meanV-2.0*v0)*exp(-2.0*lmbda*T)+
        meanV*(6.0*exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda**2*(v0-meanV)*(1-exp(-lmbda*T)))
    a = c1-12.0*sqrt(abs(c2))
    b = c1+12.0*sqrt(abs(c2))
    x = log(float(S)/K)
    k = np.arange(0,N)
    if (otype == 'C'):
        U = 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
    else:
        U = 2.0/(b-a)*(-xi(k,a,b,a,0) + psi(k,a,b,a,0))
    unit = [1.0] * N
    unit[0] = 0.5
    ret = 0
    # Note that BatesCF is independent of the strike
    BCF = BatesCF(k*pi/(b-a),T,r,sigma,lmbda,meanV,v0,rho,a_prime,b_prime,lmbda_prime)
    for i in range(N):
        ret += unit[i]*BCF[i]*exp(1j*float(k[i])*pi*(x-a)/(b-a))*U[i]
    return K*exp(-r*T)*ret.real

def HestonCF(u,T,r,sigma,lmbda,meanV,v0,rho):
    a = lmbda*meanV
    b = lmbda
    d = sqrt((1j*rho*sigma*u-b)**2+(u**2+1j*u)*sigma**2)
    g = (b-1j*rho*sigma*u-d)/(b-1j*rho*sigma*u+d)
    ret = exp(1j*u*r*T)
    ret = ret*exp((a/sigma**2)*((b - rho*1j*sigma*u - d)*T - 2.0*log((1-g*exp(-d*T))/(1-g))))
    return ret*exp((v0/sigma**2)*(b - rho*1j*sigma*u - d)*(1-exp(-d*T))/(1-g*exp(-d*T)))

def BatesCF(u,T,r,sigma,lmbda,meanV,v0,rho,a,b,lmbda_prime):
    ret = HestonCF(u,T,r,sigma,lmbda,meanV,v0,rho)# Heston is needed here - do not change it
    ret *= exp(lmbda_prime*T*(-a*u*1j + (exp(u*1j*log(1.0+a)+0.5*b*b*u*1j*(u*1j-1.0))-1.0)))
    return ret

def xi(k,a,b,c,d):
    ret = 1.0/(1+(k*pi/(b-a))**2)*(cos(k*pi*(d-a)/(b-a))*exp(d)-cos(k*pi*(c-a)/(b-a))*exp(c)
        +k*pi/(b-a)*sin(k*pi*(d-a)/(b-a))*exp(d)-k*pi/(b-a)*sin(k*pi*(c-a)/(b-a))*exp(c))
    return ret

def psi(k,a,b,c,d):
    N = len(k)
    idx = np.arange(1, N)
    ret = np.array([0.0]*N)
    ret[0] = d-c
    ret[idx] =(sin(k[idx]*pi*(d-a)/(b-a))-sin(k[idx]*pi*(c-a)/(b-a)))*(b-a)/(k[idx]*pi)
    return ret

def Error_Function(p0):
    global RMSE
    global t_cos
    kappa, theta, sigma, rho, v0, a, b, lmbda = p0
    if bIncorporateNLContraint:
       kappa = (kappa + sigma**2)/(2.0*theta)
    nk = len(chain)
    HOV_e = zeros(nk)
    HG    = zeros(nk)
    j = 0
    for o in chain:
        HOV_e[j] = BatesCOS(o.underlying,o.strike,o.maturity,r0,sigma,kappa,theta,v0,rho,a,b,lmbda, o.typ,nInt)
        HG[j]    = o.getMidPrice()
        j+=1
    r = W* (HG - HOV_e)
    RMSE = sqrt(dot(r,r)/len(r))
    return RMSE

def Error_FunctionD(candidates, args):
        global RMSE
        global t_cos

        fitness = []
        nk = len(chain)
        for p0 in candidates:
                kappa, theta, sigma, rho, v0, a, b, lmbda = p0
                if bIncorporateNLContraint:
                     kappa = (kappa + sigma**2)/(2.0*theta)
                HOV_e = zeros(nk)
                HG = zeros(nk)
                j = 0
                for o in chain:
                    HOV_e[j] = BatesCOS(o.underlying,o.strike,o.maturity,r0,sigma,kappa,theta,v0,rho,a, b, lmbda, o.typ,nInt)
                    HG[j]    = o.getMidPrice()
                    j+=1
                r = W*(HG - HOV_e)
                RMSE = sqrt(dot(r,r)/len(r))
                fitness.append(RMSE)
        return fitness
    
def my_terminator(population, num_generations, num_evaluations, args):
    best=max(population)
    tol = args['tolerance']
    if (best.fitness <=tol):
     print ("Exiting after: " + str(num_evaluations) + " evalutions")
    return ((num_evaluations >=  args['max_evaluations']) or (best.fitness <= tol))

def Generate_Function(random, args):# set a,b,lmbda as per Professor
    return [random.uniform(0.0, 5.0),random.uniform(0.0, 1.0), 
            random.uniform(0.0, 1.0),random.uniform(-1.0, 1.0),
            random.uniform(0.0, 1.0),random.uniform(-1.0, 5.0),#a=(-1,5)
            random.uniform(0.0, 5.0),random.uniform(0.0, 5.0)]#b=(0,5),lmbda=(0,5)

def defaultTimes(hazardRates, correlation, samples):

    if (hazardRates < 0).any() or (hazardRates > 1).any():
        raise ValueError('hazard rates can not be negative or larger than 1')
    if (samples < 0).any() or (samples > 1).any():
        raise ValueError('samples can not be negative or larger than 1')
    if (np.abs(correlation) > 1).any():
        raise ValueError('absolute value of correlation can not be larger than 1')
    if (correlation.diagonal() != 1).any():
        raise ValueError('correlation matrix must have all 1s in diagnal')
        
    ch = np.linalg.cholesky(correlation).T
    z = np.dot( samples, ch )
    x = norm(z)
    
    default_times = -np.log(1-x)/hazardRates

    return default_times

def portfolioLosses(notionals, maturities, recoveryRates, hazardRates, correlation, samples, r):
    
    if (hazardRates < 0).any() or (hazardRates > 1).any():
        raise ValueError('hazard rates can not be negative or larger than 1')
    if (samples < 0).any() or (samples > 1).any():
        raise ValueError('samples can not be negative or larger than 1')
    if (np.abs(correlation) > 1).any():
        raise ValueError('absolute value of correlation can not be larger than 1')
    if (correlation.diagonal() != 1).any():
        raise ValueError('correlation matrix must have all 1s in diagnal')
    if (notionals < 0).any():
        raise ValueError('notionals can not be negative')
    if (maturities < 0).any():
        raise ValueError('maturities can not be negative')
    if (recoveryRates < 0).any() or (recoveryRates > 1).any():
        raise ValueError('hazard rates can not be negative or larger than 1')        

    tau = defaultTimes(hazardRates, correlation, samples)


    V = (1-recoveryRates)*np.exp(-r*np.where(tau < maturities, tau, maturities))*notionals
    losses = V
    losses[ tau > maturities] = 0
    
    summ_losses =[sum(scenario) for scenario in losses]
    expected_loss = np.mean(summ_losses)
    return expected_loss


if __name__ == '__main__':
    np.set_printoptions(precision=3)    

# Code to test Question 2
    print('Question 2:')
    M = 10
    N = 5
    hazardRates = np.random.rand(N)
    corr = 0.6
    correlation = corr*np.ones((N,N))
    np.fill_diagonal(correlation, 1)
    samples = np.random.rand(M,N)

    print('Inputs:')
    print('Hazard rates: {:}'.format(hazardRates))
    print('Correlation matrix: {:}'.format(correlation))
    print('samples:\n{:}\n'.format(samples))
   
    
    def_times = defaultTimes(hazardRates, correlation, samples)
    print('default times:\n{:}\n'.format(def_times))
# Code to test Question 2 - end

# Code to test Question 3

    print('Question 3:')
    recoveryRates  = np.array([0.1, 0.3, 0.15, 0.4, 0.4, 0.4])
    hazardRates  = np.array([0.3, 0.03, 0.05, 0.1, 0.0010, 0.9])
    notionals = np.array([1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0])
    rho = 0.4
    r = 0.01
    M = 10000 #simulation paths
    N = len(recoveryRates)#6 #number of instruments
    T = 3
    maturities = np.array([T, T, T, T, T, T])

    correlation = rho*np.ones((N,N))
    np.fill_diagonal(correlation, 1)
    samples = np.random.rand(M,N)
    
    print('Inputs:')
    print('Notionals: {:}'.format(notionals))
    print('Recovery rates: {:}'.format(recoveryRates))
    print('Hazard rates: {:}'.format(hazardRates))
    print('Maturities: {:}'.format(maturities))
    print('rho: {:}'.format(rho))
    print('r: {:}\n'.format(r))

    
    expected_loss = portfolioLosses(notionals, maturities, recoveryRates, hazardRates, correlation, samples, r)
    print('Expected loss: $ {:.2f}'.format(expected_loss))

    # this is to find out how correlation affects loss
    exp_lss = []
    rhos = np.linspace(0.0, 1.0, num=100, endpoint=False)

    for rh in rhos:
        correlation = rh*np.ones((N,N))
        np.fill_diagonal(correlation, 1)
        exp_lss.append(portfolioLosses(notionals, maturities, recoveryRates, hazardRates, correlation, samples, r))
    
    plt.plot(rhos,exp_lss,)    
    plt.grid(True)    
    plt.xlabel('Correlation')
    plt.ylabel('Expected Loss, $M')
    plt.show()
    
# Code to test Question 3 - end    

# Implementation code for Question 1
    print('Question 1:')
    
    nInt = 32
    iMethod = 'DEA'
    bIncorporateNLContraint = False
    r0 = 0.0002
    chain = readTDAmeritradeData('data/SPX_08082013.csv', nInt)
    
    eps = 1e-8
    # kappa, ... ,a,b, lmbda
    l = [eps,eps,eps,-1.0 + eps, eps, -1+eps, eps, eps]
    u = [5.0-eps,1.0-eps,1.0-eps,1.0-eps,1.0-eps, 5.0-eps, 5.0-eps, 5.0-eps]
    
    if bIncorporateNLContraint:
        u[0] = 2.0 * u[0]* u[1] - l[2]*l[2]
       
    prng = Random()
    
    if iMethod == 'DEA':
      ea = inspyred.ec.DEA(prng)
    elif iMethod == 'SA':
      ea = inspyred.ec.SA(prng)
    else:
      print ('Error: must specify inspyred method as DEA or SA')
      sys.exit()
    
    ea.terminator = my_terminator
    final_pop = ea.evolve(generator=Generate_Function,
                              evaluator=Error_FunctionD,
                              pop_size=100,
                              bounder=ec.Bounder(l, u),
                              maximize=False,
                              max_evaluations=500,
                              tolerance = 0.0001)
    
    
    best = max(final_pop)
    print('Best Solution: \n{0}'.format(str(best)))
    xmin = best.candidate
    
    bounds_0 = [(l[0],u[0]),(l[1],u[1]),(l[2],u[2]),(l[3],u[3]),(l[4],u[4]),(l[5],u[5]),(l[6],u[6]),(l[7],u[7])]
    
    opt=fmin_tnc(Error_Function, xmin, approx_grad=True, bounds=bounds_0, epsilon=1e-04, scale=None, \
                  offset=None, messages=15, maxCGit=-1, maxfun=1000, eta=-1, stepmx=0, accuracy=0, \
                  fmin=0, ftol=2.22e-16, xtol=2.22e-16, pgtol=-1, rescale=-1, disp=5)
    print (opt[0])
    print (opt[1])
    print ("Number of options", len(chain))
    print ("Final RMSE Value: ", RMSE)
    if bIncorporateNLContraint:
        opt[0][0] = (opt[0][0] + opt[0][2]*opt[0][2])/(2.0*opt[0][1]) 
    
    print ("kappa=" + str(opt[0][0]))
    print ("theta=" + str(opt[0][1]))
    print ("sigma=" + str(opt[0][2]))
    print ("rho=" + str(opt[0][3]))
    print ("v0=" + str(opt[0][4]))        
    # print 3 additional parameters that were introduced
    print ("a=" + str(opt[0][5]))
    print ("b=" + str(opt[0][6]))
    print ("lmbda=" + str(opt[0][7]))
    print('\n')
    
# Code for Question 1 - end