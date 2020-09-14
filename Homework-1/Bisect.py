# -*- coding: utf-8 -*-

# Bisect.py - Root findin implementations of one-dimensional equations using
# Newton's and bisection methods
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import math
#from math import log, sqrt, exp
#from typing import Callable, Iterator, Union, Optional, List

def newton(target, function, derivfun, start, bounds=None,
           tols=[0.001, 0.010], maxiter=1000):
    """ Tries to find a root of one-dimensional equation using Newton's method

    Parameters
    ----------
    target : float
        the target value for the function f
    function : function
        the name of the Python function f
    derivfun : float
        the name of the Python function representing df/dx
    start : float
        the x-value to start looking at
    bounds : [float], optional
        the upper and lower bounds beyond which x shall not exceed
        the default is None
    tols : toleranse, stopping criteria, the distance between successive 
        x-values that indicate success and the difference between target 
        and the y-value that indicate success, optional
        the default is [0.001, 0.010]
    maxiter : integer, optional
        maximum iterations the solver will be allowed, the default is 1000

    Returns
    -------
    xvals[] and fdiffs[] upon success or error otherwise

    xvals : [float,float,...]
        the set of x values tried by the solver, xvals(1) represents start
    fdiffs : [float,float,...]
        the set of (target - y) values at the xvals[]
    """
    # check input parameters
    if math.isnan(target):
        raise ValueError("target argument can not be NaN")
    if math.isnan(start):
        raise  ValueError("start can not be NaN")
    if maxiter <= 0:
        raise ValueError('maxiter must be positive integer')
    
    xval = start
    yval = function(xval)
    
    xvals = [xval]
    fdiffs = [target - yval]    
    
    n = 1
    while n <= maxiter:
        try:
            xval = xvals[-1] - function(xvals[-1])/derivfun(xvals[-1])
        except ValueError as error:
            raise ValueError('Unable to devide: {0}'.format(error))
        
        xvals.append(xval)
        fdiffs.append(target - function(xvals[-1]))

        if bounds is not None:#newton.__defaults__[0]:
            if xvals[-1] < min(bounds) or xvals[-1] > max(bounds):
                raise  ValueError('value x = {0} is Out of bounds {1}'.format(
                    xvals[-1],bounds))
                
        if abs(xvals[-1]-xvals[-2])<tols[0] or abs(fdiffs[-1])<tols[1]:
            return xvals, fdiffs
      
        n = n + 1
        
    raise ValueError('Number of iterations exceeded limit:{0}'.format(maxiter))


def bisect(target, function, start=None, bounds=None,
           tols=[0.001, 0.010], maxiter=1000):
    """ Tries to find a root of one-dimensional equation using bisect method
        
    Parameters
    ----------
    target : float
        the target value for the function f
    function : function
        the name of the Python function f
    start : float
        the x-value to start looking at. If None, the mean of the upper and 
        lower bounds shall be used. Subsequent steps shall always use the mean
        of the active bounds. This input is used only when the initial bounds
        have not been supplied.
    bounds : [float, float], optional
        the upper and lower bounds beyond which x shall not exceed
        the default is None
    tols : toleranse, stopping criteria, the distance between successive 
        x-values that indicate success and the difference between target 
        and the y-value that indicate success, optional
        the default is [0.001, 0.010]
    maxiter : integer, optional
        maximum iterations the solver will be allowed, the default is 1000

    Returns
    -------
    xvals[] and fdiffs[] upon success or error otherwise

    xvals : [float,float,...]
        the set of x values tried by the solver, xvals(1) represents start
    fdiffs : [float,float,...]
        the set of (target - y) values at the xvals[]
    """        
    # check input parameters
    if math.isnan(target):
        raise ValueError("target argument can not be NaN")
    if maxiter <= 0:
        raise ValueError('maxiter must be positive integer')

    if math.isnan(start) and bounds == None:
        raise  ValueError("start and bounds can not be both NaN")        
    
    xvals = []
    fdiffs = []
    
    #xval = start
    #yval = function(xval)
    
    #xvals = [xval]
    #fdiffs = [target - yval]    
    
    n = 1
    while n <= maxiter:
        #xval = xvals[-1] - function(xvals[-1])/derivfun(xvals[-1])
        
        xvals.append(xval)
        fdiffs.append(target - function(xvals[-1]))

        #if bounds is not None:#newton.__defaults__[0]:
            #if xvals[-1] < min(bounds) or xvals[-1] > max(bounds):
                #raise  ValueError('value x = {0} is Out of bounds {1}'.format(
                    #xvals[-1],bounds))
                
        #if abs(xvals[-1]-xvals[-2])<tols[0] or abs(fdiffs[-1])<tols[1]:
            #return xvals, fdiffs
      
        n = n + 1
        
    raise ValueError('Number of iterations exceeded limit:{0}'.format(maxiter))    
    
    
    
    return xvals, fdiffs    

target = 1.24
y = lambda x: x**3 + 2*x**2 - 5
dy = lambda x: 3*x**2 + 4*x
start = 5
tols = [0.00001,0.010]
maxiter = 10
bounds = None

xvals, fdiffs = newton(target,y,dy,start,tols=tols,maxiter=maxiter)
print(xvals, fdiffs)
print(xvals[-1])
print(fdiffs[-1])

start = math.nan

xvals, fdiffs = bisect(target,y,start,tols=tols,maxiter=maxiter)
print(xvals, fdiffs)