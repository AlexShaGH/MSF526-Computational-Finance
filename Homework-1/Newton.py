# -*- coding: utf-8 -*-

# Newton.py - Root findin implementations of one-dimensional equations using
# Newton's methods
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import math

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
        the set of x values tried by the solver, xvals[0] represents start
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

        if bounds is not None:
            if xvals[-1] < min(bounds) or xvals[-1] > max(bounds):
                raise  ValueError('value x = {0} is Out of bounds {1}'.format(
                    xvals[-1],bounds))
                
        if abs(xvals[-1]-xvals[-2])<tols[0] or abs(fdiffs[-1])<tols[1]:
            return xvals, fdiffs
      
        n = n + 1
        
    raise ValueError('Number of iterations exceeded limit:{0}'.format(maxiter))