# -*- coding: utf-8 -*-

# Newton.py - Root findin implementations of one-dimensional equations using
# bisection methods
# MSF 526
# Illinois Institute of Technology
# Homework 1
# Author: Oleksandr Shashkov
# ID: A20229995
# Email: oshashko@hawk.iit.edu

__author__ = "oshashkov"

import math

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
    if math.isnan(target):
        raise ValueError("target argument can not be NaN")
    if maxiter <= 0:
        raise ValueError('maxiter must be positive integer')

    if math.isnan(start) and bounds == None:
        raise  ValueError("start and bounds can not be both NaN")
        
    if bounds is not None:
        upper_bound = max(bounds)
        lower_bound = min(bounds)
    else:
        found = False
        i = 0
        epsilon = tols[0]
        upper_bound = start
        lower_bound = start
        
        while not found and i < maxiter:
            upper_bound += epsilon*i
            lower_bound -= epsilon*i
            if function(upper_bound)*function(lower_bound) < 0:
                found = True
            i+=1
        
        if not found:
            raise ValueError(
                'Unable to expand boundaries with start={0}, tol={1},\
                maxiter={2}'.format(start,tols[0],maxiter))
    
    xvals=[]
    fdiffs=[]
    
    n = 1
    while n <= maxiter:
        midpoint = (upper_bound + lower_bound) / 2
        xvals.append(midpoint)
        fdiffs.append(target - function(midpoint))
        
        if abs(upper_bound-lower_bound)/2 < tols[0] or abs(fdiffs[-1])<tols[1]:
            return xvals, fdiffs
        
        if 0 < fdiffs[-1]:
            lower_bound = midpoint
        else:
            upper_bound = midpoint
      
        n = n + 1
        
    raise ValueError('Number of iterations exceeded limit:{0}'.format(maxiter))