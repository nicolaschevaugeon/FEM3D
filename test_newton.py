# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 10:39:29 2018

@author: chevaugeon
"""

#import math
import numpy as np
import matplotlib.pyplot as plt
import nonlinear_solver

plt.close('all')
def  r0(x) :
    return np.array([-1.+x[0]+x[0]*x[0]])

def drdx0(x) :
    return np.array([[1. +2.*x[0]]])
    
solver =nonlinear_solver.newton(r0, drdx0)
x= solver.solve(np.array([0.]))
print ("x ", x)
#x = np.linspace(-2.,2.,1000)
#plt.plot(x, r0(x))

def  r1(x) :
    return np.array([-1.+x[0]+x[0]*x[0]+x[1]*x[1], -1.+x[1]+x[1]*x[1]+x[0]*x[0]])

def drdx1(x) :
    return np.array( [[ 1.+2*x[0], 2*x[1] ], [2*x[0] , 1+2*x[1] ]])
    
solver = nonlinear_solver.newton( r1, drdx1) 
x0 = np.array([0.,0.])
solver.solve( x0 )

solver = nonlinear_solver.broyden(r1)

solver.solve(x0, drdx1(x0) )