# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:46:22 2018

@author: chevaugeon
"""
import numpy as np
import numpy.linalg as la

class newton :
    def __init__(self, r, drdx):
        self.r = r
        self.drdx = drdx
    def solve(self, x0):
        k=0
        xp = x0
        res = self.r(xp)
        #print(" k = ",k , " res = ", np.abs(res), "x = ", xp) 
        while ( (la.norm(res) > 1.e-6) & (k < 10)):
          xp= xp - la.solve(self.drdx(xp), res)
          res = self.r(xp)
          k = k+1 
          #print(" k = ",k , " res = ", np.abs(res), "x = ", xp, "slope ", self.drdx(xp)) 
        if (la.norm (res) > 1.e-6) :
            print("not converged k = ",k , " res = ", la.norm(res)) 
        return xp

#taken from broyden 1965        
class broyden:
    def __init__(self, r, drdx):
        self.r = r


class fixed_point :
    def __init__(self, r):
        self.r=r
        
    def solve(self,x0):
        k=0;
        xp = x0
        res = self.r(xp)
        while ( (np.abs(res) > 1.e-6) & (k < 100)):
            xp = xp -res;
            res = self.r(xp)
            k=k+1
            print(" k = ",k , " res = ", np.abs(res)) 
        if (np.abs (res) > 1.e-6) :
            print("not converged k = ",k , " res = ", np.abs(res)) 
        return xp
          
class newton_linesearch :
    def solve(self, x0, r, drdx):
        k=0
        xp = x0
        res = r(xp)
        
        while ( (np.abs(res) > 1.e-6) & (k < 10)):
          alpha = 1.
          slope = drdx(xp)
          xp= xp - alpha*res/slope
          resp = r(xp)
          lk = 0
          while ((np.abs(resp) > np.abs(res)) & (lk < 20)):
              alpha*=0.5
              print("alpha, ", alpha)
              xp= xp - alpha*res/slope
              resp = r(xp)
              lk = lk+1
              print("lk " , lk , " alpha ", alpha, "\n", res, " ", resp, " ", slope  )
          res = resp
          k = k+1
        if (np.abs (res) > 1.e-6) :
            print("not converged k = ",k , " res = ", np.abs(res)) 
        return xp
