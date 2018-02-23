# -*- coding: utf-8 -*-
"""

"""

import math
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

class verlet: 
    def __init__(self, afun):
        self.afun = afun
        self.__name__='verlet'
    def step(self, x,v,a,t,dt):
        x = x+v*dt+0.5*a*dt*dt
        t= t+dt
        ap =self.afun(x, t)
        v = v+(ap+a)*0.5*dt
        a = ap
        return x,v,a
    
class euler_imp :
    def __init__(self, afun, dafun):
       self.afun  = afun
       self.dafun = dafun 
       self.__name__='euler_implicit' 
    def step(self, x,v,a,t,dt):
        k =0
        xp = x+v*dt+0.5*a*dt*dt;
        ap = self.afun(xp,t+dt)
        r = xp-x-dt*v-dt*dt*ap;
  #      print(" k " , k, " ' r ",  np.fabs(r))
        while ((np.abs (r) > 1.e-6) & (k < 10)):
            drdx = 1-dt*dt*self.dafun(xp,t+dt)
            xp= xp -r/drdx
            ap = self.afun(xp, t+dt)
            r = xp-x-dt*v-dt*dt*ap;
            k = k+1
 #           print(" k " , k, "  r ",  np.fabs(r))
            #if (k ==12) return x, v ,a 
        if (np.abs (r) > 1.e-6) :
            print("not converged k = ",k , " res = ", np.fabs(r))
        t = t+dt
        x =xp
        a =ap
        v = v + dt*ap
        return x,v, a
    
class newmark :
    def __init__(self, afun, dafun, gamma, beta):
       self.afun  = afun
       self.dafun = dafun 
       self.__name__='newmark' 
       self.gamma = gamma
       self.beta = beta
    def stepold(self, x,v,a,t,dt):
        k =0
        xp = x+v*dt + 0.5*a*dt*dt;
        ap = self.afun(xp,t+dt)
        r = x+dt*v+dt*dt*0.5*((1.-2.*self.beta)*a+2*self.beta*ap) - xp;
        print(" k " , k, " ' r ",  np.fabs(r), " x " , x)
        while ((np.abs (r) > 1.e-6) & (k < 10)):
            drdx = dt*dt*self.beta*self.dafun(xp,t+dt) -1
            xp= xp -r/drdx
            ap = self.afun(xp, t+dt)
            r = x+dt*v+dt*dt*0.5*((1.-2.*self.beta)*a+2*self.beta*ap) - xp;
            k = k+1
            print(" k " , k, "  r ",  np.fabs(r), " x " , x,  " drdx ",drdx)
            #if (k ==12) return x, v ,a 
        if (np.abs (r) > 1.e-6) :
            print("not converged k = ",k , " res = ", np.fabs(r))
        t = t+dt
        x =xp
        v = v + dt*((1.-self.gamma)*a +self.gamma*ap)
        a =ap
       
        return x,v, a    
    def step(self, x,v,a,t,dt):
      
        r0 = x+dt*v+dt*dt*0.5*((1.-2.*self.beta)*a)
        ap =0.
        def r(xp):
            ap = self.afun(xp,t+dt)
            return r0 + dt*dt*self.beta*ap - xp;
        def drdx(xp):
            return dt*dt*self.beta*self.dafun(xp,t+dt) -1.      
        x0 = x+v*dt + 0.5*a*dt*dt;
        #solver = newton_linesearch();
        solver = newton( r, drdx);
        solver = fixed_point( r);
        xp = solver.solve(x0)
        
        x =xp
        ap = self.afun(xp,t+dt) 
        # !!!!!!!!!i have to recompu a here ...
        # but already computed during newton ...
        v = v + dt*((1.-self.gamma)*a +self.gamma*ap)
        a =ap
        t = t+dt
       
        return x,v, a
    
class newton :
    def __init__(self, r, drdx):
        self.r = r
        self.drdx = drdx
    def solve(self, x0):
        k=0
        xp = x0
        res = self.r(xp)
        print(" k = ",k , " res = ", np.abs(res)) 
        while ( (np.abs(res) > 1.e-6) & (k < 10)):
          xp= xp - res/self.drdx(xp)
          res = self.r(xp)
          k = k+1 
          print(" k = ",k , " res = ", np.abs(res)) 
        if (np.abs (res) > 1.e-6) :
            print("not converged k = ",k , " res = ", np.abs(res)) 
        return xp

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
    
def verlet_h( x0, v0, dt, nstep, afun, dafun):
    x = x0
    v = v0
    a = afun(x)
    tilt =dafun(x)*v;
    for i in range(nstep):
        x = x+v*dt+0.5*a*dt*dt+ tilt*dt*dt*dt/6
        ap =afun(x)
        v = v + (ap+a)*0.5*dt#+ tilt*dt*dt*0.5
        dap =dafun(x);
        tiltp = dap*v
        v = v + (tilt-tiltp)*dt*dt/12.
        a = ap
        tilt = dap*v
    return x,v

def verlet_hh( x0, v0, dt, nstep, afun, dafun):
    x = x0
    v = v0
    a = afun(x)
    tilt =dafun(x)*v;
    bling = dafun(x)*a;
    for i in range(nstep):
        x = x+v*dt+0.5*a*dt*dt+ tilt*dt*dt*dt/6 + bling*dt*dt*dt*dt/24
        ap =afun(x)
        v = v + (ap+a)*0.5*dt#+ tilt*dt*dt*0.5
        dap =dafun(x);
        tiltp = dap*v
        v = v + (tilt-tiltp)*dt*dt/12.
        a = ap
        tilt = dap*v
        bling = dafun(x)*a;
    return x,v

def euler( x0, v0, a0, t, dt, nstep, afun):
    x = x0
    v = v0
    a=  a0
    for i in range(nstep):
        x = x+dt*v;
        v = v+dt*a; 
        t=t+dt
        a = afun(x, t);
    return x,v,a



def euler_h( x0, v0, dt, nstep, afun, dafun):
    x = x0
    v = v0
    for i in range(nstep):
        a  = afun(x);
        tilt = dafun(x)*v;
        x = x+dt*v+0.5*a*dt*dt+dt*dt*dt/6.*tilt;
        v = v+dt*a+0.5*tilt*dt*dt;
    return x,v
def euler_cromer( x0, v0, dt, nstep, afun):
    x = x0
    v = v0
    for i in range(nstep):
        a = afun(x);
        v = v+dt*a;
        x = x+dt*v;
    return x,v

def second_order_centered( x0, v0, x1, v1, dt, nstep, afun):
   # print ("start" , x0, v0, x1, v1)
    for i in range(nstep):
      v2 = 2*dt*afun(x1) + v0
      x2 = 2*dt*v1 + x0
      x0 = x1
      v0 = v1
      x1 = x2
      v1 = v2 
    #print ("end", x0, v0, x1, v1)
    return (x0, v0, x1, v1)
    
