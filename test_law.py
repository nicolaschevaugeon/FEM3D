#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:42:19 2018

@author: chevaugeon
"""

import math
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import law


plt.close("all")
sig0 = 1.
E = 0.1
#mat =  law.perfect_plastic(E, sig0)
mat =  law.hardening_plastic2(E, sig0, 0.1)
#mat =  law.perfect_plastic(E, sig0)
#mat = law.linear_elastic(E)

ifig = 0
static = 0
if static:
    mat.clear()
    N= 5000
    t= np.linspace(0.,1.,N)
    epsmax = (sig0/E)*2
    eps = epsmax*np.sin(50*np.pi*t)
    sig = np.zeros(eps.shape)
    epsp = np.zeros(eps.shape)
    for i in range(eps.size):
        sig[i] = mat.stress(eps[i])
        mat.update()
        epsp[i] = mat.ep;
    plt.figure(ifig)
    plt.plot(t, eps, label = "eps")
    plt.plot(t , sig, label ="sig")
    plt.plot(t, epsp, label = "epsp")#, ttab, alpha*0.5*vtab.*vtab, ttab, alpha*0.5*vtab.*vtab+alpha*0.5+xtab.*xtab ,'-gx')
    plt.xlabel("time")
    plt.legend()
    plt.title("time")
    ifig = ifig+1
    plt.figure(ifig)
    plt.plot( eps, sig) #, ttab, alpha*0.5*vtab.*vtab, ttab, alpha*0.5*vtab.*vtab+alpha*0.5+xtab.*xtab ,'-gx')
    plt.title("sig_eps")
    ifig = ifig+1

dynamic = 1
if dynamic :
    import time_integrator 
    masse = 1.
    S = 1.
    L =1.
    ifig
    mat.clear()
   
    x0 = 0.;
    v0 = 0.;
    t0 = 0.;
    f0 = 1.4;
    ramp = 2*np.pi*25.12;
    def fe(t) :
         #return min (f0, f0*t/ramp)
         return f0*np.sin(t/25.12)
         return f0
     
    def afun (x,t):
        eps = x/L;
        sig = mat.stress(eps)
        return fe(t)-sig*S/masse
    
    def dadxfun(x,t):
        eps = x/L;
        dsig = mat.dstressdeps(eps)
        return -dsig*S/masse/L
        
   
   
    sig0 = mat.stress(x0/L)
    stored_energie0 = mat.stored_energie()*S*L
    Ec0  = 0.5*v0*v0*masse;
    a0 = afun(x0, t0);
    fe0 = fe(t0);
    W0 = 0.;
    #dt = 0.25;
   
    Nstep =2.;
    period = 2*np.pi;
    dt = period/(Nstep)
    N=int(np.floor(Nstep*50));
   
   
   
   
    

 #   integ = time_integrator.verlet(afun)
 #   integ = time_integrator.euler_imp(afun, dadxfun)
 # verlet =difference centrée = newmark gamma =0.5, beta = 0.
 # explicit, few dissipation  
    gamma = 0.5
    beta = 0.
#   Fox Goodwin 
#   implicit, few dissipation
    gamma = 0.5
    beta = 1./12
#   Accélération linéaire
    gamma = 0.5
    beta = 1./6.
#   Accélération moyenne
    gamma = 0.5
    beta = 1./4.
    gamma = 0.5
    beta = 1./4.
    integ = time_integrator.newmark(afun, dadxfun, gamma, beta)
    ttab = [];
    xtab = [];
    vtab = [];
   
    fetab=[];
    Wetab=[];
    sigtab =[];
    stored_energietab =[];
    dissipated_energietab =[];
    Ectab = [];
    epsptab=[];
    Eetab = [];
    Edtab = [];
    Wetab = [];
    
    ttab.append(t0)
    xtab.append(x0)
    vtab.append(v0)
    fetab.append(fe0)
    sigtab.append(sig0)
    stored_energietab.append(stored_energie0)
    dissipated_energietab.append(0.)
    Ectab.append(Ec0)
    Wetab.append(W0)
   # wtab.append(x0*f)
   
    t = t0
    x = x0
    v = v0
    a = a0
    W = W0
    for i in range(N):
 #     print (i+1)
      [x,v,a] = integ.step(x,  v, a, t, dt) #, dadxfun)
      t =t +dt;
      ttab.append(t);
      xtab.append(x);
      vtab.append(v);
      sigtab.append( mat.sig_trial)
      stored_energietab.append(mat.stored_energie()*S*L)
      Edp = dissipated_energietab.pop()
      dissipated_energietab.append(Edp)
      dissipated_energietab.append(Edp+mat.dissipated_increment()) 
      Wep= Wetab.pop()
      Wetab.append(Wep)
      Wetab.append(Wep+dt*v*fe(t))
      mat.update();
      Ectab.append(0.5*v*v*masse)
      fetab.append(fe(t));
   # print(xtab)  
    #Ep = 0.5*alpha*np.array(xtab)*np.array(xtab)
    plt.figure(ifig)
    #plt.hold(0)
    plt.plot(ttab, xtab, label ='x')
    #plt.hold(1)
    plt.plot(ttab, vtab, label ='v')
    plt.plot(ttab, np.array(sigtab),  label ='sig')
    #plt.plot(ttab, mat.hep, label ='Ep')
    plt.plot(ttab, np.array(fetab), label ='fe')
    plt.title(integ.__name__)
    plt.legend()
    ifig = ifig+1
    plt.figure(ifig)
  #   plt.plot(ttab, mat.hEd , label ='Ed')
    
    plt.plot(ttab, stored_energietab, label ='Ee')
    plt.plot(ttab, Ectab, label ='Ec')
    plt.plot(ttab, dissipated_energietab, label ='Ed')
    plt.plot(ttab, Wetab, label ='We')
    plt.plot(ttab, np.array(Ectab)+np.array(stored_energietab)+np.array(dissipated_energietab)-np.array(Wetab), label ='Ec+Ee+Ed-We')
#    plt.plot(ttab, Ec+mat.hEE+mat.hEd-W, label ='Ec+Ee+Ed-W')
    plt.title(integ.__name__+"_energie")
    plt.legend()
    ifig = ifig+1
    
   

