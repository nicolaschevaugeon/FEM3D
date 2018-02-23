# -*- coding: utf-8 -*-
"""

"""

import math
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import time_integrator


plt.close("all")

integ_list = [ time_integrator.verlet]
ifg =1; 
alpha =1.;
x0 = 1.;
v0 = 0.;
t0 = 0.;
#dt = 0.25;
N = 100;
Nstep = 5;
period = 2*np.pi*alpha;
dt = period/(Nstep)
N=Nstep*1;



te_tab  =  np.linspace(0., dt*N, N+1)
xe_tab  =  np.cos(alpha*te_tab);
ve_tab  = -alpha*np.sin(alpha*te_tab);
Ece_tab =  0.5*ve_tab*ve_tab;
Epe_tab =  0.5*alpha*xe_tab*xe_tab;

tef_tab  =  np.linspace(0., dt*N, (100*N)+1)
xef_tab  =  np.cos(alpha*tef_tab);
vef_tab  = -alpha*np.sin(alpha*tef_tab);
Ecef_tab =  0.5*vef_tab*vef_tab;
Epef_tab =  0.5*alpha*xef_tab*xef_tab;
def gfun (x):
    return -alpha*x

def dgfun(x):
    return -alpha



    

for integ in integ_list:
    xtab = [];
    vtab = [];
    ttab = [];
    xtab.append(x0)
    vtab.append(v0)
    ttab.append(t0)
   
    t = t0
    x = x0
    v = v0
    for i in range(N):
      [x,v] = integ(x,  v, dt, 1, gfun)
      t =t +dt;
      xtab.append(x);
      vtab.append(v);
      ttab.append(t);
    Ec = 0.5*np.array(vtab)*np.array(vtab)
    Ep = 0.5*alpha*np.array(xtab)*np.array(xtab)
    plt.figure(ifg)
    plt.plot(ttab, xtab, ttab, vtab,'gx', te_tab, xe_tab)
    plt.title(integ.__name__)
    ifg = ifg+1
    plt.figure(ifg)
    plt.plot(ttab, Ec, ttab , Ep, ttab, Ep+Ec) #, ttab, alpha*0.5*vtab.*vtab, ttab, alpha*0.5*vtab.*vtab+alpha*0.5+xtab.*xtab ,'-gx')
    plt.title(integ.__name__+"_energie")
    print(integ.__name__)
    print (" error x = ", max(abs(xe_tab-xtab)), ", error v ", max(abs(ve_tab-vtab)), ", error e ", max(abs(Ec+Ep - Ece_tab -Epe_tab)) )
    
    ifg = ifg+1
    
    
integ_tilt_list = [time_integrator.verlet_h, time_integrator.verlet_hh]  
for integ_tilt in integ_tilt_list:
    xtab = [];
    vtab = [];
    ttab = [];
    xtab.append(x0)
    vtab.append(v0)
    ttab.append(t0)
   
    t = t0
    x = x0
    v = v0
    for i in range(N):
      [x,v] = integ_tilt(x,  v, dt, 1, gfun, dgfun)
      t =t +dt;
      xtab.append(x);
      vtab.append(v);
      ttab.append(t);
    Ec = 0.5*np.array(vtab)*np.array(vtab)
    Ep = 0.5*alpha*np.array(xtab)*np.array(xtab)
    plt.figure(ifg)
    plt.plot(ttab, xtab, ttab, vtab,'-gx', ttab, mat.hsig )
    plt.title(integ_tilt.__name__)
    ifg = ifg+1
    plt.figure(ifg)
    plt.plot(ttab, Ec, ttab , Ep, ttab, Ep+Ec) #, ttab, alpha*0.5*vtab.*vtab, ttab, alpha*0.5*vtab.*vtab+alpha*0.5+xtab.*xtab ,'-gx')
    plt.title(integ_tilt.__name__+"_energie")
    print(integ_tilt.__name__)
    print (" error x = ", max(abs(xe_tab-xtab)), ", error v ", max(abs(ve_tab-vtab)), ", error e ", max(abs(Ec+Ep - Ece_tab -Epe_tab)) )
    
    ifg = ifg+1
xtab = [];
vtab = [];
ttab = [];
xtab.append(x0)
vtab.append(v0)
ttab.append(t0)
#verlet for first step 


#==============================================================================
# [x1,v1] = verlet(x0, v0, dt,1, gfun)
# t1 = t0+dt
# xtab.append(x1)
# vtab.append(v1)
# ttab.append(t1)
# integ = second_order_centered
# for i in range(N-1):
#      [x0, v0, x1, v1] = integ( x0, v0, x1, v1, dt, 1, gfun)
#      t1= t1+dt
#      xtab.append(x1);
#      vtab.append(v1);
#      ttab.append(t1);
# Ec = 0.5*np.array(vtab)*np.array(vtab)
# Ep = 0.5*alpha*np.array(xtab)*np.array(xtab)
# plt.figure(ifg)
# plt.plot(ttab, xtab, ttab, vtab,'-gx')
# plt.title(integ.__name__)  
# ifg = ifg+1
# plt.figure(ifg)
# plt.plot(ttab, Ec, ttab , Ep, ttab, Ep+Ec) #, ttab, alpha*0.5*vtab.*vtab, ttab, alpha*0.5*vtab.*vtab+alpha*0.5+xtab.*xtab ,'-gx')
# plt.title(integ.__name__+"_energie")
# ifg = ifg+1  
# print(integ.__name__)
# print (" error x = ", max(abs(xe_tab-xtab)), ", error v ", max(abs(ve_tab-vtab)), ", error e ", max(abs(Ec+Ep - Ece_tab -Epe_tab)) )
# 
# 
#==============================================================================
     
 