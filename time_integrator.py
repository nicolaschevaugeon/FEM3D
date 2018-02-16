# -*- coding: utf-8 -*-
"""

"""

import math
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


def verlet( x0, v0, dt, nstep, afun):
    x = x0
    v = v0
    a = afun(x)
    for i in range(nstep):
        x = x+v*dt+0.5*a*dt*dt
        ap =afun(x)
        v = v+(ap+a)*0.5*dt
        a = ap
    return x,v

def verlet_h( x0, v0, dt, nstep, afun, dafun):
    x = x0
    v = v0
    a = afun(x)
    tilt =dafun(x)*v;
    for i in range(nstep):
        x = x+v*dt+0.5*a*dt*dt+ tilt*dt*dt*dt/6
        ap =afun(x)
        dap =dafun(x);
        tiltp = dap*(v+(ap+a)*0.5*dt)
       # v = v+(ap+a)*0.5*dt+0.5*dap*v*dt*dt
#        v = vpred+dt*dt*0.25*(tilt+dap*vpred)
#        v = v+(ap+a)*0.5*dt+(tilt-tiltp)*dt/12.
#        v = (v+(ap+a)*0.5*dt+(tilt)*dt/12.)/(1+dap*dt/12)
        v = v+(ap+a)*0.5*dt
 #       v = (v+(ap+a)*0.5*dt+0.25*(tilt)*dt*dt)/(1..-0.25*dap*dt*dt)
        
        a = ap
    return x,v

def euler( x0, v0, dt, nstep, afun):
    x = x0
    v = v0
    for i in range(nstep):
        a = afun(x);
        x = x+dt*v;
        v = v+dt*a;
    return x,v
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
    

plt.close("all")

integ_list = [euler, euler_cromer, verlet]
ifg =1; 
alpha =1.;
x0 = 1.;
v0 = 0.;
t0 = 0.;
#dt = 0.25;
#N = 80;
Nstep = 400;
period = 2*np.pi*alpha;
dt = period/(Nstep)
N=Nstep*1;



te_tab  =  np.linspace(0., dt*N, N+1)
xe_tab  =  np.cos(alpha*te_tab);
ve_tab  = -alpha*np.sin(alpha*te_tab);
Ece_tab =  0.5*ve_tab*ve_tab;
Epe_tab =  0.5*alpha*xe_tab*xe_tab;

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
    plt.plot(ttab, xtab, ttab, vtab,'-gx')
    plt.title(integ.__name__)
    ifg = ifg+1
    plt.figure(ifg)
    plt.plot(ttab, Ec, ttab , Ep, ttab, Ep+Ec) #, ttab, alpha*0.5*vtab.*vtab, ttab, alpha*0.5*vtab.*vtab+alpha*0.5+xtab.*xtab ,'-gx')
    plt.title(integ.__name__+"_energie")
    print(integ.__name__)
    print (" error x = ", max(abs(xe_tab-xtab)), ", error v ", max(abs(ve_tab-vtab)), ", error e ", max(abs(Ec+Ep - Ece_tab -Epe_tab)) )
    
    ifg = ifg+1
    
    
integ_tilt_list = [euler_h,verlet_h]  
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
    plt.plot(ttab, xtab, ttab, vtab,'-gx')
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
[x1,v1] = verlet(x0, v0, dt,1, gfun)
t1 = t0+dt
xtab.append(x1)
vtab.append(v1)
ttab.append(t1)
integ = second_order_centered
for i in range(N-1):
     [x0, v0, x1, v1] = integ( x0, v0, x1, v1, dt, 1, gfun)
     t1= t1+dt
     xtab.append(x1);
     vtab.append(v1);
     ttab.append(t1);
Ec = 0.5*np.array(vtab)*np.array(vtab)
Ep = 0.5*alpha*np.array(xtab)*np.array(xtab)
plt.figure(ifg)
plt.plot(ttab, xtab, ttab, vtab,'-gx')
plt.title(integ.__name__)  
ifg = ifg+1
plt.figure(ifg)
plt.plot(ttab, Ec, ttab , Ep, ttab, Ep+Ec) #, ttab, alpha*0.5*vtab.*vtab, ttab, alpha*0.5*vtab.*vtab+alpha*0.5+xtab.*xtab ,'-gx')
plt.title(integ.__name__+"_energie")
ifg = ifg+1  
print(integ.__name__)
print (" error x = ", max(abs(xe_tab-xtab)), ", error v ", max(abs(ve_tab-vtab)), ", error e ", max(abs(Ec+Ep - Ece_tab -Epe_tab)) )


     
 