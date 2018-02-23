#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:41:06 2018

@author: chevaugeon
"""
import numpy as np
        
class perfect_plastic:
    def __init__(self, E, sig0, ep0 = 0.):
        self.sig0 = sig0;
        self.E = E;
       
        self.ep_trial  = ep0
        self.eps_trial = 0.
        self.sig_trial = 0.
        self.ep = ep0
        self.eps = 0.
        self.sig = 0.
        
    def stress(self, eps):
        self.eps_trial = eps
        self.ep_trial = self.ep
        self.sig_trial = self.E*(eps-self.ep)
        if (np.abs(self.sig_trial) - self.sig0 > 0):
            deltaEp = (eps - self.eps)
            self.sig_trial -= self.E*deltaEp
            self.ep_trial += deltaEp
        return self.sig_trial;
    
    def dstressdeps(self, eps):
        self.eps_trial = eps
        sigt = self.E*(eps-self.ep)
        dsigt = self.E
        if (sigt > self.sig0):
            self.ep_trial = eps - self.sig0/self.E
            sigt = self.sig0 
            dsigt = 0.
        elif(sigt < (-self.sig0)):
            self.ep_trial = eps+self.sig0/self.E
            sigt =-self.sig0 
            dsigt = 0.
        self.sig_trial =sigt
        return dsigt;
    
    def stored_energie(self):
        return 0.5*self.sig_trial*(self.eps_trial-self.ep_trial)
    
    def dissipated_increment(self):
        return self.sig_trial*(self.ep_trial-self.ep)
    
    def update( self):
        self.ep =    self.ep_trial
        self.eps =   self.eps_trial
        self.sig = self.sig_trial
        
    def clear(self):
       self.ep_trial = 0.
       self.ep  = 0.
#ecrouissage lineaire isotrope
class hardening_plastic:
    def __init__(self, E, sig0, h, ep0 = 0.):
        self.sig0 = sig0;
        self.E = E;
        self.h = h
        self.eps = 0.
        self.sig = 0
        self.ep =  0.
        self.gam = 0.
        self.ep_trial = 0.
        self.gam_trial = 0. 
        self.eps_trial = 0.
        self.sig_trial = 0
    def yieldfunction(self, sig, gam):
        return np.abs(sig) - (self.sig0+ self.h*gam) 
    def stress(self, eps):
        sigt = self.E*(eps-self.ep)
        deltagamma = 0.
        ft = self.yieldfunction( sigt, self.gam)
        if (ft > 0 ):
            deltagamma = ft/(self.E+self.h)
        deltaep = deltagamma*np.sign(sigt)
        sigt = sigt - self.E*deltaep
        self.ep_trial  = self.ep+deltaep
        self.gam_trial = self.gam+deltagamma
        self.eps_trial = eps
        self.sig_trial =sigt       
        return sigt;
    
    def dstressdeps(self, eps):
        self.eps_trial = eps
        sigt = self.E*(eps-self.ep)
        dsigt = self.E
        if (sigt > self.sig0):
            self.ep_trial = eps - self.sig0/self.E
            sigt = self.sig0 
            dsigt = 0.
        elif(sigt < (-self.sig0)):
            self.ep_trial = eps+self.sig0/self.E
            sigt =-self.sig0 
            dsigt = 0.
        self.sig_trial =sigt
        return dsigt;
    
    def stored_energie(self):
        return 0.5*self.sig_trial*(self.eps_trial-self.ep_trial)
    
    def dissipated_increment(self):
        return self.sig_trial*(self.ep_trial-self.ep)
    
    def update( self):
        self.ep =    self.ep_trial
        self.gam =    self.gam_trial
        
        self.eps =   self.eps_trial
        self.sig = self.sig_trial
        
    def clear(self):
       self.ep  = 0.
       self.ep_trial = 0.
       
#kinematic hardening (prager)
class hardening_plastic2:
    def __init__(self, E, sig0, h, ep0 = 0.):
        self.sig0 = sig0;
        self.E = E;
        self.h = h
        self.eps = 0.
        self.sig = 0
        self.ep =  0.
        self.ep_trial = 0.
        self.eps_trial = 0.
        self.sig_trial = 0
    def yieldfunction(self, sig, X):
        return np.abs(sig-X) - self.sig0
    def stress(self, eps):
        sigt = self.E*(eps-self.ep)
        ft = self.yieldfunction( sigt, self.h*self.ep)
        deltaep = 0.
        if (ft > 0 ):
            epsp = self.eps
            if (self.yieldfunction( self.E*(self.eps-self.ep) , self.h*self.ep) < 0) :
                epsp = ( np.sign( sigt - self.h*self.ep)*self.sig0 + self.h*self.ep)/self.E + self.ep
            deltaep = (eps-epsp)*self.E/(self.E+self.h)
            
        sigt = sigt - self.E*deltaep
        self.ep_trial  = self.ep+deltaep
        self.eps_trial = eps
        self.sig_trial = sigt       
        return sigt;
    
    def dstressdeps(self, eps):
        self.eps_trial = eps
        sigt = self.E*(eps-self.ep)
        ft = self.yieldfunction( sigt, self.h*self.ep)
        if (ft > 0):
            return (self.E*self.h)/(self.h+ self.E)
        return self.E;
    
    def stored_energie(self):
        return 0.5*self.sig_trial*(self.eps_trial-self.ep_trial)
    
    def dissipated_increment(self):
        return self.sig_trial*(self.ep_trial-self.ep)
    
    def update( self):
        self.ep =    self.ep_trial
        self.eps =   self.eps_trial
        self.sig = self.sig_trial
        
    def clear(self):
       self.ep  = 0.
       self.ep_trial = 0.
       self.eps =0.
       self.eps_trial = 0.
       self.sig = 0
       self.sigtrial = 0.
       
class linear_elastic:
    def __init__(self, E):
        self.E = E;
        self.eps_trial = 0.;
        self.sig_trial =0.;
        self.eps = 0.;
        self.sigma = 0.;
        
    def stress(self, eps):
        sigt = self.E*(eps)
        self.eps_trial = eps;
        self.sig_trial = sigt;
        return sigt;
    
    def dstressdeps(self, eps):
        return self.E
    
    def update(self):
        self.eps =   self.eps_trial
        self.sig = self.sig_trial
        
    def clear(self):
        self.eps_trial = 0.
        self.sig_trial = 0. 
        self.eps = 0.
        self.sig = 0.