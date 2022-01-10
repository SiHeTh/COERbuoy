#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 08:55:59 2020

@author: heiko
"""

import numpy as np;
import json;
from COERbuoy.floater import Floater;
from COERbuoy.floater_BEM_LUT import Floater_BEM;
#from scipy.fftpack import hilbert as KramerKronig;
from scipy.interpolate import interp1d;

pi=np.pi;
#wave=wavefield.wavefield(np.zeros(20),np.linspace(1,20),2)

class Floater_LIN(Floater_BEM):   
    eq_force=None;
    res=None;
    def get_forces(self, t, wave, p, v, a):
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));
        
        if self.eq_force is None:
            self.eq_force=self.Calculate(0,0,0,0)[0][1];
            self.res=self.Calculate(0*z0, 0*x0, 0*delta0, 0*eta);
            
        #res=self.Calculate(0*z0, 0*x0, 0*delta0, eta);#Calculate coefficents
        
        ret=[0,0,0];#return array
        exc1 = np.array(self.res[1]);#Exitation force
        rad = np.real(self.res[2]);#radiation over omega
        amw = np.imag(self.res[2]);#added mass over omega
        am1 = np.array(self.res[3]);#added mass @ inf
        rad1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1)*(exc1[1]*np.conjugate(exc1[1]))
        #rad1=(self.omega*(self.omega**2)/9.81)/(4*1000*3.14*9.81**a2)*(exc1[1]*np.conjugate(exc1[1]));
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        #print("Rad")
        #print([rad[1][1],np.real(rad1)])
        #print([exc1[1],rad[1][1],self.omega])
        #exit();
        rad[1][1]=np.real(rad1);
        #Generate wave from movement
        if (np.sum(np.abs(exc1))>0):
            #r1=am_omom[1][1]/(exc1[1])*v[1]+rad[0][1]/(exc1[1])*v[0];
            #r2=am_omom[0][0]/(exc1[0])*v[0]+rad[1][0]/(exc1[0])*v[1];
            wave.add_diracWave(-2/np.pi*(rad[1][1]*(v[1])),t,True);
            wave.add_diracWave2(-2/np.pi*(rad[0][0]*(v[0])),t,True);
        
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            #FK=np.sum(np.real(exc1[i]*Awave[0]));
            FK=np.sum(m(exc1[i],Awave[0])).real;
            ret[i]=FK;#buoyance + FK force
            if i==1:
                #print(np.sum(np.real(exc1[i])*Awave[0]))
                ret[i]=ret[i]-self.Area(0)*self.g*self.rho*(z0)+self.eq_force;#self.Volume(0)*self.rho*self.g;
        #print([np.real(np.sum(wave.get_rad2(t,x0))),np.real(np.sum(wave.get_rad(t,x0)*(exc1[1])))])
        Frad=[np.real(np.sum(wave.get_rad2(t,x0))),np.real(np.sum(wave.get_rad(t,x0))),0];#radiation force
        #print(exc1/(8*rad))
        #exit();
        ret=np.array(ret)+np.array(Frad);
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf
