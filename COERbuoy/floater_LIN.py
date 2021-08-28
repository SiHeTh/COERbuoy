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
    
    def get_forces(self, t, wave, p, v, a):
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        Awave=wave.get(t,x0);
        eta=np.sum(Awave[0]);
        
        res=self.Calculate(0*z0, 0*x0, 0*delta0, eta);#Calculate coefficents
        
        ret=[0,0,0];#return array
        exc1 = np.array(res[1]);#Exitation force
        am_omega = np.real(res[2]);#added mass over omega
        am1 = np.array(res[3]);#added mass @ inf
        am1=[0,Floater.added_mass(self,0),0];
        
        am_omom=am_omega#*self.omega;
        #am_omom[0][0]=am_omom[0][0]-am_omom[0][0][-1];
        #am_omom[0][1]=am_omom[0][1]-am_omom[0][1][-1];
        #am_omom[1][0]=am_omom[1][0]-am_omom[1][0][-1];
        #am_omom[1][1]=am_omom[1][1]-am_omom[1][1][-1];
        #Generate wave from movement
        if (np.sum(np.abs(exc1))>0):
            r1=am_omom[1][1]/np.conjugate(exc1[1])*(v[1]-Awave[2])+am_omom[0][1]/np.conjugate(exc1[1])*(v[0]-Awave[3]);
            r2=am_omom[0][0]/np.conjugate(exc1[0])*(v[0]-Awave[3])+am_omom[1][0]/np.conjugate(exc1[0])*(v[1]-Awave[2]);
            #r1=self.omega*am_omega[1][1]/np.conjugate(exc1[1])*(v[1]-Awave[2])+self.omega*am_omega[0][1]/np.conjugate(exc1[1])*(v[0]-Awave[3]);
            #r2=self.omega*am_omega[0][0]/np.conjugate(exc1[0])*(v[0]-Awave[3])+self.omega*am_omega[1][0]/np.conjugate(exc1[0])*(v[1]-Awave[2]);
            wave.add_diracWave(-2/np.pi*r1,t,True);
            wave.add_diracWave2(-2/np.pi*r2,t,True);
        
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(np.real(exc1[i])*Awave[0]+np.imag(exc1[i])*Awave[1]);
            ret[i]=FK;#buoyance + FK force
            if i==1:
                ret[i]=ret[i]-self.Area(0)*self.g*self.rho*(z0)+self.Volume(0)*self.rho*self.g;
        Frad=[np.real(np.sum(wave.get_rad2(t,x0)*np.conjugate(exc1[0]))),np.real(np.sum(wave.get_rad(t,x0)*np.conjugate(exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf