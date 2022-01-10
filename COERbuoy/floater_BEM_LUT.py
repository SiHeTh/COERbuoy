#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Floater_BEM_LUT.py - Calculating hydrostatic- and dyanmic forces based on BEM data
# Initally designed to handle three degree of freedom (3D) heave, surge and sway, but only two
# (heave, surge) are used.
# 2020/2021 COER Laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie

"""
Created on Fri Aug 28 08:55:59 2020

@author: heiko
"""

#Dependencies
from COERbuoy.floater import Floater;
import numpy as np;
from COERbuoy import LUT_hydro as LUT;
import os;

pi=np.pi;

class Floater_BEM(Floater):    
    BEMexc=np.array([[],[]]);
    BEMrad=np.array([[],[]]);
    BEMam=np.array([[],[]]);
    rad_set=False;
    rad_old=None;
    p_old=np.NaN;
    t_old=-0.01;
    
    def __init__ (self, xi, g, depth, CoG, *args):
        super().__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81),os.path.join(os.path.dirname(args[0]),"BEM"));
        
    #calculate hydrodnymaic parmeters from heave, surge, pitch, surface elevation
    def Calculate(self, z0, x0, delta0, eta):
        fb=np.array([0,0,0]);
        
        z0=-z0;
        
        #Get buyoancy force for each element
        for e in self.elements:
            fb = fb + e.Calculate(z0, x0, delta0, eta)[0];
        if np.sum(np.abs(fb))==0:#out of water
            en=[self.omega*0,self.omega*0,self.omega*0];
            return [fb,en,np.array([en,en,en]),[0,0,0]]
        
        draft=eta+z0;#current submergence
        #print(draft)
        res = LUT.get_fromLUT(draft,0);
        exc1 = res[0];
        #exc1 = res[0]*np.cos(res[1])+1j*res[0]*np.sin(res[1]);
        exc1 = res[0]*np.exp(1j*res[1]);

        rad1 = [[res[3][0][0]+1j*res[2][0][0],res[3][0][1]+1j*res[2][0][1],res[3][0][2]+1j*res[2][0][2]],
                [res[3][1][0]+1j*res[2][1][0],res[3][1][1]+1j*res[2][1][1],res[3][1][2]+1j*res[2][1][2]],
                [res[3][2][0]+1j*res[2][2][0],res[3][2][1]+1j*res[2][2][1],res[3][2][2]+1j*res[2][2][2]]];
        am1 = [res[4][0],res[4][1],res[4][2]];

        return [fb,exc1,rad1,am1];#F_buoynacy, F_excitation, F_radiation, F_added_mass
    
    #calculate hydro forces from time, wave, heave, surge, pitch, velocity, acceleration
    def get_forces(self, t, wave, p, v, a):
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        
        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));
        
        res=self.Calculate(z0, x0, 0*delta0, eta);#Calculate coefficents
        dam=np.array([0,0,0]);#(np.array(res[3])-np.array(self.Calculate(z0+0.01, x0, 0*delta0, eta)[3]))/0.01;#Calculate coefficents
        if self.rad_old is None:
            self.rad_old=np.real(res[2])/[np.array(res[1]),np.array(res[1]),np.array(res[1])];
        ret=[0,0,0];#return array
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        #exc1 = np.conjugate(np.array(res[1]));#Exitation force
        exc1 = np.array(res[1]);#Exitation force
        am_omega = np.real(res[2]);#added mass over omega
        am1 = np.array(res[3]);#added mass @ inf
        
      
        am_omom=am_omega/np.array([exc1,exc1,exc1]);#np.matmul(am_omega,vv)+np.diag((am_omega-self.rad_old)/np.max([t-self.t_old,1e-5]));
        #print([v[1],np.sum(np.real(Awave[1]))])
        #Generate wave from movement
        #print([np.sum(np.real(Awave[0])),np.sum(np.real(Awave[1]))])
        if (np.sum(np.abs(exc1))>0):
            #dx=v[1]*0.01;
            #res2=self.Calculate(z0, x0+dx, 0*delta0, eta);
            r1=1*am_omom[1][1]*(v[1])+am_omom[0][1]*(v[0]);
            r2=1*am_omom[0][0]*(v[0])+am_omom[1][0]*(v[1]);
            #r1=1*am_omom[1][1];
            #r2=1*am_omom[0][0]*(v[0]-np.sum(np.imag(Awave[1])))+am_omom[1][0]*(v[1]-np.sum(np.imag(Awave[1])));
            #print([v[1],np.sum(Awave[2])])
            
            #wave.add_diracWave(-2/np.pi*(am_omom[1][1]*(v[1]-0*np.sum(Awave[2]))),t,True);
            wave.add_diracWave(-2/np.pi*r1,t,True);#-0*((am_omega-self.rad_old)/np.max([t-self.t_old,1e-5]))[1][1],t,True);
            wave.add_diracWave2(-2/np.pi*r2,t,True);
            self.rad_old=am_omega/[exc1,exc1,exc1];
            self.t_old=t;
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(exc1[i],Awave[0])).real;#np.sum(np.real(exc1[i])*np.real(Awave[0])+np.imag(exc1[i]*np.imag(Awave[0])));
            ret[i]=res[0][i]+FK;#buoyance + FK force
        Frad=[np.real(np.sum(wave.get_rad2(t,x0)*np.abs(exc1[0]))),np.real(np.sum(wave.get_rad(t,x0)*(exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);
        
        self.t_old=t;
        self.rad_old=am_omega;
        self.z0_old=z0;
            
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf
