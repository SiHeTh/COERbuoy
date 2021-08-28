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

pi=np.pi;

class Floater_BEM(Floater):    
    BEMexc=np.array([[],[]]);
    BEMrad=np.array([[],[]]);
    BEMam=np.array([[],[]]);
    
    def __init__ (self, xi, g, depth, CoG, *args):
        super().__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81));
        
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
        eta=np.sum(Awave[0]);
        
        res=self.Calculate(z0, x0, 0*delta0, eta);#Calculate coefficents
        
        ret=[0,0,0];#return array
        exc1 = np.array(res[1]);#Exitation force
        am_omega = np.real(res[2]);#added mass over omega
        am1 = np.array(res[3]);#added mass @ inf
        #am1=[0,Floater.added_mass(self,z0),0];
        
        am_omom=am_omega#*self.omega;
        #print(am_omom[1][1])
        #am_omom[0][0]=am_omom[0][0]-am_omom[0][0][-1];
        #am_omom[0][1]=am_omom[0][1]-am_omom[0][1][-1];
        #am_omom[1][0]=am_omom[1][0]-am_omom[1][0][-1];
        #am_omom[1][1]=am_omom[1][1]-am_omom[1][1][-1];
        #Generate wave from movement
        if (np.sum(np.abs(exc1))>0):
            r1=am_omom[1][1]/np.conjugate(exc1[1])*(v[1]-Awave[2])+am_omom[0][1]/np.conjugate(exc1[1])*(v[0]-Awave[3]);
            r2=am_omom[0][0]/np.conjugate(exc1[0])*(v[0]+Awave[3])+am_omom[1][0]/np.conjugate(exc1[0])*(v[1]-Awave[2]);
            #r1=self.omega*am_omega[1][1]/np.conjugate(exc1[1])*(v[1]-Awave[2])+self.omega*am_omega[0][1]/np.conjugate(exc1[1])*(v[0]-Awave[3]);
            #r2=self.omega*am_omega[0][0]/np.conjugate(exc1[0])*(v[0]-Awave[3])+self.omega*am_omega[1][0]/np.conjugate(exc1[0])*(v[1]-Awave[2]);
            wave.add_diracWave(-2/np.pi*r1,t,True);
            wave.add_diracWave2(-2/np.pi*r2,t,True);
        
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(np.real(exc1[i])*Awave[0]+np.imag(exc1[i])*Awave[1]);
            ret[i]=res[0][i]+FK;#buoyance + FK force
        Frad=[np.real(np.sum(wave.get_rad2(t,x0)*np.conjugate(exc1[0]))),np.real(np.sum(wave.get_rad(t,x0)*np.conjugate(exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf