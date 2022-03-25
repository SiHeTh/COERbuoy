#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Floater_BEM_LUT.py - Calculating body-exact hydrostatic- and dyanmic forces based on BEM data calculated at different body positions
# Initally designed to handle three degree of freedom (3D) heave, surge and sway, but only two
# (heave, surge) are used.
# 2020/2021/2022 COER laboratory, Maynooth University
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
    BEMexc=np.array([[],[]]);#initialising hydrodynamic look-up-tables
    BEMrad=np.array([[],[]]);
    BEMam=np.array([[],[]]);
    t_old=-0.01; #time variable that is set to log only in intervals
    
    def __init__ (self, xi, g, depth, CoG, *args):
        super().__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81),os.path.join(os.path.dirname(args[0]),"BEM"));
        
    #calculate hydrodynamic parmeters from heave, surge, pitch, surface elevation
    def Calculate(self, z0, x0, delta0, eta):
        #1) buoyancy calculation
        fb=np.array([0,0,0]);#buoyancy array (surge, heave, pitch)
        z0=-z0;#For some reasons the sign is wrong; complain with the coder...
        
        #Get buoyancy force for each element
        for e in self.elements:
            fb = fb + e.Calculate(z0, x0, delta0, eta)[0];
        if np.sum(np.abs(fb))==0:#out of water
            en=[self.omega*0,self.omega*0,self.omega*0];
            return [fb,en,np.array([en,en,en]),[0,0,0]]
        
        #2) get hydrodynamic parameters from LUT
        draft=eta+z0;#get current submergence
        res = LUT.get_fromLUT(draft,0);
        exc1 = res[0]*np.cos(res[1])+1j*res[0]*np.sin(res[1]);
       
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
        dam=np.array([0,0,0]);
        ret=[0,0,0];#return array
        
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        
        exc1 = np.array(res[1]);#Exitation force
        rad1 = np.real(res[2]);#Radiation coefficents
        am1 = np.array(res[3]);#added mass @ inf
        
      
        #Calculate the instantanious radiated wave caused by the body's velocity
        if (np.sum(np.abs(exc1))>0):
            #Using Haskind relation to get radiation from excitation
            c_1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1);
            r0=(c_1*exc1[0].real+1j*c_1*exc1[0].imag)*v[0];#Using Hashkind to get radiation (numerically best solution)
            r1=(c_1*exc1[1].real+1j*c_1*exc1[1].imag)*v[1];
            #axisymetric devices do not really have cross terms for radiation
            wave.add_diracWave(-2/np.pi*r0,t,True);
            wave.add_diracWave2(-2/np.pi*r1,t,True);
            #TODO: Add pitch
            
        #Calculate hydro forces for each DOF
        dP=0;
        FK=0;
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(exc1[i]),Awave[0])).real;#np.sum(np.real(exc1[i])*np.real(Awave[0])+np.imag(exc1[i]*np.imag(Awave[0])));
            ret[i]=res[0][i]+FK;#buoyance + FK force
            
            if i==1:#added mass only implemented for heave
                v=v[i];
                if abs(v) > 0:#slamming force
                    dP=(am1[i]-self.Calculate(z0+v/abs(v)*0.01,x0,0*delta0,eta)[3][i])/0.01*v;
                    ret[i]=ret[i]+dP;
                if self.file and (t-self.t_old)>0.09:#code to log the forces (for debugging and in-detail analysis)
                    self.file.write(str(t)+","+str(eta)+","+str(abs(res[0][1]))+","+str(abs(FK))+","+str(abs(np.real(np.sum(m(wave.get_rad2(t,x0),exc1[1])))))+","+str(abs(dP))+"\r\n");
                    self.t_old=t;
                
        Frad=[np.real(np.sum(m(wave.get_rad(t,x0),exc1[0]))),np.real(np.sum(m(wave.get_rad2(t,x0),exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);
            
        
        self.z0_old=z0;
            
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf
