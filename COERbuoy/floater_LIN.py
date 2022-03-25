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
 
        ret=[0,0,0];#return array
        exc1 = np.array(self.res[1]);#Exitation force
        rad = np.real(self.res[2]);#radiation over omega
        amw = np.imag(self.res[2]);#added mass over omega
        am1 = np.array(self.res[3]);#added mass @ inf
        rad1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1)*(exc1[1]*np.conjugate(exc1[1]))
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;

        rad[1][1]=np.real(rad1);
        #Generate wave from movement
        if (np.sum(np.abs(exc1))>0):
            wave.add_diracWave(-2/np.pi*(rad[1][1]*(v[1])),t,True);
            wave.add_diracWave2(-2/np.pi*(rad[0][0]*(v[0])),t,True);
        
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(exc1[i]),Awave[0])).real;
            ret[i]=FK;#buoyance + FK force
            if i==1:
                ret[i]=ret[i]-self.Area(0)*self.g*self.rho*(z0)+self.eq_force;
        Frad=[np.real(np.sum(wave.get_rad2(t,x0))),np.real(np.sum(wave.get_rad(t,x0))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf


from COERbuoy import LUT_hydro as LUT;
import os;

class Floater_FK(Floater_BEM):
    def __init__ (self, xi, g, depth, CoG, *args):
        super(Floater_BEM,self).__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81),os.path.join(os.path.dirname(args[0]),"BEM"),file_exc="/results/DiffractionForce.tec");
        
    eq_force=None;
    res=None;
    def get_forces(self, t, wave, p, v, a):
        FrKr=np.array([self.omega*0,self.omega*0,self.omega*0]);
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        buoyancy=np.array([0,0,0]);
        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));
        
        if self.eq_force is None:
            self.eq_force=self.Calculate(0,0,0,0)[0][1];
            self.res=self.Calculate(0*z0, 0*x0, 0*delta0, 0);
            self.res[1]=np.array(self.res[1]);
        
        for e in self.elements:
            forces = e.Calculate(-1*z0, x0, delta0, eta);
            FrKr = FrKr-np.array(forces[1]); 
            buoyancy=buoyancy+np.array(forces[0]);
            
        ret=[0,0,0];#return array
        diff1 = np.array(self.res[1]);#Exitation force
        rad = np.real(self.res[2]);#radiation over omega
        amw = np.imag(self.res[2]);#added mass over omega
        am1 = np.array(self.res[3]);#added mass @ inf
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        #Generate wave from movement
        if (np.sum(np.abs(diff1))>0):
            wave.add_diracWave(-2/np.pi*(rad[0][0]*(v[0])),t,True);
            wave.add_diracWave2(-2/np.pi*(rad[1][1]*(v[1])),t,True);
        
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(diff1[i]+FrKr[i]),Awave[0])).real;
            ret[i]=buoyancy[i]+FK;#buoyance + FK force
        Frad=[np.real(np.sum(wave.get_rad(t,x0))),np.real(np.sum(wave.get_rad2(t,x0))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf
class Floater_FKRAD0(Floater_BEM):
    eq_force=None;
    res=None;
    def __init__ (self, xi, g, depth, CoG, *args):
        super(Floater_BEM,self).__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81),os.path.join(os.path.dirname(args[0]),"BEM"),file_exc="/results/DiffractionForce.tec");
        
    eq_force=None;
    res=None;
    def get_forces(self, t, wave, p, v, a):
        FrKr=np.array([self.omega*0,self.omega*0,self.omega*0]);
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        buoyancy=np.array([0,0,0]);
        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));
        
        if self.eq_force is None:
            self.eq_force=self.Calculate(0,0,0,0)[0][1];
            self.res=self.Calculate(0*z0, 0*x0, 0*delta0, 0);
            self.res[1]=np.array(self.res[1]);
        
        for e in self.elements:
            forces = e.Calculate(-1*z0, x0, delta0, eta);
            FrKr = FrKr-np.array(forces[1]); 
            buoyancy=buoyancy+np.array(forces[0]);
        ret=[0,0,0];#return array
        diff1 = np.array(self.res[1]);#Exitation force
        amw = np.imag(self.res[2]);#added mass over omega
        am1 = np.array(self.res[3]);#added mass @ inf
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        
        exc1=diff1;
        exc1[1]=FrKr[1];
        if (np.sum(np.abs(exc1))>0):
            c_1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1);
            r1=(c_1*exc1[1].real+1j*c_1*exc1[1].imag)*v[1];#imag is sometimes zero...
            r2=r1;
            wave.add_diracWave(-2/np.pi*r1,t,True);#-0*((am_omega-self.rad_old)/np.max([t-self.t_old,1e-5]))[1][1],t,True);
            
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(exc1[i]),Awave[0])).real;#np.sum(np.real(exc1[i])*np.real(Awave[0])+np.imag(exc1[i]*np.imag(Awave[0])));
            ret[i]=buoyancy[i]+FK;#buoyance + FK force
            
                
        Frad=[0,np.real(np.sum(m(wave.get_rad(t,x0),exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);

            
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf


class Floater_FKRAD(Floater_BEM):
    eq_force=None;
    res=None;
    def __init__ (self, xi, g, depth, CoG, *args):
        super(Floater_BEM,self).__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81),os.path.join(os.path.dirname(args[0]),"BEM"),file_exc="/results/DiffractionForce.tec");
        
    eq_force=None;
    res=None;
    def get_forces(self, t, wave, p, v, a):
        FrKr=np.array([self.omega*0,self.omega*0,self.omega*0]);
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        buoyancy=np.array([0,0,0]);
        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));
        
        if self.eq_force is None:
            self.eq_force=self.Calculate(0,0,0,0)[0][1];
            self.res=self.Calculate(0*z0, 0*x0, 0*delta0, 0);
            self.res[1]=np.array(self.res[1]);
        
        for e in self.elements:
            forces = e.Calculate(-1*z0, x0, delta0, eta);
            FrKr = FrKr-np.array(forces[1]); 
            buoyancy=buoyancy+np.array(forces[0]);
        ret=[0,0,0];#return array
        diff1 = np.array(self.res[1]);#Exitation force
        amw = np.imag(self.res[2]);#added mass over omega
        am1 = np.array(self.res[3]);#added mass @ inf
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        
        exc1=diff1;
        exc1[1]=FrKr[1]+diff1[1];
        if (np.sum(np.abs(exc1))>0):
            c_1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1);
            r1=(c_1*exc1[1].real+1j*c_1*exc1[1].imag)*v[1];#imag is sometimes zero...
            r2=r1;
            wave.add_diracWave(-2/np.pi*r1,t,True);#-0*((am_omega-self.rad_old)/np.max([t-self.t_old,1e-5]))[1][1],t,True);
            
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(exc1[i]),Awave[0])).real;#np.sum(np.real(exc1[i])*np.real(Awave[0])+np.imag(exc1[i]*np.imag(Awave[0])));
            ret[i]=buoyancy[i]+FK;#buoyance + FK force
            
                
        Frad=[0,np.real(np.sum(m(wave.get_rad(t,x0),exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);

            
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf

class Floater_FKRADAM(Floater_BEM):
    def __init__ (self, xi, g, depth, CoG, *args):
        super(Floater_BEM,self).__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81),os.path.join(os.path.dirname(args[0]),"BEM"),file_exc="/results/DiffractionForce.tec");
        
    eq_force=None;
    res=None;
    def get_forces(self, t, wave, p, v, a):
        FrKr=np.array([self.omega*0,self.omega*0,self.omega*0]);
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        buoyancy=np.array([0,0,0]);
        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));
        
        if self.eq_force is None:
            self.eq_force=self.Calculate(0,0,0,0)[0][1];
            self.res=self.Calculate(0*z0, 0*x0, 0*delta0, 0);
            self.res[1]=np.array(self.res[1]);
            
            
        
        for e in self.elements:
            forces = e.Calculate(-1*z0, x0, delta0, eta);
            FrKr = FrKr-np.array(forces[1]); 
            buoyancy=buoyancy+np.array(forces[0]);
        ret=[0,0,0];#return array
        diff1 = np.array(self.res[1]);#Exitation force
        amw = np.imag(self.res[2]);#added mass over omega
        am1 = np.array([0,self.added_mass(z0,eta),0]);#added mass @ inf
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        
        exc1=diff1;
        exc1[1]=FrKr[1]+diff1[1];
        #Generate wave from movement
        if (np.sum(np.abs(exc1))>0):
            c_1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1);
            r1=(c_1*exc1[1].real+1j*c_1*exc1[1].imag)*v[1];#imag is sometimes zero...
            
            r2=r1;
            wave.add_diracWave(-2/np.pi*r1,t,True);#-0*((am_omega-self.rad_old)/np.max([t-self.t_old,1e-5]))[1][1],t,True);
        am1=[0,self.added_mass(z0,eta),0];
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(exc1[i]),Awave[0])).real;#np.sum(np.real(exc1[i])*np.real(Awave[0])+np.imag(exc1[i]*np.imag(Awave[0])));
            ret[i]=buoyancy[i]+FK;#buoyance + FK force
            
            if i==1:#added mass only implemented for heave
                #v=v[i];
                if abs(v[i]) > 0:
                    dP=(am1[i]-self.added_mass(z0+v[i]/abs(v[i])*0.01,eta))/0.01*v[i];
                    ret[i]=ret[i]-dP;
                
        Frad=[0,np.real(np.sum(m(wave.get_rad(t,x0),exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);

            
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf

class Floater_FKRADAMDIFF(Floater_BEM):
    def __init__ (self, xi, g, depth, CoG, *args):
        super(Floater_BEM,self).__init__(xi, g, depth, CoG, *args);
        LUT.load_LUT(np.sqrt(xi*9.81),os.path.join(os.path.dirname(args[0]),"BEM"),file_exc="/results/DiffractionForce.tec");
        
    eq_force=None;
    res=None;
    def get_forces(self, t, wave, p, v, a):
        FrKr=np.array([self.omega*0,self.omega*0,self.omega*0]);
        z0=p[1];
        x0=p[0];
        delta0=p[2];
        buoyancy=np.array([0,0,0]);
        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));
        
        if self.eq_force is None:
            self.eq_force=self.Calculate(0,0,0,0)[0][1];
            self.res=self.Calculate(0*z0, 0*x0, 0*delta0, 0);
            self.res[1]=np.array(self.res[1]);
            
            
            
        self.res1=self.Calculate(z0, 0*x0, 0*delta0, eta);
        self.res[1]=np.array(self.res1[1]);
        self.res[3]=np.array(self.res1[3]);
        
        for e in self.elements:
            forces = e.Calculate(-1*z0, x0, delta0, eta);
            FrKr = FrKr-np.array(forces[1]); 
            buoyancy=buoyancy+np.array(forces[0]);
        ret=[0,0,0];#return array
        diff1 = np.array(self.res[1]);#Exitation force
        amw = np.imag(self.res[2]);#added mass over omega
        am1 = np.array(self.res[3]);#added mass @ inf
        def m(a,b):
            return a.real*b.real+a.imag*b.imag;
        
        exc1=diff1;
        exc1[1]=FrKr[1]+diff1[1];
        #Generate wave from movement
        if (np.sum(np.abs(exc1))>0):
            c_1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1);
            r1=(c_1*exc1[1].real+1j*c_1*exc1[1].imag)*v[1];#imag is sometimes zero...
            
            r2=r1;
            wave.add_diracWave(-2/np.pi*r1,t,True);#-0*((am_omega-self.rad_old)/np.max([t-self.t_old,1e-5]))[1][1],t,True);
        am1=[0,self.added_mass(z0,eta),0];
        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(exc1[i]),Awave[0])).real;#np.sum(np.real(exc1[i])*np.real(Awave[0])+np.imag(exc1[i]*np.imag(Awave[0])));
            ret[i]=buoyancy[i]+FK;#buoyance + FK force
            
            if i==1:#added mass only implemented for heave
                #v=v[i];
                if abs(v[i]) > 0:
                    dP=(am1[i]-self.added_mass(z0+v[i]/abs(v[i])*0.01,eta))/0.01*v[i];
                    ret[i]=ret[i]-dP;
                
        Frad=[0,np.real(np.sum(m(wave.get_rad(t,x0),exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);

            
        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf
