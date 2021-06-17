#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 08:55:59 2020

@author: heiko
"""

import numpy as np;
import json;
from scipy.interpolate import interp1d;

pi=np.pi;

class Cone:
  m = -1;
  q = -1;
  z1 = -1;
  z2 = 0;
  I = 0;
  g=9.81;
  zCoG = -1;
  depth=30;
  
  def __init__ (self, z1, r1, z2, r2, zCoG, xi, g, h):
      #Frequency parameters
      self.xi=xi;
      
      #geometric parameters
      self.m=(r1-r2)/(z1-z2);
      self.z1=z1;
      self.z2=z2;
      self.g=g;
      self.rho=1000;
      self.q=-r1+self.m*self.z1;
      self.omega=np.sqrt(self.xi*g);
      
      self.zCoG=zCoG;
      self.z0=0;
      self.x0=0;
      self.delta0=0;
      self.depth=h;
  
  def Calculate(self, z0, x0, delta0, eta):
      #Static FK force calculation
      
      z0=-z0;
      
      #Calculate integartion limits (dependent on buoy position)
      z1=self.z1+z0;
      z2=np.min([self.z2+z0,eta]);
      q=-self.q-(z0)*self.m;
      
      empty=[0,0,0]
      empty2=[self.omega*0,self.omega*0,self.omega*0]
      #If segment not in water, return zero
      if (z1>z2):
         return [empty,empty2];
      if (z1>eta):
         return [empty,empty2];
      
      
      #buoyancy force (only heave)
      F_st_h=-self.g*self.rho*(2*pi*self.m*(((self.m*np.power(z2,3)/3+q*np.power(z2,2)/2))-(self.m*np.power(z1,3)/3+q*np.power(z1,2)/2)))
      F_st=[0,F_st_h,0]
          
      return np.array([F_st,[0,0,0]]);
  
  def Radius(self,z0):
      if(z0>self.z1) and (z0<=self.z2):      
          r=-self.q+(z0)*self.m;
          return r;
      else:
          return 0;
      
  def max_Radius(self,z0):
      if (z0>self.z2):
          return np.max([self.z1,self.z2]);
      if(z0>self.z1) and (z0<=self.z2):      
          r=np.max([-self.q+(z0)*self.m,self.z1]);
          return r;
      else:
          return 0;
 
  def Area (self, z0):
          r=self.Radius(z0);
          return r*r*pi;
  
  def AreaSurge (self, z0):
          ru=self.Radius(z0);
          rb=self.Radius(self.z1);
          h=np.max([z0-self.z1,0]);
          return (ru+rb)*h;
      
  def Volume(self, z):
      z1l=np.min([self.z1,z])
      z2l=np.min([self.z2,z])
      z1l1=(z1l-self.q/self.m)
      z2l2=(z2l-self.q/self.m)
      r2=z2l2*self.m;
      r1=z1l1*self.m;
      return 1/3*pi*(r2*r2*z2l2-r1*r1*z1l1);


idname="test.csv";
class Floater:
    volume=0;
    mode=0;
    g=9.81;
    d=300;
    Cog=0;
    rho=1000;
    xi=np.array([]);
    
    def __init__ (self, xi, g, depth, CoG, *args):
        self.xi=xi;
        self.g=g;
        self.omega=np.sqrt(self.xi*self.g)
        self.elements=[];
        self.d=depth;
        self.CoG=CoG;
        if len(args)>0:
            with open(args[0]) as file:
                geo=json.load(file);
                for g in geo["geo"]:
                    if g["type"] == "cone":
                        self.addCone(g["coord"][0],g["coord"][1],g["coord"][2],g["coord"][3])
        
    def addCone(self, z1, r1, z2, r2):
        self.elements.append(Cone(z1, r1, z2, r2,self.CoG,self.xi,self.g,self.d))
        self.volume=self.volume+self.elements[-1].Volume(z2);
    
    def set_mode(self, mode):
        self.mode=mode;
        
    def Calc_CoG(self):
        mi=0;
        ms=0;
        rge=np.linspace(np.min([e.z1 for e in self.elements]),np.max([e.z2 for e in self.elements]),10);

        for i in rge:
            ai=self.Area(i);
            mi=mi+ai*i;
            ms=ms+ai;
        return mi/ms;
        
    def get_parameters(self, z0, x0, delta0, eta):
        a=self.Calculate();
        a=[a[9],a[1],a[2]+1j*a[3],a[4]]
    def Calculate(self, z0, x0, delta0, eta):
        #just a placeholder
        forces=np.array([[0,0,0],[0,0,0]]);
        
        m0=[0,0,0];
        
        for e in self.elements:
            forces = forces + e.Calculate(z0, x0, delta0, eta);
        if np.sum(np.abs(forces[0]))==0:
            return [forces[0],forces[1],forces[1],forces[1],forces[0]]
        
        c_dy=forces[1];
        c_ex=c_dy.copy()*0.8;
      
        #Initialization of parameters
        c_added_mass=[np.zeros(len(self.omega))+m0[0],np.zeros(len(self.omega))+m0[1],np.zeros(len(self.omega))+m0[2]];
        c_rad=[np.zeros(len(self.omega)),np.zeros(len(self.omega)),np.zeros(len(self.omega))];
          
        c_rad_cmplx=np.real(c_rad)+1j*np.real(c_rad);
        return [forces[0],c_ex,c_rad_cmplx,m0];
    
    def get_forces(self, t, wave, z0, x0, delta0, v, a):
        Awave=wave.get(t,x0);
        eta=np.sum(Awave[0]);
        res=self.Calculate(z0, x0, delta0, eta);
        ret=[0,0,0];
        if (np.sum(np.abs(res[1][1]))>0):
            wave.add_diracWave(-2/np.pi*(res[2][1])/((res[1][1]))*(v-Awave[2]),t,True);
            
        for i in range(len(res[0])):
            Fb=res[0][1]-self.Volume(0)*self.rho*self.g;
            FK=np.sum(np.real(res[1][i])*Awave[0]+np.imag(res[1][i])*Awave[1])
            ret[i]=res[0][i]+FK;
            FK=np.sum(np.real(res[1][1])*Awave[0]+np.imag(res[1][1])*Awave[1])
        Frad=np.real(np.sum(wave.get_rad(t,x0)*(res[1][1])));
        ret[1]=ret[1]+Frad;
        
        return [np.real(ret),res[4]];
        
    def get_force_lin(self, t, wave, z0, x0, delta0, v, a):
        ret=self.get_forces(t,wave,0,0,0,v,a);
        ret[1]=ret[1]-self.Area(0)*1000*z0*9.81;
        return ret;        
        
    def Area(self, z0):
        area=0;
        for e in self.elements:
            area = area + e.Area(z0);
        return area;
        
    def AreaSurge(self, z0):
        area=0;
        for e in self.elements:
            area = area + e.AreaSurge(z0);
        return area;
    
    def added_mass(self, z0):
        return self.Volume(z0)*1000*0.5;
    
    def max_radius(self,z0):
        return np.max([e.max_Radius(z0) for e in self.elements]);
    
    def Volume(self, z0):
        vol=0;
        for e in self.elements:
            vol = vol + e.Volume(z0);
        return vol;
    
    def clear(self):
        self.elements.clear();
