#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 08:55:59 2020

@author: heiko
"""

import numpy as np;
import json;
from scipy.integrate import dblquad as int2d;
from scipy.special import jv as bessel;
from scipy.signal import hilbert as KramerKronig;
from scipy.optimize import fsolve;
#from scipy.fftpack import hilbert as KramerKronig;
from scipy.interpolate import interp1d;

pi=np.pi;
#wave=wavefield.wavefield(np.zeros(20),np.linspace(1,20),2)

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
      return self.Calculate_xi(z0, x0, delta0, eta, self.xi)
  def Calculate_xi(self, z0, x0, delta0, eta, xi0):
      #Acos=0.5*Acos
      #Asin=0.5*Asin
      z0=-z0;
      #eta=np.sum(Acos);
      
      #Calculate integartion limits (dependent on buoy position)
      z1=self.z1+z0;
      z2=np.min([self.z2+z0,eta]);
      q=-self.q-(z0)*self.m;
      
      empty=[0,0,0]
      empty2=[xi0*0,xi0*0*0,xi0*0*0]
      #If segment not in water, return zero
      if (z1>z2):
         return [empty,empty2];
      if (z1>eta):
         return [empty,empty2];
      
      
      #Using wheeler strechting to make linear wave theory less linear ;-)
      xis=np.array(xi0*self.depth/(eta+self.depth))#xi dived by water depth (h) and position over the surface eta.
      
      
      #buoyancy force (only heave)
      F_st_h=-self.g*self.rho*(2*pi*self.m*(((self.m*np.power(z2,3)/3+q*np.power(z2,2)/2))-(self.m*np.power(z1,3)/3+q*np.power(z1,2)/2)))
      F_st=[0,F_st_h,0]
      
      #dynamic Froude-Krylov force
      #surge
      c_dy_s=np.array(xis);
      #TODO: implement pitch
      c_dy_p=np.array(xis);
      c_dy_h=np.array(xis);
      c_dy_h=np.array(xis);
      
      if (len(xi0)>1):
          #Avoid numerical errors
          i=np.argmax(np.diff(np.abs(c_dy_h))>0)
          if i>0:
              c_dy_h[i:]=0;
          i=np.argmax(c_dy_h<0)
          if i>0:
              c_dy_h[i:]=0;
              
      c_dy=[c_dy_s.copy(),c_dy_h.copy(),c_dy_p.copy()*0]
          
      
      #print(F_st)
      #print(c_dy)
      return np.array([F_st,F_st]);
  
  def Radius(self,z0):
      if(z0>self.z1) and (z0<=self.z2):      
          r=-self.q+(z0)*self.m;
          return r;
      else:
          return 0;
      
  def max_Radius(self,z0):
      if (z0>self.z2):
          return np.max([-self.q+(self.z1)*self.m,-self.q+(self.z2)*self.m]);
      if(z0>=self.z1) and (z0<=self.z2):      
          r=np.max([np.abs(-self.q+(z0)*self.m),-self.q+(self.z1)*self.m]);
          return r;
      else:
          return 0;
 
  def Area (self, z0):
      #get the area at z0
          r=self.Radius(z0);#r = z0*self.m-q;
          return r*r*pi;
  
  def AreaSurge (self, z0):
          ru=self.Radius(z0);#r = z0*self.m-q;
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
        #self.file = open(idname+"_f.csv", "w")
        #self.file.write("time,wave,buoyancy,FK,rad\r\n");
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
        a=[a[9],a[1],a[2],a[3]]
    def Calculate(self, z0, x0, delta0, eta):
        forces=np.array([[0,0,0],[0,0,0]]);
        
        return [forces[0],0,0,0];
    
    def get_forces(self, t, wave, p, v, a):
        
        #self.file.write(str(t)+","+str(eta)+","+str(Fb)+","+str(FK)+","+str(Frad)+"\r\n");
        return [0,0];
        
    def get_force_lin(self, t, wave, z0, x0, delta0, v, a):
        return 0;        
        
    def Area(self, z0):
        area=0;
        for e in self.elements:
            area = area + e.Area(z0);
        return area;
    
    def AreaProjectedHeave(self, z0):
        area=0;
        res=self.getGeoBox();
        zrmax=res[2];
        rmax=res[3];
        
        if (z0>zrmax):
            return rmax**2*np.pi;
        for e in self.elements:
            area = area + e.Area(z0);
        return area;
        
    def AreaProjectedSurge(self, z0):
        area=0;
        for e in self.elements:
            area = area + e.AreaSurge(z0);
        return area;
    
    def getGeoBox(self):
        z_min=0;
        z_max=0;
        r_max=0;
        z_r_max=0;
        for e in self.elements:
            if  z_min>e.z1:
                z_min=e.z1;
            if  z_max<e.z2:
                z_max=e.z2;
            if e.Radius(e.z1)>r_max:
                r_max=e.Radius(e.z1);
                z_r_max=e.z1;
            if e.Radius(e.z2)>r_max:
                r_max=e.Radius(e.z2);
                z_r_max=e.z2;
        return (z_min,z_max,z_r_max,r_max)
        
    
    def added_mass(self, z0):
        return 0;
        
    
    def max_radius(self,z0):
        return np.max([e.max_Radius(z0) for e in self.elements]);
    
    def Volume(self, z0):
        vol=0;
        for e in self.elements:
            vol = vol + e.Volume(z0);
        return vol;
    
    def clear(self):
        #self.file.close();
        self.elements.clear();
