#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 09:14:13 2020

@author: heiko
"""

import json;
import numpy as np;
from scipy.optimize import minimize_scalar;
import COERbuoy.utils as utils;
#from Floater import Floater as fl;
from floater_BEM_LUT import Floater_BEM as fl;
import os;


g=9.81;
rho=1000;
class WEC():
    omega_cut_off=20;
    states=10;
    am_old=0;
    amlow_old=0;
    v_old=0;
    t_old=-0.0001;
    file=os.path.join(utils.wec_dir,"floater.txt");
    acc=0;
    def __init__(self):
        self.force_sensor=0;
        self.cD=0.09;
        self.damping=0;
        return;
        
    def load_buoy(self,floater_class,xi,depth,cog):
        self.buoy=floater_class(xi,g,depth,cog,self.file);
        self.mass=self.buoy.Volume(0)*1000;
        #self.mass=self.buoy.Volume(0)*1000;
    
    def load_param(self):
        return 0;
    
    def Calc_fs_spring(self,x1,x2):
      return 0;
            
    def pto_mdc (self):
        m=0;
        d=0;
        c=0;
        return [m,d,c];
    
    def Calc_drag(self,x,dx):
        #print([x,dx,self.buoy.max_radius(-x)**2*3.14,-0.5*self.cD*rho*dx*np.abs(dx)*self.buoy.max_radius(-x)**2*3.14])
        return -0.5*self.cD*rho*dx*np.abs(dx)*self.buoy.max_radius(-x)**2*3.14;
    def get_translator_speed(self,x):
        return x[1]-x[7];
    def get_translator_position(self,x):
        return x[0]-x[6];
    def get_surge(self,x):
        return (x[0])*np.sin(x[2]);
    def get_force(self,x):
        #print(self.force_sensor)
        return 0;#self.force_sensor;
    def Calc(self,t,wave,x,PTO_force,valve,ulast):
      
      #m0=self.buoy.added_mass(x[0]);
      alpha=x[2];
      stroke=x[0];
      m_rot=np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
      #m_dotrot=np.array([[-np.sin(dalpha), -np.cos(dalpha)], [np.cos(dalpha), -np.sin(dalpha)]])
      
      heave=stroke;#(stroke+self.l)*np.cos(alpha)-self.l;
      surge=self.get_surge(x);
      
      
      #F_FKS = self.buoy.Calculate(-x[0],x[2],x[4],Asin,Acos)
      #F_FKS_h=F_FKS[0][1]+np.sum(np.real(F_FKS[1][1])*Acos-Asin*np.imag(F_FKS[1][1]))
      #radForce=np.zeros(3);
      
      f_hy = self.buoy.get_forces(t,wave,[surge*0,heave,alpha],[0,x[1],0],self.acc)
      
      F_h=f_hy[0];
      m0=f_hy[1];
     # print(heave,f_hy[0])
      Fdrag=0#self.Calc_drag(x[0],x[1]);
      
      
      dx=np.zeros(10);
      #m0[1]=self.buoy.Volume(heave*-1)*1000*0.5;
      #print(t,heave,m0[1])
      dP=0;
      if False:#(t-self.t_old)>0:
          print(t-self.t_old)
          #param=self.buoy.Calculate(heave, 0, 0, 0);
          if False:#x[1]>0:
              dP=0#0.5*(np.imag(param[2][1][0])-self.amlow_old)/(t-self.t_old)*x[1];
          else:
              self.am_old=self.buoy.get_forces(t,wave,heave+x[1]*0.01,0*surge,alpha,[0,x[1],0],self.acc)[1][1];
              dP=-0.5*(m0[1]-self.am_old)/(0.01)*x[1];
              #dP=-1*(m0[1]-self.am_old)/(t-self.t_old)*x[1];
          
              self.am_old=m0[1];
              #print(np.imag(param[2][1][0]))
              #self.amlow_old=np.imag(param[2][1][0]);
              self.t_old=t;
              print([self.am_old,m0[1],dP/(x[1]+0.0001)])
      
      self.v_old=x[1];
      F_sum_x=F_h[1]-self.mass*g+Fdrag-dP;
      
      F_sum_y=F_h[0];
      
      
      dx[1]=(F_sum_x)/(self.mass+m0[1]);
      #print(self.mass)
      dx[0]=x[1];
      dx[3]=0;#F_sum_y/(self.mass+m0[0]);
      dx[2]=x[3]*0;
      dx[5]=0;
      dx[4]=0;
      dx[7]=0;
      dx[6]=0;
      
      dx[8]=0;
      self.acc=dx[1]
      return dx;
      
      
    def get_data(self):
        return 
        
    def release(self):
        self.buoy.clear();
        del self.buoy;
