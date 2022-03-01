#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dynamics_coerbuoy.py - The WEC model for the COERbuoy1 WEC
# 2020/2021 COER Laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie

"""
Created on Fri Oct 30 09:14:13 2020

@author: Simon H. Thomas
"""

#The WEC class specified in this model defines the behaviour of the COERbuoy1
#module. It takes the hydrostatic and -dynamic data calculated by the Floater
#class and adds the machinery and mooring forces

#Dependencies
import json;
import numpy as np;#BSD 3-clause license
#from .floater_BEM_LUT import Floater_BEM as fl;
import os;
import COERbuoy.utils as utils;
#floaterfile=pkg_resources.resource_filename(__name__,"floater.txt");

x_lim=3.5; #stroke limit
g=9.81; #gravity acceleration
rho=1000; #density of water


class WEC():
    #number of states

    omega_cut_off=3.2; #highest frequency used for calculation    
    
    states=10;
    try:
        file=os.path.join(utils.wec_dir,"floater.txt");
    except(AttributeError): #if directory is not set in utils, use own directory
        file=os.path.join(os.path.dirname(__file__),"floater.txt");
    acc=0;
    force_sensor=0;
    
    def __init__(self):
        self.force_sensor=0;
        self.Ar=[[],[],[]];
        utils.get();
        
        return;
        
    def load_buoy(self,floater_class,xi,depth,cog):
        self.buoy=floater_class(xi,g,depth,cog,self.file);
        self.mass=self.buoy.Volume(0)*1000;
    
    def load_param(self):
        with open(self.file) as file:
            data=json.load(file);
            #generator damping (only for testing)
            self.damping =data.get("PTO_damping",100000);
            
    #Get linearised mass, damping and spring coefficent        
    def pto_mdc (self,z):
        m=(self.mass);
        d=0;#self.d_add#;+self.Calc_drag(0,1);
        c=0#-self.c_fs;
        return [m,d,c];
    
    #negative spring force + pre-tension
    def Calc_fs_spring(self,x1,x2):
      return self.c_fs*np.arctan((x1-x2)/self.l_fs)+self.Calc_fs_pretension();
    #pre-tension (constant)
    def Calc_fs_pretension(self):
      return -self.mass*(1-self.mb)*g
    #drag force heave
    def Calc_drag(self,x,dx):
      return -self.cD*rho*dx*np.abs(dx)*self.buoy.AreaProjectedHeave(-x)
    #drag force surge   
    def Calc_drag_surge(self,x,dx):
      return -self.cDs*rho*dx*np.abs(dx)*self.buoy.AreaProjectedSurge(-x)
    #additional (machinery) damping force
    def Calc_dadd(self,dx1,dx2):
      return -self.d_add*(dx1);
    
    def get_surge(self,x):
        return -1*(x[0])*np.sin(x[2]);
    def get_translator_speed(self,x):
        return x[1];
    def get_translator_speed_surge(self,x):
        return x[3];
    def get_translator_position(self,x):
        return x[0];
    def get_force(self,x):
        return self.force_sensor;
    
    #calculate generator force and absorbed energy.
    def calc_PTO(self,F_pto,dx):
        E=self.gen_clambda*dx
        X=self.gen_cL*dx
        
        if (F_pto == 0):
            return [0, 0];
        
        if (F_pto*dx>=0):
            dx=dx+1;
            
        p=E**2/(F_pto*dx);
        q=X**2;
        R=self.gen_Rc;
        if (p**2-q>0):
            R=-(p+np.sqrt(p**2-q))
        else:
            R=-p+0;
        R=R;
        Rl=R-self.gen_Rc;
        
        
        I=E/np.sqrt(R**2+X**2);
        
        gamma=np.abs(self.gen_I_lim/I);
        if gamma < 1:
            D=0.5*np.pi*(np.arctan(gamma/(1-gamma))+np.sqrt(1-gamma**2));
            X=(self.gen_cL*D)*dx;
            I=E/np.sqrt(R**2+X**2);
        
        Pabs=-Rl*I**2;
        return [-E**2*R/(R**2+X**2)*1/(dx+0.0001), Pabs];
      
    #main calculation function, called by ODE solver from main programm
    def Calc(self,t,wave,x,PTO_force,brake,ulast):
      #State vector:
      #x[0] - stroke position
      #x[1] - stroke speed
      #x[2] - anchor joint pitch angle
      #x[3] - anchor joint pitch angular velocity
      #x[4]-x[6] - not used
      #x[7] - time integrated generator force
      #x[8] - absorbed energy
      
      stroke=x[0];
            
      heave=x[0];
      surge=self.get_surge(x);
      surge_v= 0;
      heave_v= x[1];
      
      #Drag forces
      #Drag force in heave
      #Fdrag=self.Calc_drag(heave,heave_v);
      #Drag force in sway
      #Fdrags=self.Calc_drag_surge(surge,surge_v);
      
      
      f_hy = self.buoy.get_forces(t,wave,[0,heave,0],[0,x[1],0],self.acc)
      mah=np.real(f_hy[1][1])
      f_exc=f_hy[0][1]-self.mass*9.81;
      
      #Test with reactive control for 10 s period
      m_c=-406458.62;
      m_d=21864.64;
      #m_c=(self.mass+mah)*(6.28/10)**2-self.buoy.Area(0)*1000*9.81;
      #PTO_force=-1*float(m_d)*x[1]-m_c*x[0];
      #PTO_force=-1*float(274637.452)*x[1];
      #print([PTO_force,PTO_force2,brake])
      #Generator force
      F_gen=PTO_force;
      Pabs=PTO_force*x[1];
          
      #Calculate all inertia (physical mass+added mass)
      mass_sum_floater=(self.mass+mah);
      #print(heave,f_exc)
      if (x[1]>=0):
          brake=-1*abs(brake);
      F_sum_floater=f_exc;
      F_sum_floater=np.sum(np.real(F_sum_floater))+F_gen;
      self.force_sensor=F_sum_floater;
      #print([np.sum(np.real(f_exc)),brake])
      dx=np.zeros(10);
      #Fill dx vector
      if np.abs(F_sum_floater)<abs(brake):
          dx[1]=-x[1]*10;
          #print(x[1])
          dx[0]=x[1];
          F_gen=0;
          Pabs=0;
      else:
          dx[1]=(F_sum_floater+brake)/mass_sum_floater;
          self.force_sensor=F_sum_floater+brake;
          dx[0]=x[1];
      #print(PTO_force/(x[0]+1e-6),Pabs,x[8])
          
      
     
      
      dx[3]=0;
          
      dx[2]=x[3];
      dx[5]=0;
      dx[4]=0;
      dx[7]=0;
      dx[6]=0;
      dx[7]=-F_gen;
      self.acc=[0,dx[1],0];#*np.cos(x[2]);
      dx[8]=-Pabs;
      return dx;
      
      
    def get_data(self):
        return 
        
    def release(self):
        self.buoy.clear();
        del self.buoy;

#Testing  
if __name__=="__main__":
    #Testing if any function throws an error
    omega=np.round(np.linspace(0.1,1.5*2+0.1,32),10);
    xi=omega*omega/9.81;
    w=WEC();
    w.load_buoy(xi,100,0);
    w.load_param();
    res=w.buoy.Calculate(0,0,0,0);
    print(w.pto_mdc());
