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
    file=os.path.join(utils.wec_dir,"floater.txt");#floaterfile;#Path(__file__).parent()/"floater.txt";
    acc=0;
    force_sensor=0;
    
    def __init__(self):
        self.force_sensor=0;
        self.Ar=[[],[],[]];
        return;
        
    def load_buoy(self,floater_class,xi,depth,cog):
        self.buoy=floater_class(xi,g,depth,cog,self.file);
        self.mass=self.buoy.Volume(0)*1000;
    
    def load_param(self):
        with open(self.file) as file:
            data=json.load(file);
            #generator damping (only for testing)
            self.damping =data.get("PTO_damping",100000);
            #negatie spring stiffnes
            self.c_fs=data.get("negative_spring_stiffness",20*1000*9.81);
            #length of negative spring
            self.l_fs=data.get("negative_spring_length",2.0);
            #Drag coefficent heave
            self.cD=data.get("viscous_drag_heave",0.2);
            #Drag coefficent surge
            self.cDs=data.get("viscous_drag_surge",0.5);
            #Percentage mass floater
            self.mb=data.get("mass_percent_floater",0.25);
            #added mass floater
            self.m0=data.get("added_mass",self.mass);
            #generator efficency
            self.losses_coil=1/data.get("eff_generator",0.8);
            #mooring line stiffness
            self.ml_c=data.get("mooring_stiffness",4e5);
            #mooring line damping
            self.ml_d=data.get("mooring_damping",4000);
            #static friction
            self.fr_s=data.get("friction_force_static",0);
            #kinetic friction
            self.fr_d=data.get("friction_force_kinetic",0);
            #friction damping
            self.d_add=data.get("friction_damping",0);
            #maximum braking power
            self.P_mbreak=data.get("brake_Pmax",1000);
            #generator copper resistance
            self.gen_Rc=data.get("generator_Rc",0.01);
            #generator inductance coefficent
            self.gen_cL=data.get("generator_c_L",1.57);
            #generator flux coefficent
            self.gen_clambda=data.get("generator_c_lambda",365*3);
            #generator magnetic saturation current
            self.gen_I_lim=data.get("generator_I_s",100000);
            #mooring line length
            self.l=data.get("l_mooring",30);
            
    #Get linearised mass, damping and spring coefficent        
    def pto_mdc (self):
        m=(self.mass);
        d=self.d_add#;+self.Calc_drag(0,1);
        c=-self.c_fs;
        return [m,d,c];
    
    #negative spring force + pre-tension
    def Calc_fs_spring(self,x1,x2):
      return self.c_fs*np.arctan((x1-x2)/self.l_fs)+self.Calc_fs_pretension();
    #pre-tension (constant)
    def Calc_fs_pretension(self):
      return -self.mass*(1-self.mb)*g
    #drag force heave
    def Calc_drag(self,x,dx):
      return -self.cD*rho*dx*np.abs(dx)*self.buoy.Area(-x)
    #drag force surge   
    def Calc_drag_surge(self,x,dx):
      return -self.cDs*rho*dx*np.abs(dx)*self.buoy.AreaSurge(-x)
    #additional (machinery) damping force
    def Calc_dadd(self,dx1,dx2):
      return -self.d_add*(dx1);
    
    def get_surge(self,x):
        return -1*(x[0]+self.l)*np.sin(x[2]);
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
            
        p=E**2/(F_pto*dx)*0.5;
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
      
      alpha=x[2];
      stroke=x[0];
      
      if PTO_force>0.1:
          PTO_force=PTO_force+1;
      
      #Rot. matrix: Global corrdinates (hydro-forces) into body coordinates (PTO forces)
      m_rot=np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
      
      heave=(stroke+self.l)*np.cos(alpha)-self.l;
      surge=self.get_surge(x);
      surge_v= x[1]*np.cos(x[2])# + x[3]*np.sin(x[2]);
      heave_v=-x[1]*np.sin(x[2])# + x[3]*np.cos(x[2]);
      
      if x[1]>0.1:
          brake=np.max([brake,self.P_mbreak/x[1]]);
      else:
          brake=np.abs(brake);
      
      #Drag forces
      #Drag force in heave
      Fdrag=self.Calc_drag(heave,heave_v);
      #Drag force in sway
      Fdrags=self.Calc_drag_surge(surge,surge_v);
      
      
      f_hy = self.buoy.get_forces(t,wave,[surge,heave,alpha],[x[1]*np.cos(x[2])+x[3]*np.sin(x[2]),-x[1]*np.sin(x[2])+x[3]*np.cos(x[2]),0],self.acc)
      f_hy[0][0]=f_hy[0][0]-Fdrags;
      f_hy[0][1]=f_hy[0][1]-Fdrag-self.mb*g*self.mass;
      F_radax = np.matmul(m_rot,f_hy[0][:2]);
      
      
      #PTO-forces:
      #negative spring, including pre-tension
      WS=self.Calc_fs_spring(x[0],0);
      #Machinery damping
      Fd_add=self.Calc_dadd(x[1],0);
      #Generator force
      [F_gen, Pabs]=self.calc_PTO(PTO_force,x[1]);
      #Do not consider power absorbed above constraints
      if (abs(x[0]>x_lim)):
          Pabs=0;
          
      #Calculate all inertia (physical mass+added mass)
      mah=np.real(f_hy[1][0]*np.sin(np.abs(alpha))+f_hy[1][1]*np.cos(alpha))
      mas=np.real(f_hy[1][0]*np.cos(np.abs(alpha))+f_hy[1][1]*np.sin(alpha))
      mass_sum_floater=(self.mass*self.mb+np.real(mah));
      
      F_sum_floater=F_radax[1]+F_gen+WS+Fd_add+brake;
      F_sum_floater=np.sum(np.real(F_sum_floater));
      self.force_sensor=F_sum_floater;
      
      
      #Fill dx vector
      
      dx=np.zeros(10);
      #if machinery or static friction too high: PTO-stuck, no speed, no force
      if np.abs(F_sum_floater-brake)<brake or (np.abs(F_sum_floater)<self.fr_s and np.abs(x[1])<0.01):
          dx[1]=0;
          dx[0]=0;
          F_gen=0;
          Pabs=0;
      else:
          dx[1]=(F_sum_floater)/mass_sum_floater;
          dx[0]=x[1];
          
      
      dx[3]=((self.l+stroke)*(np.sum(F_radax[0])))/((self.mass+mas)*(self.l+stroke)**2+np.real(f_hy[1][2]));
          
      dx[2]=x[3];
      dx[5]=0;
      dx[4]=0;
      dx[7]=0;
      dx[6]=0;
      dx[7]=-F_gen;
      self.acc=[dx[1],dx[3],0];#*np.cos(x[2]);
      dx[8]=-Pabs;
      return dx;
      
      
    def get_data(self):
        return 
        
    def release(self):
        self.buoy.clear();
        del self.buoy;

#Testing  
if __name__=="__main__":
    #Testing if any function throughs an error
    omega=np.round(np.linspace(0.1,1.5*2+0.1,32),10);
    xi=omega*omega/9.81;
    w=WEC();
    w.load_buoy(xi,100,0);
    w.load_param();
    res=w.buoy.Calculate(0,0,0,0);
    print(w.pto_mdc());
