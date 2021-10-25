#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dynamics.py - A generic WEC model

#The WEC class specifies a generic WEC model

#Dependencies
import json;
import numpy as np;
import os;
import COERbuoy.utils as utils;

g=9.81; #gravity acceleration
rho=1000; #density of water


class WEC():
    omega_cut_off=3.2; #highest frequency used for calculation    
    
    states=10; #number of states the model needs fo its calculation
    #For compatibility with the COERbuoy platfrom the following conventions is
    #recommended:
    #state id, Description
    #0, stroke length
    #1, stroke speed
    #2, angular position
    #3, angular speed
    #4, free
    #5, free
    #6, free
    #7, Generator force
    #8, Absorbed power
    #9, free
    
    
    #load parameters (do not change)
    file=os.path.join(utils.wec_dir,"floater.txt");
    
    acc=0;
    
    def __init__(self): #initialise
        return;
        
    def load_buoy(self,floater_class,xi,depth,cog): #load buoy (do not change)
        self.buoy=floater_class(xi,g,depth,cog,self.file);
        self.mass=self.buoy.Volume(0)*1000; # mass is usually buoyancy at mean
    
    def load_param(self):
        #specify the parameters loaded from floater.txt
        with open(self.file) as file:
            data=json.load(file);
            #generator damping (only for testing)
            self.damping =data.get("PTO_damping",100000); #("parameter_name",standard_value)
    
    #Get linearised mass, damping and spring coefficent (only for calculation of eigenperiod)    
    #Includes only non-hydrodynamic parameters
    def pto_mdc (self,z):
        m=(self.mass);
        d=0;
        c=0;
        return [m,d,c];
    
    #the following values need to be specified for the
    #status message for the controller
    def get_surge(self,x): #specify the surge component
        return  0;
    def get_translator_speed(self,x): #the translator speed
        return x[1];
    def get_translator_position(self,x): # the translator position
        return x[0];
    def get_force(self,x): # the force of an integrated force sensor
        return 0;
    
    ## This function includes the main logic of the controller
    #it is called at every iteration of the ODE loop
    #Please write your WEC ogic here!
    def Calc(self,t,wave,x,PTO_force,brake,ulast):
      #Input:
      #t         - current time since start in seconds
      #wave      - the wave as an object
      #x         - the state vector with lenth n_states (see above)
      #PTO_force - force set by the controller
      #brake     - brake force set by the controller
      #ulast     - velocity during the last iteration
      
      #Example here for a simple heave only PAWEC.
      heave       = x[0];
      surge       = 0;
      pitch       = 0;
      surge_speed = 0;
      heave_speed = x[1];
      pitch_speed = 0;
      
      #get hydrodynamic forces: time, wave object, position vector, velocity vector, accleration vector
      #position, velcoity, accleration vector: [surge, heave, picth]
      F_hy = self.buoy.get_forces(t,wave,[surge,heave,pitch],[surge_speed,heave_speed,pitch_speed],self.acc)
      
      #return value of F_hy:
      #F_hy[0]-> hydrostatic/-dynamic forces
      #    F_hy[0][0]-> force surge [N]
      #    F_hy[0][1]-> force heave [N]
      #    F_hy[0][2]-> moment pitch [N/m]
      #F_hy[1]-> added inertia
      #    F_hy[1][0]-> added mass surge [kg]
      #    F_hy[2][1]-> added mass heave [kg]
      #    F_hy[3][2]-> added inertia pitch [kg m^2]
      
      #Fill dx vector
      dx=np.zeros(self.states);
      dx[1]=(F_hy[0][1]-self.mass*g+PTO_force)/(self.mass+F_hy[1][1]); #calculate heave acceleration
      dx[0]=x[1];
          
      
      dx[3]=0;
      dx[2]=x[3];
      
      dx[5]=0;
      dx[4]=0;
      dx[7]=0;
      dx[6]=0;
      
      dx[7]=-PTO_force;# store PTO_force
      self.acc=[0,dx[1],0];#set acceeration vector
      dx[8]=-PTO_force*x[1];#calculate absorbed power
      return dx;
      
      
    def get_data(self):
        return 
        
    def release(self):
        self.buoy.clear();
        del self.buoy;

#Testing  
if __name__=="__main__":
    from COERbuoy import floater_BEM_LUT as Floater;
    from COERbuoy import wavefield;
    #Testing if any function threws an error
    w=WEC();
    t=np.linspace(1,10,100);
    wave=wavefield.wavefield.set_wave(t,np.sin(t/10),w.omega_cut_off);
    w.load_buoy(getattr(Floater,utils.class_hydro),wave.xi,300,0);
    w.load_param();
    res=w.buoy.Calculate(0,0,0,0);
    print(w.pto_mdc());
