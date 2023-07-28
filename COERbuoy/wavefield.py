#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# wavefield.py - Handles the wave: Read the wave as time-domain data,
# transfers it into frequency domain (window it, FFT, cut-off unwanted frequencies),
# and then calculate the instantanious frequency domain response for any time and 
# horizontal position
# Furthermore, wave sources can be added, which is used for the radiation calculations
# 2020/2021 COER Laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie
"""
Created on Fri Feb  5 16:52:55 2021

@author: heiko
"""
import numpy as np;

def planck_window(N, epsilon):
    w0=[0];
    n=np.array(range(1,int(np.around(epsilon*N))))
    w1=1/(1+np.exp(epsilon*N/n-epsilon*N/(epsilon*N-n)));

    w2=np.ones([1,int(np.around(N/2-epsilon*N))]).flatten()

    w=np.concatenate([w0,w1,w2]);
    w=np.concatenate([w,np.flip(w,axis=0)]);
    return w;

#class for a wave created by a dirac impulse; needed for radiation calculation
class diracWave():
    A=np.array([]);
    t0=0;
    te=1;
    phase=np.array([]);
    omega=0;
    bcos=True;
    def __init__(self,omega,A,t0,bc):
        self.omega=omega;
        self.A=A*(self.omega[1]-self.omega[0]); 
        self.phase=0;
        self.t1=0.5/np.pi/np.max(omega);
        self.t0=t0#Current time; time when radiation started
        self.te=t0+0.25*np.pi/np.min(omega);#-t1;
        self.bcos=bc;
        if False:
            if t0>3:
                irf=[]; 
                for ts in np.linspace(self.t0,self.te,30):
                    irf.append((np.sum(self.evaluate(ts))));
                print(np.real(irf).tolist());
                print(np.imag(irf).tolist());
                print(self.te-self.t0);
                exit();
    def evaluate(self, t):
        if t<self.t0 or t>self.te:
            return self.A*0;
        #print(np.sum(self.A*np.cos(self.omega*(t-self.t0)-self.phase)))
        if self.bcos:
            return self.A*np.exp(self.omega*1j*(t-self.t0))#*np.cos(self.omega*(t-self.t0))#+1j*np.imag(self.A)*np.sin(self.omega*(t-self.t0))
        else:
            return -1j*self.A*np.exp(self.omega*1j*(t-self.t0))
            #return np.real(self.A)*np.sin(self.omega*(t-self.t0))+1j*np.imag(self.A)*np.cos(self.omega*(t-self.t0))
        return self.A*np.exp(self.omega*1j*(t-self.t0)-self.phase)#*(1-1/(self.te-self.t0)*(t-self.t0));

#The class can right now only handle a 1 dimensional wave, however, in the future
#it might be etended to planar wavefields
class wavefield:
    dw=[];
    dw2=[];
    g=9.81;
    rho=1000;
    
    @classmethod
    def set_wave(cls,y, t, cut_off):
      #Window the wave to avoid hard cuts, leading to high frequencies after FFT
      win = planck_window(len(y)+1,0.10)[0:len(y)]
      
      if abs(y[0]-y[-1])<0.1*np.mean(np.abs(y)):#no windowing if signal is periodic
          win=1;
      # FFT
      Y=(np.fft.fft(y*win))/y.size;
      omega=np.array(range(0,int(Y.size/2)-1))+1;
      
      #get amplitude and phase
      A=np.abs(Y);
      alpha=np.angle(Y);
      A=2*A[1:int(Y.size/2)]
      phase=alpha[1:int(Y.size/2)]
      omega=omega/t[-1]*(2*np.pi);
    
      A[np.where(A<0.001)]=0;
      #cut off high frequencies
      if np.argmax(omega>cut_off)>0:
          oi1 = 0;
          oi2 = np.argmax(omega>cut_off)
          omega = omega[oi1:oi2];
          A = A[oi1:oi2];
          phase=phase[oi1:oi2];
          
          import matplotlib.pyplot as plt
          #plt.figure();
          #plt.plot(t,y);
          #plt.plot(t,y*win);
          #plt.plot(t,np.sum(A*np.exp(1j*(omega*t.reshape(t.size,1)+phase)),1).real);
          #plt.show();
      print("Using the following angular wave frequencies: [rad/s],[m]")
      for (a,o) in zip(A,omega):
          print(str(round(o,2))+", "+str(a));
      
      return cls(A,phase,omega);
    
    def __init__(self, A, p, omega):
      self.A=A;
      self.phase=p;
      self.omega=omega;
      self.xi=omega*omega/self.g;
      
    def get_period(self):
        return 6.28/self.omega[np.argmax(np.abs(self.A))];
    
    def get_height(self):
        return sum(np.abs(self.A)+(0*self.omega[1]-0*self.omega[0]));
    
    #add the radiation wave caused by a current impulse #surge
    def add_diracWave(self,A,t0,bc):
        if (len(self.dw)>0):                        
            if t0-self.dw[-1].t0<=0.05:
                #we are too close to the previous time step
                #do not add the wave and save computational time
                return;
        self.dw.append(diracWave(self.omega,A,t0,bc));
    
    #add the radiation wave caused by a current impulse #heave
    def add_diracWave2(self,A,t0,bc):
        if (len(self.dw2)>0):
            if t0-self.dw2[-1].t0<=0.05:
                #we are too close to the previous time step
                #do not add the wave and save computational time
                return;
        self.dw2.append(diracWave(self.omega,A,t0,bc));
    
    #get the radiated wave, the sum of all waves created by impulses
    def get_rad(self,t,x):
        j=0;
        Adw1=np.zeros(len(self.omega));
        l=len(self.dw);
        t1=0;
        if isinstance(t,float):
          while j < l:
              if (j==0):
                  t1=self.dw[j].t0;
                  if t1>t:
                      j=len(self.dw)+1;
                      break;
              if self.dw[j].te<t:
                  self.dw.pop(j);
                  l=len(self.dw);
              else:
                  if t1<self.dw[j].t0:
                      Adw1=Adw1+self.dw[j].evaluate(t)*(self.dw[j].t0-t1);
                      t1=self.dw[j].t0;
                  j=j+1;
                  
        return Adw1;

    #get the radiated wave, the sum of all waves created by impulses
    def get_rad2(self,t,x):
        j=0;
        Adw1=np.zeros(len(self.omega));
        l=len(self.dw2);
        t1=0;
        if isinstance(t,float):
          while j < l:
              if (j==0):
                  t1=self.dw2[j].t0;
                  if t1>t:
                      j=len(self.dw2)+1;
                      break;
              if self.dw2[j].te<t:
                  #del self.dw[j];
                  self.dw2.pop(j);
                  l=len(self.dw2);
              else:
                  if t1<self.dw2[j].t0:
                      try:
                        Adw1=Adw1+self.dw2[j].evaluate(t)*(self.dw2[j].t0-t1);
                        t1=self.dw2[j].t0;
                      except:
                        t1=3;  
                  j=j+1;
                  
        return Adw1;
        
    #get the currrent surface elevation at time t and horizontal position x
    def get(self,t,x):
        a=self.A*np.exp(1j*(self.omega*t+self.phase-self.xi*x))#np.cos(self.omega*t+self.phase-self.xi*x);
        b=1j*self.omega*self.A*np.exp(1j*(self.omega*t+self.phase-self.xi*x))
        c=1j*1j*self.omega*self.omega*self.A*np.exp(1j*(self.omega*t+self.phase-self.xi*x))
       
        return[a,b,c];#wave_elevation, imag_wave_elevation, wave_speed, wave_speed_imag, wave_acceleration, wave_acceleration_imag
        
    def clear(self):
        self.dw.clear();
        self.dw2.clear();