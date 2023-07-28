#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple linear damping
@author: Simon H. Thomas, COER laboratory, Maynooth University
"""
import numpy as np;
from scipy.interpolate import interp1d;
import COERbuoy.connection as connection;
import COERbuoy.Parameters;
import sys;
from scipy.signal import find_peaks;

period = 10;
H=2;
if len(sys.argv)>2:
    period=float(sys.argv[1]);
    H=abs(float(sys.argv[2]));

omega = 6.28/period;
param = COERbuoy.Parameters.parameters(None,"Floater_LIN");#get hydrodynamic data based on the current buoy
param.init_hydro([omega]); #select the frequencies

[fe,R]=param.hydro(0, 1, 1)[1:3];

A=H/2;
Fe=1.4*A*fe.real[0];
u1=Fe/(2*R.real[0]);
a1=float(u1/omega);
print(a1,A,omega)
#hyparams=param.hydro(0,1,1);
#rad1=np.real(hyparams[2]);
#mad1=-np.real(hyparams[2]);
  
#print("Reactive control optimized for wave period "+str(period)+" s.")

conn_model=connection.connection();#Initialize connection
conn_model.openC();#Use client mode
while msg:=conn_model.get_control():

  ##Read incoming message
  time  =msg["time"]          #1x100 array with time series
  wave  =msg["wave"];         #1x100 array with wave data (related to time series)*
  wave_f=msg["wave_forecast"];#1x100 array with wave forecast*
  x     =msg["stroke_pos"];   #1x100 array with stroke position (related to time series)
  dx    =msg["stroke_speed"]; #1x100 array with stroke speed (related to time series)
  alpha =msg["angular_pos"];  #1x100 array with pitch angle (related to time series)
  dalpha=msg["angular_speed"];#1x100 array with pitch angular speed (related to time series)
  force =msg["force"];        #1x100 array with force sensor data (related to time series)
  #* not available during COERbuoy1 benchmark
  now=time[-1]; #The last element contains the most recent data (except the wave forecast)
  
  
  t_peak=find_peaks(wave_f)[0];
  t_peak_m=np.mean(np.diff(t_peak));
  #print(t_peak)
  #print(t_peak_m)
  pos=wave[-int(t_peak_m/4)];
  
  #ddx=(dx[-1]-dx[-2])/(time[-1]-time[-2])
  omega=6.28/(t_peak_m*(time[1]-time[0]));
  a1=np.min([abs(a1),3]);
  F_pto=400000*1.2*(a1*pos/A-x[-1])+1400000*1.2*(a1*omega*wave[-1]/A-dx[-1]);
  #print([a1*pos-x[-1],400000*(a1*pos-x[-1]),1000000*(a1*omega*wave[-1]-dx[-1])])
  #In this example we don't use the WECs brake
  brake=0;#kNs/m
  
  ##Write control message 
  answer={
          "time":np.linspace(now,now+1,9),#1x9 array with time steps in the future
          "pto":np.array([F_pto]*9),      #1x9 array with PTO force
          "brake":np.zeros(9),            #1x9 array with brake force
          "test":np.zeros(9),             #1x9 array; not used
          }
  #Send control message to model
  conn_model.set_control(answer["time"],answer["pto"],answer["brake"],answer["test"]);

conn_model.close();

