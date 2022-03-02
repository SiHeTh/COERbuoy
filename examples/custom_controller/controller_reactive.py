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
m_c=None;
d_c=None;
def calc_control(period):
    omega = 6.28/period;
    param = COERbuoy.Parameters.parameters(None,None);#get hydrodynamic data based on the current buoy
    param.init_hydro([omega]); #select the frequencies

    [m1,d1,c1]=param.mdc(0);
    m2=-1*float(c1)+omega**2*(float(m1));
    print("Reactive control optimized for wave period "+str(period)+"s.")
    return [float(d1),float(m2)]
    
if len(sys.argv)>1:
    period=float(sys.argv[1]);
    [d_c,m_c]=calc_control(period);
    
  

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
  
  if not d_c:
      t_peak=find_peaks(wave_f)[0];
      t_peak_m=np.mean(np.diff(t_peak));
      pos=wave[-int(t_peak_m/4)];
      period=(t_peak_m*(time[1]-time[0]))
      [d_c,m_c]=calc_control(period);
  
  F_pto=-1*d_c/2*dx[-1]-m_c*x[-1];
  
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

