#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple linear damping
@author: Simon H. Thomas, COER laboratory, Maynooth University
"""
import numpy as np;
from scipy.interpolate import interp1d;
import COERbuoy.connection as connection;

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
  
  ##This part contains the controller's logic
  ##TODO: replace logic with own controller idea
  
  #(here we use a simple linear damping)
  #Calculate PTO force
  gamma=100000;#kNs/m
  F_pto=-gamma*dx[-1];
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
