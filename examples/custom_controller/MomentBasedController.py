#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 23:16:53 2020
Moment-based control with TCP interface; no bugfixing yet

@author: Simon H. Thomas
"""

#Linear receeding horizon WEC control in the moment domain
#With stroke limit
#Self-parametrising using the COERbuoy params interface

#Based on the following literature:
#[1] N. Faedo et al.: Energy-maximising control of wave energy converters using a moment-domain representation, Control Engineering Practice 81, 2018
#[2] N. Faedo et al.: Receding-Horizon Energy-Maximising Optimal Control of Wave Energy Systems Based on Moments, IEEE Transaction on Sustainable Energy, 2020
#missing: additional condition to include actual position
#velocity and force limit

import COERbuoy.Parameters;
import COERbuoy.connection as connection;
import quadprog; #BSD 2-clause license
import numpy as np; #BSD 3-clause license
from scipy.interpolate import interp1d; #BSD 3-clause license
import sys;

##A) Basic set-up
#set constants
rho=1000; #water density
pi=np.pi; #pi
g=9.81;   #gravity acceleration
h=60;     #water depth

#define planck window for filtering the input data
def planck_window(N, epsilon):
    w0=[0];
    n=np.array(range(1,int(np.around(epsilon*N))))
    w1=1/(1+np.exp(epsilon*N/n-epsilon*N/(epsilon*N-n)));

    w2=np.ones([1,int(np.around(N/2-epsilon*N))]).flatten()
    w=np.concatenate([w0,w1,w1,w2]);
    if N%2==0:
        w=np.concatenate([w0,w1,w2]);
    w=np.concatenate([w,np.flip(w,axis=0)]);
    return w;


##B) Moment domain set-up
#select frequencies; frequencies must be a n*base_frequency
omega0=np.array([0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3,3.3,3.6,3.9]);# comprise all relevant angular wvae frequencies

#count of frequencies (k)
k=omega0.size;

#Create the time series so that it gets valid data for the highest and lowest frequency
ts=np.arange(-1*pi/np.min(omega0),pi/np.min(omega0),1/(2*np.max(omega0)));
dt=ts[1]-ts[0]; #time step

#set the S (deviation) matrix
S=np.zeros([omega0.size*2,omega0.size*2]);
for idx, w in enumerate(omega0):
    S[idx*2+1,idx*2]=-w;
    S[idx*2,idx*2+1]=w;

#Setting transformation matrix Se, which transfers between moment and time domain
Se=np.array([]);
for t in ts:
    a=np.array([])
    for w in omega0:
        a=np.append(a,[np.cos(w*t),-np.sin(w*t)]);
    if Se.size == 0:
        Se=a;
    else:
        Se=np.vstack([Se,a])

##C) Set-up model
#This controller example is not adapted to a specific WEC, it will get the
#basic parameter directly from the "parameter" class.
param = COERbuoy.Parameters.parameters(None,None);
param.init_hydro(omega0); #select the frequencies

if False: #dynamic parametrisation for a specific WEC model
    #To get the maximum out of a specific WEC, controller have tp be
    #specifically designed for this device. In this case the parameter class can be
    #used to do the parameterisation dynamically.
    #If the controller would be for the "COERbuoy1" WEC model, the class would be called:
    param = COERbuoy.Parameters.parameters("[data.COERbuoy1]","Floater_BEM");
    param.init_hydro(omega0); #select the frequencies
    data=param.dic_param(); #get the parameter list
    line_length=data["viscous_drag_coefficient_heave"];#The parameterlist gives access to all parameters defined for this device

#load parameters
[m,d,c]=param.mdc(0);

#calculate systemk matrix from eq. (59) in [1]
Tr=np.zeros([m.size*2,m.size*2]);
for idx, (m0, w0, d0, c0) in enumerate(zip(m,omega0,d,c)):
    X=w0*(m0*w0**2-c0);# eq. (59) [1]
    ZZ=((m0*w0**2-c0)**2+(d0*w0)**2);
    Tr[idx*2,idx*2]=w0**2*d0/ZZ;
    Tr[idx*2+1,idx*2]=X/ZZ;
    Tr[idx*2,idx*2+1]=-X/ZZ;
    Tr[idx*2+1,idx*2+1]=(w0**2*d0)/ZZ;
Tr_eye=Tr*np.eye(Tr[:,1].size) # eq. (60) in [1]

##D) Set-up constraint handling
#set heave limit
X_Max=3;
if len(sys.argv)>1: #set firstv argument as stroke limit
    X_Max=float(sys.argv[1]);
    
#set constraints matrix
nabla=np.vstack([Se, -Se]).transpose()
c1=np.matmul(np.matmul(-Tr,np.linalg.inv(S)),nabla); #eq. (66c), left hand side in [1]


##E) Controller
conn_model=connection.connection();#Open connection to simulation
host=False; #Standard: client mode
if host:
    conn_model.openH();
else:
    conn_model.openC();
    
exc=param.hydro(0,1,1)[1];#read excitation coefficients  
    
def Calc_fe(wave):
    W=np.linalg.lstsq(Se,wave,rcond=None)[0]; #wave into moment domain

    E=W.copy();
    E[::2]=W[::2]*exc.real-W[1::2]*exc.imag;
    E[1::2]=W[1::2]*exc.real-W[::2]*exc.imag;
    return E;

while True: ##main loop
  buf_l=conn_model.get_control();
  if not buf_l: break # no more data? Then knock off!
  
  #1) read the sensor reading of the wave senor and map into onto ts
  now=buf_l["time"][-1]; # last value gives current time
  #set time horizon of data (past wave data and forecast)
  t_lim=buf_l["time"][-1]-buf_l["time"][0];
  buf_time=np.linspace(-t_lim,t_lim,200)
      
  wave=(buf_l["wave"][::1]+buf_l["wave_forecast"][::1])[::1];#combine past and forecasted wave
  f=interp1d(buf_time.tolist(),wave);#Set the incoming wave into the right format
  win = planck_window(ts.size,0.3)[0:ts.size]#Windowing for receeding horizon
  windowed_sample=f(ts)*win;
  
  
  
  #2) Get the signal generator and constraint matrix ready
  Le = Calc_fe(windowed_sample.transpose())#wave*excitation into time domain -> signal generator  
  c2=(X_Max*np.ones([1,(nabla[1,:].size)])-np.matmul(np.matmul(np.matmul(Le,Tr),np.linalg.inv(S)),nabla));#vector for constraint handling, eq. (66c), right hand side in [1]
  #Use quadratic solver 'quadprog' to solve the problem
  #Minimize     1/2 x^T G x - a^T x
  #Subject to   C.T x >= b
  Lu=quadprog.solve_qp(Tr_eye,0.5*np.matmul(Tr,Le).transpose(),c1,-1*c2[0])[0]
  Vl=np.matmul(Le-np.array(Lu).transpose(),Tr);# eq. (56) in [1]
  
  #4) obtain control signal from solution
  Vl[0]=0;#remove first moment (numerical errors)
  Vl[1]=0;#remove first moment (numerical errors)
  v1=np.round(np.matmul(Vl.transpose(),Se.transpose()),2)
  p1=np.round(np.matmul(np.matmul(Vl,np.linalg.inv(S)),Se.transpose()),4)
  ei=-p1-buf_l["stroke_pos"][-1]; #position tracking
  e=-v1-buf_l["stroke_speed"][-1]; #velocity tracking
  force=-1/Vl.size*np.matmul(np.array(Lu),Se.transpose())#get force data

  force=force+1400000*e+400000*ei;
  
  #5) write and send answer
  start=int(v1.size/2-4);
  
  answer={
          "time":ts[start:start+9].flatten()+now,#1x9 array with time steps in the future
          "pto":force[start:start+9].flatten(),  #1x9 array with PTO force
          "brake":np.zeros(9),                   #1x9 array with brake force
          "test":np.zeros(9),                    #1x9 array; not used
          }
  
  conn_model.set_control(answer["time"],answer["pto"],answer["brake"],answer["test"]); #Send control message to model
  
conn_model.close();
