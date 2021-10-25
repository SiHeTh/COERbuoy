#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 23:16:53 2020
Moment-based control with TCP interface; no bugfixing yet

@author: Simon H. Thomas
"""

#Linear receeding horizon WEC control in the moment domain
#Based on the following literature:
#[1] N. Faedo et al.: Energy-maximising control of wave energy converters using a moment-domain representation, Control Engineering Practice 81, 2018
#[2] N. Faedo et al.: Receding-Horizon Energy-Maximising Optimal Control of Wave Energy Systems Based on Moments, IEEE Transaction on Sustainable Energy, 2020

import COERbuoy.Parameters;
import COERbuoy.connection as connection;
import numpy as np; #numpy under BSD 3-clause license
from scipy.interpolate import interp1d; #scipy under BSD 3-clause license
#The cvxopt library is needed to solve the qp problem, please uncomment the following line and install the CVXOPT library
#from cvxopt import matrix, solvers #CVXOPT, GPL license (copyleft);

##A) Basic set-up
#set constants
rho=1000;#water density
pi=np.pi;#pi
g=9.81;#gravity acceleration
h=60;#water depth

#do not show the qp solver progress
solvers.options["show_progress"]=False;

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
omega0=np.array([0.4,0.8,1.2,1.6]);

#number of frequencies (k)
k=omega0.size;

#Create the time series so that it gets valid data for the highest and lowest frequency
ts=np.arange(-1*pi/np.min(omega0),pi/np.min(omega0),1/(2*np.max(omega0)));
dt=ts[1]-ts[0]; #time step

#set the S matrix
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
    #used to do the paramerisation dynamically.
    #If the controller would be for the "COERbuoy1" WEC model, the class would be called:
    param = COERbuoy.Parameters.parameters("[data.COERbuoy1]","Floater_BEM");
    param.init_hydro(omega0); #select the frequencies
    data=param.dic_param(); #get the parameter list
    line_length=data["viscous_drag_coefficient_heave"];#The parameterlist gives access to all parameters defined for this device


#get pto parameter (mass, damping, stiffness) for system state-space model
[m,d,c]=param.pto_mdc(0);

#add buoyancy stiffness
c=c-param.area_vol(0)[0]*1000;

#load hydrodynamic parameters
hyparams=param.hydro(0,1,1);#usage: params.hydro(submergence=0,mode1=1(heave,mode2=1(heave))
                            #returns: [buoyancy, excitation, rad. impedande, added mass at inf.]
m0=np.imag(hyparams[3]);
rad1=np.real(hyparams[2]);
mad1=-omega0*np.real(-hyparams[2]);
exc1=np.real(hyparams[1]);


#calculate Tau_R matrix eq. (59) in [1]
Tr=np.zeros([rad1.size*2,rad1.size*2]);
for idx, (M_pw, pw,r_pw) in enumerate(zip(mad1,omega0,rad1)):
    m_pw=-1*pw*(M_pw-m0); # eq. (46) [1]
    a=pw*((m+m0)*pw**2+m_pw*pw-c);# eq. (59) [1]
    b=pw**2*r_pw**2+(a/pw)**2;
    Tr[idx*2,idx*2]=pw**2*r_pw/b;
    Tr[idx*2+1,idx*2]=-a/b;
    Tr[idx*2,idx*2+1]=a/b;
    Tr[idx*2+1,idx*2+1]=(pw**2*r_pw)/b;
Tr_eye=np.eye(Tr[:,1].size) # eq. (60) in [1]


##D) Set-up constraint handling
#set heave limit
X_Max=5;
#set constraints matrix
nabla=np.vstack([Se, -Se]).transpose()
c1=np.matmul(np.matmul(-Tr,np.linalg.inv(S)),nabla); #eq. (66c), left hand side in [1]


##E) Controller
conn_model=connection.connection();#Open connection to simulation
host=False;    
if host:
    conn_model.openH();
else:
    conn_model.openC();


def Calc_fe(wave):
    a12=np.linalg.lstsq(Se,wave,rcond=None)[0]; #wave into moment domain
    #ab=-np.sqrt(a12[::2]**2+a12[1::2]**2); #absolute value
    #phase=np.arctan2(a12[1::2],a12[::2]); #phase
    res=param.hydro(0,1,1);#get hydroparameter
    
    a12[::2]=a12[::2]*np.real(res[1]);
    a12[1::2]=a12[1::2]*np.imag(res[1]);
    
    return a12;

while True: ##main loop
  buf_l=conn_model.get_control();
  if not buf_l: break # no more data? Then knock off!
  
  #1) read the sensor reading of the wave senor and map into onto ts
  now=buf_l["time"][-1]; # last value gives current time
  buf_time=np.array(buf_l["time"])-buf_l["time"][50];#combine forecast and past data
  wave=(buf_l["wave_forecast"][50:0:-1]+buf_l["wave"][:49:-1])[::1];#combine past and upcoming wave
  f=interp1d(buf_time.tolist(),wave);#Set the incoming wave into the right format
  win = planck_window(ts.size,0.3)[0:ts.size]#Windowing for receeding horizon
  windowed_sample=f(ts)*win;
  
  
  #2) Get the signal generator and constraint matrix ready
  Le = Calc_fe(windowed_sample.transpose())#wave*excitation into time domain -> signal generator  
  c2=(X_Max*np.ones([1,(nabla[1,:].size)])-np.matmul(np.matmul(np.matmul(Le,Tr),np.linalg.inv(S)),nabla));#vector for constraint handling, eq. (66c), right hand side in [1]
  
  
  #3) Solve qp, get signal generator for PTO force (Lu)
  Lu=solvers.qp(matrix(Tr_eye),matrix(-0.5*np.matmul(Le,Tr)),matrix(c1.transpose()),matrix(c2.transpose()));#eq (66) in [1]
  Vl=np.matmul(Le-np.array(Lu['x']).transpose(),Tr);# eq. (56) in [1]
  
  
  #4) obtain control signal from solution
  #optional: position tracking: ei=np.matmul(Se,np.matmul(Vl,np.linalg.inv(S)).transpose())-buf_l["stroke_pos"][-1];
  e=np.matmul(Se,Vl.transpose())-buf_l["stroke_speed"][-1]; #velocity tracking
  f=1*np.array(np.matmul(Se,Lu['x']))#get force data
  force=e*15000*10; #velocity tracking
  
  
  #5) write and send answer
  start=int(ts.size/2-4);
  answer={
          "time":ts[start:start+9].flatten()+now,#1x9 array with time steps in the future
          "pto":force[start:start+9].flatten(),  #1x9 array with PTO force
          "brake":np.zeros(9),                   #1x9 array with brake force
          "test":np.zeros(9),                    #1x9 array; not used
          }
  
  conn_model.set_control(answer["time"],answer["pto"],answer["brake"],answer["test"]); #Send control message to model
  
conn_model.close();