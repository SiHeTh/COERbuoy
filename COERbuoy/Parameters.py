#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 13:34:05 2020

@author: heiko
"""

import numpy as np;
import pandas;
import COERbuoy.utils as utils;
from COERbuoy import floater_LIN as Floater;
from COERbuoy import wavefield;
import importlib;
import os;
 
spec=importlib.util.spec_from_file_location("dynamics.py",os.path.join(utils.wec_dir,"dynamics.py"));
dynamics=importlib.util.module_from_spec(spec);
spec.loader.exec_module(dynamics);
#import dynamics_sphere as dynamics
def run():
    wec=dynamics.WEC();
    omega=np.round(np.linspace(0.1,wec.omega_cut_off,10),10);#select frequency used
    xi=omega*omega/9.81;#
    
    
    #Initialize WEC
    wec.load_buoy(getattr(Floater,utils.class_hydro),xi,300,0);
    wec.load_param();
    
    #limit=7;#heave limit inbetween which to calculate parameters
    delta=0.1;
    geoBox=wec.buoy.getGeoBox();#get limits inbetween which to calculate parameters
    zs=np.linspace(geoBox[0],geoBox[1],18+1);#submergence levels at which to evaluate
    n_mode=3;
    area=np.zeros(zs.size);
    vol=np.zeros(zs.size);
    Fstat=np.copy([np.zeros(zs.size)]*n_mode);
    Fdyn=np.copy([np.zeros([zs.size,xi.size])]*n_mode);
    Frad=np.copy([np.zeros([zs.size,xi.size])]*n_mode);
    Amass=np.copy([np.zeros([zs.size,xi.size])]*n_mode);
    Fdiff=np.copy([np.zeros([zs.size,xi.size])]*n_mode);
    Am8=np.copy([np.zeros([zs.size])]*n_mode);
    floater_slider_spring=np.zeros(zs.size);
    
    #Madd-Damper-Spring parameters
    mdc=np.array([0,0,0]);
    
    #mode=1;#0-surge, 1-heave, 2-pitch
    
    
    wec.buoy.Calculate(0,0,0,0);
        
            
    for idx,z in enumerate(zs):#calculate parameters for each submergence level
        area[idx]=wec.buoy.Area(z);
        vol[idx]=wec.buoy.Volume(z);
        results=wec.buoy.Calculate(z,0,0,0)
        results2=wec.buoy.Calculate(z-delta,0,0,0);
        for mode in [0,1,2]:
            Fstat[mode][idx]=(results[0][mode]-results2[0][mode])*1/delta;
            Fdyn[mode][idx]=np.real(results[1][mode]);
            Frad[mode][idx]=np.real(results[2][mode][mode]);
            Amass[mode][idx]=np.imag(results[2][mode][mode]);
            Am8[mode][idx]=results[3][mode];
            Fdiff[mode][idx]=np.imag(results[1][mode]);
        floater_slider_spring[idx]=(wec.Calc_fs_spring(z,0)-wec.Calc_fs_spring(z-delta,0))*1/delta;
        
        #Get linear parameter only at zero position
        if z==0:
            mdc[1]=Frad[1][idx][0];
            mdc[2]=-Fstat[1][idx]#-floater_slider_spring[idx]);
    
    print("calculation finished")
    omega2=[];
    
    #Write data
    folder0=utils.pkg_dir+"/param/";
    mod_name={0:"surge/",1:"heave/",2:"pitch/"}
    folder=folder0;
    pandas.DataFrame(np.vstack((zs,area,vol)).transpose(),columns=["z-offset","cross-sec.area","volume"]).round(2).to_csv(folder+"HydroParam1.csv",index=False)
    pandas.DataFrame(np.vstack((zs,Fstat[1],floater_slider_spring,Fstat[1]+floater_slider_spring)).transpose(),columns=["z-offset","c_hydrostatic","c_spring","resulting"]).round(2).to_csv(folder+"Stiffness.csv",index=False)
        
    for o in omega:
        omega2.append("&omega;="+str(o.round(2))+" rad/m");
    
    omega8=omega2.copy();
    omega8.append("&omega;=inf");
 
    for mode in [0,1,2]:
        folder=folder0+mod_name[mode];
        pandas.DataFrame(np.vstack((zs,Fdyn[mode].round(0).transpose())).transpose(),columns=["Froude-Krylov force"]+omega2).to_csv(folder+"fk_force.csv",index=False)
        pandas.DataFrame(np.vstack((zs,Frad[mode].round(0).transpose())).transpose(),columns=["Radiation force"]+omega2).to_csv(folder+"radiation_force.csv",index=False)
        pandas.DataFrame(np.vstack((zs,np.concatenate((Amass[mode],Am8[mode].reshape(len(Am8[mode]),1)),axis=1).round(0).transpose())).transpose(),columns=["Added mass"]+omega8).to_csv(folder+"added_mass.csv",index=False)
        pandas.DataFrame(np.vstack((zs,Fdiff[mode].round(0).transpose())).transpose(),columns=["diffraction force"]+omega2).to_csv(folder+"diff_force.csv",index=False)
    
    #Calculate generator efficancy
    wave=wavefield.wavefield(omega*0,omega*0,omega);
    v_max=2.5;
    p_max=500000;
    n=8;
    eff=np.zeros([n,n]);
    #print(eff)
    
    speed_s=np.linspace(v_max*0.5,v_max*1.5,n);
    power_s=np.linspace(p_max*0+100,p_max*1.2,n);
    for idx,s in enumerate(speed_s):
        for jdx,p in enumerate(power_s):
            x=np.zeros(wec.states);
            x[1]=s;
            #print(1*wec.Calc(0,wave,x,-1*p/s,0,[1,1,1])[8]/(p))
            eff[jdx][idx]=1*wec.Calc(0,wave,x,-1*p/s,0,[1,1,1])[8]/(p);
            #print(eff[jdx][idx])
            
            
    print(eff)
    axis2=[];
    for s in speed_s:
        axis2.append("velocity="+str(s.round(2))+"m/s");
    pandas.DataFrame(np.vstack(((power_s/1000).round(2),eff.round(2).transpose())).transpose(),columns=["Generator_efficency_over_Power[kW]"]+axis2).to_csv(folder0+"gen_eff.csv",index=False)
    
    mdc=mdc+wec.pto_mdc();#Get lienarized data from WEC

    #... now the linearized eigenfrequencys can be calculated
    deigfreq=2*np.pi/np.sqrt(mdc[2]/mdc[0]-(0.5*mdc[1]/mdc[0])**2)
    eigfreq=2*np.pi/np.sqrt(mdc[2]/mdc[0])
    pandas.DataFrame(np.vstack((mdc[0],mdc[1],mdc[2],eigfreq,deigfreq,wec.buoy.Calc_CoG())).transpose(),columns=["mass[kg]","damping [Ns/m]","stiffness [N/m]","eigenperiod [s]","damped eigenperiod [s]","center of gravity (heave) [m]"]).round(2).to_csv(folder0+"info.csv",index=False)

    #clearning up
    wec.release();

if __name__ == "__main__":
    run();    