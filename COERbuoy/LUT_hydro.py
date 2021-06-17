#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Creating the hydrodynamic Look-Up-Table (LUT) from the NEMOH output files
# The code was initally written for three degree-of-freedom (heave, surge and pitch), however,
# only two degrees are used (heave and surge), thus only heave parameters have to be read
#
# 2020/2021 COER Laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie
#

# Import modules
import os;
import pandas;#BSD 3-clause license
import numpy as np;#BSD 3-clause license
from scipy.interpolate import interp1d;#BSD 3-clause license
import pkg_resources
BEMdir=pkg_resources.resource_filename(__name__,"data/BEM");

# Values for whitch BEM data is available min:step:max
# p-pitch, h-heave
p_max=5;
p_min=-5;
h_max=4;
h_min=-4;
p_step=2.5;
h_step=0.5;

def interpolate(x,y,o):
    fi=interp1d(x,y,fill_value='extrapolate',bounds_error=False);
    return fi(o);          

LUT={};
    
# Load look up table
def load_LUT(omegatest):  
    for entry in os.scandir(BEMdir):
        if entry.is_dir():
            a=entry.name.split("_");
            if len(a)==4:
                h=(int)((float)(a[1])*10);
                p=(int)((float)(a[3])*10);
                #print("Reading data for heave h and pitch p");
                
                #Excitation force
                b=np.array(pandas.read_csv(entry.path+"/results/ExcitationForce.tec",header=4,delimiter="\s+"));
                omega=b[:,0];
                fe_a=np.array([[0]*len(omegatest)]*3).copy();
                fe_p=np.array([[0]*len(omegatest)]*3).copy();
                
                
                
                for i in np.arange(3):
                    fe_a[i]=interpolate(omega,b[:,1+i*2],omegatest);
                    fe_p[i]=interpolate(omega,b[:,2+i*2],omegatest);
                    
                    
                #Radiation force
                fr_m=np.array([[[0]*len(omegatest)]*3]*3).copy();
                fr_r=np.array([[[0]*len(omegatest)]*3]*3).copy();#[[[0]*omega]*3]*3;
                
                o=omegatest
                for i in np.arange(3):
                    b=np.array(pandas.read_csv(entry.path+"/results/RadiationCoefficients.tec",header=4+(1+len(omega))*i,nrows=len(omega),delimiter="\s+"));
                    fr_m[0][i]=interpolate(omega,b[:,1].copy(),o)
                    fr_r[0][i]=interpolate(omega,b[:,2].copy(),o)
                    fr_m[1][i]=interpolate(omega,b[:,3].copy(),o)
                    fr_r[1][i]=interpolate(omega,b[:,4].copy(),o)
                    fr_m[2][i]=interpolate(omega,b[:,5].copy(),o)
                    fr_r[2][i]=interpolate(omega,b[:,6].copy(),o)
                LUT[str((int)(h))+"_"+str((int)(p))]=[fe_a,fe_p,fr_m,fr_r]            
                
# get parameters for body pose specified by heave and pitch (only heave working)
# interpolate (linearly) between values
def get_fromLUT(h,p):
    p=np.min([p_max,p]);
    p=np.max([p_min,p]);
    h=np.min([h_max,h]);
    h=np.max([h_min,h]);
    p1=0#(int)(np.floor(p/p_step))*0;
    h1=(int)(np.floor(h/h_step));
    p2=0#np.min([(int)(np.ceil(p/p_step)),0])*0;
    h2=(int)(np.ceil(h/h_step));
    
    d1=LUT[str((int)((h1*h_step)*10))+"_"+str((int)((p1*p_step)*10))]
    d2=LUT[str((int)((h2*h_step)*10))+"_"+str((int)((p1*p_step)*10))]
    d3=LUT[str((int)((h1*h_step)*10))+"_"+str((int)((p2*p_step)*10))]
    d4=LUT[str((int)((h2*h_step)*10))+"_"+str((int)((p2*p_step)*10))]
    
    d=[[],[],[],[]]
    for idx, (e1,e2,e3,e4) in enumerate(zip(d1,d2,d3,d4)):
        pc=np.mod(h,h_step)/h_step;
        pcp=np.mod(p,p_step)/p_step;
        e5=e1*(1-pc)+e2*pc;
        e6=e3*(1-pc)+e4*pc;
        d[idx]=e5*(1-pcp)+e6*(pcp);
        
    return d;