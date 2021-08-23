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
import COERbuoy.utils as utils;
import pandas;#BSD 3-clause license
import numpy as np;#BSD 3-clause license
from scipy.interpolate import interp1d;#BSD 3-clause license
#BEMdir=pkg_resources.resource_filename(__name__,"data/BEM");
bem_dir=os.path.join(utils.wec_dir,"BEM");

# Values for whitch BEM data is available min:step:max
# p-pitch, h-heave
p_max=0;
p_min=-0;
h_max=0;
h_min=-0;
p_step=1000;
h_step=1000;
prec=100;

def interpolate(x,y,o):
    fi=interp1d(x,y,fill_value='extrapolate',bounds_error=False);
    return fi(o);          

LUT={};
    
# Load look up table
def load_LUT(omegatest):  
    global h_min, h_max, p_min, p_max, h_step, p_step;
    p_max=0;
    p_min=-0;
    h_max=0;
    h_min=-0;
    p_step=1000;
    h_step=1000;
    prec=100;
    for entry in os.scandir(bem_dir):
        if entry.is_dir():
            a=entry.name.split("_");
            if len(a)==4:
                h0=((float)(a[1]));
                p0=((float)(a[3]));
                if np.abs(h_max-h0)<h_step:
                    h_step=np.round(np.abs(h_max-h0),4);
                if h0>h_max:
                    h_max=h0;
                if h0<h_min:
                    h_min=h0;
                if np.abs(p_max-p0)<p_step:
                    p_step=np.round(np.abs(p_max-p0),4);
                if p0>p_max:
                    p_max=p0;
                if p0<p_min:
                    p_min=p0;
                    
                h=(int)(h0*prec);
                p=(int)(p0*prec);
                #print("Reading data for heave h and pitch p");
                
                #Excitation force
                b=np.array(pandas.read_csv(entry.path+"/results/ExcitationForce.tec",header=4,delimiter="\s+"));
                omega=b[:,0];
                fe_a=np.array([[0.0]*len(omegatest)]*3).copy();
                fe_p=np.array([[0.0]*len(omegatest)]*3).copy();
                
                
                
                for i in np.arange(3):
                    fe_a[i]=interpolate(omega,b[:,1+i*2],omegatest);
                    fe_p[i]=interpolate(omega,b[:,2+i*2],omegatest);
                    
                    
                #Radiation force
                fr_m=np.array([[[0.0]*len(omegatest)]*3]*3).copy();
                fr_r=np.array([[[0.0]*len(omegatest)]*3]*3).copy();#[[[0]*omega]*3]*3;
                fr_inf=np.array([0.0]*3).copy();#[[[0]*omega]*3]*3;
                
                o=omegatest
                for i in np.arange(3):
                    #b=np.array(pandas.read_csv(entry.path+"/results/RadiationCoefficients.tec",header=4+(1+len(omega))*i,nrows=len(omega),delimiter="\s+"));
                    #b=np.array(pandas.read_csv(entry.path+"/results/RadiationCoefficients.tec",header=4+(1+len(omega))*i,nrows=len(omega),delimiter="\s+"));
                    b=np.array(pandas.read_csv(entry.path+"/results/RadiationCoefficients.tec",header=4+(1+len(omega))*i,nrows=len(omega),sep="\s+|,+"));
                    if np.isnan(b[0,0]):
                        b=b[:,1:];
                    fr_m[0][i]=interpolate(omega,b[:,1].copy(),o)
                    fr_r[0][i]=interpolate(omega,b[:,2].copy(),o)
                    fr_m[1][i]=interpolate(omega,b[:,3].copy(),o)
                    fr_r[1][i]=interpolate(omega,b[:,4].copy(),o)
                    fr_m[2][i]=interpolate(omega,b[:,5].copy(),o)
                    fr_r[2][i]=interpolate(omega,b[:,6].copy(),o)
                
                
                b=np.array(pandas.read_csv(entry.path+"/results/IRF.tec",header=4+(1+len(omega))*i,nrows=len(omega),sep="\s+|,+"));
                for i in np.arange(3):
                    #fr_inf[i]=b[1,1+i*2];#Not working wth current data
                    fr_inf[i]=fr_m[i][i][-1];
                LUT[str((int)(h))+"_"+str((int)(p))]=[fe_a,fe_p,fr_m,fr_r, fr_inf]            
    print("LUT table h spacing: "+str(h_min)+ ": "+str(h_step)+" : "+str(h_max));
    print("LUT table p spacing: "+str(p_min)+ ": "+str(p_step)+" : "+str(p_max)+"\n");
# get parameters for body pose specified by heave and pitch (only heave working)
# interpolate (linearly) between values
def get_fromLUT(h,p):
    global p_step, h_step;
    p=np.min([p_max,p]);
    p=np.max([p_min,p]);
    h=np.min([h_max,h]);
    h=np.max([h_min,h]);
    p1=0#(int)(np.floor(p/p_step))*0;
    h1=(int)(np.floor(h/h_step));
    p2=0#np.min([(int)(np.ceil(p/p_step)),0])*0;
    h2=(int)(np.ceil(h/h_step));
    
    #p2=np.min([(int)(p_max/p_step),p2]);
    #p1=np.max([(int)(p_min/p_step),p1]);
    h2=np.min([(int)(h_max/h_step),h2]);
    h1=np.max([(int)(h_min/h_step),h1]);
    
    d1=LUT[str((int)((h1*h_step)*prec))+"_"+str((int)((p1*p_step)*prec))]
    d2=LUT[str((int)((h2*h_step)*prec))+"_"+str((int)((p1*p_step)*prec))]
    d3=LUT[str((int)((h1*h_step)*prec))+"_"+str((int)((p2*p_step)*prec))]
    d4=LUT[str((int)((h2*h_step)*prec))+"_"+str((int)((p2*p_step)*prec))]
    d5=LUT[str((int)((h2*h_step)*prec))+"_"+str((int)((p2*p_step)*prec))]

    
    d=[[],[],[],[],[]]
    for idx, (e1,e2,e3,e4,e5) in enumerate(zip(d1,d2,d3,d4,d5)):
        if h_step==0:
            pc=0;
        else:
            pc=np.mod(h,h_step)/h_step;
        if p_step==0:
            pcp=0;
        else:
            pcp=np.mod(p,p_step)/p_step;
        e5=e1*(1-pc)+e2*pc;
        e6=e3*(1-pc)+e4*pc;
        d[idx]=e5*(1-pcp)+e6*(pcp);  
    return d;

if __name__=="__main__":
        
    load_LUT(np.array([0.1,0.2,0.3,0.4]));
    print("Result 1:");
    print(get_fromLUT(h_max/2, p_min/2));
    print("Result 2:");
    print(get_fromLUT(-2, p_min/2)[4]);
    print("Result 3:")
    print(get_fromLUT(0.001, p_min/2)[4]);