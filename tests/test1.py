#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# test1.py - A simple test of the COERbuoy platform
#Checks overall functionality:
    # Parameter class
    # linear and look-up-body-exact hydrodynmaic models
    # and finally the overall system checking the absorbed power
# created 2022
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie

import COERbuoy.Parameters;
import COERbuoy.utils as u;
import numpy as np;

u.get();
if u.resolution==0.101:
    print("Correct value read!");
else:
    raise Exception("Settings are not correct!");

#Test to get hydrodynamic parameters
param1 = COERbuoy.Parameters.parameters("[data.COERsimple]","Floater_BEM");
param2 = COERbuoy.Parameters.parameters("[data.COERsimple]","Floater_LIN");
param1.init_hydro([0.1,0.2,0.3]); #select the frequencies
param2.init_hydro([0.1,0.2,0.3]); #select the frequencies
hy1=param1.hydro(0,1,1);#test at equilibrium for heave
hy2=param2.hydro(0,1,1);#test at equilibrium for heave

error=0;
for (e1,e2) in zip(hy1,hy2):
    error=error+np.sum(np.abs(e1-e2)/np.abs(e1))
if error<10:
    print("Hydrodynamic parameter read at equilibrium and no errors detected")
else:
    raise Exception("Hydro-Coefficents are not matching");

for z in [-0.2,0.2]:
    param1.hydro(0,1,1)
    param1.hydro(0,1,2)

print("Hydrodynamic parameter could be read!")

p1=param1.mdc(0)
p2=param2.mdc(0)
error=0;
for (e1,e2) in zip(p1,p2):
    error=error+np.sum(np.abs(e1-e2)/np.abs(e1))
if error<1:
    print("Could read mdc values; values are matching")
else:
    raise Exception("Hydro-Coefficents are not matching");
  
#Check absorbed power with theoretical limit
H=0.5;T=8;

omega=6.28/T;
k=omega**2/9.81;
lambda1=9.81/(6.28)*T**2;
vg=(omega/k);
vg=0.5*lambda1/T;
Pmax=1000*9.81*vg/(k)*1/2*(H/2)**2
from COERbuoy.simulation import start_simu, reg_wave, bretschneider_wave, decay_test;


power=start_simu(wave=reg_wave(H,T),control="controller_reactive.py "+str(T))[2];
p_ratio=power/Pmax;
print("Absorbed power: "+str(np.round(p_ratio*100))+" % of the theoretical limit");
if np.abs(p_ratio)<1.1 and np.abs(p_ratio)>0.9:
    print("Absorbed power within 10% range of the theoretical limit");
else:
    raise Exception("Power absorption is not correct!");

print("Test run passed!")