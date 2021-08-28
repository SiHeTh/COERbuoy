#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:16:12 2021

@author: heiko
"""

import json;
import os;

pkg_dir=os.path.dirname(__file__);

fp=os.path.join(pkg_dir,"settings.txt");

if os.path.exists("coerbuoy_settings.txt"):
    fp="coerbuoy_settings.txt";
with open(fp) as file:
    data=json.load(file);
    class_hydro=data.get("hydro","Floater_BEM");
    wec_dir=data.get("WECfolder","[COERbuoy.data.COERbuoy]");
    resolution=data.get("resolution",0.1);
    msg_status=data.get("status_message",[1,1,1,1,1,1,1,1,1,1]);
 
if wec_dir[0] == '[' and wec_dir[-1]==']':
    wec_dir=wec_dir[1:-1].replace('.','/');
    wec_dir=os.path.join(pkg_dir,wec_dir);
a=0;