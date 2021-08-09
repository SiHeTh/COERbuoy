#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:16:12 2021

@author: heiko
"""

import json;
import os;
import pkg_resources;

pkg_dir=os.path.dirname(__file__);

with open(os.path.join(pkg_dir,"settings.txt")) as file:
    data=json.load(file);
    class_hydro=data.get("hydro","Floater");
    wec_dir=data.get("WECfolder","[COERbuoy.data.COERbuoy]");
    resolution=data.get("resolution",0.1);
 
if wec_dir[0] == '[' and wec_dir[-1]==']':
    wec_dir=wec_dir[1:-1].replace('.','/');
    wec_dir=os.path.join(pkg_dir,wec_dir);
a=0;