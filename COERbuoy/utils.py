#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:16:12 2021

@author: heiko
"""

import json;
import os;

def get():
    global pkg_dir;
    pkg_dir=os.path.dirname(__file__);
    
    fp=os.path.join(pkg_dir,"settings.txt");
    print(os.getcwd())
    
    if os.path.exists("coerbuoy_settings.txt"):
        fp="coerbuoy_settings.txt";
    with open(fp) as file:
        data=json.load(file);
        print(data)
        global class_hydro, wec_dir, conn_ip, conn_port, resolution, dt_controller, msg_status
        class_hydro=data.get("hydro","Floater_BEM");
        wec_dir=data.get("WECfolder","[COERbuoy.data.COERbuoy]");
        conn_ip=data.get("conn_ip","localhost");
        conn_port=data.get("conn_port",5050);
        resolution=data.get("resolution",0.1);
        dt_controller=data.get("dt_controller",0.1);
        msg_status=data.get("status_message",[1,1,1,1,1,1,1,1,1,1]);
     
    if wec_dir[0] == '[' and wec_dir[-1]==']':
        wec_dir=wec_dir[1:-1].replace('.','/');
        wec_dir=os.path.join(pkg_dir,wec_dir);
    a=0;
get();