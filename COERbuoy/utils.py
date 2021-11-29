#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:16:12 2021

@author: heiko
"""

import json;
import os;
import subprocess;
import importlib;
import numpy as np;

global pkg_dir;
pkg_dir=os.path.dirname(__file__);

def get():
    global pkg_dir;#Package location

    global cmddict;#Program <-> filename
    cmddict={};

    
    #read which command to use for which file type
    try:
        with open(os.path.join(pkg_dir,"stdcmds.txt")) as file:
            cmddict=json.load(file);
    except (FileNotFoundError):
        get_command_for_extension();
        
        
    
    def read_settings(fp): #read settings from file
        with open(fp) as file:
            data=json.load(file);
            #print(data)
            global class_hydro, ode_time_step, wec_dir, wec_dir0, conn_ip, conn_port, resolution, dt_controller, msg_status, user_dir;
            class_hydro=data.get("hydro",class_hydro);
            wec_dir=data.get("WECfolder",wec_dir);
            wec_dir0=data.get("WECfolder_ideal",wec_dir);
            conn_ip=data.get("conn_ip",conn_ip);
            conn_port=data.get("conn_port",conn_port);
            user_dir=data.get("user_dir",user_dir);
            resolution=data.get("resolution",resolution);
            ode_time_step=data.get("ODE_time_step",ode_time_step);
            dt_controller=data.get("dt_controller",dt_controller);
            msg_status=data.get("status_message",msg_status);
         
    #set default values
    global class_hydro, ode_time_step, wec_dir, wec_dir0,conn_ip, conn_port, resolution, dt_controller, msg_status, controller0, results0, user_dir;
    user_dir=os.path.join(os.path.expanduser("~"),"COERbuoy_data")
    class_hydro="floater_bem";
    wec_dir="[data.COERbuoy1]";
    wec_dir0=wec_dir;
    conn_ip="localhost";
    conn_port=5050;
    ode_time_step=0.05;
    resolution=0.1;
    dt_controller=0.1;
    msg_status=[1,1,1,1,1,1,1,1,1,1];
    
    #read settings from package dir
    read_settings(os.path.join(pkg_dir,"settings.txt"));
    controller0=os.path.join(user_dir,"controller")
    results0=os.path.join(user_dir,"results")
    
    #if available, read settings from working dir
    if os.path.exists("coerbuoy_settings.txt"):
        read_settings("coerbuoy_settings.txt");
        
    #transfer in WEC dir into absolute path, if necessary
    wec_dir=WECpath(wec_dir);
    wec_dir0=WECpath(wec_dir0);
    
def WECpath(wecdir):
    #convention: [...] indicates folder in package dir
    if wecdir[0] == '[' and wecdir[-1]==']':
        wecdir=wecdir[1:-1].replace('.','/');
        wecdir=os.path.join(pkg_dir,wecdir);
    return wecdir;

        
#set a key value pair in settings
def set_settings(key,value):
        
    with open(pkg_dir+"/settings.txt",'r') as jfile:
        sets=json.load(jfile);
        sets[key]=value;
        with open(pkg_dir+"/settings.txt",'w') as jfile:
            jfile.write(json.dumps(sets));

def get_controller(ctrl_name=None):
    
    if ctrl_name!=None:
        for path in ["",controller0]:
            try:
                l=os.listdir(path);
                if ctrl_name in l:
                    return os.path.join(path,ctrl_name);
            except(FileNotFoundError):
                    ...
        raise Exception("ctrl file not found!");

    else:        
        lctrl=["linear","TCP","none"];
        for path in [controller0]:
            try:
                l=os.listdir(path);
                lctrl=lctrl+l;
            except(FileNotFoundError):
                os.makedirs(controller0);
                os.makedirs(results0);
        return lctrl;

    
def get_command_for_extension():
    
    #Obtain commands
    
    #test python cmd
    file="div/python.py"
    for cmd in ["py","python3","python"]:
        try:
            p=subprocess.Popen([cmd,os.path.join(pkg_dir,file)])
            if p.wait()==42:
                #correct answer
                cmddict["py"]=cmd;
                break;
        except:
            ...
            

    #test octave cmd
    file="div/octave.m"
    for cmd in ["octave","matlab"]:
        try:
            p=subprocess.Popen([cmd,os.path.join(pkg_dir,file)])
            if p.wait()==42:
                #correct answer
                cmddict["m"]=cmd;
                break;
        except:
            ...
            
    
    with open(os.path.join(pkg_dir,"stdcmds.txt"),'w') as file:
        file.write(json.dumps(cmddict));
    
        
