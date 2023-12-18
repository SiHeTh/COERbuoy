#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# COERbuoy1 WEC model
# A realistic, nonlinear Wave Energy Converter model
# 2020/2021/2022 COER laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie
#

# Main file

# Dependencies
import numpy as np;
import json;
import re;
import warnings;
import COERbuoy.utils as utils;
from COERbuoy import connection;
from COERbuoy import floater_LIN as Floater;
import pandas;
from COERbuoy import wavefield;
from scipy.interpolate import interp1d;
omega_cut_off=4;
import os;
pkg_dir=os.path.dirname(__file__);

# Debugging feature
debug=not True;

class wave_series():#power is measured starting from t=0
                    #The simulation runtime te is the last element of the time series
                    #Unless the last element, the time series has to be monotonary increasing 
                    #Having a runtime shorter than the wave length is recommended when windowing
    t=[];
    y=[];
    t0=0;
    te=0;
    name="";
    def __init__(self,t,y,n="wave"):
        self.t=t;
        self.y=y;
        self.t0=np.min([t[0],0]);
        self.te=t[-1];
        self.name=n;
        if t[-2]>t[-1]:
            self.t=t[:-1];
            self.y=y[:-1];
        self.t=self.t-self.t0;
        self.te=self.te-self.t0;
       
    @classmethod
    def fromLists (cls,t,y,n=""):
        return cls(t,y,n);
    @classmethod
    def fromFile (cls,filename):
        a=np.array(pandas.read_csv(filename));#"TestFull.csv", header=None))
        return cls(a[:,0],a[:,1],filename)
    def to_file(self,name):
        t=self.t+self.t0;
        y=self.y;
        ye=np.where(np.array(self.t)>self.te*0.98)[0];
        if len(ye) > 0:
            ye=self.y[ye[0]];
            t=np.hstack([t,[self.te+self.t0]]);
            y=np.hstack([self.y,[ye]]);
        pandas.DataFrame(np.vstack((t,y)).transpose(),columns=["time","wave-elevation"]).round(6).to_csv(name,index=False)
       

#new main function        
def run (cWEC, cHyModel, wave_data, connection, omega_cut_off, filename):
    
    import time
    from scipy.integrate import solve_ivp
    
    
    init_condition=np.zeros(cWEC.states).tolist()
    
    # Setting the wave (processed in a seperate module)
    wave=wavefield.wavefield.set_wave(wave_data.y,wave_data.t,omega_cut_off);
    
    # Calculating buoy data
    cWEC.load_buoy(getattr(Floater,cHyModel),wave.xi,3000000,0);
    # Set the time steps for which we want the solution from ODE solver, equally spaced with dt
    dt=utils.resolution;
    steps=np.array(range(0,int(1/dt*wave_data.te)))*dt;
    
    # WEC reads in its parameters
    cWEC.load_param();
    
    class history():
        def __init__(self,x):
            self.t=[0];
            self.x=[];
            for (i,e) in enumerate(x):
                self.x.append([e]);
        def add(self,t,x):
            #while len(self.t)>0 and self.t[-1]+0.05>t:
            #    self.t.pop();
            #    for (i,e) in enumerate(x):
            #        self.x[i].pop();
            if (t>self.t[-1]+0.0000000125):
                self.t.append(t);
                for (i,e) in enumerate(x):
                    self.x[i].append(e);
        def get(self,i,t):
            if (len(self.t)<2):
                return t*0;
            fin=interp1d(self.t,self.x[i],fill_value=0,bounds_error=False);
            res=fin(t);
            res[-1]=self.x[i][-1];
            return res;
          
    
    # This is the ODE function
    def dynamics(t,x):

      # to provide past sensor measures to the control interface, we need previous values of the ODE solution
      # Herefore DDE (delayed differntial euqation) exists. 
      # However, the implementations tested showed to be several orders of magnitude slower than the ODE-solver
      # Needing the past values only at specific time steps, the values are stored and
      # interpolated manually:
      
      dynamics.xlast.add(t,x.tolist()+[cWEC.get_force(x)]);
      
      
      # While the ODE solver can jump forth and back in time; the control is only executed
      # the first time a specific time is passed
      # This is not _exact_, but a good compramise between computational complexity and accuracy
      if utils.dt_controller<0 or t-dynamics.tcontrol>=utils.dt_controller-0.05*utils.dt_controller or t-dynamics.tcontrol<0:
          tw=100;
          tw1=tw-1;
          tseq=np.linspace(t-tw1/4,t,int(tw));
          dynamics.tcontrol=t;
          
          # In case the TCP/IP control interface is used (normal operation)
          if connection:
              msg_status=utils.msg_status;
              # Prepare message sent to controller
              msg={"time":tseq*msg_status[0],
                   "wave":(np.sum(np.real(wave.get(tseq.reshape(tseq.size,1),0)[0]),1))*msg_status[1],
                   "wave_forecast":(np.sum(np.real(wave.get(2*t-np.flip(tseq.reshape(tseq.size,1)),0)[0]),1))*msg_status[2],
                   ##Wave dependent on the position not time; can be used for visualisation
                   #"wave":(np.sum(wave1.get(t,np.linspace(60,0,int(tw)).reshape(tseq.size,1))[0],1))*msg_status[1],
                   #"wave_forecast":(np.sum(wave1.get(t,np.linspace(0,-60,int(tw)).reshape(tseq.size,1))[0],1))*msg_status[2],
                   "stroke_pos":dynamics.xlast.get(0,tseq)*msg_status[3],
                   "stroke_speed":dynamics.xlast.get(1,tseq)*msg_status[4],
                   "angular_pos":dynamics.xlast.get(2,tseq)*msg_status[5],
                   "angular_speed":dynamics.xlast.get(3,tseq)*msg_status[6],
                   "force":dynamics.xlast.get(-1,tseq)*msg_status[7],
                   "test":(np.real(np.sum(wave.get(t,np.linspace(-14,14,int(tw)).reshape(tseq.size,1))[0],1)))*msg_status[8]#(np.sum(wave1.get(t,np.linspace(30,-30,int(tw)).reshape(tseq.size,1))[0],1))*msg_status[8]
                   
                   }
              # Exchange data with controller
              data=connection.exchange_model(msg["time"],msg["wave"],msg["wave_forecast"],msg["stroke_pos"],msg["stroke_speed"],msg["angular_pos"],msg["angular_speed"],msg["force"],msg["test"])
              
              dynamics.PTOt=data["time"]
              dynamics.PTO=data["pto"]
              dynamics.brake=data["brake"]
              
          else:
              # if the TCP/IP control interface is not used, apply a velocity dependent damping
              # (for testing)
              dynamics.PTOt=[t]
              dynamics.PTO=[-cWEC.get_translator_speed(x)*dynamics.damping];
              dynamics.brake=[0];
      #Print progress
      print("{:.0f}".format(100*t/dynamics.duration,0)+"% completed", end="\r");
      i=np.abs(np.array(dynamics.PTOt)-t).argmin();
      # Send the data to the WEC
      out=cWEC.Calc(t,wave,x,dynamics.PTO[i],dynamics.brake[i],dynamics.ulast);
      dynamics.ulast=cWEC.get_translator_speed(x);
      return out;    

    dynamics.ts=np.linspace(0,wave_data.t[-1]+4,int((wave_data.t[-1]+4)*cWEC.buoy.omega[-1]*2));
    dynamics.ulast=0;
    dynamics.latch_counter=0;
    dynamics.tcontrol=0;
    dynamics.PTO=[0];
    dynamics.damping=cWEC.damping;#TODO: No damping not implemented yet
    dynamics.PTOt=[0];
    dynamics.brake=[0];
    dynamics.duration=wave_data.te;
    dynamics.xlast=history([0]*(len(init_condition)+1));
    
   
    
    print("\nRunning the simulation...")
    t_start=time.time();

    # Start the ODE-solver
    sol = solve_ivp(dynamics,[0,wave_data.te],init_condition,t_eval=steps.tolist(),max_step=utils.ode_time_step,rtol=100,atol=100)#state vecor[z, dz, x, dx, delta, ddelta, slidex, dslidex]
    print("Elapsed time :"+str(time.time()-t_start)+"\n");

    
    #For debugging: plot simulation results
    #import matplotlib.pyplot as plt
    #plt.figure();
    #plt.plot(sol.t[:],np.transpose([np.sum(wave1.get(sol.t.reshape(sol.t.size,1),0)[0],1),sol.y[0,:]]));
    #plt.show();
    
    sol.t=sol.t+wave_data.t0;

    
    # Write the solution of the data frame
    pandas.DataFrame(np.array([sol.t[:],np.sum(np.real(wave.get(sol.t.reshape(sol.t.size,1),0)[0]),1),sol.y[0,:],sol.y[1,:],sol.y[2,:]*180/np.pi,sol.y[3,:]*180/np.pi,sol.y[7,:],sol.y[8,:]]).transpose(),columns=["time [s]","wave [m]","stroke [m]","stroke speed [m/s]","angle [deg]","angular_speed [deg/s]","F_PTO [N]","Energy [J]"]).round(3).to_csv(filename,index=False)
    
    
    
    # free data in memory
    wave.clear();
    del wave;
    cWEC.release();
    
    return sol;
         
# Main program
def start_simu (**kwargs):
    host=True;
    interface=False;
    utils.get();
    conn_ctrl=None;

    import subprocess;
    import time
    import importlib
    
    
    damping=0;
    pi=np.pi;
    process=0;
    ctrl="";
    ctrlcmd="";
    wavedata=None;
    
    #load WEC module    
    if "wec" in kwargs:
        utils.wec_dir=utils.WECpath(kwargs["wec"]);
    spec=importlib.util.spec_from_file_location("dynamics.py",os.path.join(utils.wec_dir,"dynamics.py"));
    dyn_wec=importlib.util.module_from_spec(spec);
    spec.loader.exec_module(dyn_wec);
    
    
    if "hydro" in kwargs:
        utils.class_hydro=utils.WECpath(kwargs["wec"]);
    
    print("Using the following WEC: "+utils.wec_dir);

    
    wec=dyn_wec.WEC();
    
    #Set host or client mode
    if "host" in kwargs:
        if kwargs["host"]==False:
            host=False;
    
    # Define control mode and eventually start control
    # "TCP": A TCP/IP connection is opened, the controller is started manually
    # "linear": A constant velocity damping is applied; The external control interface is not used; Implemented for debugging/testing
    # "none": No control/damping is applied
    # [string]: String specifying control command (f.ex.: "python3 MBC.py"), a TCP/IP socket is opened and the control program started as a subprocess
    ctrlname=kwargs["control"];
    if "control" in kwargs:
        if kwargs["control"]=="TCP":
            interface=True;
        elif kwargs["control"]=="" or kwargs["control"]=="linear" or kwargs["control"]=="none":
            interface=False;
            damping=0;
        elif kwargs["control"]!="":
            ctrl0=kwargs["control"].split(".",1);
            extension=ctrl0[1].split(" ")[0];
            args=ctrl0[0].split(" ")[1:];
            ctrlcmd=ctrl0[0].split(" ");
            ctrlname=ctrlcmd[-1]+"."+extension;
            try:
                fullpath=utils.get_controller(ctrlname);
            except:
                fullpath=ctrlname;
            
            if len(ctrlcmd)>1: #ctrl command already provided
                ctrl=kwargs["control"].replace(ctrlname,fullpath).split(" ");
    
            if len(ctrlcmd)==1: #get control command from extension
                ctrlcmd=utils.cmddict.get(extension,"");
                ctrl=[ctrlcmd]+kwargs["control"].replace(ctrlname,fullpath).split(" ");
                if ctrlcmd=="":
                    ctrl=kwargs["control"].replace(ctrlname,fullpath).split(" ");
                
            if not host:
                print("Start "+str(" ".join(ctrl)))
                process=subprocess.Popen(ctrl)
                ctrl="";
            time.sleep(3);
            interface=True;
    teval=0;
    
    if interface:
        ip=utils.conn_ip;
        port=utils.conn_port;
        if "ip" in kwargs:
            ip=kwargs["ip"];
        if "port" in kwargs:
            port=int(kwargs["port"]);
        conn_ctrl=connection.connection(ip = ip, port = port);
        
    # t0 specify the transient time to get the device in steady state
    # for t<t0 the data is not logged
    if "t0" in kwargs:
        warnings.warn("parameter t0 is ignored; support ended in version 0.2.0")
        #wavedata.t0=kwargs["t0"];
        
    # The "buoy_file" contains the parameters of the WEC
    if "buoy_file" in kwargs:
        wec.file=kwargs["buoy_file"];
        
    # Initial condition (position, velocity); used for decay test
    if "init" in kwargs:
        init_condition=kwargs["init"]+np.zeros(wec.states-len(kwargs["init"])).tolist();
        
    # Wave data (wave elevation y over time t)...
    if "time" in kwargs and "wave" in kwargs:
        wavedata=wave_series.fromLists(kwargs["time"],kwargs["wave"])
        
    # ... get a wave series object directly ...
    elif "wave" in kwargs:
        wavedata=kwargs["wave"];
            
        
    # ... or read wave data from file.
    elif "file" in kwargs:
        wavedata=wave_series.fromFile(kwargs["file"]);
        
    # if not specify, use default settings
    else:
        wavedata.fromList(t=np.linspace(0,200,1000),y=np.sin(1*t))   
        
    
    if wavedata.te<=wavedata.t0:
        wavedata.te=wavedata.t[-1];
    
    
    omega_cut_off=wec.omega_cut_off;
    
    #Set filename
    wave1=wavefield.wavefield.set_wave(wavedata.y,wavedata.t,omega_cut_off);
    filename=wavedata.name.replace("_","")+"_p_"+f'{wave1.get_period():.2f}'+"_h_"+f'{wave1.get_height():.2f}'+"_"+re.split("[\\,/,\s]",ctrlname.replace("_","").replace(".",""))[-1]+"_"+re.split("[\\,/,.,-,\s]",utils.wec_dir.replace("_",""))[-1]+"_"+utils.class_hydro.replace("_","")+".csv";#Standard format
    wave1.clear();
    if "name" in kwargs: #use custom file name if specified
        if os.path.isdir(kwargs["name"]):
            filename=os.path.join(kwargs["name"],filename);
        elif kwargs["name"]=="":
            filename=filename;
        else:
            filename=kwargs["name"];

    if connection:
        print("Using control interface with ip "+conn_ctrl.ip+" at port " + str(conn_ctrl.port) +".")
        if host:
            if (ctrl!=""):
                if ctrlcmd=="": #if no extension, or extension unknown, we assume it is an executaböe
                    print("Start "+str(" ".join(ctrl)))
                    process=subprocess.Popen(ctrl)
                    ctrl="";
            conn_ctrl.openH();
        else:
            conn_ctrl.openC();

    #####----------------
    sol=run(wec,utils.class_hydro,wavedata,conn_ctrl,omega_cut_off,filename);
    
    
    #Cut data to exclude transient data (if applicable) for power calculation
    s1=np.argmax(sol.t>=0)
    # Print the absorbed power
    power=(sol.y[8,-1]-sol.y[8,s1])/(sol.t[-1])
    print("Mean absorbed power: "+str((power/1000).round(2))+" kW")
    
    
    # CLose the TCP/IP control interface connection
    if conn_ctrl:
        conn_ctrl.close()
        if process != 0:
            process.wait()
        time.sleep(5);

        
    # return the generated electrical energy    
    return [sol,filename,power.round(2)];

# Decay test (Init state, name, duration, control)
def decay_test(x0, n, t, ctrl):
    t2=np.arange(0,t,1/(omega_cut_off*np.pi))
    y=0*t2;
    if isinstance(x0,list):
        return start_simu(time=t2,wave=y,name=n, t0=0, init=x0, control=ctrl )#/(9.81*9.81*1000/(32*np.pi)*A**2*p)    
    return start_simu(time=t2,wave=y,name=n, t0=0, init=[x0, 0, 0], control=ctrl )#/(9.81*9.81*1000/(32*np.pi)*A**2*p)

 
# Regular wave (Height, period)
def reg_wave(H=1,p=10,n0=8,n=8,ne=1):
    print("Regular wave")
    H=float(H);
    p=float(p);
    t0=-p*n0;
    t2=p*n;
    t=np.arange(t0,t2,1/(omega_cut_off*np.pi))
    y=H/2*np.cos(2*np.pi/p*(t-t0))
    return wave_series.fromLists(np.append(t,[t2-ne*p]),np.append(y,y[-1]),"regular");

# Bretschneider wave (significant wave height, energy period, name, control)
def bretschneider_wave(Hs=1,p=6,n0=4,n=10,ne=1):
    print("Bretschneider wave")
    Hs=float(Hs);
    p=float(p);
    omega=np.linspace(0.001,4,200);
    omega_m=2*np.pi/(p/0.856);
    S=5/16*(omega_m**4)/(omega**5)*(Hs**2)*np.exp(-5*(omega_m**4)/(4*omega**4))
    t0=-p*n0;
    t2=p*n;
    t=np.arange(t0,t2,1/(omega_cut_off*np.pi))
    np.random.seed(6)#Maybe replace by fixed phase vector; has to guaranteed that random sequecne is always the same
    phase=np.random.rand(omega.size)*2*np.pi;
    y=np.sum(np.sqrt(2*S*(omega[1]-omega[0]))*np.sin(omega*t.reshape(t.size,1)+phase),1)
    return wave_series.fromLists(np.append(t,[t2-ne*p]),np.append(y,y[-1]),"bretschneider");

def get_WEC_data():
    f = open(utils.wec_dir+"/floater.txt",'r')
    data=json.load(f);
    f.close();
    return data; 

if __name__=="__main__":
        
    #Few examples how to run different tests:
    t=np.linspace(0,10,100);
    #start_simu(time=t, wave=np.sin(t/10), name="test", t0=0, control="TCP", host=True);
    #decay_test([1, 0, 0],"decay1.csv",20,"linear")
    start_simu(time=t,wave=np.sin(t/10),name="test123", t0=0, init=[0, 0, 0], control="linear" , wec="[data.COERsimple]", hydro="Floater_BEM")
    #bretschneider_wave(1,3).to_file("testwave1.csv")
    #start_simu(wave=reg_wave(1,5),control="controller_damping.py 5")
    #start_simu(file="testwave1.csv",name="output.csv",control="linear")
    #bretschneider_wave(1.5,12,"bretschneider_wave.csv","python3 controller.py")
     
   
