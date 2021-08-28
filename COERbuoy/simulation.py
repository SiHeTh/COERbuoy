# COERbuoy1 WEC model
# A realistic, nonlinear Wave Energy Converter model
# 2020/2021 COER Laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie
#

# Main file

# Dependencies
import sys;
import numpy as np;
import json;
import time;
import COERbuoy.utils as utils;
from COERbuoy import connection;
#from COERbuoy import dynamics_coerbuoy as dyn_wec;
import pandas;
from COERbuoy import floater_LIN as Floater;
#from COERbuoy import floater_BEM_LUT as Floater;
from COERbuoy import wavefield;
from scipy.interpolate import interp1d;
#omega_cut_off=dyn_wec.omega_cut_off;
omega_cut_off=4;
import os;
pkg_dir=os.path.dirname(__file__);
# Debugging feature
debug=not True;

# Main program
def start_simu (**kwargs):

    host=True;
    interface=False;
    
    conn_ctrl=connection.connection();
    
    
    import subprocess;
    import time
    import importlib
    from scipy.integrate import solve_ivp
    
    
    pi=np.pi;
    dt=utils.resolution;
    process=0;
    ctrl="";
    #read settings
    #with open(os.path.join(pkg_dir,"settings.txt")) as file:
    #    data=json.load(file);
    #    class_hydro=data.get("hydro","Floater");
    #    WECfolder=data.get("WECfolder","COERbuoy.data.COERbuoy");
   
    #dyn_wec=importlib.import_module(utils.wec_dir+".dynamics","COERbuoy");
    spec=importlib.util.spec_from_file_location("dynamics.py",os.path.join(utils.wec_dir,"dynamics.py"));
    dyn_wec=importlib.util.module_from_spec(spec);
    spec.loader.exec_module(dyn_wec);
    
    print("Using the following WEC: "+utils.wec_dir);
    #Set filename    
    filename="output.csv"
    if "name" in kwargs:
        filename=kwargs["name"];
    #Floater.idname=filename;
    
    wec=dyn_wec.WEC();
    init_condition=np.zeros(wec.states).tolist()
    
    #Set host or client mode
    if "host" in kwargs:
        if kwargs["host"]==False:
            host=False;
    
    # Define control mode and eventually start control
    # "TCP": A TCP/IP connection is opened, the controller is started manually
    # "linear": A constant velocity damping is applied; The external control interface is not used; Implemented for debugging/testing
    # [string]: String specifying control command (f.ex.: "python3 MBC.py"), a TCP/IP socket is opened and the control program started as a subprocess
    if "control" in kwargs:
        if kwargs["control"]=="TCP":
            interface=True;
        elif kwargs["control"]=="linear":
            interface=False;
        elif kwargs["control"]!="":
            ctrl=kwargs["control"].split(" ");
            #if True:#(ctrl[0]=="octave"):
                #host=True;
                #print("host mode");
                #if ("host" in kwargs) and  (kwargs["host"]=="False"):
                    #host=False;
            print(host)
            if not host:
                print("Start "+ctrl[0]+" "+ctrl[1])
                process=subprocess.Popen(ctrl)
                ctrl="";
                time.sleep(2);
            interface=True;
    teval=0;
    
    # t0 specify the transient time to get the device in steady state
    # for t>t0 the absorbed power is measured
    if "t0" in kwargs:
        teval=kwargs["t0"];
        
    # The "buoy_file" contains the parameters of the WEC
    if "buoy_file" in kwargs:
        wec.file=kwargs["buoy_file"];
        
    # Initial condition (position, velocity); used for decay test
    if "init" in kwargs:
        init_condition=kwargs["init"]+np.zeros(wec.states-len(kwargs["init"])).tolist();
        
    # Wave data (wave elevation y over time t)...
    if "time" in kwargs and "wave" in kwargs:
        t=kwargs["time"];
        y=kwargs["wave"];
        
    # ... or read wave data from file.
    elif "file" in kwargs:
        a=np.array(pandas.read_csv(kwargs["file"]));#"TestFull.csv", header=None))
        t=a[:,0];
        y=a[:,1];
        
    # if not specify, use default settings
    else:   
        t=np.linspace(0,200,1000)
        y=np.sin(1*t)
    
       
    #import matplotlib.pyplot as pyplt;
    
    omega_cut_off=wec.omega_cut_off;
    # Setting the wave (pocessed in a seperate module)
    wave1=wavefield.wavefield.set_wave(y,t,omega_cut_off);
    
    # Calculating buoy data
    wec.load_buoy(getattr(Floater,utils.class_hydro),wave1.xi,300,0);
    # Set the time steps for which we want the solution from ODE solver, equally spaced with dt
    steps=np.array(range(0,int(1/dt*t[-1])))*dt;
    
    # WEC reads in its parameters
    wec.load_param();
    
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
            if (t>self.t[-1]+0.05):
                self.t.append(t);
                for (i,e) in enumerate(x):
                    self.x[i].append(e);
        def get(self,i,t):
            if (len(self.t)<2):
                return t*0;
            fin=interp1d(self.t,self.x[i],fill_value=0,bounds_error=False);
            return fin(t)
          
    
    # This is the ODE function
    def dynamics(t,x):

      # to provide past sensor measures to the control interface, we need previous values of the ODE solution
      # Herefore DDE (delayed differntial euqation) exists. 
      # However, the implementations tested showed to be several orders of magnitude slower than the ODE-solver
      # Needing the past values only at specific time steps, the values are stored and
      # interpolated manually:
      
      dynamics.xlast.add(t,x.tolist()+[wec.get_force(x)]);
      
      
      # While the ODE solver can jump forth and back in time; the control is only executed
      # the first time a specific time is passed
      # This is not _exact_, but a good compramise between computational complexity and accuracy
      if t-dynamics.tcontrol>0.19 or t-dynamics.tcontrol<0:
          
          tw=100;
          tw1=tw-1;
          tseq=np.linspace(t-tw1/4,t,int(tw));
          dynamics.tcontrol=t;
          
          # In case the TCP/IP control interface is used (normal operation)
          if interface:
              msg_status=utils.msg_status;
              # Prepare message sent to controller
              msg={"time":tseq*msg_status[0],
                   "wave":(np.sum(wave1.get(tseq.reshape(tseq.size,1),0)[0],1))*msg_status[1],
                   "wave_forecast":(np.sum(wave1.get(-1*tseq.reshape(tseq.size,1),0)[0],1))*msg_status[2],
                   "stroke_pos":dynamics.xlast.get(0,tseq)*msg_status[3],
                   "stroke_speed":dynamics.xlast.get(1,tseq)*msg_status[4],
                   "angular_pos":dynamics.xlast.get(2,tseq)*msg_status[5],
                   "angular_speed":dynamics.xlast.get(3,tseq)*msg_status[6],
                   "force":dynamics.xlast.get(-1,tseq),
                   "test":np.zeros(100)
                   
                   }
              
              # Exchange data with controller
              data=conn_ctrl.exchange_model(msg["time"],msg["wave"],msg["wave_forecast"],msg["stroke_pos"],msg["stroke_speed"],msg["angular_pos"],msg["angular_speed"],msg["force"],msg["test"])
              
              dynamics.PTOt=data["time"]
              dynamics.PTO=data["pto"]
              dynamics.brake=data["brake"]
          else:
              # if the TCP/IP control interface is not used, apply a velocity dependent damping
              # (for testing)
              dynamics.PTOt=[t]
              dynamics.PTO=[-wec.get_translator_speed(x)*wec.damping];
              dynamics.brake=[0];
      
      i=np.abs(np.array(dynamics.PTOt)-t).argmin();
      # Send the data to the WEC
      out=wec.Calc(t,wave1,x,dynamics.PTO[i],dynamics.brake[i],dynamics.ulast);
      dynamics.ulast=wec.get_translator_speed(x);
      return out;    
    dynamics.ts=np.linspace(0,t[-1]+4,int((t[-1]+4)*wec.buoy.omega[-1]*2));
    dynamics.ulast=0;
    dynamics.latch_counter=0;
    dynamics.tcontrol=0;
    dynamics.PTO=[0];
    dynamics.PTOt=[0];
    dynamics.brake=[0];
    dynamics.xlast=history([0]*(len(init_condition)+1));
    
    if interface:
        if host:
            if (ctrl!=""):
                print("Start "+ctrl[0]+" "+ctrl[1]);
                process=subprocess.Popen(ctrl);            
                ctrl="";
            conn_ctrl.openH();
        else:
            conn_ctrl.openC();
    
    print("Start solver")
    t_start=time.time();

    # Start the ODE-solver
    sol = solve_ivp(dynamics,[0,t[-1]],init_condition,t_eval=steps.tolist(),max_step=0.5*1/(wec.buoy.omega[-1]*2),rtol=0.8,atol=0.8)#state vecor[z, dz, x, dx, delta, ddelta, slidex, dslidex]
    print("Elapsed time :"+str(time.time()-t_start)+"\n");

    # Write the solution of the data frame
    pandas.DataFrame(np.array([sol.t[:],np.sum(wave1.get(sol.t.reshape(sol.t.size,1),0)[0],1),sol.y[0,:],sol.y[1,:],sol.y[2,:]*180/pi,sol.y[3,:]*180/pi,sol.y[7,:],sol.y[8,:]]).transpose(),columns=["time","wave [m]","stroke [m]","stroke speed [m/s]","angle [deg]","angular_speed [deg/s]","F_PTO [N]","Energy [J]"]).round(3).to_csv(filename,index=False)
    s=np.argmax(sol.t>teval)

    # free data from the wave and teh WEC
    wave1.clear();
    del wave1;
    wec.release();
    
    # Output the absorbed power
    power=(sol.y[8,-1]-sol.y[8,s])/(sol.t[-1]-sol.t[s])
    print("Absorbed power: "+str((power/1000).round(2))+" kW")
    
    
    # CLose the TCP/IP control interface connection
    if interface:
        conn_ctrl.close()
        if process != 0:
            process.wait()
        time.sleep(5);
        
    # return the generated electrical energy    
    return power.round(2);

# Decay test (Init state, name, duration, control)
def decay_test(x0, n, t, ctrl):
    t2=np.arange(0,t,1/(omega_cut_off))
    y=0*t2;
    return start_simu(time=t2,wave=y,name=n, t0=0, init=[x0, 0, 0.5], control=ctrl )#/(9.81*9.81*1000/(32*np.pi)*A**2*p)

# Regular wave (Height, period, name, control)
def reg_wave(H,p,n,ctrl):
    print("Regular wave")
    H=float(H);
    p=float(p);
    t0=np.max([p*2,0]);#4
    if (p<6):
        t0=p*4;
    t2=np.max([t0+p*3,0]);#6
    t=np.arange(0,t2,1/(omega_cut_off))
    y=H/2*np.sin(2*np.pi/p*t)
    return start_simu(time=t,wave=y,name=n, t0=t0, control=ctrl )#/(9.81*9.81*1000/(32*np.pi)*H**2*p)


# Brettschneider wave (significant wave height, energy period, name, control)
def bretschneider_wave(Hs,p,n,ctrl):
    print("Brettschneider wave")
    Hs=float(Hs);
    p=float(p);
    omega=np.linspace(0.001,4,200);
    omega_m=2*np.pi/(p/0.856);
    S=5/16*(omega_m**4)/(omega**5)*(Hs**2)*np.exp(-5*(omega_m**4)/4/(omega**4))
    t0=np.max([p*3,120]);
    t2=np.max([t0+p*6]);
    t=np.arange(0,t2,1/(omega_cut_off/np.pi))
    np.random.seed(6)#Maybe replace by fixed phase vector; has to guaranteed that rnadom sequecne is always the same
    phase=np.random.rand(omega.size)*2*np.pi;
    y=np.sum(np.sqrt(2*S*(omega[1]-omega[0]))*np.sin(omega*t.reshape(t.size,1)+phase),1)
    
        
    return start_simu(time=t,wave=y,name=n, t0=t0, control=ctrl )


def quartil (a, p):
    a=np.abs(a);
    b=np.sort(a);
    aa=np.sum(a);
    if (aa==0):
        return 0;
    return b[np.argmax(np.cumsum(b)/aa>p)];


if __name__=="__main__":
        
    #Few examples how to run different tests:
    t=np.linspace(0,10,100);
    #start_simu(time=t, wave=np.sin(t/10), name="test", t0=0, control="TCP", host=False);
    reg_wave(4,3.5,"test.csv","TCP");#"python3 /media/heiko/Windows/Control_Maynooth/COERbuoy/COERbuoy/examples/custom_controller/controller.py")
    #reg_wave(4,3.5,"test.csv","linear")
    #decay_test(0.15,"decay1.csv",10,"linear")
    #reg_wave(1,4,"output.csv","linear")
    #bretschneider_wave(1.5,12,"bretschneider_wave.csv","python3 TestController1.py")
     
       
