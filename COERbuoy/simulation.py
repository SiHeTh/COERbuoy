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
import numpy as np;
import json;
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

class wave_series:
    t=[];
    y=[];
    t0=0;
    te=0;
    def fromLists (self,t,y,t0=0,te=-1):
        if te<0:
            te=t[-1];
        self.t=t;
        self.y=y;
        self.t0=0;
        self.te=te;
    def fromFile (self,filename,t0=0,te=-1):
        a=np.array(pandas.read_csv(filename));#"TestFull.csv", header=None))
        self.t=a[:,0];
        self.y=a[:,1];
    def to_file(self,name):
        pandas.DataFrame(np.vstack((self.t,self.y)).transpose(),columns=["time","wave-elevation"]).round(6).to_csv(name,index=False)
       
# Main program
def start_simu (**kwargs):
    host=True;
    interface=False;
    utils.get();
    conn_ctrl=connection.connection(ip = utils.conn_ip, port = utils.conn_port);


    import subprocess;
    import time
    import importlib
    from scipy.integrate import solve_ivp
    
    
    damping=0;
    pi=np.pi;
    dt=utils.resolution;
    process=0;
    ctrl="";
    ctrlcmd="";
    wavedata=wave_series();
    #read settings
    #with open(os.path.join(pkg_dir,"settings.txt")) as file:
    #    data=json.load(file);
    #    class_hydro=data.get("hydro","Floater");
    #    WECfolder=data.get("WECfolder","COERbuoy.data.COERbuoy");
   
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
    # "none": No control/damping is applied
    # [string]: String specifying control command (f.ex.: "python3 MBC.py"), a TCP/IP socket is opened and the control program started as a subprocess
    if "control" in kwargs:
        if kwargs["control"]=="TCP":
            interface=True;
        elif kwargs["control"]=="linear":
            interface=False;
            damping=0;
        elif kwargs["control"]=="none":
            interface=False;
            damping=0;
        elif kwargs["control"]!="":
            ctrl0=kwargs["control"].split(" ");
            ctrlcmd=="";
            if len(ctrl0)>1: #ctrl command already provided
                ctrlcmd=ctrl0[0];
                ctrl=utils.get_controller(ctrl0[1]);
    
            if len(ctrl0)==1: #get control command from extension
                ctrl=utils.get_controller(ctrl0[0]);
                ext=ctrl0[0].split(".")
                if len(ext)>1: #check if there is file extension at all
                    ext=ext[-1];
                    ctrlcmd=utils.cmddict.get(ext,"");
                
            if not host:
                if ctrlcmd=="": #if no extension, or extension unknown, we assume it is an executaböe
                    print("Start "+ctrl)
                    process=subprocess.Popen(ctrl)
                    ctrl="";
                else:
                    print("Start "+ctrl[0]+" "+ctrl[1])
                    process=subprocess.Popen([ctrlcmd,ctrl])
                    ctrl="";
            time.sleep(2);
            interface=True;
    teval=0;
    
    # t0 specify the transient time to get the device in steady state
    # for t<t0 the data is not logged
    if "t0" in kwargs:
        wavedata.t0=kwargs["t0"];
    # te specify the transient time at the end (for a windoweds wave)
    if "te" in kwargs:
        wavedata.te=kwargs["te"];
        
    # The "buoy_file" contains the parameters of the WEC
    if "buoy_file" in kwargs:
        wec.file=kwargs["buoy_file"];
        
    # Initial condition (position, velocity); used for decay test
    if "init" in kwargs:
        init_condition=kwargs["init"]+np.zeros(wec.states-len(kwargs["init"])).tolist();
        
    # Wave data (wave elevation y over time t)...
    if "time" in kwargs and "wave" in kwargs:
        wavedata.fromList(kwargs["time"],kwargs["wave"])
        
    # ... get a wave series object directly ...
    elif "wave" in kwargs:
        wavedata=kwargs["wave"];
            
        
    # ... or read wave data from file.
    elif "file" in kwargs:
        wavedata.fromFile(kwargs["file"]);
        #a=np.array(pandas.read_csv(kwargs["file"]));#"TestFull.csv", header=None))
        #t=a[:,0];
        #y=a[:,1];
        
    # if not specify, use default settings
    else:
        wavedata.fromList(t=np.linspace(0,200,1000),y=np.sin(1*t))   
        
    
    if wavedata.te<=wavedata.t0:
        wavedata.te=wavedata.t[-1];
    
       
    #import matplotlib.pyplot as pyplt;
    
    omega_cut_off=wec.omega_cut_off;
    # Setting the wave (processed in a seperate module)
    wave1=wavefield.wavefield.set_wave(wavedata.y,wavedata.t,omega_cut_off);
    
    # Calculating buoy data
    wec.load_buoy(getattr(Floater,utils.class_hydro),wave1.xi,300,0);
    # Set the time steps for which we want the solution from ODE solver, equally spaced with dt
    steps=np.array(range(0,int(1/dt*wavedata.te)))*dt;
    
    # WEC reads in its parameters
    wec.load_param();
    
    if kwargs["control"]=="linear":
        damping=wec.damping;
            
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
      
      dynamics.xlast.add(t,x.tolist()+[wec.get_force(x)]);
      
      
      # While the ODE solver can jump forth and back in time; the control is only executed
      # the first time a specific time is passed
      # This is not _exact_, but a good compramise between computational complexity and accuracy
      if utils.dt_controller<0 or t-dynamics.tcontrol>=utils.dt_controller-0.05*utils.dt_controller or t-dynamics.tcontrol<0:
          #print(t-dynamics.tcontrol)
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
                   "wave_forecast":(np.sum(wave1.get(2*t-np.flip(tseq.reshape(tseq.size,1)),0)[0],1))*msg_status[2],
                   ##Wave dependent on the position not time; can be used for visualisation
                   #"wave":(np.sum(wave1.get(t,np.linspace(60,0,int(tw)).reshape(tseq.size,1))[0],1))*msg_status[1],
                   #"wave_forecast":(np.sum(wave1.get(t,np.linspace(0,-60,int(tw)).reshape(tseq.size,1))[0],1))*msg_status[2],
                   "stroke_pos":dynamics.xlast.get(0,tseq)*msg_status[3],
                   "stroke_speed":dynamics.xlast.get(1,tseq)*msg_status[4],
                   "angular_pos":dynamics.xlast.get(2,tseq)*msg_status[5],
                   "angular_speed":dynamics.xlast.get(3,tseq)*msg_status[6],
                   "force":dynamics.xlast.get(-1,tseq),
                   "test":np.zeros(100)
                   
                   }
              #print([msg["stroke_pos"][-1],x[0]])
              # Exchange data with controller
              data=conn_ctrl.exchange_model(msg["time"],msg["wave"],msg["wave_forecast"],msg["stroke_pos"],msg["stroke_speed"],msg["angular_pos"],msg["angular_speed"],msg["force"],msg["test"])
              
              dynamics.PTOt=data["time"]
              dynamics.PTO=data["pto"]
              dynamics.brake=data["brake"]
              
          else:
              # if the TCP/IP control interface is not used, apply a velocity dependent damping
              # (for testing)
              dynamics.PTOt=[t]
              dynamics.PTO=[-wec.get_translator_speed(x)*dynamics.damping];
              dynamics.brake=[0];
      print(str(np.round(100*t/dynamics.duration,0))+"% completed", end="\r");
      i=np.abs(np.array(dynamics.PTOt)-t).argmin();
      # Send the data to the WEC
      out=wec.Calc(t,wave1,x,dynamics.PTO[i],dynamics.brake[i],dynamics.ulast);
      dynamics.ulast=wec.get_translator_speed(x);
      return out;    

    dynamics.ts=np.linspace(0,wavedata.t[-1]+4,int((wavedata.t[-1]+4)*wec.buoy.omega[-1]*2));
    dynamics.ulast=0;
    dynamics.latch_counter=0;
    dynamics.tcontrol=0;
    dynamics.PTO=[0];
    dynamics.damping=damping;
    dynamics.PTOt=[0];
    dynamics.brake=[0];
    dynamics.duration=wavedata.t[-1];
    dynamics.xlast=history([0]*(len(init_condition)+1));
    
    if interface:
        print("Using control interface with ip "+conn_ctrl.ip+" at port " + str(conn_ctrl.port) +".")
        if host:
            if (ctrl!=""):
                if ctrlcmd=="": #if no extension, or extension unknown, we assume it is an executaböe
                    print("Start "+ctrl)
                    process=subprocess.Popen(ctrl)
                    ctrl="";
                else:
                    print("Start "+ctrlcmd+" "+ctrl)
                    process=subprocess.Popen([ctrlcmd,ctrl])
                    ctrl="";
            conn_ctrl.openH();
        else:
            conn_ctrl.openC();
    
    print("\nRunning the simulation...")
    t_start=time.time();

    # Start the ODE-solver
    sol = solve_ivp(dynamics,[0,wavedata.te],init_condition,t_eval=steps.tolist(),max_step=utils.ode_time_step,rtol=100,atol=100)#state vecor[z, dz, x, dx, delta, ddelta, slidex, dslidex]
    print("Elapsed time :"+str(time.time()-t_start)+"\n");

    #Cut data so that only data after transient times (if applicable) is considered
    s1=np.argmax(sol.t>wavedata.t0)
    ydata=sol.y[:,s1:];
    tdata=sol.t[s1:];
    # Write the solution of the data frame
    pandas.DataFrame(np.array([tdata[:],np.sum(wave1.get(tdata.reshape(tdata.size,1),0)[0],1),ydata[0,:],ydata[1,:],ydata[2,:]*180/pi,ydata[3,:]*180/pi,ydata[7,:],ydata[8,:]]).transpose(),columns=["time [s]","wave [m]","stroke [m]","stroke speed [m/s]","angle [deg]","angular_speed [deg/s]","F_PTO [N]","Energy [J]"]).round(3).to_csv(filename,index=False)
    
    
    #import matplotlib.pyplot as plt
    #plt.figure();
    #plt.plot(sol.t[:],np.transpose([np.sum(wave1.get(sol.t.reshape(sol.t.size,1),0)[0],1),sol.y[0,:]]));
    #plt.show();
    
    # free data from the wave and the WEC
    wave1.clear();
    del wave1;
    wec.release();
    
    # Output the absorbed power
    power=(ydata[8,-1])/(tdata[-1])
    print("Mean absorbed power: "+str((power/1000).round(2))+" kW")
    
    
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
    t2=np.arange(0,t,1/(omega_cut_off*np.pi))
    y=0*t2;
    if isinstance(x0,list):
        return start_simu(time=t2,wave=y,name=n, t0=0, init=x0, control=ctrl )#/(9.81*9.81*1000/(32*np.pi)*A**2*p)    
    return start_simu(time=t2,wave=y,name=n, t0=0, init=[x0, 0, 0], control=ctrl )#/(9.81*9.81*1000/(32*np.pi)*A**2*p)

 
# Regular wave (Height, period)
def reg_wave(H=1,p=10,n0=8,n=8):
    print("Regular wave")
    w=wave_series();
    H=float(H);
    p0=abs(float(p));
    t0=p0*n0;
    t2=t0+p0*n;
    w.t=np.arange(0,t2,1/(omega_cut_off*np.pi))
    w.y=H/2*np.cos(2*np.pi/p0*w.t)
    w.t0=t0;
    w.te=t2-p0;
    if float(p)>0:
        w.t0=0;
    return w;

# Brettschneider wave (significant wave height, energy period, name, control)
def bretschneider_wave(Hs=1,p=6,n0=4,n=6):
    print("Brettschneider wave")
    w=wave_series();
    Hs=float(Hs);
    p0=abs(float(p));
    omega=np.linspace(0.001,4,200);
    omega_m=2*np.pi/(p0/0.856);
    S=5/16*(omega_m**4)/(omega**5)*(Hs**2)*np.exp(-5*(omega_m**4)/4/(omega**4))
    t0=p0*n0;
    t2=t0+p0*n;
    w.t=np.arange(0,t2,1/(omega_cut_off*np.pi))
    np.random.seed(6)#Maybe replace by fixed phase vector; has to guaranteed that random sequecne is always the same
    phase=np.random.rand(omega.size)*2*np.pi;
    w.te=t2-p0;
    if float(p)>0:
        w.t0=0;
    w.y=np.sum(np.sqrt(2*S*(omega[1]-omega[0]))*np.sin(omega*w.t.reshape(w.t.size,1)+phase),1)
    return w;

def get_WEC_data():
    f = open(utils.wec_dir+"/floater.txt",'r')
    data=json.load(f);
    f.close();
    return data; 
    
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
    #start_simu(time=t, wave=np.sin(t/10), name="test", t0=0, control="TCP", host=True);
    #reg_wave(4,3.5,"test.csv","Controller1.py");
    #reg_wave(1,10,"test.csv","controller_reactive.py")
    #decay_test(0.15,"decay1.csv",20,"linear")
    #start_simu(wave=reg_wave(1,3),name="output.csv",control="linear")
    bretschneider_wave(1,3).to_file("testwave1.csv")
    #start_simu(wave=bretschneider_wave(1,3),name="output.csv",control="linear")
    start_simu(file="testwave1.csv",name="output.csv",control="linear")
    #bretschneider_wave(1.5,12,"bretschneider_wave.csv","python3 controller.py")
     
   
