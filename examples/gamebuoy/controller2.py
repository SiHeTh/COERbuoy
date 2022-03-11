#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 18:56:19 2021

@author: heiko
"""

import asyncio
import websockets
import numpy as np;
import COERbuoy;
import threading;
import json;
import tempfile;
import os;

latch=0;
def threadfkt(path):
    os.chdir(path)
    print("start server")
    #COERbuoy.reg_wave(4,12,"test.csv","TCP");
    #COERbuoy.bretschneider_wave(3,8,"bretschneider_wave.csv","TCP")
    COERbuoy.start_simu(wave=COERbuoy.bretschneider_wave(4,6,n=23,n0=0,ne=0),control="TCP")
    #COERbuoy.bretschneider_wave(5,12,"bretschneider_wave.csv","TCP")
    
#async def input_handler(ws, path):
#    global latch;
#    #print("bla1")
#    #ms= await ws.recv();
#    #print(ms)
#    async for msg1 in ws:
#        ms=json.loads(msg1);
#        #print(ms)
#        if ms["action"]=="latch":
#            latch=1;
#            break;
#    print("finished")
               
async def echo(ws, path):
    try:
        from random import randint;
        #tdir=tempfile.TemporaryDirectory(dir=os.path.join(os.path.realpath(__file__),"tmp"));
        tdir=tempfile.TemporaryDirectory();
        conn_ip="localhost";
        conn_port=randint(50000,60000);
        with open(os.path.join(tdir.name,"coerbuoy_settings.txt"),"w+") as f:
            f.write("""{"hydro":"Floater_BEM",\n"WECfolder":"[data.COERsimple]",
"dt_controller": 0.1,
"resolution":0.01,
"conn_ip":"""+'"'+conn_ip+'"'+""",
"conn_port":"""+str(conn_port)+""",
"status_message":[1,1,1,1,1,1,1,1,1]}""");
    
        global fscore;
        #Start controller
        import COERbuoy.connection as connection;
        import time as timem;
        score=0;
        tlast=0;
        t=threading.Thread(target=threadfkt,args=(tdir.name,))
        t.start();
        timem.sleep(4);
        print("start")
        print("Using ip "+conn_ip+" and port " + str(conn_port) +".")
        conn_model=connection.connection(ip = conn_ip, port = conn_port);
        conn_model.openC();#Use client mode
        
        import numpy as np;
        from scipy.interpolate import interp1d;
        print("0")
        timei=-1;
        input_task = None;#asyncio.ensure_future(
        #input_handler(ws, path))
            
            
            
            
            
        #Send buoy geometry
        data=COERbuoy.get_WEC_data();
        xy=[];
        for d in data["geo"]:
            xy=xy+d["coord"];
        cy=xy[:1]+xy[2::4];
        cx=xy[1:2]+xy[3::4];
        msgdata={"type":"geo","cx":cx,"cy":cy};
        await ws.send(json.dumps(msgdata));
        latch=0;
            
        while msg:=conn_model.get_control():
              #async ms= await ws.recv();
               #   print(ms);
              #for msg1 in ms:
              #    if msg1=="latch":
              #        latch=True;
              
              if timei<0:
                  timei=timem.time();
              ##Read incoming message
              time  =msg["time"]          #1x100 array with time series
              wave  =msg["wave"];         #1x100 array with wave data (related to time series)*
              wave_f=msg["wave_forecast"];#1x100 array with wave forecast*
              x     =msg["stroke_pos"];   #1x100 array with stroke position (related to time series)
              dx    =msg["stroke_speed"]; #1x100 array with stroke speed (related to time series)
              alpha =msg["angular_pos"];  #1x100 array with pitch angle (related to time series)
              dalpha=msg["angular_speed"];#1x100 array with pitch angular speed (related to time series)
              force =msg["force"];        #1x100 array with force sensor data (related to time series)
              wavex =msg["test"];        #1x100 array with force sensor data (related to time series)
              #* not available during COERbuoy1 benchmark
              now=time[-1]; #The last element contains the most recent data (except the wave forecast)
              
              #msgdata={"type":"status","wave":(-1*wave_f[14*2:0:-1]+-1*wave[:100-14*2:-1])[::1],"y":(x[-1]+np.cos(alpha[-1])*20)-20,"x":x[-1]*np.sin(alpha[-1])*20,"score":str(int(score/1000))};
              l=60;
              #msgdata={"type":"status","wave":(wave_f[50:0:-1]+wave[:50:-1])[::1],"y":1*(x[-1]+np.cos(alpha[-1])*l)-l,"x":x[-1]*np.sin(alpha[-1])*l,"score":str(int(score/1000))};
              msgdata={"type":"status","wave":(wavex)[::1],"y":1*(x[-1]+np.cos(alpha[-1])*l)-l,"x":x[-1]*np.sin(alpha[-1])*l,"score":str(int(score/1000))};
              await ws.send(json.dumps(msgdata));
              #print("test");
    
              ##This part contains the controller's logic
              ##TODO: replace logic with own controller idea
              
              #(here we use a simple linear damping)
              #Calculate PTO force
              gamma=100000;#kNs/m
              v=(x[-1]-x[-2])/(time[-1]-time[-2]);
              F_pto=-gamma*dx[-1];
              brake=0;#kNs/m
              #global latch;
              
              if latch<1:
                  score=score+0.5*gamma*dx[-1]**2*(now-tlast);
              tlast=now;
     
        
              latch=0;
              #if input_task!=None:
              #    try:
              #        await asyncio.wait_for(input_task,0.01);
              #        print("latch")
              #        latch=1;
              #    except asyncio.TimeoutError:
              #        latch=0;
              msg1=await ws.recv();
              ms=json.loads(msg1);
              if ms["type"]=="brake":
                latch=ms["brake"];
              else:
                latch=ms["brake"];
            
              ##Write control message 
              answer={
                      "time":np.linspace(now,now+1,9),#1x9 array with time steps in the future
                      "pto":np.array([F_pto]*9),      #1x9 array with PTO force
                      "brake":np.ones(9)*latch*2000000,            #1x9 array with brake force
                      "test":np.zeros(9),             #1x9 array; not used
                      }
              #print(now-(timem.time()-timei))
              timem.sleep(np.max([0,now-(timem.time()-timei)]));
              #Send control message to model
              #print(latch)
              conn_model.set_control(answer["time"],answer["pto"],answer["brake"],answer["test"]);
              latch=0;
              #input_task = asyncio.ensure_future(
              #input_handler(ws, path))
        
    finally:        
        conn_model.close();
        tdir.cleanup();
        
    
        score=int(score/1000);
        with open(fscore,'r') as f:
            slist=json.load(f);
            print(slist)
        for idx,s in enumerate(slist["scores"]):
            if s<score:
                print("hurray")
                #hurray, in highscore list
                #reorder table and insert new score
                for i in range(len(slist["scores"])-1,idx,-1):
                    slist["scores"][i]=slist["scores"][i-1];
                    slist["names"][i]=slist["names"][i-1];
                slist["scores"][idx]=score;
                #send the good news
                msgdata={"type":"highscore","names":slist["names"],"scores":slist["scores"],"position":idx};
                await ws.send(json.dumps(msgdata));
                #And wait for name to be graved into the list
                name="noname"
                ms={"type":"none"};
                while ms["type"]!="name":
                    msg1=await ws.recv();
                    ms=json.loads(msg1); 
                    print(ms)   
                name=ms["name"];
                slist["names"][idx]=name;
                with open(fscore,'w') as f:
                    f.write(json.dumps(slist));
                break;
        msgdata={"type":"finish","wave":(wavex)[::1],"heave":x[-1],"score":str(int(score/1000))};
        await ws.send(json.dumps(msgdata));   
            
        #async for msg in ws:
        #    print(msg)
        #    await ws.send(msg);
            
async def main():
    async with websockets.serve (echo, "localhost",8888):
        await asyncio.Future();
        
global fscore;
fscore=os.path.realpath(os.path.join(os.path.dirname(__file__),"highscore.txt"));
asyncio.run(main());
