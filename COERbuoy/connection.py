#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# connection.py - Implementation of the WEC control interface
# 2020/2021 COER Laboratory, Maynooth University
# in cooperation with CorPower Ocean AB

"""
Created on Mon Oct 12 16:05:19 2020

@author: Simon H. Thomas
"""

#import modules
import socket
import numpy as np;
import threading;
import time;
import struct

#set default address
TCP_IP = 'localhost'
TCP_PORT = 5050
BUFFER_SIZE = 1024*8;
BUFFER_SIZE_IN=9*4*8;

#model message
msg_model={"time":np.zeros(100),
           "wave":np.zeros(100),
           "wave_forecast":np.zeros(100),
           "stroke_pos":np.zeros(100),
           "stroke_speed":np.zeros(100),
           "angular_pos":np.zeros(100),
           "angular_speed":np.zeros(100),
           "force":np.zeros(100),
           "test":np.zeros(100)}

# control message
msg_ctrl={"time":np.zeros(9),
           "pto":np.zeros(9),
           "brake":np.zeros(9),
           "test":np.zeros(9)}
    
#class connection
class connection():
    socket=[];
    conn=[];
    
    #connect socket as client
    def openC(self):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.socket.settimeout(10);
        self.socket.connect((TCP_IP, TCP_PORT))
    def close(self):
        print("close socket")
        if self.conn:
            self.conn.close();
        self.socket.close()
        
    #to be called by the model: Send model message and return control answer
    def exchange_model(self,time,wave,wave_forecast,stroke_pos,stroke_speed,angular_pos,angular_speed,force,test):
        #Create model messaga array
        msg=time.tolist()+wave.tolist()+wave_forecast.tolist()+stroke_pos.tolist()+stroke_speed.tolist()+angular_pos.tolist()+angular_speed.tolist()+force.tolist()+test.tolist();
        msg=struct.pack(">{}d".format(len(msg)),*msg)
        
        #send (as host or client)
        if self.conn:
            self.conn.send(msg)
            #...and receive control message
            msg2=self.conn.recv(BUFFER_SIZE_IN,socket.MSG_WAITALL);  
        else:
            self.socket.send(msg)
            #...and receive control message
            msg2=self.socket.recv(BUFFER_SIZE_IN,socket.MSG_WAITALL);
        if not msg2:
            print("failed receiving controller data")
            self.close()
            return;
        msg2=struct.unpack(">{}d".format(int(len(msg2)/8)),msg2)
        i=0;
        msg_ctrl["time"]=np.array(msg2[i*9:i*9+9])#/1000
        i=1;
        msg_ctrl["pto"]=np.array(msg2[i*9:i*9+9])#*1000000
        i=2;
        msg_ctrl["brake"]=np.array(msg2[i*9:i*9+9])#/1000
        i=3;
        msg_ctrl["test"]=np.array(msg2[i*9:i*9+9])#/1000
        return msg_ctrl;
    
    #to be called by the model: get model data
    def get_control(self):
        msg=[];
        #time.sleep(3)
        if self.conn:
            msg=self.conn.recv(BUFFER_SIZE);
        else:
            msg=self.socket.recv(BUFFER_SIZE);
        if not msg:
            print("failed receiving model data")
            self.close()
            return;
        msg=struct.unpack(">{}d".format(int(len(msg)/8)),msg);
        msg=list(msg)
        i=0;
        msg_model["time"]=msg[i*100:i*100+100]
        i=1;
        msg_model["wave"]=msg[i*100:i*100+100]
        i=2;
        msg_model["wave_forecast"]=msg[i*100:i*100+100]
        i=3;
        msg_model["stroke_pos"]=msg[i*100:i*100+100]
        i=4;
        msg_model["stroke_speed"]=msg[i*100:i*100+100]
        i=5;
        msg_model["angular_pos"]=msg[i*100:i*100+100]
        i=6;
        msg_model["angular_speed"]=msg[i*100:i*100+100]
        i=7;
        msg_model["force"]=msg[i*100:i*100+100]
        i=8;
        msg_model["test"]=msg[i*100:i*100+100]
        #print("Model Message: "+str(msg_model))
        return msg_model;
    
    #to be called by the controller: Set (and send) the model control answer
    def set_control(self,time,pto,brake,test):
        msg=time.tolist()+pto.tolist()+brake.tolist()+test.tolist();
        msg=struct.pack(">{}d".format(len(msg)),*msg)
        if self.conn:
            self.conn.send(msg)
        else:
            self.socket.send(msg); 
    #open connection as host
    def openH(self):
        self.socket = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.socket.settimeout(10);
        self.socket.bind((TCP_IP, TCP_PORT))
        self.socket.listen()
        self.conn, addr = self.socket.accept()
        
#------------------------------------------------------------------------------
#Only for testing
def client():
    c=connection();
    c.openC();
    print("opened client")
    t=np.linspace(-1,0,50);
    for idx in range(0,10):
        print(c.exchange_model(t,np.sin(t),-np.sin(t),np.cos(t),np.sin(t)*0.1,np.cos(t)*0.1,t*2,t*3,t*0))
    c.close();
    time.sleep(2);
def host():
    h=connection();
    print("opened connection")
    h.openH();
    while True:
        t=np.linspace(-1,0,50);
        v=np.linspace(0,49,50);
    
        print(h.get_control())
        t=np.linspace(-0.5,0.5,9);
        h.set_control(t,np.sin(t),np.cos(t),t*0);
    h.close();
    time.sleep(2);
    
#test the connection class
if __name__ == "__main__":
    t1=threading.Thread(target=host,args=())
    t2=threading.Thread(target=client,args=())
    t1.start();
    time.sleep(0.5);
    t2.start();
    
    while (t1.is_alive()):# or t2.is_alive()):
        time.sleep(3);
    