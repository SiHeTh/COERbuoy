#!/Hellusr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 10:45:57 2020

@author: heiko
"""

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import time
from COERbuoy.simulation import start_simu, reg_wave, bretschneider_wave, decay_test;
import numpy as np;
import os.path
import json;
import threading;
import COERbuoy.Parameters as Parameters;
import webbrowser;
import COERbuoy.utils;
#import svgwrite
#from svgwrite import mm, deg
#import concurrent.futures

folder=COERbuoy.utils.pkg_dir+"/web/";
def about():
    print("WEC Demonstrator, COER")

fktdict={"regular_wave":reg_wave,
         "brettschneider_wave":bretschneider_wave,
         "decay_test":decay_test}
busy=False;

def threadfkt():
    global busy
    idx=0;
    #busy=True;
    while len(jobs)>idx:
        job=jobs[idx];
        if job["status"] != "Done":
            jobs[idx]["status"]="Running";
            print("Started job "+job["name"])
            jobs[idx]["power"]=fktdict[job["fkt"]](*job["args"])
            jobs[idx]["status"]="Done";
            print("Finished job "+job["name"])#+"; Absorbed: "+str(jobs[jobidx]["status"])+" %")
        idx=idx+1;
    busy=False;
    print("Work done!")
jobs=[];


        

##Application is ready for multithreading; before max_worker is set to > 1,
##please ensure the intferface simulation<>controller can handle multiple threads
#pool = concurrent.futures.ThreadPoolExecutor(max_workers=4);
#busy=False;
class GUIServer(BaseHTTPRequestHandler):
    #t=threading.Thread(target=threadfkt,args=())    
    print(threading.active_count())
    def send_html_file(self, file):
        self.send_response(200)
        self.send_header("Content-type","text/html")
        self.end_headers()
        f = open(folder+str(file))
        #if f.read(6)=="navbar":
        #    print("enter navbar")
        #    fhead = open("head.html");
        #    self.wfile.write(bytes(fhead.read(),"utf-8"))
        self.wfile.write(bytes(f.read(),"utf-8"))
        f.close();
    
    def do_GET(self):
        print(self.path)
        global busy;
        p=self.path.split("?");  
        if self.path[-3:]=="svg":
            if os.path.isfile(folder+self.path[1:]):
                print("delivering "+self.path[1:])
                f = open(folder+self.path[1:],'rb')
                self.send_response(200)
                self.send_header("Content-type","image/svg+xml")
                self.end_headers()
                self.wfile.write(f.read())
                f.close();
        if self.path[-3:]=="png":
            if os.path.isfile(folder+self.path[1:]):
                print("delivering "+self.path[1:])
                f = open(folder+self.path[1:],'rb')
                self.send_response(200)
                self.send_header("Content-type","image/png")
                self.end_headers()
                self.wfile.write(f.read())
                f.close();
        if self.path[-3:]=="jpg":
            if os.path.isfile(folder+self.path[1:]):
                print("delivering "+self.path[1:])
                f = open(folder+self.path[1:],'rb')
                self.send_response(200)
                self.send_header("Content-type","image/jpg")
                self.end_headers()
                self.wfile.write(f.read())
                f.close();
        if self.path[-3:]=="css":
            if os.path.isfile(folder+self.path[1:]):
                print("delivering "+self.path[1:])
                f = open(folder+self.path[1:],'r')
                self.send_response(200)
                self.send_header("Content-type","text/css")
                self.end_headers()
                self.wfile.write(bytes(f.read(),"utf-8"))
                f.close();
        if p[0][-3:]=="csv":
            if os.path.isfile(COERbuoy.utils.pkg_dir+"/"+p[0][1:]):
                print("delivering "+p[0][1:])
                f = open(COERbuoy.utils.pkg_dir+"/"+p[0][1:],'r')
                self.send_response(200)
                self.send_header("Content-type","text/csv")
                self.end_headers()
                self.wfile.write(bytes(f.read(),"utf-8"))
                f.close();
        elif self.path=="csv":
            if os.path.isfile(COERbuoy.utils.pkg_dir+"/"+self.path[1:]):
                print("delivering "+self.path[1:])
                f = open(COERbuoy.utils.pkg_dir+"/"+self.path[1:],'r')
                self.send_response(200)
                self.send_header("Content-type","text/csv")
                self.end_headers()
                self.wfile.write(bytes(f.read(),"utf-8"))
                f.close();
                
        elif p[0][-3:]=="pdf":
            if os.path.isfile(self.path[1:]):
                print("delivering "+self.path[1:])
                f = open(self.path[1:],'rb')
                self.send_response(200)
                self.end_headers()
                self.wfile.write(bytes(f.read()))
                f.close();
                
        elif p[0][-3:]=="zip":
            if os.path.isfile(self.path[1:]):
                print("delivering "+self.path[1:])
                f = open(COERbuoy.utils.pkg_dir+"/"+self.path[1:],'rb')
                self.send_response(200)
                self.end_headers()
                self.wfile.write(bytes(f.read()))
                f.close();
                
                
        elif self.path=="/visuals.html":
            self.send_html_file("testchart.html");
            
        elif self.path=="/test.html":
            self.send_html_file("test.html");
            
        elif self.path=="/start.html":
            self.send_html_file("start.html");
            
            
        elif self.path=="/jobs.json":
            print("refreshing job list")
            self.send_response(200)
            self.send_header("Content-type","text/json")
            self.end_headers()

            self.wfile.write(bytes(json.dumps({"status":busy,"jobs":tuple(jobs)}),"utf-8"))
        elif self.path=="/params.json":
            print("delivering parameter list")
            #print(COERbuoy.utils.wec_dir)
            f = open(COERbuoy.utils.wec_dir+"/floater.txt")
            
            self.send_response(200)
            self.send_header("Content-type","text/text")
            self.end_headers()
            self.wfile.write(bytes(f.read(),"utf-8"))
            f.close();

        elif self.path=="/settings.json":
            f = open(COERbuoy.utils.pkg_dir+"/settings.txt",'r')
            self.send_response(200)
            self.send_header("Content-type","text/json")
            self.end_headers()
            self.wfile.write(bytes(f.read(),"utf-8"))
            f.close();

            
        elif self.path[-3:]==".js":
            print("delivering "+self.path[1:])
            f = open(folder+self.path[1:],'r')
            self.send_response(200)
            self.send_header("Content-type","text/javascript")
            self.end_headers()
            self.wfile.write(bytes(f.read(),"utf-8"))
            #f.close();
        elif p[0]=="/results.html":
            self.send_html_file("results.html");
        elif p[0]=="/about.html":
            self.send_html_file("about.html");
        elif p[0]=="/doc.html":
            self.send_html_file("params.html");                
        else:
            self.send_html_file("start.html");
            
    def do_POST(self):
        global busy;
        print(self.path)
        if self.path=="/set_settings.json":
            txt=self.rfile.read(int(self.headers["Content-Length"]))
            print(txt)
            txt=json.loads(txt);
            print(txt)
            with open(COERbuoy.utils.pkg_dir+"/settings.txt",'r') as jfile:
                print(jfile);
                sets=json.load(jfile);
                print(sets);
            for k,v in txt.items():
                sets[k]=v;
            Parameters.run();
            with open(COERbuoy.utils.pkg_dir+"/settings.txt",'w') as jfile:
                jfile.write(json.dumps(sets));
            #f.write(txt);
            #f.close();
            #Parameters.run();
            self.send_response(200);
            #self.send_header('Content-type','text/html');
            #self.send_header('Location','/start.html');
            self.end_headers();
            
        if self.path=="/new_param.json":
            txt=self.rfile.read(int(self.headers["Content-Length"]))
            print("Write params to "+COERbuoy.utils.wec_dir+"/floater.txt.")
            f = open(COERbuoy.utils.wec_dir+"/floater.txt",'wb+')
            f.write(txt);
            f.close();
            Parameters.run();
            self.send_response(200);
            #self.send_header('Content-type','text/html');
            #self.send_header('Location','/start.html');
            self.end_headers();
            
        elif self.path=="/run.html" :   
            json_txt=self.rfile.read(int(self.headers["Content-Length"]))
            print(json_txt)
            data1=json.loads(json_txt);
            self.send_html_file("results.html");
            
            if not busy:# not busy:
                print(busy)
                jobs.clear();
                #For security reasons only allow known command to be executed
                if not (data1["ctrl"]=="python3 Controller_NL_TCP.py" 
                        or data1["ctrl"]=="python3 Controller_LIN_TCP.py"
                        or data1["ctrl"]=="python3 LatchingForce_TCP.py"
                        or data1["ctrl"]=="linear"
                        or data1["ctrl"]=="python3 Conroller_brake_1.py"
                        or data1["ctrl"]=="ocatve Extremiumcontroller2.m"):
                    data1["ctrl"]="linear";
                for data in data1["sea_states"]:
                            print(data)
                            strdata="";
                            for e in data.items():
                                strdata=strdata+e[0]+str(e[1])+"_";
                            strdata=strdata[:-1];
                            name="results/"+strdata+data1["ctrl"].replace(' ','_')+".csv";
                            desc=strdata+". Using "+data1["ctrl"]+".";
                            print("Created job for "+desc+".")
                            if data["wave"]!="decay_test":
                                job={"fkt":data["wave"], "power":0, "status":"Not running","desc":desc, "args":(float(data["H"]), float(data["P"]), COERbuoy.utils.pkg_dir+"/"+name, data1["ctrl"]), "name":name, "id":len(jobs)};
                            else:
                                job={"fkt":data["wave"], "power":0, "status":"Not running","desc":desc, "args":(float(data["x0"]), name, float(data["t"]), data1["ctrl"]), "name":name, "id":len(jobs)};
                            
                            jobs.append(job);
                busy=True;
                self.send_response(200);
                self.end_headers();
            
                self.t=threading.Thread(target=threadfkt,args=())
                self.t.start();
          
            
                
    
def run():
    hostName="localhost"
    serverPort = 8080;
   
    with open(COERbuoy.utils.pkg_dir+"/settings.txt") as file:
        data=json.load(file);
        hostName=data.get("host","localhost");
        serverPort=data.get("port",8080);
    print("Bind to "+hostName+" at port "+str(serverPort)+".")
    #HTTPServer.socket.setsockopt(socket.SOL_SOCKeT, socket.SO_REUSEADDR,1);
    webServer = ThreadingHTTPServer((hostName, serverPort), GUIServer)
    webbrowser.open("http://localhost:8080")
    try:
        webServer.serve_forever()
    except KeyboardInterrupt:
        print("Interrupt");
        #pool.shutdown();
        pass
    
    webServer.server_close()
    
if __name__ == "__main__":
    run();