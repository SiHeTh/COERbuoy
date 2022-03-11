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
import COERbuoy.analyzer;
from shutil import copyfile;
#import svgwrite
#from svgwrite import mm, deg
#import concurrent.futures

def about():
    print("WEC Demonstrator, COER")

fktdict={"regular_wave":reg_wave,
         "bretschneider_wave":bretschneider_wave,
         "decay_test":decay_test}
busy=False;
idx=0;


def threadfkt():
    global busy
    global idx;
    #busy=True;
    while len(jobs)>idx:
        print("Running job no. "+str(idx))
        job=jobs[idx];
        job["status"] = "Running";
        if job["status"] != "Done":
            print("Started job "+job["name"])
            filename=os.path.join(COERbuoy.utils.pkg_dir,"results");
            if job["fkt"]== "regular_wave":
                [a,b,c]=start_simu(wave=reg_wave(job["H"],job["P"]),name=filename,control=job["ctrl"]);
            elif job["fkt"]== "bretschneider_wave":
                [a,b,c]=start_simu(wave=bretschneider_wave(job["H"],job["P"]),name=filename,control=job["ctrl"]);
            elif job["fkt"]=="decay_test":
                [a,b,c]=decay_test(job["x0"],job["name"],job["t"],job["ctrl"]);
            job["name"]=os.path.join("results",os.path.basename(b));
            print(b)
            job["power"]=c;
            jobs[idx]["status"]="Done";
            print("Finished job "+job["name"])#+"; Absorbed: "+str(jobs[jobidx]["status"])+" %")
            try:
                copyfile(os.path.join(COERbuoy.utils.pkg_dir,job["name"]),os.path.join(COERbuoy.utils.user_dir,job["name"]))
            except(FileNotFoundError):
                print("COERbuoy_data/results does not exist!")
        idx=idx+1;
    busy=False;
    print("Work done!")
jobs=[];


        
class GUIServer(BaseHTTPRequestHandler):
    print(threading.active_count())
    def send_html_file(self, file):
        self.send_response(200)
        self.send_header("Content-type","text/html")
        self.end_headers()
        f = open(folder+str(file))
        self.wfile.write(bytes(f.read(),"utf-8"))
        f.close();
    
    def do_GET(self):
        #print(self.path)
        global busy;
        p=self.path.split("?")[0];  
        if p[-3:]=="svg":
            if os.path.isfile(folder+p[1:]):
                print("delivering "+p[1:])
                f = open(folder+p[1:],'rb')
                self.send_response(200)
                self.send_header("Content-type","image/svg+xml")
                self.end_headers()
                self.wfile.write(f.read())
                f.close();
        if p[-3:]=="png":
            if os.path.isfile(folder+p[1:]):
                print("delivering "+p[1:])
                f = open(folder+p[1:],'rb')
                self.send_response(200)
                self.send_header("Content-type","image/png")
                self.end_headers()
                self.wfile.write(f.read())
                f.close();
        if p[-3:]=="jpg":
            if os.path.isfile(folder+p[1:]):
                print("delivering "+p[1:])
                f = open(folder+p[1:],'rb')
                self.send_response(200)
                self.send_header("Content-type","image/jpg")
                self.end_headers()
                self.wfile.write(f.read())
                f.close();
        if p[-3:]=="css":
            if os.path.isfile(folder+p[1:]):
                print("delivering "+p[1:])
                f = open(folder+p[1:],'r')
                self.send_response(200)
                self.send_header("Content-type","text/css")
                self.end_headers()
                self.wfile.write(bytes(f.read(),"utf-8"))
                f.close();
        if p[-3:]=="csv":
            file=None;
            print("requested: "+p)
            if os.path.isfile(os.path.join(COERbuoy.utils.pkg_dir,p[1:])):#local dir
                file=os.path.join(COERbuoy.utils.pkg_dir,p[1:]);
            elif os.path.isfile(p) and not secureGUI:# absolute dir
                file=p;
            elif os.path.isfile(os.path.join(os.path.expanduser("~"),"COERbuoy_data/results",p.split("/")[-1])):#local dir
                file=os.path.join(os.path.expanduser("~"),"COERbuoy_data/results",p.split("/")[-1]);
            if file:
                f = open(file,'r');
                print("delivering "+p);
                self.send_response(200)
                self.send_header("Content-type","text/csv")
                self.end_headers()
                self.wfile.write(bytes(f.read(),"utf-8"))
                f.close();
        elif p[-3:]=="csv":
            print("requested: "+[COERbuoy.utils.pkg_dir,p[1:],"+s"])
            print([COERbuoy.utils.pkg_dir,p[1:],"+s"])
            file=None;
            if os.path.isfile(os.path.join(COERbuoy.utils.pkg_dir,p[1:])):
                print("delivering "+p[1:]);
                file=os.path.join(COERbuoy.utils.pkg_dir,p[1:]);
            elif os.path.isfile(p[1:]):
                print("delivering "+p[1:]);
                file=p[1:];
            if file:
                f = open(file,'r');
                self.send_response(200)
                self.send_header("Content-type","text/csv")
                self.end_headers()
                self.wfile.write(bytes(f.read(),"utf-8"))
                f.close();
                
        elif p[-3:]=="pdf":
            if os.path.isfile(p[1:]):
                print("delivering "+p[1:])
                f = open(p[1:],'rb')
                self.send_response(200)
                self.end_headers()
                self.wfile.write(bytes(f.read()))
                f.close();
                
        elif p[-3:]=="zip":
            if os.path.isfile(p[1:]):
                print("delivering "+p[1:])
                f = open(COERbuoy.utils.pkg_dir+"/"+p[1:],'rb')
                self.send_response(200)
                self.end_headers()
                self.wfile.write(bytes(f.read()))
                f.close();
                
                
            
        elif p=="/jobs.json":
            self.send_response(200)
            self.send_header("Content-type","text/json")
            self.end_headers()
            self.wfile.write(bytes(json.dumps({"status":busy,"jobs":tuple(jobs)}),"utf-8"))
            
        elif p=="/files.json":
            self.send_response(200)
            self.send_header("Content-type","text/json")
            self.end_headers()
            a=COERbuoy.analyzer.Analyzer();
            files=[];
            path=os.path.join(COERbuoy.utils.user_dir,"results");
            #a.read_folder(path);
            #abcd=a.table.to_json(orient="records");
            #self.wfile.write(bytes(a.table.to_json(orient="records"),"utf-8"))
            for f in os.listdir(path):#look for files in a folder
                if f[-4:] == ".csv":
                    files.append({"file":os.path.join(path,f),"name":f[:-4]});
            self.wfile.write(bytes(json.dumps({"list":tuple(files)}),"utf-8"))
            
        elif p=="/params.json":
            print("delivering parameter list")
            #print(COERbuoy.utils.wec_dir)
            f = open(os.path.join(COERbuoy.utils.wec_dir,"floater.txt"))
            
            self.send_response(200)
            self.send_header("Content-type","text/text")
            self.end_headers()
            self.wfile.write(bytes(f.read(),"utf-8"))
            f.close();
        elif p=="/controllers.json":
            print("delivering controller list")
            self.send_response(200)
            self.send_header("Content-type","text/json")
            self.end_headers()
            self.wfile.write(bytes(json.dumps(COERbuoy.utils.get_controller()),"utf-8"));
            
        elif p=="/settings.json":
            f = open(COERbuoy.utils.pkg_dir+"/settings.txt",'r')
            self.send_response(200)
            self.send_header("Content-type","text/json")
            self.end_headers()
            self.wfile.write(bytes(f.read(),"utf-8"))
            f.close();

            
        elif p[-3:]==".js":
            print("delivering "+p[1:])
            f = open(folder+p[1:],'r')
            self.send_response(200)
            self.send_header("Content-type","text/javascript")
            self.end_headers()
            self.wfile.write(bytes(f.read(),"utf-8"))
        elif p[-5:]==".html":    
            self.send_html_file(p[1:]);           
        else:
            self.send_html_file("start.html");
            
    def do_POST(self):
        global busy;
        print(self.path)
        
        #set settings with key value pairs
        if self.path=="/set_settings.json":
            txt=self.rfile.read(int(self.headers["Content-Length"]))
            txt=json.loads(txt);
            with open(COERbuoy.utils.pkg_dir+"/settings.txt",'r') as jfile:
                print(jfile);
                sets=json.load(jfile);
            for k,v in txt.items():
                sets[k]=v;
            Parameters.run();
            with open(COERbuoy.utils.pkg_dir+"/settings.txt",'w') as jfile:
                jfile.write(json.dumps(sets));
            self.send_response(200);
            self.end_headers();
            COERbuoy.utils.get();
            Parameters.run();
        
        #replace settings file
        if self.path=="/new_settings.json" and not secureGUI:
            txt=self.rfile.read(int(self.headers["Content-Length"]))
            with open(COERbuoy.utils.pkg_dir+"/settings.txt",'w') as jfile:
                jfile.write(txt);
            self.send_response(200);
            self.end_headers();
            COERbuoy.utils.get();
            
        #replace WEC parameters with new values
        if self.path=="/new_param.json":
            txt=self.rfile.read(int(self.headers["Content-Length"]))
            print("Write params to "+COERbuoy.utils.wec_dir+"/floater.txt.")
            f = open(COERbuoy.utils.wec_dir+"/floater.txt",'wb+')
            f.write(txt);
            f.close();
            Parameters.run();
            self.send_response(200);
            self.end_headers();
            
        elif self.path=="/run.html" :   
            json_txt=self.rfile.read(int(self.headers["Content-Length"]))
            print(json_txt)
            data1=json.loads(json_txt);
            self.send_html_file("results.html");
            
            if not busy:# not busy:
                print(busy)
                #jobs.clear();
                #For security reasons only allow known commands to be executed;
                #can be removed in local installation
                for data in data1["sea_states"]:
                            if secureGUI:
                                if not (data1["ctrl"]=="controller_damping.py" 
                                        or data1["ctrl"]=="controller_reactive.py"
                                        or data1["ctrl"]=="MomentBasedController.py"
                                        or data1["ctrl"]=="linear"
                                        or data1["ctrl"]=="none"):
                                    data1["ctrl"]="linear";
                            print(data)
                            strdata="";
                            for e in data.items():
                                strdata=strdata+e[0]+str(e[1])+"_";
                            strdata=strdata[:-1];
                            name=os.path.join("results",strdata+data1["ctrl"].replace(' ','_')+".csv");
                            desc=strdata+". Using "+data1["ctrl"]+".";
                            print("Created job for "+desc+".")
                            if data["wave"]!="decay_test":
                                job={"fkt":data["wave"], "power":0, "status":"Not running","desc":desc, "H":float(data["H"]), "P":float(data["P"]), "name":os.path.join(COERbuoy.utils.pkg_dir,name), "ctrl":data1["ctrl"], "name":name, "id":len(jobs)};
                            else:
                                job={"fkt":data["wave"], "power":0, "status":"Not running","desc":desc, "x0":float(data["x0"]), "name":os.path.join(COERbuoy.utils.pkg_dir,name), "t":float(data["t"]), "ctrl":data1["ctrl"], "name":name, "id":len(jobs)};
                            
                            jobs.append(job);
                busy=True;
                self.send_response(200);
                self.end_headers();
            
                self.t=threading.Thread(target=threadfkt,args=())
                self.t.start();
          
            
                
    
def run():
    COERbuoy.utils.get();
    global secureGUI
    secureGUI=COERbuoy.utils.secureGUI;
    hostName="localhost"
    serverPort = 8080;
    global folder;
    folder=COERbuoy.utils.pkg_dir+"/web/";
    
    resdir=os.path.join(COERbuoy.utils.pkg_dir,"results");
    reslist=os.listdir(resdir);
    
    for e in reslist:
        if e[0]!=".":
            os.remove(os.path.join(resdir,e));
   
    with open(COERbuoy.utils.pkg_dir+"/settings.txt") as file:
        data=json.load(file);
        hostName=data.get("host","localhost");
        serverPort=data.get("port",8080);
    print("Bind to "+hostName+" at port "+str(serverPort)+".")
    webServer = ThreadingHTTPServer((hostName, serverPort), GUIServer)
    webbrowser.open("http://"+str(hostName)+":"+str(serverPort))
    try:
        webServer.serve_forever()
    except KeyboardInterrupt:
        print("Interrupt");
        #pool.shutdown();
        pass
    
    webServer.server_close()
    
if __name__ == "__main__":
    run();