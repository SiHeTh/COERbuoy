#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# COERbuoy - Data Analyzer
# A small library to simplify the analyses of dta prodiced with the COERbuo platform;
# It is assued that the standard naming sheme is used

# 2021 COER Laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie
#

import os;
import pandas;
import numpy as np;
import warnings;
try:
    import matplotlib.pyplot as plt;
except:
    plt=None;

def quartil (a, p):
    a=np.abs(a);
    b=np.sort(a);
    aa=np.sum(a);
    if (aa==0):
        return 0;
    return b[np.argmax(np.cumsum(b)/aa>p)];
    
class Analyzer():
    d={"p":"wave period [s]","h":"wave height [m]","k_rel":"ka",
       "wave":"wave_type","RAO":"RAO","ctrl":"control method",
       "P":"Mean absorbed power [W]"};
    table=None;
    def read_file(self,file):
        data0=pandas.read_csv(file).values.transpose();
        s1=np.argmax(data0[:1,:]>=0)
        return [data0[:,:s1],data0[:,s1:]];#return transient data - setady state data
                        
    def read_folder(self,folder,a=1):#get a datatable from a folder
        self.table=None;
        for f in os.listdir(folder):#look for files in a folder
            if f[-4:] == ".csv":#check if csv
                info=f[:-4].split("_");
                if len(info) == 8 and info[1]=="p" and info[3]=="h":#Check if it follows the naming convention
                    try:
                        data = self.read_file(os.path.join(folder,f))[1];
                        
                        per=float(info[2]);
                        omega=6.28/per;
                        k=omega**2/9.18;
                        k_rel=k*a;
                        
                        row={"P":[(data[-1,-1]-data[-1,0])/(data[0,-1])],
                        "z075":[quartil(data[2,:],0.75)],
                        "z095":[np.max(data[2,:])-np.min(data[2,:])],
                        "dz075":[quartil(data[3,:],0.75)],
                        "alpha075":[quartil(data[4,:],0.75)],
                        "dalpha075":[quartil(data[5,:],0.75)],
                        "RAO":[quartil(data[2,:],0.95)/float(info[4])],
                        "wtype":[info[0]],
                        "p":[per],
                        "h":[float(info[4])*2],
                        "wave":[info[0]],
                        "ctrl":[info[5]],
                        "WEC":[info[6]],
                        "omega":[omega],
                        "k":[k],
                        "k_rel":[k_rel],
                        "model":[info[7]],
                        "data":[data],
                        "name":[f],
                        "file":[os.path.join(folder,f)]};
                             
                        if abs(data[2,-1])>1e2:#Stuff got instable
                            raise  	ValueError;
                        if not isinstance(self.table,pandas.DataFrame):
                            self.table=pandas.DataFrame(row);
                        else:
                            self.table=self.table.append(pandas.DataFrame(row));
                    except:
                        warnings.warn("File "+str(f)+" is malformatted.");
    def get_set(self,wtype=None,p=None,H=None,ctrl=None,WEC=None,model=None):
        t0=self.table;
        for s,e in zip(["p","h","ctrl","WEC","model"],[p,H,ctrl,WEC,model]):
            if e is not None:
                if isinstance(e,float):
                    t0=t0[(t0[s]>e-0.4) & (t0[s]<e+0.4)]; # bitwise & is equal to logical-and for booleans
                elif isinstance(e,list):
                    t0=t0[t0[s].isin(e)];
                else:
                    t0=t0[t0[s]==e];
        return t0;
    def plot(self,xaxis,yaxis,group,t0,**kwargs):
        if plt is None:
            print("Matplotlib package required to plot graphs!");
            raise ModuleNotFoundError;
            return 0;
        plt.figure();
        scaley=1;
        if "colors" in kwargs:
            c=kwargs["colors"];
        else:
            c=['r','g','b','y','k']*5;
        if "scaley" in kwargs:
            scaley=kwargs["scaley"];
        t0=t0.sort_values(by=[xaxis,group]);
        gps=t0[group].unique().tolist();
        for i,e in t0.iterrows():
            i=gps.index(e[group]);
            plt.plot(e[xaxis],e[yaxis]*scaley,c[i]+'o');
        plt.xlabel(self.d.get(xaxis,xaxis));
        plt.ylabel(self.d.get(yaxis,yaxis));
        plt.title(kwargs.get("title",""));
        plt.legend(gps)
        plt.show()
    def plot_time(self,i,name,t0,**kwargs):
        t0=t0.sort_values(by=["model"]);
        sy=None;
        scaley=1;
        if "axes" in kwargs:
            axes=kwargs["axes"]
        else:
            fig,axes=plt.subplots(1,1)
        if "colors" in kwargs:
            c=kwargs["colors"];
        else:
            c=['r','g','b','y','k']*5;
        if "scaley" in kwargs:
            scaley=kwargs["scaley"];
        if "symbols" in kwargs:
            sy=kwargs["symbols"];
        
        k=0;
        for j,t in t0.iterrows():
            s2=len(t.data[1,:])-1;
            s1=0;
            if "time" in kwargs:
                s1=np.argmax(t.data[:1,:]>=kwargs["time"][0])
                s2=np.argmax(t.data[:1,:]>=kwargs["time"][1])
                if (s2<=s1):
                    s2=len(t.data[:1,:]);
                
                
            if k==0:
                axes.plot(t.data[0,s1:s2],t.data[1,s1:s2]*scaley,"k--");
            value0=0;
            if i>6:
                value0=t.data[i,s1]*scaley;
            if sy is None:
                axes.plot(t.data[0,s1:s2],t.data[i,s1:s2]*scaley-value0,c[k],linewidth=(len(t0)-k)*1.25+0.5);
            else:
                nk=int((s2-s1)/4);
                st=int(nk/len(t0));
                axes.plot(t.data[0,s1:s2],t.data[i,s1:s2]*scaley-value0,c[k],linewidth=1);
                axes.plot(t.data[0,s1+k*st:s2:nk],t.data[i,s1+k*st:s2:nk]*scaley-value0,c[k]+sy[k]);
            k=k+1;
        axes.legend(["Wave"]+t0[name].tolist())
        axes.set_xlabel("time [s]");
        axes.set_title(kwargs.get("title",""));
        #if not "axes" in kwargs:
           # axes.show()
        

if __name__=="__main__":
    a=Analyzer();
    a.read_folder("/home/gast/Dokumente/Validation/data");
    
    abcd=a.table.to_json(orient="records");
    #b=a.get_set(None,None,0.5,"none","COERsimple",None);
    #a.plot("p","RAO","model",b,title="Optimal damping")
    #b=a.get_set("regular",None,0.5,"controllerreactivepy","COERsimple",None);
    #a.plot("p","P","model",b,title="Reactive control")
    #b=a.get_set("regular",None,3.1,"controllerreactivepy","COERsimple",None);
    #a.plot("p","P","model",b,title="Reactive control")
    #b=a.get_set(None,None,0.5,"controllerdampingpy","COERsimple",None);
    #a.plot("p","P","model",b,title="Optimal damping")
    #b=a.get_set(None,None,3.1,"controllerdampingpy","COERsimple",None);
    #a.plot("p","P","model",b,title="Optimal damping")
    
    b=a.get_set("regular",4,0.5,"controllerdampingpy","COERsimple",None);
    a.plot_time(2,"model",b,time=[0,22],title="Stroke [m] for 0.5m, 4s, opt. damping", symbols=["-","*","h"])