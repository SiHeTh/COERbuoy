#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# COERbuoy - Data Analyzer
# A small library to simlify the analyses of dta prodiced with the COERbuo platform;
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
import matplotlib.pyplot as plt;


def quartil (a, p):
    a=np.abs(a);
    b=np.sort(a);
    aa=np.sum(a);
    if (aa==0):
        return 0;
    return b[np.argmax(np.cumsum(b)/aa>p)];
    
class Analyzer():
    table=None;
    def read_folder(self,folder):#gat a datatable from a folder
        for f in os.listdir(folder):#look for files in a folder
            if f[-4:] == ".csv":#check if csv
                info=f[:-4].split("_");
                if len(info) == 8 and info[1]=="p" and info[3]=="h":#Check if it follows the naming convention
                    try:
                        data0=pandas.read_csv(f).values.transpose();
                        s1=np.argmax(data0[:1,:]>=0)
                        data=data0[:,s1:];
                        
                        row={"P":[(data[-1,-1]-data[-1,0])/(data[0,-1])],
                        "z075":[quartil(data[1,:],0.75)],
                        "dz075":[quartil(data[2,:],0.75)],
                        "alpha075":[quartil(data[3,:],0.75)],
                        "dalpha075":[quartil(data[4,:],0.75)],
                        "RAO":[quartil(data[1,:],0.75)/float(info[4])],
                        "p":[float(info[2])],
                        "h":[float(info[4])],
                        "wave":[info[0]],
                        "ctrl":[info[5]],
                        "WEC":[info[6]],
                        "model":[info[7]]};
                             
                        if not isinstance(self.table,pandas.DataFrame):
                            self.table=pandas.DataFrame(row);
                        else:
                            self.table=self.table.append(pandas.DataFrame(row));
                            #print(self.table);
                    
                    except:
                        print("File "+str(f)+" is malformatted.")
    def get_set(self,p=None,H=None,ctrl=None,WEC=None,model=None):
        t0=self.table;
        for s,e in zip(["p","h","ctrl","WEC","model"],[p,H,ctrl,WEC,model]):
            if e!=None:
                t0=t0[t0[s]==e];
        print(t0);
        return t0;
    def plot(self,xaxis,yaxis,group,t0):
        plt.figure();
        c=['r','g','b','y','k']*5;
        gps=t0[group].unique().tolist()
        for i,e in t0.iterrows():
            i=gps.index(e[group]);
            plt.plot(e[xaxis],e[yaxis],c[i]+'o');
            print(e[xaxis],e[yaxis])
        plt.legend(gps)
        plt.show()
        

if __name__=="__main__":
    a=Analyzer();
    a.read_folder(None);
    b=a.get_set(None,-0.75,None,"COERbuoy1","FloaterBEM");
    a.plot("p","P","ctrl",b)