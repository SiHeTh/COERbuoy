#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 08:55:59 2020

@author: heiko
"""
# An analytic method to calculate hydrodynamic forces for point-absorber wave energy converter
# Base class for hydrodynamic modules for the COERbuoy wave energy platform
# 2020/2021/2022 COER laboratory, Maynooth University
# in cooperation with CorPower Ocean AB
#
# Author:
# Simon H. Thomas, simon.thomas.2021@mumail.ie
#

import numpy as np;
import json;
from scipy.signal import hilbert as KramerKronig;
from scipy.optimize import fsolve;

pi=np.pi;

#calculates the Froude-Krylov/buoyancy forces for a cone segment
class Cone:
  m = -1;
  q = -1;
  z1 = -1;
  z2 = 0;
  I = 0;
  g=9.81;
  zCoG = -1;
  depth=30;

  def __init__ (self, z1, r1, z2, r2, zCoG, xi, g, h):
      #Frequency parameters
      self.xi=xi;

      #geometric parameters
      self.m=(r1-r2)/(z1-z2);
      if self.m==0:
          self.m=0.1;#avoid divisin by zero
      self.z1=z1;
      self.z2=z2;
      self.g=g;
      self.rho=1000;
      self.q=-r1+self.m*self.z1;
      self.omega=np.sqrt(self.xi*g);

      self.zCoG=zCoG;
      self.z0=0;
      self.x0=0;
      self.delta0=0;
      self.depth=h;

  def Ap(self, z, xi0, xis, eta, q):
      return np.exp(xis*(z-eta))/np.power(xis,4)*(3*xis*self.m**3*q*(xis**2*z**2-2*xis*z+2)+xis*self.m*q*(xis**2*(2*z**2+q**2)-4*xis*z+4)+xis**2*q**2*(xis*z-1)+self.m**4*(xis**3*z**3-3*xis**2*z**2+6*xis*z-6)+self.m**2*(xis**3*(z**3+3*z*q**2)-3*xis**2*(z**2+q**2)+6*xis*z-6));

  def Bp(self, z, xi0, xis, eta, q):
      return np.exp(xis*(z-eta))/np.power(xis,6)*(xis**4*q**4*(xis*z-1)+xis**2*self.m**2*q**2*(xis**3*(6*z**3+5*z*q**2)-xis**2*(18*z**2+5*q**2)+36*xis*z-36)+xis**3*self.m*q**3*(xis**2*(4*z**2+q**2)-8*xis*z+8)+5*xis*self.m**5*q*(xis**4*z**4-4*xis**3*z**3+12*xis**2*z**2-24*xis*z+24)+2*xis*self.m**3*q*(xis**4*(2*z**4+5*z**2*q**2)-2*xis**3*(4*z**3+5*z*q**2)+2*xis**2*(12*z**2+5*q**2)-48*xis+z+48)+self.m**4*(xis**5*(z**5+10*z**3*q**2)-5*xis**4*(z**4+6*z**2*q**2)+20*xis**3*(z**3+3*z*q**2)-60*xis**2*(z**2+q**2)+120*xis*z-120)+self.m**6*(xis**5*z**5-5*xis**4*z**4+20*xis**3*z**3-60*xis**2*z**2+120*xis*z-120))

  def As(self, z, xi0, xis, eta, q):
      return np.exp(xis*(z-eta))/np.power(xis,3)*(np.power(self.m,2)*(np.power(xis,2)*np.power(z,2)-2*xi0*z+2)+2*self.m*q*xis*(xis*z-1)+np.power(q,2)*np.power(xis,2));

  def Bs(self, z, xi0, xis, eta, q):
      return np.power(xi0,2)*np.exp(xis*(z-eta))/(8*np.power(xis,5))*(np.power(self.m,4)*(np.power(xis,4)*np.power(z,4)-4*np.power(xis,3)*np.power(z,3)+12*np.power(xis,2)*np.power(z,2)-24*xis*z+24)+4*np.power(self.m,3)*q*xis*(np.power(xis,3)*np.power(z,3)-3*np.power(xis,2)*np.power(z,2)+6*xis*z-6)+6*np.power(self.m,2)*np.power(q,2)*np.power(xis,2)*(np.power(xis,2)*np.power(z,2)-2*xis*z+2)+4*self.m*np.power(q,3)*np.power(xis,3)*(xis*z-1)+np.power(q,4)*np.power(xis,4))

  def Ah(self, z, xi0, xis, eta, q):
      return np.exp(xis*(z-eta))/np.power(xis,2)*(self.m*(xis*z-1)+xis*q)

  def Bh(self, z, xi0, xis, eta, q):
      return np.exp(xis*(z-eta))/np.power(xis,4)*(np.power(xis,3)*np.power(q,3)+3*xis*self.m*self.m*q*(np.power(xis,2)*np.power(z,2)-2*xis*z+2)+3*np.power(xis,2)*self.m*q*q*(xis*z-1)+np.power(self.m,3)*(np.power(xis,3)*np.power(z,3)-3*np.power(xis,2)*np.power(z,2)+6*xis*z-6))

  def Ch(self, z, xi0, xis, eta, q):
      return np.exp(xis*(z-eta))/np.power(xis,6)*(np.power(xis,5)*np.power(q,5)+5*np.power(xis,4)*self.m*np.power(q,4)*(xis*z-1)+10*np.power(xis,2)*self.m*self.m*self.m*q*q*(np.power(xis,3)*np.power(z,3)-3*xis*xis*z*z+6*xis*z-6)+10*np.power(xis,3)*self.m*self.m*np.power(q,3)*(xis*xis*z*z-2*xis*z+2)+5*xis*np.power(self.m,4)*q*(np.power(xis,4)*np.power(z,4)-4*np.power(xis,3)*np.power(z,3)+12*np.power(xis,2)*np.power(z,2)-24*xis*z+24)+np.power(self.m,5)*(np.power(xis,5)*np.power(z,5)-5*np.power(xis,4)*np.power(z,4)+20*np.power(xis,3)*np.power(z,3)-60*np.power(xis,2)*np.power(z,2)+120*xis*z-120))

  def Calculate(self, z0, x0, delta0, eta):
      return self.Calculate_xi(z0, x0, delta0, eta, self.xi)
  def Calculate_xi(self, z0, x0, delta0, eta, xi0):
      z0=-z0;

      #Calculate integration limits (dependent on buoy position)
      z1=self.z1+z0;
      z2=np.min([self.z2+z0,eta]);
      q=-self.q-(z0)*self.m;

      empty=[0,0,0]
      empty2=[np.array(xi0)*0,np.array(xi0)*0,np.array(xi0)*0]

      #If segment not in water, return zero
      if (z1>z2):
         return [empty,empty2];
      if (z1>eta):
         return [empty,empty2];


      #Using wheeler strechting to make linear wave theory less linear ;-)
      xis=np.array(xi0*self.depth/(eta+self.depth))#xi dived by water depth (h) and position over the surface eta.


      #buoyancy force (only heave)
      F_st_h=-self.g*self.rho*(2*pi*self.m*(((self.m*np.power(z2,3)/3+q*np.power(z2,2)/2))-(self.m*np.power(z1,3)/3+q*np.power(z1,2)/2)))
      F_st=[0,F_st_h,0]

      #dynamic Froude-Krylov force
      #surge
      c_dy_s=-self.g*self.rho*pi*xi0*(self.As(z2, xi0, xis, eta, q)-self.Bs(z2, xi0, xis, eta, q)-(self.As(z1, xi0, xis, eta, q)-self.Bs(z1, xi0, xis, eta, q)))

      #pitch
      c_dy_p=-self.g*self.rho*pi*xi0*(self.Ap(z2, xi0, xis, eta, q)-xi0**2/8*(self.Bp(z2, xi0, xis, eta, q))-(self.Ap(z1, xi0, xis, eta, q)-xi0**2/8*(self.Bp(z1, xi0, xis, eta, q))));

      #heave
      c_dy_h=(2*self.Ah(z2, xi0, xis, eta, q)-np.power(xi0,2)/2*self.Bh(z2, xi0, xis, eta, q)+np.power(xi0,4)/32*self.Ch(z2, xi0, xis, eta, q)-(2*self.Ah(z1, xi0, xis, eta, q)-np.power(xi0,2)/2*self.Bh(z1, xi0, xis, eta, q)+np.power(xi0,4)/32*self.Ch(z1, xi0, xis, eta, q)))
      c_dy_h=-pi*self.m*self.g*self.rho*c_dy_h;


      if (len(xi0)>1):
          #Avoid numerical errors
          i=np.argmax(np.diff(np.abs(c_dy_h))>0)
          if i>0:
              c_dy_h[i:]=0;
          i=np.argmax(c_dy_h<0)
          if i>0:
              c_dy_h[i:]=0;

      c_dy=[c_dy_s.copy(),c_dy_h.copy(),c_dy_p.copy()*0]

      return [F_st,c_dy];

  def Radius(self,z0): #return radius for this section at z0
      if(z0>=self.z1) and (z0<=self.z2):
          r=-self.q+(z0)*self.m;
          return r;
      else:
          return 0; #return 0 if z0 is not in this section

  def max_Radius(self,z0):
      #get maximum radius for this section if the lowest part of this section is below z0
      if (z0>self.z2):
          return np.max([-self.q+(self.z1)*self.m,-self.q+(self.z2)*self.m]);
      if(z0>=self.z1) and (z0<=self.z2):
          r=np.max([np.abs(-self.q+(z0)*self.m),-self.q+(self.z1)*self.m]);
          return r;
      else:
          return 0;


  def Area (self, z0):
      #get the area at z0
          r=self.Radius(z0);
          return r*r*pi;

  def AreaSurge (self, z0):
          ru=self.Radius(z0);
          rb=self.Radius(self.z1);
          h=np.max([z0-self.z1,0]);
          return (ru+rb)*h;

  def Volume(self, z):
      z1l=np.min([self.z1,z])
      z2l=np.min([self.z2,z])
      z1l1=(z1l-self.q/self.m)
      z2l2=(z2l-self.q/self.m)
      r2=z2l2*self.m;
      r1=z1l1*self.m;
      return 1/3*pi*(r2*r2*z2l2-r1*r1*z1l1);


idname="test.csv";
class Floater:
    volume=0;
    mode=0;
    g=9.81;
    d=300;
    Cog=0;
    rho=1000;
    xi=np.array([]);

    def __init__ (self, xi, g, depth, CoG, *args):
        self.xi=xi;
        self.g=g;
        self.omega=np.sqrt(self.xi*self.g)
        self.elements=[];
        self.d=depth;
        self.CoG=CoG;
        self.file = open(idname+"_f.txt", "w")
        self.file.write("time,wave,buoyancy,exc,rad,slamming\r\n");
        if len(args)>0:
            with open(args[0]) as file:
                geo=json.load(file);
                for g in geo["geo"]:
                    if g["type"] == "cone":
                        self.addCone(g["coord"][0],g["coord"][1],g["coord"][2],g["coord"][3])



    def addCone(self, z1, r1, z2, r2):
        self.elements.append(Cone(z1, r1, z2, r2,self.CoG,self.xi,self.g,self.d))
        self.volume=self.volume+self.elements[-1].Volume(z2);

    def set_mode(self, mode):
        self.mode=mode;

    def Calc_CoG(self):
        mi=0;
        ms=0;
        rge=np.linspace(np.min([e.z1 for e in self.elements]),np.max([e.z2 for e in self.elements]),10);

        for i in rge:
            ai=self.Area(i);
            mi=mi+ai*i;
            ms=ms+ai;
        return mi/ms;

    def get_parameters(self, z0, x0, delta0, eta):
        a=self.Calculate();
        a=[a[9],a[1],a[2],a[3]]
    def Calculate(self, z0, x0, delta0, eta):
        forces=[np.array([0j,0j,0j]),np.array([self.omega*0j,self.omega*0j,self.omega*0j])];

        #Calc added mass
        m_tmp=self.added_mass(z0,eta);
        m0=np.array([m_tmp,m_tmp,m_tmp]);
        z0=-z0;


        #Get buoyancy force and Froude-Krylov forces for each element
        for e in self.elements:
            e1=e.Calculate(z0, x0, delta0, eta);
            for i in range(2):
                forces[i] = forces[i] + e1[i];
        if np.sum(np.abs(forces[0]))==0:#out of water
            en=[self.omega*0,self.omega*0,self.omega*0];
            return [forces[0],en,np.array([en,en,en]),[0,0,0]]

        #First assumption for the excitation force parameters is the the dyn. Froude-Krylov force, neglecting diffraction effects
        c_dy=np.array(forces[1]);
        c_ex=c_dy.copy()*0.8;

        #Initialization of parameters
        c_added_mass=[np.zeros(len(self.omega))+m0[0],np.zeros(len(self.omega))+m0[1],np.zeros(len(self.omega))+m0[2]];
        c_rad=[np.zeros(len(self.omega)),np.zeros(len(self.omega)),np.zeros(len(self.omega))];

        #infinity depth approximation
        D=1 #finity depth: (1+(4*self.xi*self.depth*np.exp(-2*self.depth*self.xi))/(1-np.exp(-4*self.xi*self.depth)))*(1-np.exp(-2*self.xi*self.depth))/(1+np.exp(-2*self.xi*self.depth))


        fdr=np.zeros(len(self.omega));
        fdr=[fdr,fdr,fdr];
        fdi=np.zeros(len(self.omega));
        fdi=[fdi,fdi,fdi];

        #Iterative process to calculate radiation, diffraction and added mass
        for i in range(1):
          for ii in range(len(c_dy)):

              m1=m0[ii]*np.ones(len(self.omega));
              fdr[ii]=1*fdr[ii];
              fdi[ii]=0.24*fdi[ii];
              c_1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*D);
              def func(x):
                  return (c_1*((np.abs(c_dy[ii])+fdr[ii])**2+(x*self.omega)**2))-x;

              #First estimation of radiation force:FrKr forces
              c_ex[ii]=c_dy[ii]#-(fdr[ii]+1j*fdi[ii]);

              #From this first estimation get the radiation forces
              c_rad[ii]=c_1*np.abs(c_ex[ii])**2;


              #get added mass via Kramer Kronig relation;
              om_omega=np.imag(KramerKronig(c_rad[ii]));
              #get diffraction forces from added mass (small-body approximation)
              #use only 75% of the calculated value
              fdr[ii]=0.75*(om_omega-m1)*self.omega**2;

              #filtering: when diffraction is more than 20% of the Froude-Krylov force
              #we are out of the area where this approach is valid, thus limit the diffraction value to this value
              for iii in range(1,len(fdr[ii])):
                  if abs(fdr[ii][iii])>0.2*np.max(abs(c_dy[ii])) or fdr[ii][iii-1]<fdr[ii][iii]:
                      fdr[ii][iii]=fdr[ii][iii-1];

              fdr[ii][c_dy[ii]-fdr[ii]>0]=c_dy[ii][c_dy[ii]-fdr[ii]>0];
              fdr[ii][c_dy[ii]>0]=0;

              #calculate excitation force: FrKr-force + Diffraction force (complex)
              c_ex[ii]=c_dy[ii]-(fdr[ii]+1j*fdi[ii]);

              #Solve equation to get radiation force
              c_rad[ii]=fsolve(func,c_rad[ii],xtol=1e-3);
              #calculate imaginary part of the diffraction force using small-body approximation
              fdi[ii]=-0.75*c_rad[ii]*self.omega;
              #calculate excitation force: FrKr-force + Diffraction force (complex)
              c_ex[ii]=(c_dy[ii]-(fdr[ii]+1j*fdi[ii]));
              #calculate radiation force from the new excitation force
              c_rad[ii]=c_1*np.abs(c_ex[ii])**2;

              c_added_mass[ii]=m1.copy();
        c_rad_cmplx=np.real(c_rad)+1j*np.real(c_added_mass);
        return [forces[0],-1*np.array(c_ex.tolist()),[c_rad_cmplx,c_rad_cmplx,c_rad_cmplx],m0];

    def get_forces(self, t, wave, p, v, a):
        z0=p[1];
        x0=p[0];
        delta0=p[2];

        Awave=wave.get(t,x0);
        eta=np.sum(np.real(Awave[0]));

        res=self.Calculate(z0, x0, 0*delta0, eta);#Calculate coefficents

        ret=[0,0,0];#return array

        def m(a,b):
            return a.real*b.real+a.imag*b.imag;

        exc1 = np.array(res[1]);#Exitation force
        am1 = np.array(res[3]);#added mass @ inf

        if (np.sum(np.abs(exc1))>0):
            c_1=2*3.14*(self.omega*self.omega*self.omega/self.g)/(4*pi*self.g**2*self.rho*1);
            r1=(c_1*exc1[1].real+1j*c_1*exc1[1].imag)*v[1];

            r0=(c_1*exc1[0].real+1j*c_1*exc1[0].imag)*v[0]#1*am_omom[0][0]*(v[0])+am_omom[1][0]*(v[1]);
            wave.add_diracWave(-2/np.pi*r1,t,True);
            wave.add_diracWave2(-2/np.pi*r0,t,True);

        #Calculate hydro forces for each DOF
        for i in range(len(ret)):
            FK=np.sum(m(np.conjugate(exc1[i]),Awave[0])).real;
            ret[i]=res[0][i]+FK;#buoyance + FK force
            if i==1 and v[i]>0:#added mass only implemented for heave
                dP=(am1[i]-self.added_mass(z0+v[i]/abs(v[i])*0.01,eta))/0.01*v[i];
                ret[i]=ret[i]-dP;
        Frad=[np.real(np.sum(wave.get_rad2(t,x0)*(exc1[0]))),np.real(np.sum(m(wave.get_rad(t,x0),exc1[1]))),0];#radiation force
        ret=np.array(ret)+np.array(Frad);

        return [np.real(ret),[am1[0],am1[1],am1[2]]];#hydro force, added mass @ inf


    def get_force_lin(self, t, wave, z0, x0, delta0, v, a):
        ret=self.get_forces(t,wave,0,0,0,v,a);
        ret[1]=ret[1]-self.Area(0)*1000*z0*9.81;
        return ret;

    def Area(self, z0):
        area=0;
        for e in self.elements:
            area = e.Area(z0);
            if area >0:
                return area;
        return area;

    def AreaProjectedHeave(self, z0):
        area=0;
        res=self.getGeoBox();
        zrmax=res[2];
        rmax=res[3];

        if (z0>zrmax):
            return rmax**2*np.pi;
        for e in self.elements:
            area = area + e.Area(z0);
        return area;

    def AreaProjectedSurge(self, z0):
        area=0;
        for e in self.elements:
            area = area + e.AreaSurge(z0);
        return area;

    def getGeoBox(self,z0=np.NaN):
        #gets the max and min heave position (z_min and z_max) and the maximal radius r_max and its vertical position z_r_max
        z_min=np.NaN;
        z_max=np.NaN;
        r_max=0;
        z_r_max=0;

        for e in self.elements:
            z1=e.z1;
            z2=e.z2;
            if (not np.isnan(z0)):
                z1=np.min([z1,z0]);
                z2=np.min([z2,z0]);

            if np.isnan(z_min):
                z_min=z2;
                z_max=z1;

            if  z_min>z1:
                z_min=z1;
            if  z_max<z2:
                z_max=z2;
            if e.Radius(z1)>r_max:
                r_max=e.Radius(z1);
                z_r_max=z1;
            if e.Radius(z2)>r_max:
                r_max=e.Radius(z2);
                z_r_max=z2;
        return (z_min,z_max,z_r_max,r_max)



    def added_mass(self, z0, eta):
        am=0;
        z0=-z0;
        [z_min,z_max,z_r_max,r_max]=self.getGeoBox(z0+eta);
        h=(z_r_max-z_min);
        a=r_max;
        pol=[1.22426941e-05,-4.11227691e-04,5.44460874e-03,-3.57578423e-02,1.20574494e-01,4.36864944e-01,-3.07493344e-02];
        return np.real(np.polyval(pol,a/h)*self.Volume(z_r_max)*self.rho);

    def max_radius(self,z0):
        #get maximal radius for the body below z0
        return np.max([e.max_Radius(z0) for e in self.elements]);


    def radius(self,z0):
        #get radius for the body at z0
        return np.max([e.Radius(z0) for e in self.elements]);

    def Volume(self, z0):
        vol=0;
        for e in self.elements:
            vol = vol + e.Volume(z0);
        return vol;

    def clear(self):
        self.file.close();
        self.elements.clear();
