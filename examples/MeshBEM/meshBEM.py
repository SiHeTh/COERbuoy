#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 15:15:41 2021

@author: heiko
"""

import pygmsh
import json
import os;
import xarray as xr;
from COERbuoy import utils;
from capytaine import FloatingBody, assemble_dataset,RadiationProblem,DiffractionProblem, HierarchicalToeplitzMatrixEngine,BEMSolver;
import numpy as np;
import quadpy
#from quadpy.quadrilateral import stroud_c2_7_2
from capytaine.io.legacy import write_dataset_as_tecplot_files

utils.get()
r_max=0;
z_min=0;
z_max=0;

def calc_params (mesh,p,h,omega_max):
    body=FloatingBody.from_meshio(mesh);
    body.rotated_y(p)
    body.translated_z(h)
    body.show()
    body.add_translation_dof(name="Heave")
    body.add_translation_dof(name="Surge")
    body.add_rotation_dof(name="Pitch")
    body.keep_immersed_part()
    #body.mesh.compute_quadrature(method=quadpy.c2._stroud.stroud_c2_7_2())

    
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, omega_max, 30),
        'wave_direction': [0],
        'radiating_dof': list(body.dofs),
    })
        
    problemR = RadiationProblem(body=body, radiating_dof="Heave", omega=1.0, g=9.81, rho=1000)
    problemD = DiffractionProblem(body=body, omega=1.0);
    dataset=BEMSolver(engine=HierarchicalToeplitzMatrixEngine()).fill_dataset(test_matrix, [body]);#foraxisymetric
    #rR=solver.solve(problemR);
    #rD=solver.solve(problemD);
    #dataset = assemble_dataset(solver.solve)
    folder="results/"+"h_"+str(h)+"_p_"+str(p);
    os.makedirs(folder)
    write_dataset_as_tecplot_files(folder, dataset)
    #print(rR)
    
with pygmsh.occ.Geometry() as geo:
    with open(utils.wec_dir+"/floater.txt",'r') as f:
        data=json.load(f);
        cone=[];
        for g in data["geo"]:
                    if g["type"] == "cone":
                        cone.append(geo.add_cone([0,0,g["coord"][0]],[0,0,-g["coord"][0]+g["coord"][2]],g["coord"][1],g["coord"][3]));
                        if np.max([g["coord"][1],g["coord"][3]])>r_max:
                            r_max=np.max([g["coord"][1],g["coord"][3]]);
                            
                        if g["coord"][0]<z_min:
                            z_min=g["coord"][0];
                        if g["coord"][2]>z_max:
                            z_max=g["coord"][2];
        geo.boolean_union(cone);
        mesh=geo.generate_mesh(dim=2);

zs=np.linspace(z_min,z_max,10);
omega_max=9.81/r_max;
for z in zs[1:-1]:
    calc_params(mesh,0,z,20);
        

