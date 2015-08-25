# -*- coding: utf-8 -*-
"""
/***************************************************************************
 xsection.rivaq.py
                                 CRIVapp

 This file gather functions which provide pre- and post-processing of the 
 numerical model which provide the "river-coefficient" and Xfar. 

 Grid of the numerical model SUTRA are provide by mesh genretor Gmsh

			      -------------------
        begin                : 2015-07-20
        copyright            : (C) 2015 by Cousquer
        email                : yohann.cousquer@ensegid.fr
 ***************************************************************************/
 This plugin use SUTRA VERSION 2.2  
     Copyright (C) 1984-2010  A.M. Provost & C.I. Voss
     WEB : http://water.usgs.gov/nrp/gwsoftware/sutra.html
 
 And Gmsh 2.8
     Copyright (C) 1997-2013 Christophe Geuzaine, Jean-Francois Remacle
     WEB : http://geuz.org/gmsh/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ****
"""


# Load python packages
from subprocess import call
import os 
from itertools import islice
import csv
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata


# ============ Read input parameter file for Gmesh and Sutra ========================================================
#  -- returns dictionary with parameter values
def read_input_file(input_parameter_file = 'param.txt'):
    # init parameter dictionary 
    param = {}
    # -- Read input file and load cross-section parameters
    with open(input_parameter_file,'r') as input_parameter:
	    e = float( input_parameter.readline().split()[0] )
	    l = float( input_parameter.readline().split()[0] )
	    L = float( input_parameter.readline().split()[0] )
	    P = float( input_parameter.readline().split()[0] )
	    a = float( input_parameter.readline().split()[0] )
	    b = float( input_parameter.readline().split()[0] )
	    s1 = float( input_parameter.readline().split()[0] )
	    s2 = float( input_parameter.readline().split()[0] )
	    anis = float( input_parameter.readline().split()[0] )
	    Kh = float( input_parameter.readline().split()[0] )
	    Khb = float( input_parameter.readline().split()[0] )
	    Kvb = float( input_parameter.readline().split()[0] )
	    h_left = float( input_parameter.readline().split()[0] )
	    h_riv = float( input_parameter.readline().split()[0] )
	    h_right = float( input_parameter.readline().split()[0] )
	    calcXfar = float(input_parameter.readline().split()[0] )


    # fill dictionary with parameter values 
    param  = { 'e' : e , 'l' : l, 'L' : L, 'P' : P, 's1' : s1, 's2': s2,
	    'a' : a, 'b' : b, 'anis' : anis, 'Kh' : Kh, 'Khb' : Khb, 'Kvb' : Kvb,
	    'h_left' : h_left, 'h_riv' : h_riv, 'h_right' : h_right, 'calcXfar' : calcXfar}

    return(param)


# ============= Create finite element mesh with Gmesh ============================================================
def build_mesh(param, geo_file = 'x_section.geo', mesh_file = 'x_section.msh', gmesh = '../bin/gmsh' ):
    # param : dictionary with 2D stream-aquifer cross-section characteristics
    # geo_file    : intermediate instruction file read by Gmesh executable.
    # output_file : finite element mesh, to be read by SUTRA.
    # gmesh       : path to Gmesh executable.

    # load parameter values from dictionary
    e = param['e'] 
    l = param['l'] 
    L = param['L'] 
    P = param['P'] 
    a = param['a'] 
    b = param['b'] 
    s1 = param['s1'] 
    s2 = param['s2'] 
    calcXfar = param['calcXfar']


    # compute additional geometry data 
    r4 = (l/2) - (L/2)  # distance between left limit and top of left river bank
    r3 = (l/2) + (L/2)  # distance between left limit and top of right river bank
    r1 = (r4 - ((math.tan(math.radians(a-90))*P)) ) # distance between left limit and bottom of left river bank
    r2 = (r3 + ((math.tan(math.radians(a-90))*P)) )  # distance between left limit and bottom of right river bank
    b1 = e - P # bottom height of left river bank
    b2 = e - P # bottom height of right river bank
    rt1 = r1 - b # distance between left limit and left riverbed top
    rt2 = r2 + b # distance between left limit and right riverbed top
    rt3 = r4 - b # distance between left limit and left riverbed bottom
    rt4 = r3 + b # distance between left limit and right riverbed bottom
    rb1 = b1 - b # distance between river bottom and riverbed bottom
    b2 = e - P # bottom height of right river bank
    o = 0 # origin

    # count number of elements (WHERE???)
    nnodesX = int(float(e)/float(s1)) + 2 #number of nodes at left boundary 
    nelemY = l / s1




    # --  Write .geo FILE 

    geo = open(geo_file, "w")

    geo.write('Point(1) = {0, 0, 0, ' +str(s1)+'};\n')
    geo.write('Point(2) = {0, ' +str(e)+', 0, ' +str(s1)+'};\n') 
    geo.write('Point(3) = {'+str(r1)+', '+str(e)+', 0, ' +str(s2)+'};\n') 
    geo.write('Point(4) = {'+str(r4)+', '+str(b1)+', 0, ' +str(s2)+'};\n') 
    geo.write('Point(5) = {'+str(r3)+', '+str(b2)+', 0, ' +str(s2)+'};\n') 
    geo.write('Point(6) = {'+str(r2)+', '+str(e)+', 0, ' +str(s2)+'};\n') 
    geo.write('Point(7) = {'+str(l)+', '+str(e)+', 0, ' +str(s1)+'};\n') 
    geo.write('Point(8) = {'+str(l)+', 0, 0, ' +str(s1)+'};\n')

    geo.write('Point(9) = {'+str(rt1)+', '+str(e)+', 0, ' +str(s2)+'};\n') 
    geo.write('Point(10) = {'+str(rt3)+', '+str(rb1)+', 0, ' +str(s2)+'};\n')
    geo.write('Point(11) = {'+str(rt4)+', '+str(rb1)+', 0, ' +str(s2)+'};\n')  
    geo.write('Point(12) = {'+str(rt2)+', '+str(e)+', 0, ' +str(s2)+'};\n') 


    geo.write('Line(1) = {1, 2};\n')
    geo.write('Line(2) = {2, 9};\n')
    geo.write('Line(3) = {9, 10};\n')
    geo.write('Line(4) = {10, 11};\n')
    geo.write('Line(5) = {11, 12};\n')
    geo.write('Line(6) = {12, 7};\n')
    geo.write('Line(7) = {7, 8};\n')
    geo.write('Line(8) = {8, 1};\n')

    geo.write('Line(12) = {9, 3};\n')
    geo.write('Line(13) = {3, 4};\n')
    geo.write('Line(14) = {4, 5};\n')
    geo.write('Line(15) = {5, 6};\n')
    geo.write('Line(16) = {6, 12};\n')


    geo.write('Line Loop(17) = {1,2,3,4,5,6,7,8};\n')
    geo.write('Line Loop(19) = {12,13,14,15,16,-5,-4,-3};\n')



    geo.write('Plane Surface(18) = {17, 19};\n')
    geo.write('Plane Surface(20) = {19};\n')

    geo.write('Recombine Surface{18};\n')
    geo.write('Recombine Surface{20};\n')

    # represent centroids of the 3 large-scale cells.
    geo.write('Point(13) = {'+str(l/2)+', '+str(e/2)+', 0, ' +str(s1)+'};\n')
    geo.write('Point{13} In Surface {18};\n')
    geo.write('Point(14) = {'+str(l/6)+', '+str(e/2)+', 0, ' +str(s1)+'};\n')
    geo.write('Point{14} In Surface {18};\n')
    geo.write('Point(15) = {'+str(l*(5./6))+', '+str(e/2)+', 0, ' +str(s1)+'};\n')
    geo.write('Point{15} In Surface {18};\n')

    if calcXfar > 0.: 
        # points to calculate xfar
        xfar_range = np.arange(0.025, 1, 0.025)
        xfar_range = xfar_range.ravel()
        point_id = np.arange(16, 133, 3)
        point_id = point_id.ravel()
        x_far = zip(xfar_range, point_id)
        for i,j in x_far:
	        geo.write('Point('+str(j)+') = {'+str((l/2)+(i)*(l/2))+', '+str(e-(b+P+3))+', 0, ' +str(s1)+'};\n')
	        geo.write('Point{'+str(j)+'} In Surface {18};\n')
	        geo.write('Point('+str(j+1)+') = {'+str((l/2)+(i)*(l/2))+', '+str(e-(e-5))+', 0, ' +str(s1)+'};\n')
	        geo.write('Point{'+str(j+1)+'} In Surface {18};\n')
	        geo.write('Point('+str(j+2)+') = {'+str((l/2)+(i)*(l/2)+1)+', '+str(e-(b+P+3))+', 0, ' +str(s1)+'};\n')
	        geo.write('Point{'+str(j+2)+'} In Surface {18};\n')


    geo.close()

    # --- Call Gmesh to build mesh 
    call([gmesh,"-2", geo_file,"-o",mesh_file])


# ======== Write files input files for SUTRA ============================================================
def sutra_preproc(param, 
	mesh_file = 'x_section.msh', nodes_file = 'nodes.csv', elements_file = 'element.csv',
	speci_nodes_file = 'speci_nodes.csv', sutra_files_prefix = 'sutra_model'):
    # param : dictionary with 2D stream-aquifer cross-section characteristics
    # mesh_file : mesh file generated by Gmesh
    # nodes_file : nodes file generated by Gmesh
    # elements_file : elements file generated by Gmesh
    # speci_nodes_file : output boundary condition file 
    # sutra_files_prefix : prefix for all sutra files 

    # name sutra files from prefix
    sutra_fil_file = 'SUTRA.FIL' 
    sutra_inp_file = sutra_files_prefix + '.inp'
    sutra_ics_file = sutra_files_prefix + '.ics'
    sutra_lst_file = sutra_files_prefix + '.lst'
    sutra_rst_file = sutra_files_prefix + '.rst'
    sutra_nod_file = sutra_files_prefix + '.nod'
    sutra_ele_file = sutra_files_prefix + '.ele'
    sutra_obs_file = sutra_files_prefix + '.obs'
    sutra_smy_file = sutra_files_prefix + '.smy'

    # Load parameter from param dictionary

    # load parameter values from dictionary
    e = param['e'] 
    l = param['l'] 
    L = param['L'] 
    P = param['P'] 
    a = param['a'] 
    s1 = param['s1'] 
    s2 = param['s2'] 
    anis = param['anis']
    Kh = param['Kh']
    Khb = param['Khb']
    Kvb = param['Kvb']
    h_left = param['h_left']
    h_riv = param['h_riv']
    h_right = h_left
    Temp = 15.0

    # compute additional geometry data 
    r1 = (l/2) - (L/2)  # distance between left limit and top of left river bank
    r2 = (l/2) + (L/2)  # distance between left limit and top of right river bank
    r4 = (r1 + ((math.tan(a-90)*P)) ) # distance between left limit and bottom of left river bank
    r3 = (r2 - ((math.tan(a-90)*P)) ) # distance between left limit and bottom of right river bank
    b1 = e - P # bottom height of left river bank
    b2 = e - P # bottom height of right river bank
    o = 0 # origin

    # Compute additional parameters
    Kv = 1./ anis * Kh
    Km = 2.e-10
    
    # WHAT ???
    BoundNodes = ( ( ((e+e+l)+(r1+(l-r2)))/s1) + 
	    ((r3-r4) + (math.sqrt(((r2-r3)**2) + ((e-b2)**2)) +  math.sqrt(((r4-r1)**2) + ((e-b1)**2)))) / s2 
	    + (100-(r2-r1))/2)  #approximate
    BoundNodes = int(BoundNodes)

    # Transport (unused, but required)
    ANGLES1 = 0.
    ALMAX = 10.
    ALMIN = 2.5
    ATMAX = 0.1
    ATMIN = 0.1

    # -- Simulation type
    simul = 'steady'
    #simul = 'transient'

    # -- Solver
    #solver = 'DIRECT'
    solver = 'ITERATIVE'

    # -- Create nodes & element files ---
    # SHOULD BE RE-ORGANIZED SO AS TO OPEN THE MESH FILE ONLY ONE TIME
    # - Nodes 
    lookup = '$EndNodes'
    mesh = open(mesh_file, "r")
    for num, line in enumerate(mesh, 1):
	    if lookup in line: 
		    EndNodes = num 

    End = EndNodes - 1

    with open(mesh_file,'r') as input: 
	    nodes_values = islice(input, 5,  End )    
	    with open('nodes.csv','w') as output:
		    output.writelines(nodes_values)

    # - Elements
    lookup1 = '1 2 0 15'
    lookup2 = '1 2 0 16'
    lookup3 = '$EndElements'

    mesh = open(mesh_file,"r")
    for num, line in enumerate(mesh, 1):
	    if lookup1 in line: 
		    BeginElement = num
	    elif lookup2 in line: 
		    BeginElement = num

    Begin = BeginElement  

    mesh = open(mesh_file,"r")
    for num, line in enumerate(mesh, 1):
	    if lookup3 in line: 
		    EndElement = num

    End = EndElement -1 

    with open(mesh_file,'r') as input: 
	    element_values = islice(input,Begin,End)    
	    with open('element.csv','w') as output:
		    output.writelines(element_values)

    # --------------  Model properties  ---------------------------

    # -- Specified Nodes 
    with open(nodes_file,'rb') as csvfile: 
	    nodes = [row for row in csv.reader(csvfile, delimiter = ' ')]
	    nodes = np.array(nodes)
	    with open(speci_nodes_file,'w') as output:
		    SpeciNodes = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in nodes:
			    pleft = 1000 * 9.81 *  ( h_left - float(row[2])) + 101325.0
			    h = ((pleft - 101325.0) /(1000*9.81)) + float(row[2]) 
			    if float(row[1]) == o and  float(row[2]) <= h_left :
				    SpeciNodes.writerow((row[0], str(pleft), str(Temp), str(h)))
		    for row in nodes:
			    pright = 1000 * 9.81 *  ( h_right - float(row[2])) + 101325.0
			    h = ((pright - 101325.0) /(1000*9.81)) + float(row[2]) 
			    if float(row[1]) == l  and  float(row[2]) <= h_right :
				    SpeciNodes.writerow((row[0], str(pright), str(Temp), str(h)))
		    for row in nodes :#[0:BoundNodes,:]:
			    priv = 1000 * 9.81 * ( h_riv - float(row[2])) + 101325.0
			    h = ((priv - 101325.0 )/(1000*9.81)) + float(row[2]) 
			    if (\
			    float(row[2]) >= (e - P) and 
			    float(row[1]) >= r1 and float(row[1]) <= r2 and
			    float(row[2]) <= h_riv ): 
				    SpeciNodes.writerow((row[0], str(priv), str(Temp), str(h)))

    with open(speci_nodes_file,'r') as input: # OPTIMIZE SO AS TO AVOID READING A FILE YOU HAVE JUST WRITTEN
	    nspnodes = sum(1 for row in input) # count specified nodes number#----- Specified Concentrations -------


    # --- WRITE INP FILE ---------------------------
    conceptual = open(sutra_inp_file, "w")

    conceptual.write('This is a SUTRA input conceptual.inp prepared with Python script\n')
    conceptual.write('Valid for 2D regular meshes\n')

    #------ Dataset 2A - Simulation type
    conceptual.write('# Dataset 2A\n')
    conceptual.write('\'SUTRA VERSION 2.2 SOLUTE TRANSPORT\'\n')

    #------ Dataset 2B - Mesh structure 
    conceptual.write('# Dataset 2B\n')
    conceptual.write('\'2D IRREGULAR MESH\'\n')

    #------ Dataset 3 : Simulation control numbers
    conceptual.write('# Dataset 3\n')

    with open(nodes_file,'r') as input: # MAY BE A WAY TO AVOID RE-OPENING THIS FILE
	    nnodes = sum(1 for row in input) # count nodes number

    with open(elements_file,'r') as input: 
	    nelement = sum(1 for row in input) # count nodes number


    conceptual = open(sutra_inp_file, "a")
    conceptual.write(\
    str(nnodes)+' '  # number of nodes 
    +str(nelement)+' '  # number of elements 
    +str(nspnodes)+' ' # number of nodes with specified pressure (right boundary)
    ' 0' # number of nodes with specified concentration
    ' 0' # number of nodes with fuild source/sink (top boundary, recharge)
    ' 0' # number of nodes with energy source/sink
    ' 0' # number of points at which observation will be made
    '\n')

    #------ Dataset 4 :Simulation mode options
    conceptual.write('# Dataset 4\n')

    if simul == 'steady':
	    conceptual.write('\'UNSATURATED\' \'STEADY FLOW\' \'STEADY TRANSPORT\' \'COLD\' ' +str(9999)+'\n')

    if simul == 'transient':
	    conceptual.write('\'SATURATED\' \'TRANSIENT FLOW\' \'TRANSIENT TRANSPORT\' \'COLD\' ' +str(9999)+'\n')

    #------ Dataset 5: Numerical control parameters
    conceptual.write('# Dataset 5\n')
    conceptual.write('0.0 0.1 1. \n')

    #------ Dataset 6: temporal control and solution cycling data
    conceptual.write('# Dataset 6\n')

    if simul == 'steady' :
	    conceptual.write('0 1 1 \n')
	    conceptual.write('- \n')

    if simul == 'transient':
	    conceptual.write(\
	    '1' # number of schedules
	    ' 1' # number of time steps in pressure solution cycle
	    ' 1' # number of time steps in concentration solution cycle
	    '\n')
	    conceptual.write(\
	    ' \'TIME_STEPS\'', # Initial schedule name, required for transient transport
	    ' \'TIME CYCLE\'', # Schedule type, "time cycle" for sequence of times generated at specified intervals
	    ' \'ABSOLUTE\'', # Simulation clock, "elapsed" for t=0 when starting simulation
	    ' 1' # Scale factor to be applied to each time value in the list
	    ' 1500' # Maximum number of times allowed
	    ' 0' # Initial time
	    ' 1.e+99' # Limiting time
	    ' 60' # Initial time increment
	    ' 1000' # Number of cycles after which the time increment is updated
	    ' 2' # Factor by which the time increment is multiplied after preceding number of cycles
	    ' 1.e-20' # Minimum time increment allowed
	    ' 1.e+99' # Maximum time increment allowed
	    '\n')
	    conceptual.write('-', '\n')

    #------ Dataset 7: solver control

    if solver == 'DIRECT':
	    conceptual.write(\
	    '# Dataset 7A\n'
	    '1 \n'
	    '# Dataset 7B\n'
	    '\'DIRECT\'\n'
	    '# Dataset 7C\n'
	    '\'DIRECT\'\n')

    if solver == 'ITERATIVE':
	    conceptual.write(\
	    '# Dataset 7A\n'
	    '1\n'
	    '# Dataset 7B\n'
	    '\'CG\' 3000 5.e-14\n'
	    '# Dataset 7C\n'
	    '\'ORTHOMIN\' 1 1.e-10\n')


    #------ Dataset 8: output controls
    conceptual.write('# Dataset 8A\n')
    conceptual.write('50 \'N\' \'N\' \'N\' \'Y\' \'Y\' \'Y\' \'Y\' \'Y\' \'N\'\n')

    conceptual.write('# Dataset 8B\n') # output control of the .nod file
    conceptual.write('-200 \'X\' \'Y\' \'P\' \'U\' \'S\' \'-\' \n')

    conceptual.write('# Dataset 8C\n') # output control of the .ele file
    conceptual.write('-200 \'X\' \'Y\' \'VX\' \'VY\' \'-\' \n') 

    conceptual.write('# Dataset 8D\n') # Observation points

    conceptual.write('# Dataset 8E\n')
    conceptual.write('9999 9999 9999 9999 \'Y\' \n')

    #------ Dataset 9: fluid properties
    conceptual.write('# Dataset 9\n')

    conceptual.write(\
    '4.47e-10 ' # Fluid compressibility
    '1 ' # Fluid specific heat [kg/m.s2]
    '1.5e-09 ' #Fluid diffusivity [m2/s]
    '1000 ' # Density of fluid at base concentration
    '0 ' # Base value of solute concentration
    '700 ' # Coefficient of fluid density change with concentration
    '0.001002 ' # Fluid viscosity [kg/m.s]
    '\n')

    #------ Dataset 10: solid matrix properties
    conceptual.write('# Dataset 10\n')

    conceptual.write(\
    '2.5e-9 ' # solid matrix compressibility
    '0. ' # solid grain specific heat
    '0. ' # solid grain diffusivity
    '2600.' # density of a grain
    '\n') 

    #------ Dataset 11: adsorption parameters
    conceptual.write('# Dataset 11\n')
    conceptual.write('\'NONE\' \n')

    #------ Dataset 12: production of energy or solut mass
    conceptual.write('# Dataset 12\n')
    conceptual.write('0. 0. 0. 0. \n')

    #------ Dataset 13: orientation of gravity vector
    conceptual.write('# Dataset 13\n')
    conceptual.write('-9.81 0. 0.\n')

    #------ Dataset 14: nodewise data
    # scale factors
    conceptual.write('# Dataset 14A\n')

    conceptual.write('\'NODE\' 1. 1. 1. 1. \n')

    # nodewise data (porosity) (one line for each of NN nodes)
    conceptual.write('# Dataset 14B\n')

    conceptual.close()

    with open(nodes_file,'rb') as csvfile: # 
	    nodes = csv.reader(csvfile, delimiter = ' ')
	    with open(sutra_inp_file,'a') as output: # open write access to conceptual file with csv
		    conceptual = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in nodes:
			    conceptual.writerow((row[0], str(0), row[2], row[1], str(1) , str(0.1)))

    #------ Dataset 15: elementwise data
    conceptual = open(sutra_inp_file, "a")

    # scale factors
    conceptual.write('# Dataset 15A\n')
    conceptual.write('\'ELEMENT\' 1. 1. 1. 1. 1. 1. 1.\n')

    # elementwise data (permeability,dispersion)
    conceptual.write('# Dataset 15B\n')

    with open(elements_file,'rb') as csvfile:
	    element = [row for row in csv.reader(csvfile, delimiter = ' ')]
	    element = np.array(element)
	    nbegin = element[0,0]
	    nbegin = int(nbegin)

    with open(elements_file,'rb') as csvfile:
	    element = csv.reader(csvfile, delimiter = ' ')
	    with open(sutra_inp_file,'a') as output:
		    conceptual = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in element:
			    ID = int(row[0]) - (nbegin -1)
			    if float(row[4]) == 20:
				    conceptual.writerow(\
				    (ID, str(0) , str(Khb) , str(Kvb),
				    str(ANGLES1), str(ALMAX), str(ALMIN) ,
				    str(ATMAX), str(ATMIN)))
			    else :
				    conceptual.writerow(\
				    (ID, str(0) , str(Kh) , str(Kv),
				    str(ANGLES1), str(ALMAX), str(ALMIN) ,
				    str(ATMAX), str(ATMIN)))


    #------ Dataset 17: fluid sources and sinks
    conceptual = open(sutra_inp_file, "a")
    conceptual.write('# Dataset 17\n')
    
    #------ Dataset 18: energy or solute mass sources and sinks
    conceptual.write('# Dataset 18\n')

    #------ Dataset 19: specified pressure nodes
    conceptual.write('# Dataset 19\n')

    with open(speci_nodes_file,'rb') as csvfile:
	    spnodes = csv.reader(csvfile, delimiter = ' ')
	    with open(sutra_inp_file,'a') as output:
		    conceptual = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in spnodes:
			    conceptual.writerow((row[0], row[1], row[2]))

    conceptual = open(sutra_inp_file, "a")

    conceptual.write(str(0)+'\n')

    #------------------------------------------------------------
    #------ Dataset 20: specified concentration or temperature nodes
    conceptual.write('# Dataset 20\n')

    #with open('./Gmsh/SpeciC.csv','rb') as csvfile:
    #	spc = csv.reader(csvfile, delimiter = ' ')
    #	with open(sutra_inp_file,'a') as output:
    #		conceptual = csv.writer(output, delimiter = ' ', lineterminator = '\n')
    #		for row in spc:
    #			conceptual.writerow((row[0], row[1]))
    #
    #conceptual = open(sutra_inp_file, "a")
    #
    #conceptual.write(str(0)+'\n')

    #------ Dataset 22: element incidence
    conceptual.write('# Dataset 22\n')

    # Element incidence data (cell connectivity)
    conceptual.write('\'INCIDENCE\'\n')

    with open(elements_file,'rb') as csvfile:
	    element = csv.reader(csvfile, delimiter = ' ')
	    with open(sutra_inp_file,'a') as output:
		    conceptual = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in element:
			    ID = int(row[0]) - (nbegin -1)
			    conceptual.writerow((ID, row[5], row[6], row[7], row[8]))

    # --------------------------------------------------------
    # --- Write file -----------------------------------------

    conceptual = open(sutra_ics_file, "w")

    #-------------------------------------------------------
    #------ Dataset 1 - Simulation  row[12]Starting Time
    conceptual.write('# Dataset 1\n')
    conceptual.write('0\n')

    #-------------------------------------------------------
    #------ Dataset 2 - Initial Pressure Value at Nodes
    conceptual.write('# Dataset 2\n')
    conceptual.write('\'NONUNIFORM\'\n')
    #conceptual.write(str(0)+'\n')


    with open(nodes_file,'rb') as csvfile:
	    nodes = csv.reader(csvfile, delimiter = ' ')
	    with open(sutra_ics_file,'a') as output:
		    conceptual = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in nodes:
			    conceptual.writerow(str(0))
    # ------ Dataset 3 - Initial Temperature or Concentration Value et Nodes

    conceptual = open(sutra_ics_file, "a")

    conceptual.write('# Dataset 3\n')
    conceptual.write('\'NONUNIFORM\'\n')

    with open(nodes_file,'rb') as csvfile:
	    nodes = csv.reader(csvfile, delimiter = ' ')
	    with open(sutra_ics_file,'a') as output:
		    conceptual = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in nodes:
			    conceptual.writerow(str(0))


    #-------------------------------------------------------
    # ------  WRITE SUTRA.FIL 

    SUTRA = open(sutra_fil_file, "w")

    SUTRA.write('INP '+ '50 ' + sutra_inp_file +'\n')
    SUTRA.write('ICS '+ '55 ' + sutra_ics_file +'\n')
    SUTRA.write('LST '+ '60 ' + sutra_lst_file +'\n')
    SUTRA.write('RST '+ '66 ' + sutra_rst_file +'\n')
    SUTRA.write('NOD '+ '30 ' + sutra_nod_file +'\n')
    SUTRA.write('ELE '+ '40 ' + sutra_ele_file +'\n')
    SUTRA.write('OBS '+ '70 ' + sutra_obs_file +'\n')
    SUTRA.write('SMY '+ '80 ' + sutra_smy_file +'\n')

    SUTRA.close()

# === RUN SUTRA =========================================================================
def sutra_run(sutra = '../bin/sutra'):
    call([sutra])


# ==== BUDGET ===========================================================================
def sutra_budget(param,
    speci_nodes_file = 'speci_nodes.csv', sutra_lst_file = 'sutra_model.lst', nodes_file = 'nodes.csv', 
    main_budget_file = 'budget.csv', river_budget_file = 'river_budget.txt', sutra_nod_file = 'sutra_model.nod', head_file = 'head.csv', head_aq = 'head_cell_aq.csv',  head_xfar = 'head_xfar.csv', calc_xfar = 'xfar.csv' ):
    # sutra_lst_file : SUTRA  lst file	
    # main_budget_file : output csv file with budget for all cells (?)
    # river_budget_file : output csv file with budget for river cells (?)

    e = param['e'] 
    P = param['P']
    b = param['b'] 
    s1 = param['s1'] 
    s2 = param['s2']
    l = param['l'] 
    Kh = param['Kh']
    calcXfar = param['calcXfar']
    nnodesX = float(e)/float(s1) + 2 #number of nodes at left boundary 

    with open(speci_nodes_file,'r') as input: 
	    nspnodes = sum(1 for row in input)



    # read sutra output lst file and search sources or sinks section and extract river nodes fluid sources
    sutra_lst = open(sutra_lst_file,"r") 
    lookup_begin = '                          4' # the first river node specified head is always node number 4
    lookup_end = 'RESULTS FOR TIME STEP        1' # identification of the end the following section
    for num, line in enumerate(sutra_lst, 1):
	    if lookup_begin in line: 
		    begin = num -1
	    if lookup_end in line: 
		    end = num - 3 #the end of the sources section is 3 lines before 'lookup_end'
		    break
    sutra_lst.close()

    main_budget = []
    with open(sutra_lst_file,'r') as input:  # AVOID MULTIPLE OPENING OF THE SAME FILE
	    river_budget_section = islice(input, begin,  end )
	    for line in river_budget_section:
		 main_budget.append(float(line.split()[1]))

    main_budget = np.array(main_budget)
    total_flow = main_budget.sum()/1000 # conversion of kg water to m3 water

    # write aquifer-to-stream flow into river budget file
    with open(river_budget_file,'w') as river_budget :
	river_budget.write(str(total_flow))

    # read sutra output nod file and search sources or sinks section and extract river nodes fluid sources

    with open(nodes_file,'r') as input: # MAY BE A WAY TO AVOID RE-OPENING THIS FILE
	    nnodes = sum(1 for row in input) # count nodes number

    with open(sutra_nod_file,'rb') as csvfile:
    	    nodes = csv.reader(csvfile, delimiter = ' ')
	    with open(head_file,'w') as output:
		    head = csv.writer(output, delimiter = ' ', lineterminator = '\n')
		    for row in islice(nodes, 18, int(nnodes)):
			    h = (((float(row[7])) - 101325.0)/(1000*9.81)) + float (row[3])
			    if float(row[7]) >= float(101325.0):
				    head.writerow((row[3], row[5], row[7], str(h))) #x,y,p,cc,h


    node_file = np.genfromtxt(head_file, delimiter=' ')
    idx1 = np.logical_and(node_file[:,1] == float(l/6), node_file[:,0] == float(e/2.))
    idx2 = np.logical_and(node_file[:,1] == float(l*(5./6)), node_file[:,0] == float(e/2.))
    idx3 =  np.logical_and(node_file[:,1] == float(l/2), node_file[:,0] == float(e/2.))
    
    if calcXfar > 0.: 
        calc = open(calc_xfar, "w+")
        calc.close

        xfar = []
        xfar_range = np.arange(0.025, 1, 0.025)
        xfar_range = xfar_range.ravel()
        point_id = np.arange(16, 133, 3)
        point_id = point_id.ravel()
        x_far = zip(xfar_range, point_id)
        for i,j in x_far:
	    idxfarb =  np.logical_and(node_file[:,1] == float((l/2)+(i)*(l/2)), node_file[:,0] == float(e-(b+P+3)))
	    idxfara =  np.logical_and(node_file[:,1] == float((l/2)+(i)*(l/2)), node_file[:,0] == float(e-(e-5.)))
	    idxfarc =  np.logical_and(node_file[:,1] == float((l/2)+(i)*(l/2)+1), node_file[:,0] == float(e-(b+P+3)))
	    afar = float(node_file[:,3][idxfara]) 
    	    bfar = float(node_file[:,3][idxfarb])
	    cfar = float(node_file[:,3][idxfarc])
	    opp = float((afar -  bfar)/(cfar - bfar))
 	    angle = math.degrees(math.atan(opp / ((e-(P+5.))-(e-(e-5.)))))
	    comp_x = (90 - angle) *(100/90.)
   	    with open(calc_xfar,'a') as csvfile : 
	        x_far = csv.writer(csvfile, delimiter = ' ', lineterminator = '\n')
                x_far.writerow((i, comp_x))



    Hl = node_file[:,3][idx1]
    Hr = node_file[:,3][idx2]
    H = node_file[:,3][idx3]
    Hr = float(Hr)
    Ha = ( ((float(total_flow) * (l/3))/float((Kh/1e-7)*e)) + float(Hl) + float(Hr)) /2

    with open(head_aq,'w') as aquifer_head :
	aquifer_head.write(str(Ha))

    with open(head_xfar,'w') as xfar_head :
	xfar_head.write(str(Hr))
 
#    with open(calc_xfar,'w') as csvfile : 
#	x_far = csv.writer(csvfile, delimiter = ' ', lineterminator = '\n')
#        x_far.writerow(calc_xfar)




# === Run model ================================================================================ 
def gen_mesh(input_parameter_file = 'param.txt', river_budget_file = 'river_budget.txt', mesh=True, gmesh = '../bin/gmsh', sutra = '../bin/sutra'):

    # Read cross-section parameter file
    param = read_input_file(input_parameter_file)  


    if mesh == True : 
	# Build mesh with Gmesh
	build_mesh(param, gmesh = gmesh)

def run(input_parameter_file = 'param.txt', river_budget_file = 'river_budget.txt', mesh=True, gmesh = '../bin/gmsh', sutra = '../bin/sutra'):

    # Read cross-section parameter file
    param = read_input_file(input_parameter_file)  

    # Pre-processing for SUTRA : write input files
    sutra_preproc(param)

    # Run SUTRA
    sutra_run()

    # Compute budget
    sutra_budget(param, river_budget_file = river_budget_file)




