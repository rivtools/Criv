# -*- coding: utf-8 -*-
"""
/***************************************************************************
 prep_file.py
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

import subprocess
import numpy as np
import math
import os 
import matplotlib.pyplot as plt
import csv
import sys
from pylab import *
import itertools
from scipy import stats
from itertools import islice
from numpy import genfromtxt

# Load xsection_rivaq.py script 
import xsection_rivaq


# ============ Read input parameter file for Gmesh and Sutra ========================================================
#  -- returns dictionary with parameter values
def read_input_file(input_parameter_file = '../user_param.txt'):
    # init parameter dictionary 
    param = {}
    # -- Read input file and load cross-section parameters
    with open(input_parameter_file,'r') as input_parameter:
	    ct = float( input_parameter.readline().split()[1] )
	    cw = float( input_parameter.readline().split()[1] )
	    w = float( input_parameter.readline().split()[1] )
	    d = float( input_parameter.readline().split()[1] )
	    a = float( input_parameter.readline().split()[1] )
	    m = float( input_parameter.readline().split()[1] )
	    anis = float( input_parameter.readline().split()[1] )
	    kh = float( input_parameter.readline().split()[1] )
	    khb = float( input_parameter.readline().split()[1] )
	    dh = float( input_parameter.readline().split()[1] )

    # fill dictionary with parameter values 
    param  = { 'ct' : ct , 'cw' : cw, 'w' : w, 'd' : d, 'a' : a, 'm': m,
	    'anis' : anis, 'kh' : kh, 'khb' : khb, 'dh' : dh}

    return(param)


#=================================================================================================
#==== CREATE GEOMETRY ====================================

def build_model(param, model_shape = '../output/model.png'):
	#---  Geometry - User input
	ct = param['ct'] # cell thickness [m]
	cw = param['cw'] # cell width [m]
	w = param['w'] # river width [m]
	d = param['d'] # river depth [m]
	a = param['a'] # bank angle [°]
	m = param['m']#riverbed thikness [m]
	#==== build geometry ============================================================================
	r4 = (cw/2) - (w/2)  # distance between left limit and top of left river bank
	r3 = (cw/2) + (w/2)  # distance between left limit and top of right river bank
	r1 = (r4 - ((math.tan(math.radians(a-90))*d)) ) # distance between left limit and bottom of left river bank
	r2 = (r3 + ((math.tan(math.radians(a-90))*d)) ) # distance between left limit and bottom of right river bank
	b1 = ct - d # bottom height of left river bank
	b2 = ct - d # bottom height of right river bank
	rt1 = r1 - m # distance between left limit and left riverbed top
	rt2 = r2 + m # distance between left limit and right riverbed top
	rt3 = r4 - m # distance between left limit and left riverbed bottom
	rt4 = r3 + m # distance between left limit and right riverbed bottom
	rb1 = b1 - m # elevation of riverbed bottom
	b2 = ct - d # bottom height of right river bank
	xg = (0,0,str(r1),str(r4),str(r3),str(r2),str(cw),str(cw), 0)
	yg = (0,str(ct),str(ct),str(b1),str(b1),str(ct),str(ct),0,0)
	xh = (str(rt1),str(rt3), str(rt4), str(rt2))
	yh = (str(ct),str(rb1), str(rb1), str(ct))
	plt.scatter(xg, yg)
	plt.scatter(xh, yh)
	if m == 0:
		plt.plot(xg, yg)
	else : 
		plt.plot(xg, yg)
		plt.plot(xh, yh)
	savefig(model_shape)
	plt.close()


#=================================================================================================
#==== Compute CRIV ===================================

def compute_CRIV(param, criv = './CRIV.csv', sutra_inp_table = 'param_table.csv',  
	sutra_inp_file = 'param.txt', q_riv = 'river_budget.txt', aq_head = 'head_cell_aq.csv'):
	#---  Geometry - User input
	ct = param['ct'] # cell thickness [m]
	cw = param['cw'] # cell width [m]
	w = param['w'] # river width [m]
	d = param['d'] # river depth [m]
	a = param['a'] # bank angle [°]
	m = param['m']#riverbed thikness [m]
	s1 = 1.# (cw)/100. # general mesh element size [m]
	s2 = 0.2#(cw)/100. # river mesh element size [m]

	#---  Geometry - Model input
	cw = cw * 3.
	d = d + 1.
	#---   Hydrodynamic properties 
	anis = param['anis']#anisotropy = Kh/Kv
	Kh = param['kh']  #Permeability in [L^2]  1m/s = 10^-7m^2      
	Khb = param['khb'] #riverbed horizontal permeability [L^2]
	Kh = Kh * 1e-7 
	Khb = Khb * 1e-7 
	if m == 0:
		Khb = Kh

	if m == 0:
		m = ct/20. 

	Kvb = Khb 
	#----- Flow Boundary conditions 
	h_riv = (ct - 1) # river head
	h_left = np.arange(h_riv - 5, h_riv + 5.5, 0.5) # left head 
	h_right = 38
	#----- Calc Xfar (need to be 0 in this operation) 
        calcXfar =  0
	#tested parameter
	Param = h_left     
	#output file
	param_output = criv
	conceptual = open(param_output, "w+")
	conceptual.close
	with open(sutra_inp_table,'w') as output:
		param_table = csv.writer(output, delimiter = ' ', lineterminator = '\n') #write table for each value of tested parameter 
		for i in Param:
			param_table.writerow((\
			float(ct), float(cw), float(w), float(d),
			float(a), float(m), float(s1), float(s2),
			float(anis), float(Kh), float(Khb), float(Kvb),
			float(i),float(h_riv),float(h_right),calcXfar )) 

	with open(sutra_inp_table) as csv_data:
		reader = csv.reader(csv_data)
		rows = [row for row in reader if row]
		for row in rows:
			row = row
			lis=[x.split() for x in row]
			with open(sutra_inp_file,'w') as output:
				for x in zip(*lis):
					for y in x:
						output.write(y+'\n')
			xsection_rivaq.gen_mesh()
			break
		for row in rows:
			row = row
			lis=[x.split() for x in row]
			with open(sutra_inp_file,'w') as output:
				for x in zip(*lis):
					for y in x:
						output.write(y+'\n')
			xsection_rivaq.run() 
			with open(param_output,'a') as output:
				with open(sutra_inp_file,'r') as input_parameter:
					ct = float( input_parameter.readline().split()[0] )
					cw = float( input_parameter.readline().split()[0] )
					w = float( input_parameter.readline().split()[0] )
					d = float( input_parameter.readline().split()[0] )
					a = float( input_parameter.readline().split()[0] )
					m = float( input_parameter.readline().split()[0] )
					s1 = float( input_parameter.readline().split()[0] )
					s2 = float( input_parameter.readline().split()[0] )
					anis = float( input_parameter.readline().split()[0] )
					Kh = float( input_parameter.readline().split()[0] )
					Khb = float( input_parameter.readline().split()[0] )
					Kvb = float( input_parameter.readline().split()[0] )
					h_left = float( input_parameter.readline().split()[0] )
					h_riv = float( input_parameter.readline().split()[0] )
					h_right = float( input_parameter.readline().split()[0] )
				parameter = csv.writer(output, delimiter = ' ', lineterminator = '\n')
				flow = np.loadtxt(q_riv)
				TOTAL_FLOW = float(flow)
				cell_head = np.loadtxt(aq_head)
				cell_head = float(cell_head)
				parameter.writerow((TOTAL_FLOW, cell_head-h_riv))



#=================================================================================================
#==== Plot regression line ===================================

def plot_CRIV(criv = './CRIV.csv', plot = '../output/regression_line.png'):
	data = np.genfromtxt(criv, dtype=[('Q',float),('param',float)], delimiter = ' ')
	Q = data['Q']
	param = data['param']
	#stat
	x = param[:,np.newaxis]
	y= Q
	a, _, _, _ = np.linalg.lstsq(x, y)
	ylabel('QRIV [m3/s]', fontsize =17 )
	xlabel(' Delta h [m]', fontsize =17)
	grid(True)
	plt.plot(x, y, 'bo')
	plt.plot(x, a*x, 'r-')
	savefig(plot)
	plt.close()


#=================================================================================================
#==== Get CRIV ===================================

def write_CRIV(param, criv = './CRIV.csv', criv_value = '../output/CRIV_value.txt', R2_value = '../output/R2_value.txt'):
	Kh = param['kh']        
	data1 = np.genfromtxt(criv, dtype=[('Q',float),('param',float)], delimiter = ' ')
	Q = data1['Q']
	param = data1['param']
	#stat
	x = param[:,np.newaxis]
	y= Q
	a, _, _, _ = np.linalg.lstsq(x, y)
	slope, intercept, r_value, p_value, std_err = stats.linregress(param,Q)
	#CRIV
	CRIV = (float(-a))
	criv_file = open(criv_value, "w")
	criv_file.write(''+str(CRIV)+'\n')
	criv_file.close()
	#R2
	r_squared = r_value**2
	R2_file = open(R2_value, "w")
	R2_file.write(''+str(r_squared)+'\n')
	R2_file.close()



#=================================================================================================
#==== Calculate Xfar ===================================

def compute_Xfar(param, xfar = './CRIV.csv', sutra_inp_table = 'param_table.csv',  
	sutra_inp_file = 'param.txt', q_riv = 'river_budget.txt', aq_head = 'head_cell_aq.csv', calc_xfar = 'xfar.csv', plot = '../output/plot_xfar.png'):
	#---  Geometry - User input
	ct = param['ct'] # cell thickness [m]
	w = param['w'] # river width [m]
	d = param['d'] # river depth [m]
	a = param['a'] # bank angle [°]
	m = param['m']#riverbed thikness [m]
	dh = param['dh']

	#---  Geometry - Model input
	cw = w*10 # cell width [m]
	cw = cw * 3.
	d = d + 1.
	#---mesh element size [m]
	s1 = 1.#(cw)/100. # general mesh element size [m]
	s2 = 0.2#(cw)/100. # river mesh element size [m]
	#---   Hydrodynamic properties 
	anis = param['anis']#anisotropy = Kh/Kv
	Kh = param['kh']  #Permeability in [L^2]  1m/s = 10^-7m^2      
	Khb = param['khb'] #riverbed horizontal permeability [L^2]
	Kh = Kh * 1e-7 
	Khb = Khb * 1e-7 
	if m == 0:
		Khb = Kh

	if m == 0:
		m = ct/20.  

	Kvb = Khb 
	#----- Flow Boundary conditions 
        r4 = (cw/2) - (w/2)
	h_riv = (ct - 1) # river head
	h_left = h_riv + (r4 * (dh/100.)) #  
	h_right = h_left
	#----- Calc Xfar (need to be 1 in this operation) 
        calcXfar =  1.


	with open(sutra_inp_table,'w') as output:
		param_table = csv.writer(output, delimiter = ' ', lineterminator = '\n') #write table for each value of tested parameter 
		param_table.writerow((\
		float(ct), float(cw), float(w), float(d),
		float(a), float(m), float(s1), float(s2),
		float(anis), float(Kh), float(Khb), float(Kvb),
		float(h_left),float(h_riv),float(h_right), float(calcXfar) ))
	with open(sutra_inp_table) as csv_data:
		reader = csv.reader(csv_data)
		rows = [row for row in reader if row]
		for row in rows:
			row = row
			lis=[x.split() for x in row]
			with open(sutra_inp_file,'w') as output:
				for x in zip(*lis):
					for y in x:
						output.write(y+'\n')
			xsection_rivaq.gen_mesh()
			xsection_rivaq.run() 

	data = np.genfromtxt(calc_xfar, dtype=[('d',float),('comp_x',float)], delimiter = ' ')
	d = data['d']
	comp_x = data['comp_x']
	x = (d)*(cw/2)
	y = comp_x
	ylabel(' x component [%]', fontsize =17 )
	xlabel(' Distance from river\'s center [m]', fontsize =17)
	grid(True)
	plt.plot(x, y, 'bo')
	savefig(plot)
	plt.close()






