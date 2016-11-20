"""
 *Tarang Post processor
 *
 * Copyright (C) 2014  Mahendra K. Verma (mkv@iitk.ac.in)
 *
 * @Author: Vijay Jain, Abhishek Kumar
 *
 * 
 *
 * Tarang Post processor is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-2; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 *
"""

from lib.fields import Csf
from lib.fields import Cvf
from lib.fluid.fluid_base import FluidSF
from lib.fluid.fluid_base import FluidVF
from lib.fluid.fluid_base import Correlation
from lib.fluid.incompressible.EnergyTr import EnergyTr
import yaml


import numpy
import os.path
import sys
import constants
import helperfunc

d = None
para = None
U = None
W = None
T = None
N = None
K = None
correlation = None
energytr = None
time = None
Rlambda = None
time_dt = None
no_spheres = None

pi = numpy.pi
def init(Path):
	global d
	global para
	global U
	global B
	global T
	global N,K, no_spheres
	global correlation
	global energytr
	global time
	global Rlambda
	global time_dt
	path  = Path
	constants.directory = path
	if not os.path.isfile(path + '/in/para.yaml'):
		print "Please enter correct absolute path of directory"
		sys.exit()
	else:
		stream = open(path + '/in/para.yaml', 'r')
		d = yaml.load(stream)		# it will convert yaml into python dictionary
		helperfunc.read_glob()
		para = Struct(d)	# para will be object containing all information of para.yaml file

		N = para.field.N
		K = para.field.kfactor
		#K = [pi,2*pi,2*pi]
		no_spheres = para.energy_transfer.flux.no_spheres

		if(para.program.kind == "FLUID_INCOMPRESS"):
			U = FluidVF.FluidVF("U", para.program.kind)
			correlation = Correlation.Correlation(para.program.kind)
			energytr = EnergyTr.EnergyTr()
			time = constants.data[:,0]
			Rlambda = constants.data[:,10]
			time_dt = constants.data[:,11]
			helperfunc.read_spectrum_time_fluid()
			helperfunc.read_flux_time_fluid()
		elif(para.program.kind == "RBC"):
			U = FluidVF.FluidVF("U", para.program.kind)
			T = FluidSF.FluidSF("T", para.program.kind)
			correlation = Correlation.Correlation(para.program.kind)
			energytr = EnergyTr.EnergyTr()
			time = constants.data[:,0]
			Rlambda = constants.data[:,16]
			time_dt = constants.data[:,17]
			if para.field.N[1] == 1:
				helperfunc.read_spectrum_time_RBC_2d()
			else:
				helperfunc.read_spectrum_time_RBC()
			helperfunc.read_flux_time_RBC()
		elif(para.program.kind == "MHD_INCOMPRESS"):
			U = FluidVF.FluidVF("U", para.program.kind)
			B = FluidVF.FluidVF("B", para.program.kind)
			correlation = Correlation.Correlation(para.program.kind)
			energytr = EnergyTr.EnergyTr()
			time = constants.data[:,0]
			Rlambda = constants.data[:,21]
			time_dt = constants.data[:,22]
			helperfunc.read_spectrum_time_MHD()
			helperfunc.read_flux_time_MHD()
			

# helperfunc.read_dir()


#code to convert nested dictionary to object
class Struct(object):
	def __init__(self, entries):
		self.__dict__.update(entries)
		for k, v in entries.items():
			if isinstance(v, dict):
				self.__dict__[k] = Struct(v)


################################## Spectrum ##################################
k = None
# displays list of time read from spectrum.d
#if para.program.kind == "RBC":
	#if 
def spectrum_showtime():
	for i in range(len(constants.spectrum_timelist)):
		print str(i+1) + ". " + constants.spectrum_timelist[i]



def spectrum_read(*arg):
	global k
	if para.program.kind == "FLUID_INCOMPRESS":
		if para.field.N[1] == 1:
			for lineno in arg:
				helperfunc.read_spectrum_data_fluid_2d(U, correlation, lineno-1)
			k = constants.K
		else:
			for lineno in arg:
				helperfunc.read_spectrum_data_fluid(U, correlation, lineno-1)
			k = constants.K
	elif para.program.kind == "RBC":
		if para.field.N[1] == 1:
			for lineno in arg:
				helperfunc.read_spectrum_data_RBC_2d(U, T, correlation, lineno-1)
			k = constants.K
		else:
			for lineno in arg:
				helperfunc.read_spectrum_data_RBC(U, T, correlation, lineno-1)
			k = constants.K
	elif para.program.kind == "MHD_INCOMPRESS":
		if para.field.N[1] == 1:
			for lineno in arg:
				helperfunc.read_spectrum_data_MHD_2d(U, B, correlation, lineno-1)
			k = constants.K
		else:
			for lineno in arg:
				helperfunc.read_spectrum_data_MHD(U, B, correlation, lineno-1)
			k = constants.K
		

################################## Spectrum ##################################

################################## functions for flux.d ######################################

# displays list of time read from spectrum.d
kf = None
def flux_showtime():
	for i in range(len(constants.flux_timelist)):
		print str(i+1) + ". " + constants.flux_timelist[i]

def flux_read(*arg):
	global kf
	if para.program.kind == "FLUID_INCOMPRESS":
		for lineno in arg:
			helperfunc.read_flux_data_fluid(energytr, lineno-1)
		kf = constants.KF
	elif para.program.kind == "RBC":
		for lineno in arg:
			helperfunc.read_flux_data_RBC(energytr, lineno-1)
		kf = constants.KF
	elif para.program.kind == "MHD_INCOMPRESS":
		for lineno in arg:
			helperfunc.read_flux_data_MHD(energytr, lineno-1)
		kf = constants.KF

################################## functions for flux.d ######################################

################################## functions for field_k_out.d ###############################

def field_k_read(no_files):
	if para.program.kind == "FLUID_INCOMPRESS":
		helperfunc.read_field_k_out_Fluid(U, no_files)
	elif para.program.kind == "RBC":
		helperfunc.read_field_k_out_RBC(U, no_files)
	elif para.program.kind == "MHD_INCOMPRESS":
		helperfunc.read_field_k_out_MHD(U, no_files)

def field_k_showmode():
	for mode in U.cvf.Modek:
		print "%s\t%s\t%s" %(int(mode.mode1), int(mode.mode2),int(mode.mode3))

def field_k_readmode(mode1, mode2, mode3):
	for mode in U.cvf.Modek:
		if mode.mode1 == mode1 and mode.mode2 == mode2 and mode.mode3 == mode3:
			return mode
################################## functions for field_k_out.d ###############################

################################## functions for field_r_out.d ###############################

def field_r_read(no_files):
	if para.program.kind == "FLUID_INCOMPRESS":
		helperfunc.read_field_r_out_Fluid(U, no_files)
	elif para.program.kind == "RBC":
		helperfunc.read_field_r_out_RBC(U, no_files)
	elif para.program.kind == "MHD_INCOMPRESS":
		helperfunc.read_field_r_out_MHD(U, no_files)

def field_r_showmode():
	for mode in U.cvf.Moder:
		print "%s\t%s\t%s" %(int(mode.mode1), int(mode.mode2), int(mode.mode3))

def field_r_readmode(mode1, mode2, mode3):
	for mode in U.cvf.Moder:
		if mode.mode1 == mode1 and mode.mode2 == mode2 and mode.mode3 == mode3:
			return mode
################################## functions for field_r_out.d ###############################

################################## functions for real_field.d ################################
def real_cutter():
	helperfunc.cutter()

def profile():
	if para.program.kind == "FLUID_INCOMPRESS":
		return helperfunc.profile_Fluid(N[0], N[1], N[2])
	elif para.program.kind == "RBC":
		return helperfunc.profile_RBC(N[0], N[1], N[2])
	elif para.program.kind == "MHD_INCOMPRESS":
		return helperfunc.profile_MHD(N[0], N[1], N[2])
		
def visual():
	if para.program.kind == "FLUID_INCOMPRESS":
		return helperfunc.visual_Fluid(N[0], N[1], N[2])
	elif para.program.kind == "RBC":
		return helperfunc.visual_RBC(N[0], N[1], N[2])
	elif para.program.kind == "MHD_INCOMPRESS":
		return helperfunc.visual_MHD(N[0], N[1], N[2])
	
################################## functions for real_field.d ################################

################################## shell_to_shell ##################################
def shell_showtime():
	if para.program.kind == "FLUID_INCOMPRESS":
		helperfunc.read_shell_time_Fluid()
		for i in range(len(constants.shell_to_shell_timelist)):
			print str(i+1) + ". " + constants.shell_to_shell_timelist[i]
	elif para.program.kind == "RBC":
		helperfunc.read_shell_time_RBC()
		for i in range(len(constants.shell_to_shell_timelist)):
			print str(i+1) + ". " + constants.shell_to_shell_timelist[i]
	elif para.program.kind == "MHD_INCOMPRESS":
		helperfunc.read_shell_time_MHD()
		for i in range(len(constants.shell_to_shell_timelist)):
			print str(i+1) + ". " + constants.shell_to_shell_timelist[i]
	

def shell_to_shell_read(time,mode):
	if para.program.kind == "FLUID_INCOMPRESS":
		return helperfunc.read_shell_data_Fluid(time, mode)
	elif para.program.kind == "RBC":
		return helperfunc.read_shell_data_RBC(time, mode)
	elif para.program.kind == "MHD_INCOMPRESS":
		return helperfunc.read_shell_data_MHD(time, mode)
			

################################## shell_to_shell ##################################

################################## ring_spectrum ##################################
def ring_spectrum_showtime():
	if para.program.kind == "FLUID_INCOMPRESS":
		helperfunc.read_ring_spectrum_time_Fluid()
		for i in range(len(constants.ring_spectrum_timelist)):
			print str(i+1) + ". " + constants.ring_spectrum_timelist[i]
	elif para.program.kind == "RBC":
		helperfunc.read_ring_spectrum_time_RBC()
		for i in range(len(constants.ring_spectrum_timelist)):
			print str(i+1) + ". " + constants.ring_spectrum_timelist[i]
	elif para.program.kind == "MHD_INCOMPRESS":
		helperfunc.read_ring_spectrum_time_MHD()
		for i in range(len(constants.ring_spectrum_timelist)):
			print str(i+1) + ". " + constants.ring_spectrum_timelist[i]
	

def ring_spectrum_read(time,mode):
	if para.program.kind == "FLUID_INCOMPRESS":
		return helperfunc.read_ring_spectrum_data_Fluid(time, mode)		
	elif para.program.kind == "RBC":
		return helperfunc.read_ring_spectrum_data_RBC(time, mode)		
	elif para.program.kind == "MHD_INCOMPRESS":
		return helperfunc.read_ring_spectrum_data_MHD(time, mode)		
	

################################## ring_spectrum ##################################

################################## Global functions and classes ###################



def Mag(real,imag):
	return ((real**2 + imag**2)**0.5)

def Energy_vector(A,B,C):
	return (A**2 + B**2 +C**2)/2.0
	
def Energy_scalar(A):
	return (A**2)/2.0
	

class FFF_SLAB:
	def Max_radius_inside(self, N):
		if N[1] > 1:
			Kmag = min( (N[0]/2) * K[0] , (N[1]/2) * K[1])
			Kmag = min(Kmag, (N[2]/2) * K[2])
		elif N[1] == 1:
			Kmag = min( (N[0]/2) * K[0] , (N[2]/2) * K[2])
		return int(Kmag)
		

class SFF_SLAB:
	def Max_radius_inside(self, N):
		if N[1]>1:
			Kmag = min( (N[0]-1) * K[0] , (N[1]/2) * K[1])
			Kmag = min(Kmag, (N[2]/2) * K[2])
		elif N[1] == 1:
			Kmag = min( (N[0]-1) * K[0] , (N[2]/2) * K[2])
		return int(Kmag)

def flux_radii(Max_radius_inside):
	radii = numpy.zeros(no_spheres)
	radii[0] = 0.0
	radii[1] = 2.0
	radii[2] = 4.0
	radii[3] = 8.0
	radii[no_spheres - 2] = Max_radius_inside/2.0
	radii[no_spheres - 1] = Max_radius_inside
	
	if no_spheres > 6:
		s = numpy.log2(Max_radius_inside/16.0) / (no_spheres - 5)
		for i in numpy.arange(4,(no_spheres - 3) + 1, 1):
			radii[i] = 8 * 2**(s*(i-3))
	
	return radii

################################## Global functions and classes ###################
