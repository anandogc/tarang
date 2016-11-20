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
from __future__ import division
import numpy
import linecache
import constants
import field_k
import field_r
import Global_vars
import tarang
import h5py

	

def read_dir():
	constants.directory = raw_input("Please give absolute path of sample dictionary : ")
	print 'Directory : ' + constants.directory

#################################### functions for glob.d ###################################################

def read_glob():
	constants.data = numpy.loadtxt(constants.directory + '/out/glob.d')
	

#################################### functions for glob.d ###################################################

#################################### functions for spectrum.d ###############################################
def read_spectrum_time_fluid():
	del constants.spectrum_timelist[:]
	path = constants.directory + "/out/spectrum.d"
	lineno = 0
	line = linecache.getline(path, 14*lineno+1)
	while (line) != '':
		constants.spectrum_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 14*lineno+1)		
		
def read_spectrum_time_RBC_2d():
	del constants.spectrum_timelist[:]
	path = constants.directory + "/out/spectrum.d"
	lineno = 0
	line = linecache.getline(path, 18*lineno+1)
	while (line) != '':
		constants.spectrum_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 18*lineno+1)

def read_spectrum_time_RBC():
	del constants.spectrum_timelist[:]
	path = constants.directory + "/out/spectrum.d"
	lineno = 0
	line = linecache.getline(path, 22*lineno+1)
	while (line) != '':
		constants.spectrum_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 22*lineno+1)

def read_spectrum_time_MHD():
	del constants.spectrum_timelist[:]
	path = constants.directory + "/out/spectrum.d"
	lineno = 0
	line = linecache.getline(path, 30*lineno+1)
	while (line) != '':
		constants.spectrum_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 30*lineno+1)		

def read_spectrum_data_fluid(U, correlation, lineno):
	path = constants.directory + "/out/spectrum.d"
	U.cvf.Uek_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+4).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+5).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+6).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk1.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+8).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk2.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+9).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk3.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+10).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+12).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+13).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 14*lineno+14).split('\t')[:-1]]))
	constants.K = numpy.arange(0,len(U.cvf.Uek_shell_ek1[0]),1) # set value of k
	
	
def read_spectrum_data_RBC_2d(U, T, correlation, lineno):
	path = constants.directory + "/out/spectrum.d"
	U.cvf.Uek_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+4).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+5).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+7).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+8).split('\t')[:-1]]))
	T.csf.Tek_shell_ek.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+10).split('\t')[:-1]]))
	correlation.U_T_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+12).split('\t')[:-1]]))
	correlation.U_T_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+13).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+15).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+16).split('\t')[:-1]]))
	T.csf.FT_T_shell_ek.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+18).split('\t')[:-1]]))
	constants.K = numpy.arange(0,len(U.cvf.Uek_shell_ek1[0]),1) # set value of k	
	
def read_spectrum_data_RBC(U, T, correlation, lineno):
	path = constants.directory + "/out/spectrum.d"
	U.cvf.Uek_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+4).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+5).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+6).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk1.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+8).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk2.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+9).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk3.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+10).split('\t')[:-1]]))
	T.csf.Tek_shell_ek.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+12).split('\t')[:-1]]))
	correlation.U_T_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+14).split('\t')[:-1]]))
	correlation.U_T_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+15).split('\t')[:-1]]))
	correlation.U_T_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+16).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+18).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+19).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+20).split('\t')[:-1]]))
	T.csf.FT_T_shell_ek.append(numpy.array([float(num) for num in linecache.getline(path, 22*lineno+22).split('\t')[:-1]]))
	constants.K = numpy.arange(0,len(U.cvf.Uek_shell_ek1[0]),1) # set value of k

def read_spectrum_data_MHD_2d(U, B, correlation, lineno):
	path = constants.directory + "/out/spectrum.d"
	U.cvf.Uek_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+4).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+5).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+7).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+8).split('\t')[:-1]]))
	T.csf.Tek_shell_ek.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+10).split('\t')[:-1]]))
	correlation.U_T_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+12).split('\t')[:-1]]))
	correlation.U_T_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+13).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+15).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+16).split('\t')[:-1]]))
	T.csf.FT_T_shell_ek.append(numpy.array([float(num) for num in linecache.getline(path, 18*lineno+18).split('\t')[:-1]]))
	constants.K = numpy.arange(0,len(U.cvf.Uek_shell_ek1[0]),1) # set value of k	
	
def read_spectrum_data_MHD(U, B, correlation, lineno):
	path = constants.directory + "/out/spectrum.d"
	U.cvf.Uek_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+4).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+5).split('\t')[:-1]]))
	U.cvf.Uek_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+6).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk1.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+8).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk2.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+9).split('\t')[:-1]]))
	U.cvf.UDk_shell_dissk3.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+10).split('\t')[:-1]]))
	
	B.cvf.Uek_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+12).split('\t')[:-1]]))
	B.cvf.Uek_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+13).split('\t')[:-1]]))
	B.cvf.Uek_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+14).split('\t')[:-1]]))
	B.cvf.UDk_shell_dissk1.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+16).split('\t')[:-1]]))
	B.cvf.UDk_shell_dissk2.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+17).split('\t')[:-1]]))
	B.cvf.UDk_shell_dissk3.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+18).split('\t')[:-1]]))
	
	
	correlation.U_B_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+20).split('\t')[:-1]]))
	correlation.U_B_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+21).split('\t')[:-1]]))
	correlation.U_B_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+22).split('\t')[:-1]]))
	
	U.cvf.Fv_v_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+24).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+25).split('\t')[:-1]]))
	U.cvf.Fv_v_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+26).split('\t')[:-1]]))
	
	B.cvf.Fv_v_shell_ek1.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+28).split('\t')[:-1]]))
	B.cvf.Fv_v_shell_ek2.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+29).split('\t')[:-1]]))
	B.cvf.Fv_v_shell_ek3.append(numpy.array([float(num) for num in linecache.getline(path, 30*lineno+30).split('\t')[:-1]]))
	
	constants.K = numpy.arange(0,len(U.cvf.Uek_shell_ek1[0]),1) # set value of k	

#################################### functions for spectrum.d ###############################################

#################################### functions for flux.d ###################################################
def read_flux_time_fluid():
	del constants.flux_timelist[:]
	path = constants.directory + "/out/flux.d"
	lineno = 0
	line = linecache.getline(path, 5*lineno+1)
	while (line) != '':
		constants.flux_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 5*lineno+1)


def read_flux_time_RBC():
	del constants.flux_timelist[:]
	path = constants.directory + "/out/flux.d"
	lineno = 0
	line = linecache.getline(path, 9*lineno+1)
	while (line) != '':
		constants.flux_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 9*lineno+1)

def read_flux_time_MHD():
	del constants.flux_timelist[:]
	path = constants.directory + "/out/flux.d"
	lineno = 0
	line = linecache.getline(path, 21*lineno+1)
	while (line) != '':
		constants.flux_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 21*lineno+1)
		
def read_flux_data_fluid(energytr, lineno):
	path = constants.directory + "/out/flux.d"
	energytr.U2U.append(numpy.array([float(num) for num in linecache.getline(path, 5*lineno+3).split('\t')[:-1]]))
	energytr.Fv_v.append(numpy.array([float(num) for num in linecache.getline(path, 5*lineno+5).split('\t')[:-1]]))

	
	
def read_flux_data_RBC(energytr, lineno):
	path = constants.directory + "/out/flux.d"
	energytr.U2U.append(numpy.array([float(num) for num in linecache.getline(path, 9*lineno+3).split('\t')[:-1]]))
	energytr.T2T.append(numpy.array([float(num) for num in linecache.getline(path, 9*lineno+5).split('\t')[:-1]]))
	energytr.Fv_v.append(numpy.array([float(num) for num in linecache.getline(path, 9*lineno+7).split('\t')[:-1]]))
	energytr.FT_T.append(numpy.array([float(num) for num in linecache.getline(path, 9*lineno+9).split('\t')[:-1]]))

def read_flux_data_MHD(energytr, lineno):
	path = constants.directory + "/out/flux.d"
	energytr.U2U.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+3).split('\t')[:-1]]))
	energytr.flux_VF_Uin_Wout.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+5).split('\t')[:-1]]))
	energytr.flux_VF_Uin_Win.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+7).split('\t')[:-1]]))
	energytr.flux_VF_Win_Wout.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+9).split('\t')[:-1]]))
	energytr.flux_VF_Win_Uout.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+11).split('\t')[:-1]]))
	energytr.flux_VF_Uout_Wout.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+13).split('\t')[:-1]]))
	energytr.flux_Elsasser_plus.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+15).split('\t')[:-1]]))
	energytr.flux_Elsasser_minus.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+17).split('\t')[:-1]]))
	energytr.Fv_v.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+19).split('\t')[:-1]]))
	energytr.FB_B.append(numpy.array([float(num) for num in linecache.getline(path, 21*lineno+21).split('\t')[:-1]]))


#################################### functions for flux.d ###################################################

#################################### functions for field_k_out.d ############################################

def read_field_k_out_Fluid(U, no_files):
	for i in range(no_files):
		path = constants.directory + "/out/field_k_out_" + str(i) + ".d"
		print path
		d = numpy.loadtxt(path)
		for i in range(len(d)):
			found = False
			for mode in U.cvf.Modek:
				if(mode.check_mode(d[i][1], d[i][2], d[i][3])):
					mode.insert_Fluid(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][0])
					found = True
					break

			if not found:
				mode = field_k.Field_k(d[i][1], d[i][2], d[i][3])			
				mode.insert_Fluid(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][0])
				U.cvf.Modek.append(mode)
				
				
def read_field_k_out_RBC(U, no_files):
	for i in range(no_files):
		path = constants.directory + "/out/field_k_out_" + str(i) + ".d"
		print path
		d = numpy.loadtxt(path)
		for i in range(len(d)):
			found = False
			for mode in U.cvf.Modek:
				if(mode.check_mode(d[i][1], d[i][2], d[i][3])):
					mode.insert_RBC(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][10], d[i][11], d[i][0])
					found = True
					break

			if not found:
				mode = field_k.Field_k(d[i][1], d[i][2], d[i][3])			
				mode.insert_RBC(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][10], d[i][11], d[i][0])
				U.cvf.Modek.append(mode)

def read_field_k_out_MHD(U, no_files):
	for i in range(no_files):
		path = constants.directory + "/out/field_k_out_" + str(i) + ".d"
		print path
		d = numpy.loadtxt(path)
		for i in range(len(d)):
			found = False
			for mode in U.cvf.Modek:
				if(mode.check_mode(d[i][1], d[i][2], d[i][3])):
					mode.insert_MHD(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][10], d[i][11], d[i][12], d[i][13], d[i][14], d[i][15], d[i][0])
					found = True
					break

			if not found:
				mode = field_k.Field_k(d[i][1], d[i][2], d[i][3])			
				mode.insert_MHD(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][10], d[i][11], d[i][12], d[i][13], d[i][14], d[i][15], d[i][0])
				U.cvf.Modek.append(mode)
#################################### functions for field_k_out.d ############################################

#################################### functions for field_r_out.d ############################################

def read_field_r_out_Fluid(U, no_files):
	for i in range(no_files):
		path = constants.directory + "/out/field_r_out_" + str(i) + ".d"
		print path
		d = numpy.loadtxt(path)	
		for i in range(len(d)):
			found = False
			for mode in U.cvf.Moder:
				if(mode.check_mode(d[i][1], d[i][2], d[i][3])):
					mode.insert_Fluid(d[i][4], d[i][5], d[i][6], d[i][0])
					found = True
					break

			if not found:
				mode = field_r.Field_r(d[i][1], d[i][2], d[i][3])			
				mode.insert_Fluid(d[i][4], d[i][5], d[i][6], d[i][0])
				U.cvf.Moder.append(mode)
				
def read_field_r_out_RBC(U, no_files):
	for i in range(no_files):
		path = constants.directory + "/out/field_r_out_" + str(i) + ".d"
		print path
		d = numpy.loadtxt(path)	
		for i in range(len(d)):
			found = False
			for mode in U.cvf.Moder:
				if(mode.check_mode(d[i][1], d[i][2], d[i][3])):
					mode.insert_RBC(d[i][4], d[i][5], d[i][6], d[i][7], d[i][0])
					found = True
					break

			if not found:
				mode = field_r.Field_r(d[i][1], d[i][2], d[i][3])			
				mode.insert_RBC(d[i][4], d[i][5], d[i][6], d[i][7], d[i][0])
				U.cvf.Moder.append(mode)
				
def read_field_r_out_MHD(U, no_files):
	for i in range(no_files):
		path = constants.directory + "/out/field_r_out_" + str(i) + ".d"
		print path
		d = numpy.loadtxt(path)	
		for i in range(len(d)):
			found = False
			for mode in U.cvf.Moder:
				if(mode.check_mode(d[i][1], d[i][2], d[i][3])):
					mode.insert_MHD(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][0])
					found = True
					break

			if not found:
				mode = field_r.Field_r(d[i][1], d[i][2], d[i][3])			
				mode.insert_MHD(d[i][4], d[i][5], d[i][6], d[i][7], d[i][8], d[i][9], d[i][0])
				U.cvf.Moder.append(mode)

#################################### functions for field_r_out.d ############################################


#################################### functions for read_field.d #############################################
def cutter():
	time = str(input('Please Enter time : '))
	path = constants.directory + "/out/real_" + time + "/"
	start_dim = int(input('Please Enter Starting diamension : '))
	end_dim = int(input("Please Enter ending diamension : "))
	read_file = h5py.File(path + "T.Fr.h5")
	dataset = read_file['/T.Fr']
	write_file = h5py.File(path + "T.Fr_new.h5", "w")
	write_file['/T.Fr'] = dataset[:, start_dim:end_dim, :]
	read_file.close()
	write_file.close()
	print 'File is created successfully'

def profile_Fluid(Nx, Ny, Nz):
	data = []
	lx = 1
	x_increment = lx/(Nx-1)
	time = str(input('Please Enter time : '))	
	path = constants.directory + "/out/real_" + time + "/"
	
	file_V1_read = h5py.File(path + "U.V1r.h5")
	dataset_V1_read = file_V1_read['/U.V1r']
	file_V2_read = h5py.File(path + "U.V2r.h5")
	dataset_V2_read = file_V2_read['/U.V2r']
	file_V3_read = h5py.File(path + "U.V3r.h5")
	dataset_V3_read = file_V3_read['/U.V3r']
	
	for i in numpy.arange(0,Nx): 
		x = x_increment*i
		
		v1_sqr =  (numpy.sum(numpy.sum(dataset_V1_read[i,:,:]**2,axis=0)))
		v1 =  numpy.sum(numpy.sum(dataset_V1_read[i,:,:],axis=0))
		v1_rms = numpy.sqrt(v1_sqr/(Ny*Nz))
		v1_avg = v1/(Ny*Nz)
		
		v2_sqr =  (numpy.sum(numpy.sum(dataset_V2_read[i,:,:]**2,axis=0)))
		v2 =  numpy.sum(numpy.sum(dataset_V2_read[i,:,:],axis=0))
		v2_rms = numpy.sqrt(v2_sqr/(Ny*Nz))
		v2_avg = v2/(Ny*Nz)
		
		v3_sqr =  (numpy.sum(numpy.sum(dataset_V3_read[i,:,:]**2,axis=0)))
		v3 =  numpy.sum(numpy.sum(dataset_V3_read[i,:,:],axis=0))
		v3_rms = numpy.sqrt(v3_sqr/(Ny*Nz))
		v3_avg = v3/(Ny*Nz)
		
		data.append([x,v1_avg, v2_avg, v3_avg, v1_rms, v2_rms, v3_rms])
	data = numpy.array(data)
	return data
	
def profile_RBC(Nx, Ny, Nz):
	data = []
	lx = 1
	x_increment = lx/(Nx-1)
	time = str(input('Please Enter time : '))	
	path = constants.directory + "/out/real_" + time + "/"
	
	file_T_read = h5py.File(path + "T.Fr.h5")
	dataset_T_read = file_T_read['/T.Fr']
	file_V1_read = h5py.File(path + "U.V1r.h5")
	dataset_V1_read = file_V1_read['/U.V1r']
	file_V2_read = h5py.File(path + "U.V2r.h5")
	dataset_V2_read = file_V2_read['/U.V2r']
	file_V3_read = h5py.File(path + "U.V3r.h5")
	dataset_V3_read = file_V3_read['/U.V3r']
	
	for i in numpy.arange(0,Nx): 
		x = x_increment*i
		
		v1_sqr =  (numpy.sum(numpy.sum(dataset_V1_read[i,:,:]**2,axis=0)))
		v1 =  numpy.sum(numpy.sum(dataset_V1_read[i,:,:],axis=0))
		v1_rms = numpy.sqrt(v1_sqr/(Ny*Nz))
		v1_avg = v1/(Ny*Nz)
		
		v2_sqr =  (numpy.sum(numpy.sum(dataset_V2_read[i,:,:]**2,axis=0)))
		v2 =  numpy.sum(numpy.sum(dataset_V2_read[i,:,:],axis=0))
		v2_rms = numpy.sqrt(v2_sqr/(Ny*Nz))
		v2_avg = v2/(Ny*Nz)
		
		v3_sqr =  (numpy.sum(numpy.sum(dataset_V3_read[i,:,:]**2,axis=0)))
		v3 =  numpy.sum(numpy.sum(dataset_V3_read[i,:,:],axis=0))
		v3_rms = numpy.sqrt(v3_sqr/(Ny*Nz))
		v3_avg = v3/(Ny*Nz)
		
		theta_sqr =  (numpy.sum(numpy.sum(dataset_T_read[i,:,:]**2,axis=0)))
		theta =  numpy.sum(numpy.sum(dataset_T_read[i,:,:],axis=0))
		theta_rms = numpy.sqrt(theta_sqr/(Ny*Nz))
		theta_avg = theta/(Ny*Nz)
		
		data.append([x,v1_avg, v2_avg, v3_avg, theta_avg, v1_rms, v2_rms, v3_rms, theta_rms])
	data = numpy.array(data)
	return data

def profile_MHD(Nx, Ny, Nz):
	data = []
	lx = 1
	x_increment = lx/(Nx-1)
	time = str(input('Please Enter time : '))	
	path = constants.directory + "/out/real_" + time + "/"
	
	file_T_read = h5py.File(path + "T.Fr.h5")
	dataset_T_read = file_T_read['/T.Fr']
	file_V1_read = h5py.File(path + "U.V1r.h5")
	dataset_V1_read = file_V1_read['/U.V1r']
	file_V2_read = h5py.File(path + "U.V2r.h5")
	dataset_V2_read = file_V2_read['/U.V2r']
	file_V3_read = h5py.File(path + "U.V3r.h5")
	dataset_V3_read = file_V3_read['/U.V3r']
	
	for i in numpy.arange(0,Nx): 
		x = x_increment*i
		
		v1_sqr =  (numpy.sum(numpy.sum(dataset_V1_read[i,:,:]**2,axis=0)))
		v1 =  numpy.sum(numpy.sum(dataset_V1_read[i,:,:],axis=0))
		v1_rms = numpy.sqrt(v1_sqr/(Ny*Nz))
		v1_avg = v1/(Ny*Nz)
		
		v2_sqr =  (numpy.sum(numpy.sum(dataset_V2_read[i,:,:]**2,axis=0)))
		v2 =  numpy.sum(numpy.sum(dataset_V2_read[i,:,:],axis=0))
		v2_rms = numpy.sqrt(v2_sqr/(Ny*Nz))
		v2_avg = v2/(Ny*Nz)
		
		v3_sqr =  (numpy.sum(numpy.sum(dataset_V3_read[i,:,:]**2,axis=0)))
		v3 =  numpy.sum(numpy.sum(dataset_V3_read[i,:,:],axis=0))
		v3_rms = numpy.sqrt(v3_sqr/(Ny*Nz))
		v3_avg = v3/(Ny*Nz)
		
		theta_sqr =  (numpy.sum(numpy.sum(dataset_T_read[i,:,:]**2,axis=0)))
		theta =  numpy.sum(numpy.sum(dataset_T_read[i,:,:],axis=0))
		theta_rms = numpy.sqrt(theta_sqr/(Ny*Nz))
		theta_avg = theta/(Ny*Nz)
		
		data.append([x,v1_avg, v2_avg, v3_avg, theta_avg, v1_rms, v2_rms, v3_rms, theta_rms])
	data = numpy.array(data)
	return data

def visual_RBC(Nx, Ny, Nz):
	
	time = str(input('Please Enter time : '))	
	path = constants.directory + "/out/real_" + time + "/"

	lx = 1
	ly = 1
	lz = 1
	
	x_increment = lx/(Nx-1)
	y_increment = ly/(Ny-1)
	z_increment = lz/(Nz-1)
	
	x = []
	y = []
	z = []
	
	for i in numpy.arange(0,Nx): 
		x.append(x_increment*i)
	for i in numpy.arange(0,Ny): 
		y.append(y_increment*i)
	for i in numpy.arange(0,Nz): 
		z.append(z_increment*i)
	X = numpy.zeros((Nx,Ny,Nz))
	Y = numpy.zeros((Nx,Ny,Nz))
	Z = numpy.zeros((Nx,Ny,Nz))
	
	
	
	for i in numpy.arange(0,Nx):
		X[i,:,:] = x[i]
	
	for i in numpy.arange(0,Ny):
		Y[:,i,:] = y[i]
	
	for i in numpy.arange(0,Nz):
		Z[:,:,i] = z[i]
	
	theta_data = h5py.File(path + "T.Fr.h5")
	theta = theta_data['/T.Fr']
	theta = numpy.array(theta[:,:,0:Nz])
	T = 1 -X +theta
	
	V1_data = h5py.File(path + "U.V1r.h5")
	V1 = V1_data['/U.V1r']
	V1 = numpy.array(V1[:,:,0:Nz])
	
	V2_data = h5py.File(path + "U.V2r.h5")
	V2 = V2_data['/U.V2r']
	V2 = numpy.array(V2[:,:,0:Nz])
	
	V3_data = h5py.File(path + "U.V3r.h5")
	V3 = V3_data['/U.V3r']
	V3 = numpy.array(V3[:,:,0:Nz])

	return (X, Y, Z, V1, V2, V3, theta, T)
#################################### functions for read_field.d #############################################

#################################### functions for shell_to_shell.d ###############################################
def read_shell_time_Fluid():
	del constants.shell_to_shell_timelist[:]
	path = constants.directory + "/out/shell_to_shell.d"
	lineno = 0
	line = linecache.getline(path, 21*lineno+1)
	while (line) != '':
		constants.shell_to_shell_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 21*lineno+1)

def read_shell_data_Fluid(time, mode):
	path = constants.directory + "/out/shell_to_shell.d"
	if(mode == "U2U"):
		for i in range(21*(time-1)+3, 21*(time-1)+22):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	else:
		print "Please enter U2U or T2T as second parameter"
	return constants.shell_data
	
			
def read_shell_time_RBC():
	del constants.shell_to_shell_timelist[:]
	path = constants.directory + "/out/shell_to_shell.d"
	lineno = 0
	line = linecache.getline(path, 41*lineno+1)
	while (line) != '':
		constants.shell_to_shell_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 41*lineno+1)

def read_shell_data_RBC(time, mode):
	path = constants.directory + "/out/shell_to_shell.d"
	if(mode == "U2U"):
		for i in range(41*(time-1)+3, 41*(time-1)+22):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "T2T"):
		for i in range(41*(time-1)+23, 41*(time-1)+42):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	else:
		print "Please enter U2U or T2T as second parameter"
	return constants.shell_data

def read_shell_time_MHD():
	del constants.shell_to_shell_timelist[:]
	path = constants.directory + "/out/shell_to_shell.d"
	lineno = 0
	line = linecache.getline(path, 103*lineno+1)
	while (line) != '':
		constants.shell_to_shell_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 103*lineno+1)

def read_shell_data_MHD(time, mode):
	path = constants.directory + "/out/shell_to_shell.d"
	if(mode == "U2U"):
		for i in range(103*(time-1)+3, 103*(time-1)+22):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "shelltoshell_VF_WtoW"):
		for i in range(103*(time-1)+23, 103*(time-1)+42):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "shelltoshell_VF_UtoW"):
		for i in range(103*(time-1)+43, 103*(time-1)+62):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "shelltoshell_Elsasser_plus"):
		for i in range(103*(time-1)+63, 103*(time-1)+82):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "shelltoshell_Elsasser_minus"):
		for i in range(103*(time-1)+83, 103*(time-1)+102):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	else:
		print "Please enter U2U or T2T as second parameter"
	return constants.shell_data
#################################### functions for shell_to_shell.d ###############################################

#################################### functions for ring_spectrum.d ###############################################
def read_ring_spectrum_time_Fluid():
	del constants.ring_spectrum_timelist[:]
	path = constants.directory + "/out/ring_spectrum.d"
	lineno = 0
	line = linecache.getline(path, 301*lineno+1)
	while (line) != '':
		constants.ring_spectrum_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 301*lineno+1)

def read_ring_spectrum_data_Fluid(time, mode):
	path = constants.directory + "/out/ring_spectrum.d"
	if(mode == "Uek"):
		for i in range(301*(time-1)+3, 301*(time-1)+102):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "UDk"):
		for i in range(301*(time-1)+103, 301*(time-1)+202):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "Fv_v"):
		for i in range(301*(time-1)+203, 301*(time-1)+302):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])			
	else:
		print "Please enter Uek or UDk or Fv_v as second parameter"
	return constants.shell_data

def read_ring_spectrum_time_RBC():
	del constants.ring_spectrum_timelist[:]
	path = constants.directory + "/out/ring_spectrum.d"
	lineno = 0
	line = linecache.getline(path, 2779*lineno+1)
	while (line) != '':
		constants.ring_spectrum_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 2779*lineno+1)

def read_ring_spectrum_data_RBC(time, mode):
	path = constants.directory + "/out/ring_spectrum.d"
	if(mode == "Uek"):
		for i in range(2779*(time-1)+3, 2779*(time-1)+597):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "UDk"):
		for i in range(2779*(time-1)+1193, 2779*(time-1)+1192):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "U.T"):
		for i in range(2779*(time-1)+1392, 2779*(time-1)+1391):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])	
	elif(mode == "Fv.v"):
		for i in range(2779*(time-1)+1987, 2779*(time-1)+1986):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])	
	elif(mode == "FT.T"):
		for i in range(2779*(time-1)+2582, 2779*(time-1)+2581):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])			
	else:
		print "Please enter Uek or UDk or Fv_v as second parameter"
	return constants.shell_data
	
def read_ring_spectrum_time_MHD():
	del constants.ring_spectrum_timelist[:]
	path = constants.directory + "/out/ring_spectrum.d"
	lineno = 0
	line = linecache.getline(path, 301*lineno+1)
	while (line) != '':
		constants.ring_spectrum_timelist.append(line[3:-1])
		lineno += 1
		line = linecache.getline(path, 301*lineno+1)

def read_ring_spectrum_data_MHD(time, mode):
	path = constants.directory + "/out/ring_spectrum.d"
	if(mode == "Uek"):
		for i in range(301*(time-1)+3, 301*(time-1)+102):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "UDk"):
		for i in range(301*(time-1)+103, 301*(time-1)+202):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])
	elif(mode == "Fv_v"):
		for i in range(301*(time-1)+203, 301*(time-1)+302):
			constants.shell_data.append([float(num) for num in linecache.getline(path, i).split('\t')[:-1]])			
	else:
		print "Please enter Uek or UDk or Fv_v as second parameter"
	return constants.shell_data
#################################### functions for ring_spectrum.d ###############################################
