/* Tarang-2
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-2 .
 *
 * Tarang-2 is free software; you can redistribute it and/or
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
 */

/*! \file SSF_pencil.cc 
 * 
 * @sa SSF_pencil.h
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "SSF_pencil.h"
#include "BasicIO_SP.h"


SSF_PENCIL::SSF_PENCIL()
{
			// kfactor, xfactor,  L, Delta_x
	if (global.field.kfactor.size() == 4) {
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		global.field.L[1] = M_PI/global.field.kfactor[1];
		global.field.L[2] = M_PI/global.field.kfactor[2];
		global.field.L[3] = 2*M_PI/global.field.kfactor[3];
	}
	
	else if (global.field.L.size() == 4) {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		global.field.kfactor[1] = M_PI/global.field.L[1];
		global.field.kfactor[2] = M_PI/global.field.L[2];
		global.field.kfactor[3] = 2*M_PI/global.field.L[3];
	}
	
	else  {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		
		if (global.program.kind == "RBC") {
			global.field.kfactor[1] =  M_PI;
			global.field.kfactor[2] =  M_PI/sqrt(2.0);
			global.field.kfactor[3] =  M_PI/sqrt(2.0);
		}
		else {
			global.field.kfactor[1] = 1.0;
			global.field.kfactor[2] = 1.0;
			global.field.kfactor[3] = 1.0;
		}
		
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		global.field.L[1] = M_PI/global.field.kfactor[1];
		global.field.L[2] = M_PI/global.field.kfactor[2];
		global.field.L[3] = 2*M_PI/global.field.kfactor[3];
	}
	
	for (int i=1; i<=3; i++) {
		global.field.Delta_x[i] =  global.field.L[i]/global.field.N[i];
		global.field.xfactor[i] = 1/global.field.kfactor[i];
	}
	
	
		// array size proper or not?	
	bool Nxproperdiv, Nyproperdiv, Nzproperdiv; 
	
	Nxproperdiv = ((global.field.Nx%global.mpi.num_p_hor == 0) && (global.field.Nx%global.mpi.num_p_vert == 0));
	Nyproperdiv = (global.field.Ny%global.mpi.num_p_vert == 0);
	Nzproperdiv = ((global.field.Nz/2+1)%global.mpi.num_p_hor == 0);
	
	if (!(Nxproperdiv && Nyproperdiv && Nzproperdiv)) {
		if (global.mpi.master) cerr << "ERROR: Array not being divided equally.  Check dimensions" << endl;
		exit(1);
	}
	
	
	//global.mpi.num_p_hor is assigned in global via mpirun argument
	spectralTransform.Init("SSF", "PENCIL", Nx, Ny, Nz, global.mpi.num_p_hor);

	global.mpi.num_p_vert = spectralTransform.num_p_vert;
	global.mpi.my_hor_pcoord = spectralTransform.my_hor_pcoord;
	global.mpi.my_vert_pcoord = spectralTransform.my_vert_pcoord;
	
	
	global.field.shape_complex_array = spectralTransform.local_Ny_hor,Nz/2+1,spectralTransform.local_Nx_vert;
	global.field.shape_real_array = spectralTransform.local_Ny_vert,2*spectralTransform.local_Nz_hor,Nx;

	//************
	
	
	global.mpi.num_x_procs = global.mpi.num_p_vert;
	global.mpi.num_y_procs = global.mpi.num_p_hor;
	global.mpi.num_z_procs = 1;
	
	global.mpi.my_x_pcoord = global.mpi.my_vert_pcoord;
	global.mpi.my_y_pcoord = global.mpi.my_hor_pcoord;
	global.mpi.my_z_pcoord = 0;
	
	global.mpi.num_x_procs_real = 1;
	global.mpi.num_y_procs_real = global.mpi.num_p_vert;
	global.mpi.num_z_procs_real = global.mpi.num_p_hor;
	
	global.mpi.my_x_pcoord_real = 0;
	global.mpi.my_y_pcoord_real = global.mpi.my_vert_pcoord;
	global.mpi.my_z_pcoord_real = global.mpi.my_hor_pcoord;
	
	
    
	global.field.maxlx = spectralTransform.local_Nx_vert-1;
	global.field.maxly = spectralTransform.local_Ny_hor-1;
	global.field.maxlz = Nz/2;
	
	
    // alias..
	kfactor.resize(4);
	kfactor[0]= global.field.kfactor[0];
	kfactor[1]= global.field.kfactor[1];
	kfactor[2]= global.field.kfactor[2];
	kfactor[3]= global.field.kfactor[3];
	
	local_Nx_vert=spectralTransform.local_Nx_vert;
	local_Ny_vert=spectralTransform.local_Ny_vert;
	local_Nz_vert=spectralTransform.local_Nz_vert;
	
	local_Nx_hor=spectralTransform.local_Nx_hor;
	local_Ny_hor=spectralTransform.local_Ny_hor;
	local_Nz_hor=spectralTransform.local_Nz_hor;
	
	local_Nx_start=spectralTransform.local_Nx_start;
	local_Ny_start=spectralTransform.local_Ny_start;
	local_Nz_start=spectralTransform.local_Nz_start;
	
	shape_complex_array = global.field.shape_complex_array;
    shape_real_array = global.field.shape_real_array;

	num_x_procs = global.mpi.num_x_procs;
	num_y_procs = global.mpi.num_y_procs;
	num_z_procs = global.mpi.num_z_procs;

	my_x_pcoord = global.mpi.my_x_pcoord;
	my_y_pcoord = global.mpi.my_y_pcoord;
	my_z_pcoord = global.mpi.my_z_pcoord;

	num_x_procs_real = global.mpi.num_x_procs_real;
	num_y_procs_real = global.mpi.num_y_procs_real;
	num_z_procs_real = global.mpi.num_z_procs_real;

	my_x_pcoord_real = global.mpi.my_x_pcoord_real;
	my_y_pcoord_real = global.mpi.my_y_pcoord_real;
	my_z_pcoord_real = global.mpi.my_z_pcoord_real;

	local_Nx=local_Nx_vert;
	local_Ny=local_Ny_hor;
	local_Nz=Nz/2+1;
	
	local_Nx_real=Nx;
	local_Ny_real=local_Ny_vert;
	local_Nz_real=2*local_Nz_hor;
	
	global.temp_array.X.resize(shape_complex_array);
    global.temp_array.X2.resize(shape_complex_array);
	
	global.temp_array.Xr.resize(shape_real_array);
    global.temp_array.Xr2.resize(shape_real_array);

	
	BasicIO::Array_properties<3> array_properties;
	array_properties.shape_full_complex_array = Ny, Nz/2+1, Nx;
	array_properties.shape_full_real_array = Ny, Nz+2, Nx;

	array_properties.id_complex_array = my_y_pcoord, my_z_pcoord, my_x_pcoord;
	array_properties.id_real_array = my_y_pcoord_real, my_z_pcoord_real, my_x_pcoord_real;

	array_properties.numprocs_complex_array = num_y_procs, num_z_procs, num_x_procs;
	array_properties.numprocs_real_array = num_y_procs_real, num_z_procs_real, num_x_procs_real;

	if (global.io.N_in_reduced.size() == 3)
		array_properties.shape_N_in_reduced = global.io.N_in_reduced[1], global.io.N_in_reduced[2]/2+1, global.io.N_in_reduced[0];
	
	if (global.io.N_out_reduced.size() == 3)
		array_properties.shape_N_out_reduced = global.io.N_out_reduced[1], global.io.N_out_reduced[2]/2+1, global.io.N_out_reduced[0];

	array_properties.Fourier_directions = 0,1,0;
	array_properties.Z = 1;

	array_properties.datatype_complex_space = BasicIO::H5T_COMPLX;
	array_properties.datatype_real_space = BasicIO::H5T_DP;

	BasicIO::Set_H5_plans(array_properties, this);
}
