/* Tarang-2
 *
 * Copyright_PENCIL( (C) 2008, 2009  Mahendra K. Verma
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

/*! \file FFF_slab.cc 
 * 
 * @sa FFF_slab.cc
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "FFF_slab.h"
#include "BasicIO_SP.h"
 

FFF_SLAB::FFF_SLAB()
{
	
	// kfactor, xfactor, L, Delta_x
	if (global.field.kfactor.size() == 4) {
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		for (int i=1; i<=3; i++) 
			global.field.L[i] = 2*M_PI/global.field.kfactor[i];
	}
	
	else if (global.field.L.size() == 4) {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		for (int i=1; i<=3; i++) 
			global.field.kfactor[i] = 2*M_PI/global.field.L[i];
	}
	
	else  {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		
		if (global.program.kind == "RBC") {
			global.field.kfactor[1] =  1.0;
			global.field.kfactor[2] =  1.0;
			global.field.kfactor[3] =  1.0;
		}
		else {
			global.field.kfactor[1] = 1.0;
			global.field.kfactor[2] = 1.0;
			global.field.kfactor[3] = 1.0;
		}
		
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		for (int i=1; i<=3; i++) 
			global.field.L[i] = 2*M_PI/global.field.kfactor[i];
	}
	
	for (int i=1; i<=3; i++) {
		global.field.Delta_x[i] =  global.field.L[i]/global.field.N[i];
		global.field.xfactor[i] = 1/global.field.kfactor[i];
	}
	
		// Proper allocation or not?
	
	bool Nxproperdiv, Nyproperdiv, Nzproperdiv; 
	
	if (global.field.Ny > 1) {
		Nxproperdiv = (global.field.Nx%global.mpi.numprocs == 0);
		Nyproperdiv = (global.field.Ny%global.mpi.numprocs == 0);
		Nzproperdiv = true;
	}
	
	else if (global.field.Ny == 1) {
		Nxproperdiv = (global.field.Nx%global.mpi.numprocs == 0);
		Nzproperdiv = ((global.field.Nz/2+1)%global.mpi.numprocs == 0);
		Nyproperdiv = true;
	}
	
	if (!(Nxproperdiv && Nyproperdiv && Nzproperdiv)) {
		cout << "N div " << global.field.Nx << " "<< global.field.Ny << " " << global.field.Nz << " " << Nxproperdiv << " " << Nyproperdiv << "  " << Nzproperdiv << endl;
		if (global.mpi.master) cerr << "ERROR in FFF_slab.cc: Array not being divided equally.  Check dimensions" << endl;
		exit(1);
	}

	spectralTransform.Init("FFF", Nx, Ny, Nz, 1);
	
	global.field.maxlx = spectralTransform.maxfx;
	global.field.maxly = spectralTransform.maxfy;
	global.field.maxlz = spectralTransform.maxfz;

	if (Ny > 1) {
		global.field.shape_complex_array = spectralTransform.FA_shape;
		global.field.shape_real_array = spectralTransform.RA_shape;
		
		BasicIO::Array_properties<3> array_properties;

		array_properties.shape_full_complex_array = Nx, Ny, Nz/2+1;
		array_properties.shape_full_real_array = Nx, Ny, Nz+2;

		//For HDF5 this is always the case, because it is more efficient to divide the data along the slowest direction.
		array_properties.id_complex_array = my_id, 0, 0;
		array_properties.id_real_array = my_id, 0, 0;

		array_properties.numprocs_complex_array = numprocs, 1, 1;
		array_properties.numprocs_real_array = numprocs, 1, 1;

		if (global.io.N_in_reduced.size() == 3)
			array_properties.shape_N_in_reduced = global.io.N_in_reduced[0], global.io.N_in_reduced[1], global.io.N_in_reduced[2]/2+1;
		
		if (global.io.N_out_reduced.size() == 3)
			array_properties.shape_N_out_reduced = global.io.N_out_reduced[0], global.io.N_out_reduced[1], global.io.N_out_reduced[2]/2+1;

		array_properties.Fourier_directions = 1,1,1;
		array_properties.Z = 2;

		array_properties.datatype_complex_space = BasicIO::H5T_Complex;
		array_properties.datatype_real_space = BasicIO::H5T_Real;

		BasicIO::Set_H5_plans(array_properties, this);
	}
	
	else if (Ny == 1) {
		if (master)
			cerr << "ERROR: 2D Not implemented for FFF basis, Please use FFTW basis (uses original FFTW functions)" << endl;
		exit(1);
	}


	// alias..
	kfactor.resize(4);
	kfactor[0]= global.field.kfactor[0];
	kfactor[1]= global.field.kfactor[1];
	kfactor[2]= global.field.kfactor[2];
	kfactor[3]= global.field.kfactor[3];
	
	local_Nx=spectralTransform.maxfx;
	local_Ny=spectralTransform.maxfy;
	local_Nz=spectralTransform.maxfz;
	
	local_Nx_start=spectralTransform.fx_start;
	local_Ny_start=spectralTransform.fy_start;
	local_Nz_start=spectralTransform.fz_start;

	shape_complex_array = global.field.shape_complex_array;
	shape_real_array = global.field.shape_real_array;
	
	global.temp_array.X.resize(shape_complex_array);
	global.temp_array.X2.resize(shape_complex_array);
	global.temp_array.X_transform.resize(shape_complex_array);
	
	global.temp_array.Xr.resize(shape_real_array);
	global.temp_array.Xr2.resize(shape_real_array);
	
	// temp arrays
	// Being used in void ArrayOps::Get_XY_plane(Array<Complex,3> A, Array<Complex,2> plane_xy, int kz, string configuration)
	global.temp_array.plane_xy.resize(Nx, Ny);

	global.temp_array.plane_xy_inproc.resize(Nx,spectralTransform.maxfy);
	

}

