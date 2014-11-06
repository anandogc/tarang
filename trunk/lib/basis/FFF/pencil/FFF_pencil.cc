/*
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

/*! \file FFF_pencil.cc 
 * 
 * @sa FFF_pencil.cc
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "FFF_pencil.h"
#include "BasicIO_SP.h"

FFF_PENCIL::FFF_PENCIL()
{

		// kfactor, xfactor,  L, Delta_x
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
		for (int i=1; i<=3; i++) 
			global.field.L[i] = 2*M_PI/global.field.kfactor[i];
	}
	
	for (int i=1; i<=3; i++) {
		global.field.Delta_x[i] =  global.field.L[i]/global.field.N[i];
		global.field.xfactor[i] = 1/global.field.kfactor[i];
	}
	
		// array size proper or not?
	
	
	bool Nxproperdiv, Nyproperdiv, Nzproperdiv; 
	
	Nxproperdiv = true;//((global.field.Nx%global.mpi.num_p_hor == 0) && (global.field.Nx%global.mpi.num_p_vert == 0));
	Nyproperdiv = true;//(global.field.Ny%global.mpi.num_p_hor == 0);
	Nzproperdiv = true;//((global.field.Nz/2+1)%global.mpi.num_p_vert == 0);
	
	if (!(Nxproperdiv && Nyproperdiv && Nzproperdiv)) {
		if (global.mpi.master) cerr << "ERROR: Array not being divided equally.  Check dimensions" << endl;
		exit(1);
	}
	
	//global.mpi.num_p_hor is assigned in global via mpirun argument

	spectralTransform.Init("FFF", Nx, Ny, Nz, global.mpi.num_p_cols);

	global.field.shape_complex_array = Nx, spectralTransform.local_Ny_row, spectralTransform.local_Nz_col;
	global.field.shape_real_array = spectralTransform.local_Nx_row, spectralTransform.local_Ny_col, Nz+2;


	//********
	
	global.mpi.num_x_procs = 1;
	global.mpi.num_y_procs = global.mpi.num_p_cols;
	global.mpi.num_z_procs = global.mpi.num_p_rows;
	
	global.mpi.my_x_pcoord = 0;
	global.mpi.my_y_pcoord = spectralTransform.my_id_row;
	global.mpi.my_z_pcoord = spectralTransform.my_id_col;
	
	global.mpi.num_x_procs_real = global.mpi.num_p_cols;
	global.mpi.num_y_procs_real = global.mpi.num_p_rows;
	global.mpi.num_z_procs_real = 1;
	
	global.mpi.my_x_pcoord_real = spectralTransform.my_id_row;
	global.mpi.my_y_pcoord_real = spectralTransform.my_id_col;
	global.mpi.my_z_pcoord_real = 0;
	
	
	global.field.maxlx = Nx-1;
	global.field.maxly = spectralTransform.local_Ny_row-1;
	global.field.maxlz = spectralTransform.local_Nz_col-1;

	global.field.lx_start = global.mpi.my_x_pcoord*global.field.maxlx;
	global.field.ly_start = global.mpi.my_y_pcoord*global.field.maxly;
	global.field.lz_start = global.mpi.my_z_pcoord*global.field.maxlz;


	global.field.maxrx = spectralTransform.local_Nx_row;
	global.field.maxry = spectralTransform.local_Ny_col;
	global.field.maxrz = spectralTransform.Nz+2;

	global.field.rx_start = global.mpi.my_x_pcoord_real*global.field.maxrx;
	global.field.ry_start = global.mpi.my_y_pcoord_real*global.field.maxry;
	global.field.rz_start = global.mpi.my_z_pcoord_real*global.field.maxrz;

    // alias..
	kfactor.resize(4);
	kfactor[0]= global.field.kfactor[0];
	kfactor[1]= global.field.kfactor[1];
	kfactor[2]= global.field.kfactor[2];
	kfactor[3]= global.field.kfactor[3];
	
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

	maxlx = global.field.maxlx;
	maxly = global.field.maxly;
	maxlz = global.field.maxlz;

	lx_start = global.field.lx_start;
	ly_start = global.field.ly_start;
	lz_start = global.field.lz_start;

	maxrx = global.field.maxrx;
	maxry = global.field.maxry;
	maxrz = global.field.maxrz;

	rx_start = global.field.rx_start;
	ry_start = global.field.ry_start;
	rz_start = global.field.rz_start;

	global.temp_array.X_transform.resize(shape_complex_array);

	global.temp_array.X.resize(shape_complex_array);
    global.temp_array.X2.resize(shape_complex_array);
	
	global.temp_array.Xr.resize(shape_real_array);
    global.temp_array.Xr2.resize(shape_real_array);

    // temp arrays
    // Being used in void ArrayOps::Get_XY_plane(Array<complx,3> A, Array<complx,2> plane_xy, int kz, string configuration)
	global.temp_array.plane_xy.resize(Nx, Ny);
	
    global.temp_array.plane_xy_inproc.resize(Nx, maxly);
	
	
	BasicIO::Array_properties<3> array_properties;
	array_properties.shape_full_complex_array = Nx, Ny, Nz/2+1;
	array_properties.shape_full_real_array = Nx, Ny, Nz+2;

	array_properties.id_complex_array = my_x_pcoord, my_y_pcoord, my_z_pcoord;
	array_properties.id_real_array = my_x_pcoord_real, my_y_pcoord_real, my_z_pcoord_real;

	array_properties.numprocs_complex_array = num_x_procs, num_y_procs, num_z_procs;
	array_properties.numprocs_real_array = num_x_procs_real, num_y_procs_real, num_z_procs_real;

	if (global.io.N_in_reduced.size() == 3)
		array_properties.shape_N_in_reduced = global.io.N_in_reduced[0], global.io.N_in_reduced[1], global.io.N_in_reduced[2]/2+1;
	
	if (global.io.N_out_reduced.size() == 3)
		array_properties.shape_N_out_reduced = global.io.N_out_reduced[0], global.io.N_out_reduced[1], global.io.N_out_reduced[2]/2+1;

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 2;

	array_properties.datatype_complex_space = BasicIO::H5T_COMPLX;
	array_properties.datatype_real_space = BasicIO::H5T_DP;


	BasicIO::Set_H5_plans(array_properties, this);

	int count=Nx;
	int blocklength=2*maxly;
	int stride=2*Ny;
	MPI_Type_vector(count,blocklength,stride,MPI_DP,&MPI_Vector_y_plane_block);
	MPI_Type_commit(&MPI_Vector_y_plane_block);

	int lower_bound=0;
	int extent=2*maxly*sizeof(DP);
	MPI_Type_create_resized(MPI_Vector_y_plane_block, lower_bound, extent, &MPI_Vector_resized_y_plane_block);
	MPI_Type_commit(&MPI_Vector_resized_y_plane_block);	
}


//*****************************************************************************************
