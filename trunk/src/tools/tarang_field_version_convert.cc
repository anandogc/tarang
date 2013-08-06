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

/*! \file tarang_field_version_convert.cc
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author A. G. Chatterjee, M. K. Verma
 * @date 18 October 2012
 * 
 * @bugs  No known bug
 */ 

#include "BasicIO.h"
#include <iostream>

using namespace std;

Global global;

int Convert_PRINFTY(string file_name, vector<H5_dataset_meta> field_meta);
int Convert_FLUID(string file_name, vector<H5_dataset_meta> field_meta);
int Convert_RBC(string file_name, vector<H5_dataset_meta> field_meta);
int Convert_MHD(string file_name, vector<H5_dataset_meta> field_meta);
int Convert_MRBC(string file_name, vector<H5_dataset_meta> field_meta);


int main(int argc, char** argv)
{
	  	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);

	global.mpi.master = (global.mpi.my_id == 0);
	
	if (argc != 2){
        cerr<< "Usage:\n    " << argv[0] << " /path/to/input/file/\nThis path must contain input.h5" << endl;
        return 1;
	}

	vector<H5_dataset_meta> field_meta;

	global.Global_init_default();
	global.Process_global_vars_basic();

	BasicIO::data_in_folder=argv[1];
	BasicIO::data_out_folder=argv[1];
	
	BasicIO::Initialize();

	int local_Nx, Ny, Nz;

	string file_name = BasicIO::data_in_folder + "/input.h5";

	field_meta = BasicIO::Get_meta(file_name);

	int num_fields = field_meta.size();

	int status;

	switch (num_fields){
		case 1:		
			status = Convert_PRINFTY(file_name, field_meta);
			break;
		case 3:		
			status = Convert_FLUID(file_name, field_meta);
			break;
		case 4:		
			status = Convert_RBC(file_name, field_meta);
			break;
		case 6:		
			status = Convert_MHD(file_name, field_meta);
			break;
		case 7:		
			status = Convert_MRBC(file_name, field_meta);
			break;
		default:
			cerr << "Not a tarang-1 field" << endl;
			break;
	}

	BasicIO::Finalize();
	MPI_Finalize();
	
	return 0;
}


int Convert_PRINFTY(string file_name, vector<H5_dataset_meta> field_meta)
{
	local_Nx = field_meta[0].dimensions[0]/global.mpi.numprocs;
	Ny = field_meta[0].dimensions[1];
	Nz = field_meta[0].dimensions[2]-2;
	
	Array<complx, 3> T(local_Nx, Ny, Nz/2+1);

	BasicIO::CV_Input_HDF5(1, &T);
	BasicIO::CV_Output_HDF5(1, &T);

	return 0;
}

int Convert_FLUID(string file_name, vector<H5_dataset_meta> field_meta)
{
	Nx = field_meta[0].dimensions[0];
	Ny = field_meta[0].dimensions[1];
	Nz = field_meta[0].dimensions[2]-2;

	local_Nx = Nx/global.mpi.numprocs;

	if (global.mpi.master)
		cout << "INC_FLUID: [" << Nx << ", " << Ny << ", " << Nz <<  "] " << endl;
	
	Array<complx, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV3(local_Nx, Ny, Nz/2+1);

	if (field_meta[2].dimensions[2] == 2){
		BasicIO::CV_Input_plane_HDF5(3, &CV1, &CV2, &CV3);
		BasicIO::CV_Output_plane_HDF5(3, &CV1, &CV2, &CV3);
	}
	else{
		BasicIO::CV_Input_HDF5(3, &CV1, &CV2, &CV3);
		BasicIO::CV_Output_HDF5(3, &CV1, &CV2, &CV3);
	}
	
	return 0;
}

int Convert_RBC(string file_name, vector<H5_dataset_meta> field_meta)
{

	local_Nx = field_meta[0].dimensions[0]/global.mpi.numprocs;
	Ny = field_meta[0].dimensions[1];
	Nz = field_meta[0].dimensions[2]-2;
	
	Array<complx, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV3(local_Nx, Ny, Nz/2+1);

	Array<complx, 3> T(local_Nx, Ny, Nz/2+1);

	if (field_meta[2].dimensions[2] == 2){
		BasicIO::CV_Input_plane_HDF5(4, &CV1, &CV2, &CV3, &T);
		BasicIO::CV_Output_plane_HDF5(4, &CV1, &CV2, &CV3, &T);
	}
	else{
		BasicIO::CV_Input_HDF5(4, &CV1, &CV2, &CV3, &T);
		BasicIO::CV_Output_HDF5(4, &CV1, &CV2, &CV3, &T);
	}

	return 0;
}

int Convert_MHD(string file_name, vector<H5_dataset_meta> field_meta)
{
	local_Nx = field_meta[0].dimensions[0]/global.mpi.numprocs;
	Ny = field_meta[0].dimensions[1];
	Nz = field_meta[0].dimensions[2]-2;
	
	Array<complx, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV3(local_Nx, Ny, Nz/2+1);

	Array<complx, 3> CW1(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CW2(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CW3(local_Nx, Ny, Nz/2+1);

	if (field_meta[2].dimensions[2] == 2){
		BasicIO::CV_Input_plane_HDF5(6, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3);
		BasicIO::CV_Output_plane_HDF5(6, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3);
	}
	else{
		BasicIO::CV_Input_HDF5(6, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3);
		BasicIO::CV_Output_HDF5(6, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3);
	}
	
	return 0;
}

int Convert_MRBC(string file_name, vector<H5_dataset_meta> field_meta)
{
	local_Nx = field_meta[0].dimensions[0]/global.mpi.numprocs;
	Ny = field_meta[0].dimensions[1];
	Nz = field_meta[0].dimensions[2]-2;
	
	Array<complx, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CV3(local_Nx, Ny, Nz/2+1);

	Array<complx, 3> CW1(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CW2(local_Nx, Ny, Nz/2+1);
	Array<complx, 3> CW3(local_Nx, Ny, Nz/2+1);

	Array<complx, 3> T(local_Nx, Ny, Nz/2+1);

	if (field_meta[2].dimensions[2] == 2){
		BasicIO::CV_Input_plane_HDF5(7, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3, &T);
		BasicIO::CV_Output_plane_HDF5(7, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3, &T);
	}
	else{
		BasicIO::CV_Input_HDF5(7, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3, &T);
		BasicIO::CV_Output_HDF5(7, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3, &T);
	}

	return 0;
}
