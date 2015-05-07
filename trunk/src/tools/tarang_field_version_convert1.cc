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
#include "BasicIO_SP.h"

Global global;

// SpectralTransform spectralTransform;

Uniform<Real> SPECrand;

struct H5_Planner{
    BasicIO::H5_plan H5_real;
    BasicIO::H5_plan H5_full;
    BasicIO::H5_plan H5_kz0_full;
    BasicIO::H5_plan H5_in_reduced;
    BasicIO::H5_plan H5_out_reduced;
    BasicIO::H5_plan H5_in_kz0_reduced;
    BasicIO::H5_plan H5_out_kz0_reduced;
};

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);


	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);

	global.Parse(argc, argv, false);

    BasicIO::Initialize(global.io.data_dir);

	master = (global.mpi.my_id == 0);

    vector<BasicIO::H5_dataset_meta> meta;

	if (!BasicIO::File_exists(BasicIO::data_in_folder  + "/input.h5")) {
		if (master)
			cout << BasicIO::data_in_folder + "/input.h5 does not exists." << endl;
		return 1;
	}

    meta = BasicIO::Get_meta(BasicIO::data_in_folder  + "/input.h5", "/");

    global.field.N[1]=meta[0].dimensions[0];
    global.field.N[2]=meta[0].dimensions[1];
    global.field.N[3]=meta[0].dimensions[2]-2;

	global.Process_basic_vars();

	int local_Ny=Ny/numprocs;

    Array<Complex,3> A(Nx,local_Ny,Nz/2+1);
    Array<Complex,3> B(local_Ny,Nz/2+1,Nx);


 	//Configure tarang 1 input
    H5_Planner tarang1;

	BasicIO::Array_properties<3> array_properties;

	array_properties.shape_full_complex_array = Nx, Ny, Nz+2;
	array_properties.shape_full_real_array = Nx, Ny, Nz+2;

	array_properties.id_complex_array = 0, my_id, 0;
	array_properties.id_real_array = 0, my_id, 0;

	array_properties.numprocs_complex_array = 1, numprocs, 1;
	array_properties.numprocs_real_array = 1, numprocs, 1;

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 2;

	array_properties.datatype_complex_space = BasicIO::H5T_Real;
	array_properties.datatype_real_space = BasicIO::H5T_Real;
	
	BasicIO::Set_H5_plans(array_properties, &tarang1);

	
	//Configure tarang 2.4 output
	H5_Planner tarang_2p4;

	array_properties.shape_full_complex_array = Nx, Ny, Nz/2+1;
	array_properties.shape_full_real_array = Nx, Ny, Nz+2;

	array_properties.id_complex_array = 0, my_id, 0;
	array_properties.id_real_array = my_id, 0, 0;

	array_properties.numprocs_complex_array = 1, numprocs, 1;
	array_properties.numprocs_real_array = numprocs, 1, 1;

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 2;

	array_properties.datatype_complex_space = BasicIO::H5T_Complex;
	array_properties.datatype_real_space = BasicIO::H5T_Real;

	BasicIO::Set_H5_plans(array_properties, &tarang_2p4);

	static map<string, string> transition_table;
	transition_table["CV1"]="U.V1";
	transition_table["CV2"]="U.V2";
	transition_table["CV3"]="U.V3";
	
	transition_table["CW1"]="B.V1";
	transition_table["CW2"]="B.V2";
	transition_table["CW3"]="B.V3";
	
	transition_table["T"]="T.F";

	string output_dataset_name;
	bool vx_vy_switch;


	for (int i=0; i<meta.size(); i++) {
		Nx=meta[i].dimensions[0];
		Ny=meta[i].dimensions[1];
		Nz=meta[i].dimensions[2]-2;

		vx_vy_switch=(Nz==0);

		output_dataset_name = transition_table[meta[i].name];

		if (vx_vy_switch)
			output_dataset_name += "kz0";

		if (master)
			cout << "Converting " << meta[i].name << " -> " << output_dataset_name << " .." << flush;

		if (!vx_vy_switch) {
			BasicIO::Read(A.data(), tarang1.H5_full, "input", "/"+meta[i].name);
			BasicIO::Write(A.data(), tarang_2p4.H5_full, "", output_dataset_name);
		}
		else {
			BasicIO::Read(A.data(), tarang1.H5_kz0_full, "input", "/"+meta[i].name);
			BasicIO::Write(A.data(), tarang_2p4.H5_kz0_full, "", output_dataset_name);
		}

		if (master) cout << "done" << endl;
	}

	
	BasicIO::Finalize();
	MPI_Finalize();
	
	return 0;
}
