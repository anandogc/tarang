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

inline bool file_exists(const std::string& name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}


int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);


	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);

	global.Parse(argc, argv, false);


	BasicIO::data_in_folder = "fortran";
	BasicIO::data_out_folder = "tarang";
	
    BasicIO::Initialize();

	master = (global.mpi.my_id == 0);

    vector<BasicIO::H5_dataset_meta> meta;

    meta = BasicIO::Get_meta(global.io.data_dir + "/" + BasicIO::data_in_folder + "/input.h5", "/");


    global.field.N[1]=meta[0].dimensions[0];
    global.field.N[2]=meta[0].dimensions[1];
    global.field.N[3]=meta[0].dimensions[2];

	//global.io.input_vx_vy_switch=file_exists(BasicIO::data_in_folder + "/U.V3kz0.h5");

    // global.io.data_dir=".";
        
	//global.io.output_vx_vy_switch=global.io.input_vx_vy_switch;

	global.Process_basic_vars();

	int local_Nx=Nx/numprocs;
    int local_Ny=Ny/numprocs;
    int local_Nz=(Nz/2+1)/numprocs;

    Array<Real,3> A(Nz,Ny,Nx);
    Array<Real,3> B(Ny,Nz+2,Nx);

    Real* A_data = A.data(); 
    Real* B_data = B.data();
    

 	//Configure fortran input

    H5_Planner fortran_output;

	BasicIO::Array_properties<3> array_properties;

	array_properties.shape_full_complex_array = Nx, Ny, Nz;
	array_properties.shape_full_real_array = Nz, Ny, Nx;

	array_properties.id_complex_array = 0, my_id, 0;
	array_properties.id_real_array = 0, my_id, 0;

	array_properties.numprocs_complex_array = 1, numprocs, 1;
	array_properties.numprocs_real_array = 1, numprocs, 1;

	if (global.io.N_in_reduced.size() == 3)
		array_properties.shape_N_in_reduced = global.io.N_in_reduced[1], global.io.N_in_reduced[2]/2+1, global.io.N_in_reduced[0];
	
	if (global.io.N_out_reduced.size() == 3)
		array_properties.shape_N_out_reduced = global.io.N_out_reduced[1], global.io.N_out_reduced[2]/2+1, global.io.N_out_reduced[0];

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 2;

	array_properties.datatype_complex_space = BasicIO::H5T_Real;
	array_properties.datatype_real_space = BasicIO::H5T_Real;

	BasicIO::Set_H5_plans(array_properties, &fortran_output);



	//Configure tarang 2.3 output
	H5_Planner tarang_2_3;

	// BasicIO::Array_properties<3> array_properties;

	array_properties.shape_full_complex_array = Ny, Nz/2+1, Nx;
	array_properties.shape_full_real_array = Ny, Nz+2, Nx;

	array_properties.id_complex_array = my_id, 0, 0;
	array_properties.id_real_array = 0, 0, my_id;

	array_properties.numprocs_complex_array = numprocs, 1, 1;
	array_properties.numprocs_real_array = 1, 1, numprocs;

	if (global.io.N_in_reduced.size() == 3)
		array_properties.shape_N_in_reduced = global.io.N_in_reduced[1], global.io.N_in_reduced[2]/2+1, global.io.N_in_reduced[0];
	
	if (global.io.N_out_reduced.size() == 3)
		array_properties.shape_N_out_reduced = global.io.N_out_reduced[1], global.io.N_out_reduced[2]/2+1, global.io.N_out_reduced[0];

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 1;

	array_properties.datatype_complex_space = BasicIO::H5T_Complex;
	array_properties.datatype_real_space = BasicIO::H5T_Real;

	BasicIO::Set_H5_plans(array_properties, &tarang_2_3);

	

	if (master) cout << "Processing PS/vx .." << flush;
	BasicIO::Read(A_data, fortran_output.H5_real, "input", "PS/vx");

	for (int x=0; x<Nx; x++){
		for (int y=0; y<Ny; y++){
			for (int z=0; z<Nz; z++){
				B(y,z,x) = A(z,y,x);
			}
		}
	}

	BasicIO::Write(B_data, tarang_2_3.H5_real, "", "U.V1r");
	if (master) cout << "done" << endl;


	
	if (master) cout << "Processing PS/vy .." << flush;
	BasicIO::Read(A_data, fortran_output.H5_real, "input", "PS/vy");

	for (int x=0; x<Nx; x++){
		for (int y=0; y<Ny; y++){
			for (int z=0; z<Nz; z++){
				B(y,z,x) = A(z,y,x);
			}
		}
	}

	BasicIO::Write(B_data, tarang_2_3.H5_real, "", "U.V2r");
	if (master) cout << "done" << endl;



	if (master) cout << "Processing PS/vz .." << flush;
	BasicIO::Read(A_data, fortran_output.H5_real, "input", "PS/vz");

	for (int x=0; x<Nx; x++){
		for (int y=0; y<Ny; y++){
			for (int z=0; z<Nz; z++){
				B(y,z,x) = A(z,y,x);
			}
		}
	}

	BasicIO::Write(B_data, tarang_2_3.H5_real, "", "U.V3r");
	if (master) cout << "done" << endl;

	
	BasicIO::Finalize();
	MPI_Finalize();
	
	return 0;
}

