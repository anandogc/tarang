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
	/*MPI_Init(&argc, &argv);


	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);

	global.Parse(argc, argv, false);

	
    BasicIO::Initialize();

	BasicIO::data_in_folder = "tarang2.2";
	BasicIO::data_out_folder = "tarang2.3";

    vector<BasicIO::H5_dataset_meta> meta;

    meta = BasicIO::Get_meta(global.io.data_dir + "/" + BasicIO::data_in_folder + "/U.V1.h5");

    global.field.N[1]=meta[0].dimensions[0];
    global.field.N[2]=meta[0].dimensions[1];
    global.field.N[3]=meta[0].dimensions[2]-2;

	global.io.input_vx_vy_switch=BasicIO::File_exists(global.io.data_dir + "/" + BasicIO::data_in_folder + "/U.V3kz0.h5");

    
    
	global.io.output_vx_vy_switch=global.io.input_vx_vy_switch;


	global.Process_basic_vars();

	if (master) {
		cout << "N = [" << Nx << "," << Ny << "," << Nz << "]" << endl;
		cout << boolalpha << "input_vx_vy_switch = " << global.io.input_vx_vy_switch << endl;
		cout << boolalpha << "output_vx_vy_switch = " << global.io.output_vx_vy_switch << endl;
		
		cout << endl;
	}

	int local_Nx=Nx/numprocs;
    int local_Ny=Ny/numprocs;
    int local_Nz=(Nz/2+1)/numprocs;

    Array<Complex,3> A(Nx,local_Ny,Nz/2+1);
    Array<Complex,3> B(local_Ny,Nz/2+1,Nx);

    Complex* A_data = A.data(); 
    Complex* B_data = B.data();
    

 	//Configure tarang 2.2 input

    H5_Planner tarang_2_2;

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

	BasicIO::Set_H5_plans(array_properties, &tarang_2_2);



	//Configure tarang 2.3 output
	H5_Planner tarang_2_3;

	array_properties.shape_full_complex_array = Ny, Nz/2+1, Nx;
	array_properties.shape_full_real_array = Ny, Nz+2, Nx;

	array_properties.id_complex_array = my_id, 0, 0;
	array_properties.id_real_array = 0, 0, my_id;

	array_properties.numprocs_complex_array = numprocs, 1, 1;
	array_properties.numprocs_real_array = 1, 1, numprocs;

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 1;

	array_properties.datatype_complex_space = BasicIO::H5T_Complex;
	array_properties.datatype_real_space = BasicIO::H5T_Real;

	BasicIO::Set_H5_plans(array_properties, &tarang_2_3);
	
	//T.F
	if (master) cout << "Processing T.F .." << flush;
	BasicIO::Read(A_data, tarang_2_2.H5_full, "T.F");

	for (int x=0; x<Nx; x++){
		for (int y=0; y<local_Ny; y++){
			for (int z=0; z<Nz/2+1; z++){
				*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
			}
		}
	}

	BasicIO::Write(B_data, tarang_2_3.H5_full, "", "T.F");
	if (master) cout << "done" << endl;

	
	//U.V1
	if (master) cout << "Processing U.V1 .." << flush;
	BasicIO::Read(A_data, tarang_2_2.H5_full, "U.V1");

	for (int x=0; x<Nx; x++){
		for (int y=0; y<local_Ny; y++){
			for (int z=0; z<Nz/2+1; z++){
				*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
			}
		}
	}

	BasicIO::Write(B_data, tarang_2_3.H5_full, "", "U.V1");
	if (master) cout << "done" << endl;

	//U.V2
	if (master) cout << "Processing U.V2 .." << flush;
	BasicIO::Read(A_data, tarang_2_2.H5_full, "U.V2");

	for (int x=0; x<Nx; x++){
		for (int y=0; y<local_Ny; y++){
			for (int z=0; z<Nz/2+1; z++){
				*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
			}
		}
	}

	BasicIO::Write(B_data, tarang_2_3.H5_full, "", "U.V2");
	if (master) cout << "done" << endl;

	//To be checked after this

		//U.V3
		if (global.io.input_vx_vy_switch){
			if (master) cout << "Processing U.V3kz0 .." << flush;
			A=0;
			B=0;
			BasicIO::Read(A_data, "U.V3kz0", tarang_2_2.H5_kz0_full);

			for (int x=0; x<Nx; x++){
				for (int y=0; y<local_Ny; y++){
					for (int z=0; z<1; z++){
						*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
					}
				}
			}

			BasicIO::Write(B_data, "U.V3kz0", tarang_2_3.H5_kz0_full);
			if (master) cout << "done" << endl;
		}
	
		else {
			if (master) cout << "Processing U.V3 .." << flush;
			BasicIO::Read(A_data, "U.V3", tarang_2_2.H5_full);

			for (int x=0; x<Nx; x++){
				for (int y=0; y<local_Ny; y++){
					for (int z=0; z<Nz/2+1; z++){
						*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
					}
				}
			}

			BasicIO::Write(B_data, "U.V3", tarang_2_3.H5_full);
			if (master) cout << "done" << endl;
		}
	}

	if (file_exists(BasicIO::data_in_folder + "/B.V1.h5")) {
		//B.V1
		if (master) cout << "Processing B.V1 .." << flush;
		BasicIO::Read(A_data, "B.V1", tarang_2_2.H5_full);

		for (int x=0; x<Nx; x++){
			for (int y=0; y<local_Ny; y++){
				for (int z=0; z<Nz/2+1; z++){
					*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
				}
			}
		}

		BasicIO::Write(B_data, "B.V1", tarang_2_3.H5_full);
		if (master) cout << "done" << endl;

		//B.V2
		if (master) cout << "Processing B.V2 .." << flush;
		BasicIO::Read(A_data, "B.V2", tarang_2_2.H5_full);

		for (int x=0; x<Nx; x++){
			for (int y=0; y<local_Ny; y++){
				for (int z=0; z<Nz/2+1; z++){
					*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
				}
			}
		}

		BasicIO::Write(B_data, "B.V2", tarang_2_3.H5_full);
		if (master) cout << "done" << endl;
	
		//B.V3
		if (global.io.input_vx_vy_switch){
			if (master) cout << "Processing B.V3kz0 .." << flush;
			A=0;
			B=0;
			BasicIO::Read(A_data, "B.V3kz0", tarang_2_2.H5_kz0_full);

			for (int x=0; x<Nx; x++){
				for (int y=0; y<local_Ny; y++){
					for (int z=0; z<1; z++){
						*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
					}
				}
			}

			BasicIO::Write(B_data, "B.V3kz0", tarang_2_3.H5_kz0_full);
			if (master) cout << "done" << endl;
		}
	
		else {
			if (master) cout << "Processing B.V3 .." << flush;
			BasicIO::Read(A_data, "B.V3", tarang_2_2.H5_full);

			for (int x=0; x<Nx; x++){
				for (int y=0; y<local_Ny; y++){
					for (int z=0; z<Nz/2+1; z++){
						*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
					}
				}
			}

			BasicIO::Write(B_data, "B.V3", tarang_2_3.H5_full);
			if (master) cout << "done" << endl;
		}
	}

	BasicIO::Finalize();
	MPI_Finalize();*/
	
	return 0;
}

