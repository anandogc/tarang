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


	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	master = (my_id == 0);

	global.Parse(argc, argv, false);

	BasicIO::LOG_filename="";
    BasicIO::Initialize(global.io.data_dir);
    
    BasicIO::data_in_folder=global.io.data_dir + "/tarang2.3/";
    BasicIO::data_out_folder=global.io.data_dir + "/tarang2.4/";

    vector<string> list = BasicIO::Get_file_list(BasicIO::data_in_folder);

    vector<BasicIO::H5_dataset_meta> meta;

	meta = BasicIO::Get_meta(BasicIO::data_in_folder + list[0], "/");

    Ny=meta[0].dimensions[0];
    Nz=meta[0].dimensions[1]*2-2;
    Nx=meta[0].dimensions[2];

	int local_Ny=Ny/numprocs;

    Array<Complex,3> A(local_Ny,Nz/2+1,Nx);
    Array<Complex,3> B(Nx,local_Ny,Nz/2+1);


 	//Configure tarang 2.3 input
    H5_Planner tarang_2p3;

	BasicIO::Array_properties<3> array_properties;

	array_properties.shape_full_complex_array = Ny, Nz/2+1, Nx;
	array_properties.shape_full_real_array = Ny, Nz+2, Nx;

	array_properties.id_complex_array = my_id, 0, 0;
	array_properties.id_real_array = my_id, 0, 0;

	array_properties.numprocs_complex_array = numprocs, 1, 1;
	array_properties.numprocs_real_array = numprocs, 1, 1;

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 1;

	array_properties.datatype_complex_space = BasicIO::H5T_Complex;
	array_properties.datatype_real_space = BasicIO::H5T_Real;
	
	BasicIO::Set_H5_plans(array_properties, &tarang_2p3);

	
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



    for (int i=0; i<list.size(); i++)
    {
	    
		meta = BasicIO::Get_meta(BasicIO::data_in_folder + list[i], "/");

	    Ny=meta[0].dimensions[0];
	    Nz=meta[0].dimensions[1]*2-2;
	    Nx=meta[0].dimensions[2];

		bool vx_vy_switch;

		vx_vy_switch=(Nz==0);

		if (master)
			cout << "Rearranging " << meta[0].name << ": (" << Nx << "x" << Ny << "x" << Nz/2+1 << ") (y,z,x) -> (x,y,z)" << " .." << flush;

		if (!vx_vy_switch)
			BasicIO::Read(A.data(), tarang_2p3.H5_full, meta[0].name, "/"+meta[0].name);
		else
			BasicIO::Read(A.data(), tarang_2p3.H5_kz0_full, meta[0].name, "/"+meta[0].name);


		for (int x=0; x<Nx; x++){
			for (int y=0; y<local_Ny; y++){
				for (int z=0; z<Nz/2+1; z++){
					B(x,y,z) = A(y,z,x);
				}
			}
		}

		if (!vx_vy_switch)
			BasicIO::Write(B.data(), tarang_2p4.H5_full, "", meta[0].name);
		else
			BasicIO::Write(B.data(), tarang_2p4.H5_kz0_full, "", meta[0].name);

		if (master) cout << "done" << endl;
	}
	
	BasicIO::Finalize();
	MPI_Finalize();
	
	return 0;
}
