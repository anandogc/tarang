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


/*! \file  test_BasicIO.cc
 *
 * @brief  Main file for test-BasicIO
 * @sa BasicIO.h
 *
 *
 *	@author  A. G. Chatterjee, M. K. Verma
 *	@version 4.0
 *	@date Aug 2012
 */

#include "fluid.h"
#include <unistd.h>

//*********************************************************************************************

SpectralTransform spectralTransform;

Universal *universal;
Global global;

Uniform<Real> SPECrand;

int main(int argc, char** argv)
{
	/*MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);
	
    global.mpi.master_id=0;
	global.Parse(argc, argv, true);
  
	master = (global.mpi.my_id == 0);

	char cwd[256];
    getcwd(cwd, sizeof(cwd));

	//Set up global
	global.Init_defaults();

	global.program.basis_type = "FFF";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "FFF";

	global.field.N[0] = 0;
	global.field.N[1] = 16;
	global.field.N[2] = 4;
	global.field.N[3] = 8;

	global.time.now=0;

	global.io.data_dir=cwd;
	global.io.N_in_reduced.resize(4);
	global.io.N_in_reduced[0]=0;
	global.io.N_in_reduced[1]=32;
	global.io.N_in_reduced[2]=32;
	global.io.N_in_reduced[3]=32;

	global.io.N_out_reduced.resize(4);
	global.io.N_out_reduced[0]=0;
	global.io.N_out_reduced[1]=8;
	global.io.N_out_reduced[2]=8;
	global.io.N_out_reduced[3]=8;

	//Set up BasicIO
	BasicIO::data_in_folder="/test/time_"+To_string(global.time.now)+"/";
	BasicIO::data_out_folder="/test/";
	
	BasicIO::Initialize();

	//Instantiate a Basis
	//global.Process_global_vars_basic();
	//universal=new FFF_SLAB;
	//global.Process_global_vars_advanced();

	//Or Perform the test directly
	Array<int,1> filter[3];
	filter[0].resize(8);
	filter[1].resize(5);
	filter[2].resize(6);

	filter[0] = 1;
	filter[1] = 1;
	filter[2] = 1;


	int my_id[3];
	my_id[0] = global.mpi.my_id;
	my_id[1] = 0;
	my_id[2] = 0;

	int numprocs[3];
	numprocs[0] = global.mpi.numprocs;
	numprocs[1] = 1;
	numprocs[2] = 1;



	BasicIO::H5_plan cube = BasicIO::Set_plan(3, my_id, numprocs, filter, filter, BasicIO::H5T_Complex);

	Array<Complex,3> A(filter[0].size()/global.mpi.numprocs, filter[1].size(), filter[2].size());
	A=0;
	int index=0;
	for (int i=0; i<A.extent(0); i++)
		for (int j=0; j<A.extent(1); j++)
			for (int k=0; k<A.extent(2); k++){
				index+=2;
				A(i,j,k) = Complex(index, index+1);
			}

	cout << A << endl;

	BasicIO::Write(A.data(), "abcd", cube);


	BasicIO::Finalize();
	MPI_Finalize();*/
	
	return 0;
}



//********************************** End of main.cc *******************************************



