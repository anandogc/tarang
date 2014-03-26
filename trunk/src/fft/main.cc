
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


/*! \file  main.cc
 *
 * @brief  Main file for Tarang (sequential). Calls various modules.
 * @sa main.h
 *
 * @note Reads which turbulence to start (Fluid, MHD, Passive scalcar, RBslip),
 *					the dimensionality of simulation (2D or 3D), and
 *					which directory to use for input and output.
 *
 *	@author  M. K. Verma
 *	@version 4.0 MPI
 *	@date Sept 2008
 *
 *  USAGE: mpirun -np numprocs test-fft num_p_hor
 *  For pencil numproc must be > 1 and multiple of hum_p_hor
 *
 */

#include "main.h"

Universal *universal;
Global global;

SpectralTransform spectralTransform;

Uniform<DP> SPECrand;

Array<DP,3> Ar, Ar_old;
Array<complx,3> A, A_old;

//*********************************************************************************************

int main(int argc, char** argv)
{
	
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);
	global.mpi.master= (global.mpi.my_id == global.mpi.master_id);
	
	global.Parse(argc, argv, false);

	int result = 0;

	global.Init_defaults();
	BasicIO::Initialize();

	SPECrand.seed((unsigned int)time(0));				// Initialize random number seed
	
	// if (master)
		cout << "************ SLAB TESTS STARTS HERE *********** " << endl;
	// 
	time_t          start;
    time_t          end;
	
	global.field.kfactor[1] = 3.14159;
	global.field.kfactor[2] = 2.22144;
	global.field.kfactor[3] = 2.22144;
	
	// for slab
	global.field.N[0] = 0;
	global.field.N[1] = 32;
	global.field.N[2] = 32;
	global.field.N[3] = 32;
	
	cout << "Nx = " << global.field.N[1] << ", Ny = " << global.field.N[2] << " , Nz = " << global.field.N[3] << endl;
	
	
    // test_CFFF_slab();
    //test_FFFW_slab();
    // test_FFF_slab();
	test_SFF_slab();
	//test_SSF_slab();
	//test_SSS_slab();
	//test_ChFF_slab();
	//test_transpose_slab();
	
	//test_FFF_pencil();
	//test_SFF_pencil();
	//test_SSF_pencil();
	//test_SSS_pencil();
	
	//test_transpose_pencil();
	
	BasicIO::Finalize();
    MPI_Finalize();
	return 0;
}

//********************************** End of main.cc *******************************************


