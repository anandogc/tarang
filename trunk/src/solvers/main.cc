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
 * Usage: mpirun -np numprocs tarang datadir num_p_hor
 *	@author  M. K. Verma
 *	@version 4.0 MPI
 *	@date Sept 2008
 */

#include "main.h"

Universal *universal = NULL;

Global global;

SpectralTransform spectralTransform;

Uniform<Real> SPECrand;

FluidIO_incompress fluidIO_incompress;

//*********************************************************************************************


int main(int argc, char** argv)
{

	MPI_Init(&argc, &argv);


	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);

	global.Parse(argc, argv);
	global.Read();

	global.Process_basic_vars();
	
    BasicIO::Initialize(global.io.data_dir, global.mpi.num_p_rows, global.mpi.num_p_cols); 

	if (global.program.basis_type == "FFFW"){
		// if (global.program.decomposition == "SLAB")
		// 	universal=new FFFW_SLAB;
	}

	else if (global.program.basis_type == "FFF"){
		// if (global.program.decomposition == "SLAB")
		// 	universal=new FFF_SLAB;
		if (global.program.decomposition == "PENCIL")
			universal=new FFF_PENCIL;
	}
	
   else if (global.program.basis_type == "SFF"){
		// if (global.program.decomposition == "SLAB")
			// universal=new SFF_SLAB;
		if (global.program.decomposition == "PENCIL")
			universal=new SFF_PENCIL;
	}
   else if (global.program.basis_type == "SSF"){
		// if (global.program.decomposition == "SLAB")
			// universal=new SSF_SLAB;
		if (global.program.decomposition == "PENCIL")
			universal=new SSF_PENCIL;
	}
	else if (global.program.basis_type == "SSS"){
		// if (global.program.decomposition == "SLAB")
			// universal=new SSS_SLAB;	
		if (global.program.decomposition == "PENCIL")
			universal=new SSS_PENCIL;
	}
/*    else if (global.program.basis_type == "ChFF"){
		if (global.program.decomposition == "SLAB")
			universal=new ChFF_SLAB;
	} */
/*    else if (global.program.basis_type == "CFFF"){
//		if (global.program.decomposition == "SLAB")
//			universal=new CFFF_SLAB;
	}*/
	
	if (universal == NULL){
			if(master) cerr << "program.basis_type must be one of FFF, SFF, SSF, SSS. And program.decomposition must be PENCIL"<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
	}

	global.Process_advanced_vars();
	global.Print();

	
	
	Correlation::Initialize();
	
	time_t start;
	time_t end;

	time (&start);

	
	SPECrand.seed((unsigned int)time(0));				// Initialize random number seed

	if (global.program.kind == "FLUID_INCOMPRESS")
		Ifluid_main();

	else if (global.program.kind == "SCALAR_INCOMPRESS")
		Iscalar_main();
	
	else if (global.program.kind == "RBC")
		Iscalar_main();
	
	else if (global.program.kind == "STRATIFIED")
		Iscalar_main();

	else if (global.program.kind == "MRBC")
		MRBC_main();

	else if (global.program.kind == "MHD_INCOMPRESS")
		IMHD_main();
	
	else if (global.program.kind == "KEPLERIAN")
		IMHD_main();
	
	else if (global.program.kind == "MHD_SCALAR_INCOMPRESS")
		IMHDscalar_main();
	
	else if (global.program.kind == "MHD_ASTRO_INCOMPRESS")
		IMHDastro_main();
	
	else if (global.program.kind == "GP")
		GP_main();

	time(&end);

	Real vm, rss;
	process_mem_usage(vm, rss);

	Real dif = difftime(end, start);
	if (my_id == master_id){
		cout << endl << "Memory usage by Process " << my_id
		<< "   VM: " << vm << "KB  RSS: " << rss << "KB" << endl;

		cout << endl << "Approximate total memory usage,  VM:"
		<< vm*numprocs << "KB   RSS:" << rss*numprocs << endl;

		cout << endl << "PROGRAM TERMINATING HERE" << endl << endl
			<< "TOTAL TIME ELAPSED: " << dif << " sec" << endl;
	}

	BasicIO::Finalize();
	MPI_Finalize();

	return 0;
}



//********************************** End of main.cc *******************************************


