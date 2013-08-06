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

Uniform<DP> SPECrand;

int main(int argc, char** argv)
{
  	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);
	
    global.mpi.master_id=0;
	global.Global_Parse(argc, argv, true);
    

	char cwd[256];
    getcwd(cwd, sizeof(cwd));

	//Set up global
	global.Global_init_default();
	
	
	global.io.data_dir=cwd;
	global.time.now=0;
	global.program.basis_type = "FFF";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "FFF";

	global.field.N[0] = 0;
	global.field.N[1] = 16;
	global.field.N[2] = 4;
	global.field.N[3] = 8;

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
	
	
	
	global.Process_global_vars_basic();
	universal=new FFF_SLAB;
	global.Process_global_vars_advanced();

	//Set up BasicIO
	BasicIO::data_in_folder="/test/time_"+To_string(global.time.now)+"/";
	BasicIO::data_out_folder="/test/";
	
	BasicIO::Initialize();


	Array<complx, 3> A(local_Ny,Nz/2+1,Nx);
	Array<complx, 3> B(2*local_Ny,2*(Nz/2+1),2*Nx);

	firstIndex i;
	secondIndex j;
	thirdIndex k;


	//Initialize A
	int dim1=A.extent(0);
	int dim2=A.extent(1);
	int dim3=A.extent(2);

	real(A) = i*2*dim2*dim3 + j*2*dim3 + k*2 + global.mpi.my_id*2*dim1*dim2*dim3;
	imag(A) = i*2*dim2*dim3 + j*2*dim3 + k*2+1 + global.mpi.my_id*2*dim1*dim2*dim3;

	bool array_matches;

	//Testing pair 'Write_complex_array', 'Read_complex_array'
	cout << "my_id: " << my_id <<": Testing pair 'Write(full)', 'Read(full)': \n";

	cout << "BasicIO::full: \n" << BasicIO::full;
	B=complx(0,0);
	BasicIO::Write(A, "test_BasicIO_full", BasicIO::full);
	BasicIO::Read(B, "test_BasicIO_full", BasicIO::reduced);
	BasicIO::Write(B, "test_BasicIO_full", BasicIO::full);

	array_matches=true;

	/*for (int i=0; i<local_Ny; i++)
		for (int j=0; j<Nz/2+1; j++)
			for (int k=0; k<Nx; k++)
				if (A(i,j,k)!=B(i,j,k)){
					array_matches = false;
					//cout << "A("<<i<<","<<j<<","<<k<<") = " << A(i,j,k) << endl;
				}

	if (array_matches)
		cout << "passed" << endl;
	else
		cout << "failed" << endl;*/
	/*
	//Testing pair 'Write_kz0plane', 'Read_kz0plane'
	cout << "my_id: " << my_id <<": Testing pair 'Write_kz0plane', 'Read_kz0plane': ";
	B=complx(0,0);
	BasicIO::Write_kz0plane(A, "test_BasicIO_plane");
	BasicIO::Read_kz0plane(B, "test_BasicIO_plane");

	array_matches=true;

	for (int i=0; i<local_Nx; i++)
		for (int j=0; j<Ny; j++)
			for (int k=0; k<1; k++)
				if (A(i,j,k)!=B(i,j,k))
					array_matches = false;

	if (array_matches)
		cout << "passed" << endl;
	else
		cout << "failed" << endl;
	
	//Testing pair 'Write_kz0plane', 'Read_kz0plane'
	//cout << "my_id: " << my_id <<": Testing pair 'Write_reduced_array', 'Read_reduced_array': ";
	
	B=complx(0,0);
		//BasicIO::Write_complex_array(A, "test_BasicIO_reduced");
	
		//return 0;
	
	BasicIO::Read_reduced_complex_array(B, "test_BasicIO_reduced");
	
	array_matches=true;
	
	cout << "reduced = " << global.io.N_in_reduced[1] << " " << global.io.N_in_reduced[2] << " " << global.io.N_in_reduced[3] << endl;
	
	if (numprocs==1) {
		for (int lx=0; lx<Nx; lx++)
			for (int ly=0; ly<Ny; ly++)
				for (int lz=0; lz<=(Nz/2); lz++) {
						//cout << "lx,ly,lz " << lx << " " << ly << " "<< lz << " " << ((lx<global.io.N_in_reduced[1]/2) && (ly<global.io.N_in_reduced[2]/2) && (lz <= global.io.N_in_reduced[3]/2)) << endl;
					
					if ( (lx<global.io.N_in_reduced[1]/2) && (ly<global.io.N_in_reduced[2]/2) && (lz <= global.io.N_in_reduced[3]/2) ) {
							//cout << "Here1" << endl;
						if (B(lx,ly,lz) != A(lx,ly,lz)) 
							cout << "ERROR: MISMATCH at lx,ly,lz = " << lx << " " << ly << " "<< lz << " "  << B(lx,ly,lz) << " " << A(lx,ly,lz) << endl;
					}
					
					else if ( (lx<global.io.N_in_reduced[1]/2) && (ly>= (Ny-global.io.N_in_reduced[2]/2)) && (lz <= global.io.N_in_reduced[3]/2) ) {
						if (B(lx,ly,lz) != A(lx,global.io.N_in_reduced[2]/2+ly-(Ny-global.io.N_in_reduced[2]/2),lz)) 
							cout << "ERROR: MISMATCH at lx,ly,lz = " << lx << " " << ly << " "<< lz << " "  << B(lx,ly,lz) << " " << A(lx,ly-(Ny-global.io.N_in_reduced[2]/2),lz)  << endl;
					}
					
					else if ( (lx>= (Nx-global.io.N_in_reduced[1]/2)) && (ly<global.io.N_in_reduced[2]/2) && (lz <= global.io.N_in_reduced[3]/2) ) {
						if (B(lx,ly,lz) != A(global.io.N_in_reduced[1]/2+lx-(Nx-global.io.N_in_reduced[1]/2),ly,lz)) 
							cout << "ERROR: MISMATCH at lx,ly,lz = " << lx << " " << ly << " "<< lz << " "  << B(lx,ly,lz) << " " << A(lx-(Nx-global.io.N_in_reduced[1]/2),ly,lz) << endl;
					}
					
					else if ( (lx>= (Nx-global.io.N_in_reduced[1]/2)) && (ly>= (Ny-global.io.N_in_reduced[2]/2)) && (lz <= global.io.N_in_reduced[3]/2) ) {
						if (B(lx,ly,lz) != A(global.io.N_in_reduced[1]/2+lx-(Nx-global.io.N_in_reduced[1]/2),global.io.N_in_reduced[2]/2+ly-(Ny-global.io.N_in_reduced[2]/2),lz)) 
							cout << "ERROR: MISMATCH at lx,ly,lz = " << lx << " " << ly << " "<< lz << " "  << B(lx,ly,lz) << " " << A(lx-(Nx-global.io.N_in_reduced[1]/2),ly-(Ny-global.io.N_in_reduced[2]/2),lz) << endl;					
					}
					
					else {
						if ( B(lx,ly,lz) != complx(0,0)) 
							cout << "ERROR: NONZERO B at lx,ly,lz = " << lx << " " << ly << " "<< lz << B(lx,ly,lz) << endl;
					}
				}
	}
	else if (numprocs==4) {
		if (my_id==0 || my_id==3){
			for (int i=0; i<local_Nx; i++)
				for (int j=0; j<Ny; j++)
					for (int k=0; k<global.io.N_in_reduced[3]/2+1; k++){
						if ( (j<global.io.N_in_reduced[2]/2 || j>=Ny-global.io.N_in_reduced[2]/2) ){
							if (A(i,j,k)!=B(i,j,k)){
								cout << "F1: (" << i << "," << j << "," << k << "): " << A(i,j,k) << " " << B(i,j,k) << endl;
								array_matches = false;
							}
						}
						else
							if (B(i,j,k) != complx(0,0)){
								cout << "F2: (" << i << "," << j << "," << k << "): " << A(i,j,k) << " " << B(i,j,k) << endl;
								array_matches = false;
							}
					}
		}
		else
			for (int i=0; i<local_Nx; i++)
				for (int j=0; j<Ny; j++)
					for (int k=0; k<Nz/2+1; k++)
						if (B(i,j,k)!=complx(0,0))
							array_matches=false;
	}
	if (array_matches)
		cout << "passed" << endl;
	else
		cout << "failed" << endl;

	if (global.program.decomposition == "PENCIL") {
		MPI_Comm_free(&global.mpi.MPI_COMM_VERT_SEGMENT);
		MPI_Comm_free(&global.mpi.MPI_COMM_HOR_SEGMENT);
	}*/

	BasicIO::Finalize();
	MPI_Finalize();
	
	return 0;
}



//********************************** End of main.cc *******************************************
