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

/*! \file Ifluid_main.cc 
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 


#include "main.h"
#include "transpose_pencil.h"

//*********************************************************************************************					
void test_FFF_pencil(){
	global.program.basis_type = "FFF";
	global.program.decomposition = "PENCIL";
	global.program.sincostr_switch = "FFF";
	
	global.field.shape_complex_array = shape(local_Ny_vert,local_Nz_hor,Nx);
	global.field.shape_real_array = shape(local_Ny_hor,Nz+2,local_Nx_vert);
	
	
	global.Process_global_vars_basic();
	universal=new FFF_PENCIL;
	global.Process_global_vars_advanced();
	
	//global.Print();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		//cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	//	result += test_1dtransform();
	
	/*
	 
	spectralTransform.Init("FFF", "PENCIL", Nx, Ny, Nz, global.mpi.num_p_hor);Array<DP,3> Ax, Ay, Az;
	
	Ax.resize(spectralTransform.local_Ny_vert,2*spectralTransform.local_Nz_hor,Nx);
	Ay.resize(Ny,2*spectralTransform.local_Nz_hor,spectralTransform.local_Nx_vert);
	Az.resize(spectralTransform.local_Ny_hor,2*(Nz/2+1),spectralTransform.local_Nx_vert);
	

	int index=0;
	for (int ly=0; ly<Ax.extent(0); ly++)
		for (int lz=0; lz<Ax.extent(1); lz++)
			for (int lx=0; lx<Ax.extent(2); lx++){
				//Ax(ly,lz,lx)=complx(index+my_id*1000, index+1+my_id*1000);
				Ax(ly,lz,lx)=index+my_id*1000;
				index+=1;
			}
	if (my_id==0)
		cout << "Before Transpose" << endl;
	ArrayOps::Print_array_all_procs(Ax);
	
	spectralTransform.Transpose_array_real_XY_PENCIL(Ax, Ay);
	
	if (my_id==0)
		cout << "After Transpose" << endl;
	
	ArrayOps::Print_array_all_procs(Ay);*/
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,0,1, A, complx(0,-1));
	universal->Assign_spectral_field(-1,0,1, A, complx(0,1));
	
	test_transform();
}


void test_SFF_pencil(){
	global.program.basis_type = "SFF";
	global.program.decomposition = "PENCIL";
	global.program.sincostr_switch = "SFF";
	
	
	global.Process_global_vars_basic();
	universal=new SFF_PENCIL;
	global.Process_global_vars_advanced();
	
	
	
	//global.Print();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		//cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	//	result += test_1dtransform();
	
	/*
	spectralTransform.Init("SFF", "PENCIL", Nx, Ny, Nz, global.mpi.num_p_hor);
	
	Array<DP,3> Ax;
	Array<complx,3> Ay;
	Array<complx,3> Az;
	
	Ax.resize(spectralTransform.shape_x_array_3d);
	Ay.resize(spectralTransform.shape_y_array_3d);
	Az.resize(spectralTransform.shape_z_array_3d);
	

	int index=0;
	for (int ly=0; ly<Ay.extent(0); ly++)
		for (int lz=0; lz<Ay.extent(1); lz++)
			for (int lx=0; lx<Ay.extent(2); lx++){
				Ay(ly,lz,lx)=complx(index+my_id*1000, index+1+my_id*1000);
				index+=2;
			}
	if (my_id==0)
		cout << "Before Transpose" << endl;
	ArrayOps::Print_array_all_procs(Ay);
	
	spectralTransform.Transpose_array_real_XY_PENCIL(Ay, Ax);
	
	if (my_id==0)
		cout << "After Transpose" << endl;
	
	ArrayOps::Print_array_all_procs(Ax);*/
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,0,1, A, complx(1,0));
	
	test_transform();
}


void test_SSF_pencil(){
	global.program.basis_type = "SSF";
	global.program.decomposition = "PENCIL";
	global.program.sincostr_switch = "SCF";
	
	
	global.Process_global_vars_basic();
	universal=new SSF_PENCIL;
	global.Process_global_vars_advanced();
	
	
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		//cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}

	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,1,1, A, complx(1,0));
	
	test_transform();
}

void test_SSS_pencil(){
	global.program.basis_type = "SSS";
	global.program.decomposition = "PENCIL";
	global.program.sincostr_switch = "SSS";
	
	
	global.Process_global_vars_basic();
	universal=new SSS_PENCIL;
	global.Process_global_vars_advanced();
	
	
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		//cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,1,1, A, 1.0);
	
	test_transform();

}

void test_transpose_pencil()
{
	 global.Process_global_vars_basic();
	
	 spectralTransform.Init("SSS", "PENCIL", Nx, Ny, Nz, global.mpi.num_p_hor);
	 
	 Array<DP,3> Ax;
	 Array<complx,3> Ay;
	 Array<complx,3> Az;
	 
	 Ax.resize(spectralTransform.shape_x_array_3d);
	 Ay.resize(spectralTransform.shape_y_array_3d);
	 Az.resize(spectralTransform.shape_z_array_3d);
	 
	 
	 int index=0;
	 for (int ly=0; ly<Ay.extent(0); ly++)
		 for (int lz=0; lz<Ay.extent(1); lz++)
			 for (int lx=0; lx<Ay.extent(2); lx++){
				 Ay(ly,lz,lx)=complx(index+my_id*1000, index+1+my_id*1000);
				 index+=2;
			 }
	
	 if (my_id==0)
		 cout << "Before Transpose" << endl;
	 ArrayOps::Print_array_all_procs(Ay);
	 
	MPI_Barrier(MPI_COMM_WORLD);
	
	spectralTransform.Transpose_array(Ay, Az, spectralTransform.YZ);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	 if (my_id==0)
		 cout << "After Transpose" << endl;
	 MPI_Barrier(MPI_COMM_WORLD);
	
	 ArrayOps::Print_array_all_procs(Az);
}

//********************************** End of Ifluid_main.cc ************************************



