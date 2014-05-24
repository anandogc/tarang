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
 * along with Tarang-2; if 2D Not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */
/*! \file fourier.cc 
 * 
 * @sa fourier.h
 * 
 * @author  M. K. Verma
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		ArrayIFFT
 */ 


#include "sss_slab_alltoall_3d.h"


	
//*********************************************************************************************
SSS_slab_Alltoall_3D::SSS_slab_Alltoall_3D(int my_id, int numprocs, int num_iter, int Nx, int Ny, int Nz): SpectralPlan_Slab_3D("SSS", my_id, numprocs, Nx, Ny, Nz)
{

	X_3d.resize(local_Nx, Ny, Nz/2);	   // temp array
	Xr_3d.resize(Nx, local_Ny, Nz);


	//Initialize plans
	int Nx_dims[]={Nx};
	int Ny_dims[]={Ny};
	int Nz_dims[]={Nz};

	fftw_r2r_kind kind[1];

	//Sin x
	kind[0]=FFTW_RODFT10;
	plan_sintr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Nz*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Nz*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);

	//Cos x
	kind[0]=FFTW_REDFT10;
	plan_costr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Nz*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);

	
	kind[0]=FFTW_REDFT01;
	plan_icostr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Nz*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);


	//Sin y
	cout << "SSS_slab_Alltoall_3D::SSS_slab_Alltoall_3D: " << X_3d.shape() << " " << Nx << " " << Ny << " " << Nz << endl;
	kind[0]=FFTW_RODFT10;
	plan_sintr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    kind, FFTW_PLAN_FLAG);
	
	//Cos y
	kind[0]=FFTW_REDFT10;
	plan_costr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	plan_icostr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    Nz, 1,
	    kind, FFTW_PLAN_FLAG);

	//Sin z
	kind[0]=FFTW_RODFT10;
	plan_sintr_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, Ny*local_Nx,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, Ny*local_Nx,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    kind, FFTW_PLAN_FLAG);
	
	//Cos z
	kind[0]=FFTW_REDFT10;
	plan_costr_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, Ny*local_Nx,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	plan_icostr_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, Ny*local_Nx,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    1, Nz,
	    kind, FFTW_PLAN_FLAG);

	Init_array();
	Evaluate_time_per_step("CSC", num_iter);

	X_3d.free();
	Xr_3d.free();
}
//*********************************************************************************************

void SSS_slab_Alltoall_3D::Init_array()
{
	DP k0 = 1;
	DP x,y,z;
	
	DP Lx=M_PI;
	DP Ly=M_PI;
	DP Lz=M_PI;

	for (int rx=0; rx<Nx; rx++)
		for (int ry=0; ry<local_Ny; ry++)
			for (int rz=0; rz<Nz; rz++)
			{
				x = (rx+0.5)*Lx/Nx;
				y = (local_Ny_start+ry+0.5)*Ly/Ny;
				z = (rz+0.5)*Lz/Nz;
				Xr_3d(rx, ry, rz) = 8*cos(k0*x)*sin(k0*y)*cos(k0*z);
			}
}

//*************************************************************************************

void SSS_slab_Alltoall_3D::Normalize(Array<complx,3> A)
{
	A /=(8.0 * DP(Nx) * DP(Ny) * DP(Nz));
}



void SSS_slab_Alltoall_3D::Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A)
{

	SinCostr_x(sincostr_option[0],Ar);

	Transpose(Ar,A);

	SinCostr_y(sincostr_option[1],A);

	SinCostr_z(sincostr_option[2],A);

	Normalize(A);
}


void SSS_slab_Alltoall_3D::Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar)
{
	ISinCostr_z(sincostr_option[2],A);

	ISinCostr_y(sincostr_option[1],A);

	Transpose(A,Ar);

	ISinCostr_x(sincostr_option[0], Ar);
}



//Transpose
void SSS_slab_Alltoall_3D::Transpose(Array<DP,3> Ar, Array<complx,3> A) {
	Alltoall(Ar,A,XZ);
}

void SSS_slab_Alltoall_3D::Transpose(Array<complx,3> A, Array<DP,3> Ar) {
	Alltoall(A,Ar,YZ);
}


void SSS_slab_Alltoall_3D::SinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
	if (sincostr_option == 'S') {
		FFTW_EXECUTE_R2R(plan_sintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'X');
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_costr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SSS_slab_Alltoall_3D::ISinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
    if (sincostr_option == 'S') {
		ArrayShiftLeft(Ar, 'X');
		FFTW_EXECUTE_R2R(plan_isintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_icostr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SSS_slab_Alltoall_3D::SinCostr_y(char sincostr_option, Array<complx,3> Ar)
{
	if (sincostr_option == 'S') {
		for (int lx=0; lx<local_Nx; lx++)
			FFTW_EXECUTE_R2R(plan_sintr_y, reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()), reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()));
		ArrayShiftRight(Ar, 'Y');
	}
	
	else if (sincostr_option == 'C')
		for (int lx=0; lx<local_Nx; lx++)
			FFTW_EXECUTE_R2R(plan_costr_y, reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()), reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()));
}

void SSS_slab_Alltoall_3D::ISinCostr_y(char sincostr_option, Array<complx,3> Ar)
{
    if (sincostr_option == 'S') {
		ArrayShiftLeft(Ar, 'Y');
		for (int lx=0; lx<local_Nx; lx++)
			FFTW_EXECUTE_R2R(plan_isintr_y, reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()), reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()));
	}
	
	else if (sincostr_option == 'C')
		for (int lx=0; lx<local_Nx; lx++)
			FFTW_EXECUTE_R2R(plan_icostr_y, reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()), reinterpret_cast<DP*>(Ar(lx,Range::all(),Range::all()).data()));
}



void SSS_slab_Alltoall_3D::SinCostr_z(char sincostr_option, Array<complx,3> Ar)
{
	if (sincostr_option == 'S') {
		FFTW_EXECUTE_R2R(plan_sintr_z, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'Z');
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_costr_z, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SSS_slab_Alltoall_3D::ISinCostr_z(char sincostr_option, Array<complx,3> Ar)
{
    if (sincostr_option == 'S') {
		ArrayShiftLeft(Ar, 'Z');
		FFTW_EXECUTE_R2R(plan_isintr_z, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_icostr_z, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

//********************************	End of four_tr.cc *****************************************
