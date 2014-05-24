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


#include "ssf_slab_alltoall_3d.h"


	
//*********************************************************************************************
SSF_slab_Alltoall_3D::SSF_slab_Alltoall_3D(int my_id, int numprocs, int num_iter, int Nx, int Ny, int Nz): SpectralPlan_Slab_3D("SSF", my_id, numprocs, Nx, Ny, Nz)
{

	X_3d.resize(local_Nx, Ny, Nz/2+1);	   // temp array
	Xr_3d.resize(Nx, local_Ny, Nz+2);


	//Initialize plans
	int Nx_dims[]={Nx};
	int Ny_dims[]={Ny};
	int Nz_dims[]={Nz};

	fftw_r2r_kind kind[1];

	//Sin x
	kind[0]=FFTW_RODFT10;
	plan_sintr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);

	//Cos x
	kind[0]=FFTW_REDFT10;
	plan_costr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);

	
	kind[0]=FFTW_REDFT01;
	plan_icostr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);


	//Sin y
	kind[0]=FFTW_RODFT10;
	plan_sintr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    kind, FFTW_PLAN_FLAG);
	
	//Cos y
	kind[0]=FFTW_REDFT10;
	plan_costr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	plan_icostr_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Nz,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2), 1,
	    kind, FFTW_PLAN_FLAG);


	//FT z
	plan_ft_r2c_z = FFTW_PLAN_MANY_DFT_R2C(1, Nz_dims, Ny*local_Nx,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, (Nz+2),
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, (Nz/2+1),
		FFTW_PLAN_FLAG);

	plan_ift_c2r_z = FFTW_PLAN_MANY_DFT_C2R(1, Nz_dims, Ny*local_Nx,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, (Nz/2+1),
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, (Nz+2),
		FFTW_PLAN_FLAG);

	Init_array();
	Evaluate_time_per_step("CCF", num_iter);

	X_3d.free();
	Xr_3d.free();
}
//*********************************************************************************************

void SSF_slab_Alltoall_3D::Init_array()
{
	DP k0 = 1;
	DP x,y,z;
	
	DP Lx=M_PI;
	DP Ly=M_PI;
	DP Lz=2*M_PI;

	for (int rx=0; rx<Nx; rx++)
		for (int ry=0; ry<local_Ny; ry++)
			for (int rz=0; rz<Nz; rz++)
			{
				x = (rx+0.5)*Lx/Nx;
				y = (local_Ny_start+ry+0.5)*Ly/Ny;
				z = rz*Lz/Nz;
				Xr_3d(rx, ry, rz) = 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
			}
}

//*************************************************************************************

void SSF_slab_Alltoall_3D::Normalize(Array<complx,3> A)
{
	A /=(4.0 * DP(Nx) * DP(Ny) * DP(Nz));
}



void SSF_slab_Alltoall_3D::Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A)
{
	Zero_pad_last_plane(Ar);	// Zero_pad the plane Nz

	SinCostr_x(sincostr_option[0],Ar);

	Transpose(Ar,A);

	SinCostr_y(sincostr_option[1],A);

	FT_r2c_z(A);

	Normalize(A);
}


void SSF_slab_Alltoall_3D::Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar)
{
	IFT_c2r_z(A);

	ISinCostr_y(sincostr_option[1],A);

	Transpose(A,Ar);

	ISinCostr_x(sincostr_option[0], Ar);

	Zero_pad_last_plane(Ar); 
}



//Transpose
void SSF_slab_Alltoall_3D::Transpose(Array<DP,3> Ar, Array<complx,3> A) {
	Alltoall(Ar,A,XZ);
}

void SSF_slab_Alltoall_3D::Transpose(Array<complx,3> A, Array<DP,3> Ar) {
	Alltoall(A,Ar,YZ);
}


void SSF_slab_Alltoall_3D::SinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
	if (sincostr_option == 'S') {
		FFTW_EXECUTE_R2R(plan_sintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'X');
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_costr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SSF_slab_Alltoall_3D::ISinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
    if (sincostr_option == 'S') {
		ArrayShiftLeft(Ar, 'X');
		FFTW_EXECUTE_R2R(plan_isintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_icostr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SSF_slab_Alltoall_3D::SinCostr_y(char sincostr_option, Array<complx,3> Ar)
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

void SSF_slab_Alltoall_3D::ISinCostr_y(char sincostr_option, Array<complx,3> Ar)
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

void SSF_slab_Alltoall_3D::FT_r2c_z(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_R2C(plan_ft_r2c_z, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}


void SSF_slab_Alltoall_3D::IFT_c2r_z(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_C2R(plan_ift_c2r_z, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

//********************************	End of four_tr.cc *****************************************
