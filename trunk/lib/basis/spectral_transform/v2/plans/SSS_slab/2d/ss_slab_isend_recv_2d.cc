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


#include "sf_slab_isend_recv_2d.h"


	
//*********************************************************************************************
SF_slab_Isend_Recv_2D::SF_slab_Isend_Recv_2D(int my_id, int numprocs, double num_iter, int Nx, int Nz): SpectralPlan_Slab_2D(my_id, numprocs, Nx, Nz)
{

	X_2d.resize(local_Nx, Nz/2+1);
	Xr_2d.resize(Nx, 2*local_Nz);

	//Initialize plans
	int Nx_dims[]={Nx};
	int Nz_dims[]={Nz};

	fftw_r2r_kind kind[1];

	kind[0]=FFTW_RODFT10;
	plan_sintr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, 2*local_Nz,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT10;
	plan_costr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, 2*local_Nz,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, 2*local_Nz,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	plan_icostr_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, 2*local_Nz,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    reinterpret_cast<DP*>(Xr_2d.data()), NULL,
	    2*local_Nz, 1,
	    kind, FFTW_PLAN_FLAG);


	plan_ft_r2c_z = FFTW_PLAN_MANY_DFT_R2C(1, Nz_dims, local_Nx,
		reinterpret_cast<DP*>(X_2d.data()), NULL,
		1, (Nz+2),
		reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), NULL,
		1, (Nz/2+1),
		FFTW_PLAN_FLAG);

	plan_ift_c2r_z = FFTW_PLAN_MANY_DFT_C2R(1, Nz_dims, local_Nx,
		reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), NULL,
		1, (Nz/2+1),
		reinterpret_cast<DP*>(X_2d.data()), NULL,
		1, (Nz+2),
		FFTW_PLAN_FLAG);

	Init_array();
	Evaluate_time_per_step("CF", num_iter);

	X_2d.free();
	Xr_2d.free();
}
//*********************************************************************************************

DP SF_slab_Isend_Recv_2D::f(string sincostr_option, int rx, int rz)
{
	DP k0 = 1;
	DP x,z;
	
	DP Lx=M_PI;
	DP Lz=2*M_PI;

	x = (rx+0.5)*Lx/Nx;
	z = rz*Lz/Nz;

	if (sincostr_option=="SF")
		return 4*sin(k0*x)*cos(k0*z);
	else if (sincostr_option=="CF")
		return 4*cos(k0*x)*cos(k0*z);
}
 
void SF_slab_Isend_Recv_2D::Init_array()
{
	Xr_2d=0;
	for (int rx=0; rx<Nx; rx++)
		for (int rz=0; rz<2*local_Nz; rz++)
			Xr_2d(rx, rz) = f("CF", rx,2*local_Nz_start + rz);
}

//*********************************************************************************************

void SF_slab_Isend_Recv_2D::Normalize(Array<complx,2> A)
{
	A /= (2.0 * DP(Nx) * DP(Nz));
}

//*************************************************************************************


void SF_slab_Isend_Recv_2D::Forward_transform(string sincostr_option, Array<DP,2> Ar, Array<complx,2> A)
{
	SinCostr_x(sincostr_option[0],Ar);

	Transpose(Ar,A);

	Zero_pad_last_plane(A);	// Zero_pad the row=Nz

	FT_r2c_z(A);

	Normalize(A);
}


//*********************************************************************************************

void SF_slab_Isend_Recv_2D::Inverse_transform(string sincostr_option, Array<complx,2> A, Array<DP,2> Ar)
{
	IFT_c2r_z(A);

	Zero_pad_last_plane(A);	// Zero_pad the row=Nz

	Transpose(A,Ar);

	ISinCostr_x(sincostr_option[0],Ar);
}

//Transpose
void SF_slab_Isend_Recv_2D::Transpose(Array<DP,2> Ar, Array<complx,2> A){
	Isend_x(Ar);
	Recv_z(A);
}

void SF_slab_Isend_Recv_2D::Transpose(Array<complx,2> A, Array<DP,2> Ar){
	Isend_z(A);
	Recv_x(Ar);
}


//Transforms
void SF_slab_Isend_Recv_2D::SinCostr_x(char sincostr_option, Array<DP,2> Ar)
{
	if (sincostr_option == 'S') {
		FFTW_EXECUTE_R2R(plan_sintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'X');
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_costr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SF_slab_Isend_Recv_2D::ISinCostr_x(char sincostr_option, Array<DP,2> Ar)
{
    if (sincostr_option == 'S') {
		ArrayShiftLeft(Ar, 'X');
		FFTW_EXECUTE_R2R(plan_isintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_icostr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SF_slab_Isend_Recv_2D::FT_r2c_z(Array<complx,2> A)
{
	FFTW_EXECUTE_DFT_R2C(plan_ft_r2c_z, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}


void SF_slab_Isend_Recv_2D::IFT_c2r_z(Array<complx,2> A)
{
	FFTW_EXECUTE_DFT_C2R(plan_ift_c2r_z, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

//********************************	End of four_tr.cc *****************************************


