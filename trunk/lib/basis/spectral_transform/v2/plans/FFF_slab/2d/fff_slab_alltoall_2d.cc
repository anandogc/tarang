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


#include "fff_slab_alltoall_2d.h"


	
//*********************************************************************************************
FFF_slab_Alltoall_2D::FFF_slab_Alltoall_2D(int my_id, int numprocs, int num_iter, int Nx, int Nz): SpectralPlan_Slab_2D("FF", my_id, numprocs, Nx, Nz)
{

	X_2d.resize(Nx, local_Nz);
	Xr_2d.resize(local_Nx, Nz+2);

	//Initialize plans
	int Nx_dims[]={Nx};
	int Nz_dims[]={Nz};

	plan_ft_c2c_x = FFTW_PLAN_MANY_DFT(1, Nx_dims, local_Nz,
		reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), NULL,
		local_Nz, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), NULL,
		local_Nz, 1,
		FFTW_FORWARD, FFTW_PLAN_FLAG);


	plan_ift_c2c_x = FFTW_PLAN_MANY_DFT(1, Nx_dims,local_Nz,
		reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), NULL,
		local_Nz, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), NULL,
		local_Nz, 1,
		FFTW_BACKWARD, FFTW_PLAN_FLAG);

	plan_ft_r2c_z = FFTW_PLAN_MANY_DFT_R2C(1, Nz_dims, local_Nx,
		reinterpret_cast<DP*>(Xr_2d.data()), NULL,
		1, (Nz+2),
		reinterpret_cast<FFTW_COMPLEX*>(Xr_2d.data()), NULL,
		1, (Nz/2+1),
		FFTW_PLAN_FLAG);

	plan_ift_c2r_z = FFTW_PLAN_MANY_DFT_C2R(1, Nz_dims, local_Nx,
		reinterpret_cast<FFTW_COMPLEX*>(Xr_2d.data()), NULL,
		1, (Nz/2+1),
		reinterpret_cast<DP*>(Xr_2d.data()), NULL,
		1, (Nz+2),
		FFTW_PLAN_FLAG);

	Init_array();
	Evaluate_time_per_step(num_iter);

	X_2d.free();
	Xr_2d.free();
}
//*********************************************************************************************

DP FFF_slab_Alltoall_2D::f(int rx, int rz)
{
	DP k0 = 1;
	DP x,z;
	DP L=2*M_PI;

	x = rx*L/Nx;
	z = rz*L/Nz;
	return 4*sin(k0*x)*cos(k0*z);
}
 

void FFF_slab_Alltoall_2D::Init_array()
{
	for (int rx=0; rx<local_Nx; rx++)
		for (int rz=0; rz<Nz; rz++){
			Xr_2d(rx, rz) = f(local_Nx_start+rx,rz);
		}
	Zero_pad_last_plane(Xr_2d);
}

//*********************************************************************************************

void FFF_slab_Alltoall_2D::Normalize(Array<complx,2> A)
{
	A /= (DP(Nx) * DP(Nz));
}

//*************************************************************************************

void FFF_slab_Alltoall_2D::Forward_transform(Array<DP,2> Ar, Array<complx,2> A)
{
	Zero_pad_last_plane(Ar);	// Zero_pad the row=Nz

	FT_r2c_z(Ar);

	Transpose(Ar,A);
	
	FT_c2c_x(A);

	Normalize(A);
}


void FFF_slab_Alltoall_2D::Inverse_transform(Array<complx,2> A, Array<DP,2> Ar)
{
	Zero_pad_last_plane(Ar); 

	IFT_c2c_x(A);

	Transpose(A,Ar);

	IFT_c2r_z(Ar);

	Zero_pad_last_plane(Ar);
}


void FFF_slab_Alltoall_2D::Transpose(Array<DP,2> Ar, Array<complx,2> A)
{
	Alltoall(Ar,A,Z);
}

void FFF_slab_Alltoall_2D::Transpose(Array<complx,2> A, Array<DP,2> Ar)
{
	Alltoall(A,Ar,X);
}


void FFF_slab_Alltoall_2D::FT_r2c_z(Array<DP,2> A)
{
	FFTW_EXECUTE_DFT_R2C(plan_ft_r2c_z, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}


void FFF_slab_Alltoall_2D::IFT_c2r_z(Array<DP,2> A)
{
	FFTW_EXECUTE_DFT_C2R(plan_ift_c2r_z, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

void FFF_slab_Alltoall_2D::FT_c2c_x(Array<complx,2> A)
{
	FFTW_EXECUTE_DFT(plan_ft_c2c_x, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}

void FFF_slab_Alltoall_2D::IFT_c2c_x(Array<complx,2> A)
{
	FFTW_EXECUTE_DFT(plan_ift_c2c_x, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}

//********************************	End of four_tr.cc *****************************************


