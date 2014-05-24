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


#include "fff_slab_alltoall_3d.h"


	
//*********************************************************************************************
FFF_slab_Alltoall_3D::FFF_slab_Alltoall_3D(int my_id, int numprocs, int num_iter, int Nx, int Ny, int Nz): SpectralPlan_Slab_3D("FFF", my_id, numprocs, Nx, Ny, Nz)
{

	X_3d.resize(Nx, local_Ny, Nz/2+1);
	Xr_3d.resize(local_Nx, Ny, Nz+2);

	
	//Initialize plans
	int Nx_plan[]={Nx};
	int Nyz_plan[]={Ny,Nz};

	plan_ft_c2c_x = FFTW_PLAN_MANY_DFT(1, Nx_plan, (Nz/2+1)*local_Ny,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		(Nz/2+1)*local_Ny, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		(Nz/2+1)*local_Ny, 1,
		FFTW_FORWARD, FFTW_PLAN_FLAG);

	plan_ift_c2c_x = FFTW_PLAN_MANY_DFT(1, Nx_plan, (Nz/2+1)*local_Ny,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		(Nz/2+1)*local_Ny, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		(Nz/2+1)*local_Ny, 1,
		FFTW_BACKWARD, FFTW_PLAN_FLAG);

	plan_ft_r2c_yz = FFTW_PLAN_MANY_DFT_R2C(2, Nyz_plan, local_Nx,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, (Nz+2)*Ny,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, (Nz/2+1)*Ny,
		FFTW_PLAN_FLAG);

	plan_ift_c2r_yz = FFTW_PLAN_MANY_DFT_C2R(2, Nyz_plan, local_Nx,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, (Nz/2+1)*Ny,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, (Nz+2)*Ny,
		FFTW_PLAN_FLAG);

	Init_array();
	Evaluate_time_per_step(num_iter);

	X_3d.free();
	Xr_3d.free();
}
//*********************************************************************************************

DP FFF_slab_Alltoall_3D::f(int rx, int ry, int rz)
{
	DP k0 = 1;
	DP x,y,z;
	DP L=2*M_PI;

	x = rx*L/Nx;
	y = ry*L/Ny;
	z = rz*L/Nz;
	return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
}
 

void FFF_slab_Alltoall_3D::Init_array()
{
	Xr_3d=0;
	for (int rx=0; rx<local_Nx; rx++)
		for (int ry=0; ry<Ny; ry++)
			for (int rz=0; rz<Nz; rz++) {
				Xr_3d(rx, ry, rz) = f(local_Nx_start+rx,ry,rz);
			}
}

//*********************************************************************************************

void FFF_slab_Alltoall_3D::Normalize(Array<complx,3> A)
{
	A /= (DP(Nx) * DP(Ny) * DP(Nz));
}

//*************************************************************************************


void FFF_slab_Alltoall_3D::Forward_transform(Array<DP,3> Ar, Array<complx,3> A)
{

	Zero_pad_last_plane(Ar);	// Zero_pad the row=Nz

	FT_r2c_yz(Ar);

	Transpose(Ar,A);
	
	FT_c2c_x(A);

	Normalize(A);
}


void FFF_slab_Alltoall_3D::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar)
{

	IFT_c2c_x(A);

	Transpose(A,Ar);

	IFT_c2r_yz(Ar);

	Zero_pad_last_plane(Ar); 
}


void FFF_slab_Alltoall_3D::Transpose(Array<DP,3> Ar, Array<complx,3> A)
{
	Alltoall(Ar,A,YZ);
}

void FFF_slab_Alltoall_3D::Transpose(Array<complx,3> A, Array<DP,3> Ar)
{
	Alltoall(A,Ar,XZ);
}


void FFF_slab_Alltoall_3D::FT_r2c_yz(Array<DP,3> A)
{
	FFTW_EXECUTE_DFT_R2C(plan_ft_r2c_yz, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}


void FFF_slab_Alltoall_3D::IFT_c2r_yz(Array<DP,3> A)
{
	FFTW_EXECUTE_DFT_C2R(plan_ift_c2r_yz, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

void FFF_slab_Alltoall_3D::FT_c2c_x(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT(plan_ft_c2c_x, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}

void FFF_slab_Alltoall_3D::IFT_c2c_x(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT(plan_ift_c2c_x, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}

//********************************	End of four_tr.cc *****************************************


