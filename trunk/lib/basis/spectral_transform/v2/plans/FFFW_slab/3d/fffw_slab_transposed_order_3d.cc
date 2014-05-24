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


#include "fffw_slab_transposed_order_3d.h"


	
//*********************************************************************************************
FFFW_slab_transposed_order_3D::FFFW_slab_transposed_order_3D(int my_id, int numprocs, int num_iter, int Nx, int Ny, int Nz): SpectralPlan_Slab_3D("FFFW", my_id, numprocs, Nx, Ny, Nz)
{

	X_3d.resize(local_Nx, Ny, Nz/2+1);
	Xr_3d.resize(local_Ny, Nx, Nz+2);

	fftw_mpi_init();

	//Initialize plans
	plan_ft_r2c_xyz = FFTW_MPI_PLAN_DFT_R2C_3D(Ny, Nx, Nz,
		reinterpret_cast<DP*>(Xr_3d.data()), reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()),
		MPI_COMM_WORLD, FFTW_PLAN_FLAG | FFTW_MPI_TRANSPOSED_OUT);

	plan_ift_c2r_xyz = FFTW_MPI_PLAN_DFT_C2R_3D(Ny, Nx, Nz,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), reinterpret_cast<DP*>(Xr_3d.data()),
		MPI_COMM_WORLD, FFTW_PLAN_FLAG | FFTW_MPI_TRANSPOSED_IN);


	Init_array();
	Evaluate_time_per_step(num_iter);

	X_3d.free();
	Xr_3d.free();
}
//*********************************************************************************************

DP FFFW_slab_transposed_order_3D::f(int rx, int ry, int rz)
{
	DP k0 = 1;
	DP x,y,z;
	DP L=2*M_PI;

	x = rx*L/Nx;
	y = ry*L/Ny;
	z = rz*L/Nz;
	return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
}
 

void FFFW_slab_transposed_order_3D::Init_array()
{
	for (int rx=0; rx<local_Nx; rx++)
		for (int ry=0; ry<Ny; ry++)
			for (int rz=0; rz<Nz; rz++) {
				Xr_3d(rx, ry, rz) = f(local_Nx_start+rx,ry,rz);
			}
}

//*********************************************************************************************

void FFFW_slab_transposed_order_3D::Normalize(Array<complx,3> A)
{
	A /= (DP(Nx) * DP(Ny) * DP(Nz));
}

//*************************************************************************************


void FFFW_slab_transposed_order_3D::Forward_transform(Array<DP,3> Ar, Array<complx,3> A)
{

	Zero_pad_last_plane(Ar);	// Zero_pad the row=Nz

	FFTW_MPI_EXECUTE_DFT_R2C(plan_ft_r2c_xyz, Ar.data(), reinterpret_cast<FFTW_COMPLEX*>(A.data()));

	Normalize(A);
}


void FFFW_slab_transposed_order_3D::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar)
{

	FFTW_MPI_EXECUTE_DFT_C2R(plan_ift_c2r_xyz, reinterpret_cast<FFTW_COMPLEX*>(A.data()), Ar.data());

	Zero_pad_last_plane(Ar); 
}

void FFFW_slab_transposed_order_3D::Transpose(Array<DP,3> Ar, Array<complx,3> A)
{
	//Alltoall(Ar,A,YZ);
}

void FFFW_slab_transposed_order_3D::Transpose(Array<complx,3> A, Array<DP,3> Ar)
{
	//Alltoall(A,Ar,XZ);
}

//********************************	End of four_tr.cc *****************************************


