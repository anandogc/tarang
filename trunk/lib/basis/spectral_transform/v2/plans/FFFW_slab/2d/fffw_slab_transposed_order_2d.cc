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


#include "fffw_slab_transposed_order_2d.h"


	
//*********************************************************************************************
FFFW_slab_transposed_order_2D::FFFW_slab_transposed_order_2D(int my_id, int numprocs, int num_iter, int Nx, int Nz): SpectralPlan_Slab_2D("FFFW", my_id, numprocs, Nx, Nz)
{

	X_2d.resize(local_Nx, Nz/2+1);
	Xr_2d.resize(local_Nx, 2*(Nz/2+1));

	//Initialize plans

	fftw_mpi_init();

	//Initialize plans
	plan_ft_r2c_xz = FFTW_MPI_PLAN_DFT_R2C_2D(Nx, Nz,
		reinterpret_cast<DP*>(Xr_2d.data()), reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()),
		MPI_COMM_WORLD, FFTW_PLAN_FLAG);

	plan_ift_c2r_xz = FFTW_MPI_PLAN_DFT_C2R_2D(Nx, Nz,
		reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), reinterpret_cast<DP*>(Xr_2d.data()),
		MPI_COMM_WORLD, FFTW_PLAN_FLAG);

	Init_array();
	Evaluate_time_per_step(num_iter);

	X_2d.free();
	Xr_2d.free();
}
//*********************************************************************************************

DP FFFW_slab_transposed_order_2D::f(int rx, int rz)
{
	DP k0 = 1;
	DP x,z;
	DP L=2*M_PI;

	x = rx*L/Nx;
	z = rz*L/Nz;
	return 4*sin(k0*x)*cos(k0*z);
}
 

void FFFW_slab_transposed_order_2D::Init_array()
{
	for (int rz=0; rz<2*local_Nz; rz++)
		for (int rx=0; rx<Nx; rx++)
			Xr_2d(rz, rx) = f(2*local_Nz_start + rz, rx);
	Zero_pad_last_plane(Xr_2d);
}

//*********************************************************************************************

void FFFW_slab_transposed_order_2D::Normalize(Array<complx,2> A)
{
	A /= (DP(Nx) * DP(Nz));
}

//*************************************************************************************

void FFFW_slab_transposed_order_2D::Forward_transform(Array<DP,2> Ar, Array<complx,2> A)
{

	Zero_pad_last_plane(Ar);	// Zero_pad the row=Nz

	FFTW_MPI_EXECUTE_DFT_R2C(plan_ft_r2c_xz, Ar.data(), reinterpret_cast<FFTW_COMPLEX*>(A.data()));

	Normalize(A);
}


void FFFW_slab_transposed_order_2D::Inverse_transform(Array<complx,2> A, Array<DP,2> Ar)
{
	FFTW_MPI_EXECUTE_DFT_C2R(plan_ift_c2r_xz, reinterpret_cast<FFTW_COMPLEX*>(A.data()), Ar.data());

	Zero_pad_last_plane(Ar); 
}


void FFFW_slab_transposed_order_2D::Transpose(Array<DP,2> Ar, Array<complx,2> A)
{
}

void FFFW_slab_transposed_order_2D::Transpose(Array<complx,2> A, Array<DP,2> Ar)
{
}


//********************************	End of four_tr.cc *****************************************


