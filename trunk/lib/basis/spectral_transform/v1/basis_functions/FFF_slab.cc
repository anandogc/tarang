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


#include "spectral_transform.h"
#include "transpose.h"
#include "set_transpose_config.h"


	
//*********************************************************************************************
 void SpectralTransform::Init_FFF_SLAB()
{
	if (Ny == 1) {
		if (my_id == 0) cerr << "ERROR: 2D Not implemented for FFF basis, Please use FFTW basis (uses original FFTW functions)" << endl;
		return;
	}

	local_Nx = Nx/numprocs;
	local_Ny = Ny/numprocs;
	local_Nz = Nz/2+1;
	
	local_Nx_start = my_id * local_Nx;
	local_Ny_start = my_id * local_Ny;
	local_Nz_start = 0;
	
	shape_vertical_array_3d = shape(local_Ny,Nz/2+1,Nx);
	shape_horizontal_array_3d = shape(Ny,Nz/2+1,local_Nx);

	X_3d.resize(shape_vertical_array_3d);	   // temp array
		
	//Initialize plans

	int Nyz_plan[]={Ny,Nz};
	int Nx_plan[]={Nx};
	
	//Strided plans
	c2c_x_forward_plan = FFTW_PLAN_MANY_DFT(1, Nx_plan, local_Ny*(Nz/2+1),
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, Nx,
		FFTW_FORWARD, FFTW_MEASURE);

	c2c_x_inverse_plan = FFTW_PLAN_MANY_DFT(1, Nx_plan, local_Ny*(Nz/2+1),
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, Nx,
		FFTW_BACKWARD, FFTW_MEASURE);

	r2c_yz_plan = FFTW_PLAN_MANY_DFT_R2C(2, Nyz_plan, local_Nx,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		local_Nx, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx, 1,
		FFTW_MEASURE);

	c2r_yz_plan = FFTW_PLAN_MANY_DFT_C2R(2, Nyz_plan, local_Nx,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		local_Nx, 1,
		FFTW_MEASURE);


	//Confugure Transpose
	Set_transpose_config<complx>(ZX_PLANE, YZ_PLANE, MPI_COMM_WORLD, XY);
	YX=XY.conj();

}
//*********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_FFF_SLAB(Array<DP,3> Ar)
{
	Ar(Range::all(),Range(Nz,Nz+1),Range::all()) = 0.0;
}

void SpectralTransform::Zero_pad_last_col_FFF_SLAB(Array<DP,2> Ar)
{
	if (my_id == 0) cerr << "ERROR: 2D Not implemented for FFF basis, Please use FFTW basis (uses original FFTW functions)" << endl;
}

//*********************************************************************************************

void SpectralTransform::Norm_FFF_SLAB(Array<complx,3> A)
{
	A = A/(DP(Nx) * DP(Ny) * DP(Nz));
}

void SpectralTransform::Norm_FFF_SLAB(Array<complx,2> A)
{
	A = A/(DP(Nx) * DP(Nz));
}

//*************************************************************************************


void SpectralTransform::Forward_transform_FFF_SLAB(Array<DP,3> Ar, Array<complx,3> A)
{

	Zero_pad_last_plane_FFF_SLAB(Ar);							 // Zero_pad the row=Nz

	FTr2c_yz(Ar);
	Transpose_array(Ar, A, YX);
	FTc2c_x(A);

	Norm_FFF_SLAB(A);
}
		

void SpectralTransform::Forward_transform_FFF_SLAB(Array<DP,2> Ar, Array<complx,2> A)
{
	if (my_id == 0) cerr << "ERROR: 2D Not implemented for FFF basis, Please use FFTW basis (uses original FFTW functions)" << endl;
}

//*********************************************************************************************

void SpectralTransform::Inverse_transform_FFF_SLAB(Array<complx,3> A, Array<DP,3> Ar)
{

	X_3d = A;

	IFTc2c_x(X_3d);
	Transpose_array(X_3d, Ar, XY);
	FTc2r_yz(Ar);

	Zero_pad_last_plane_FFF_SLAB(Ar); 
}


//**************************************************************************************


void SpectralTransform::Inverse_transform_FFF_SLAB(Array<complx,2> A, Array<DP,2> Ar)
{
	
	if (my_id == 0) cerr << "ERROR: 2D Not implemented for FFF basis, Please use FFTW basis (uses original FFTW functions)" << endl;
}


//********************************	End of four_tr.cc *****************************************


