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
 void SpectralTransform::Init_FFF_PENCIL()
{

	local_Nx_vert = Nx/num_p_vert;
	local_Ny_vert = Ny/num_p_vert;
	local_Nz_vert = (Nz/2+1)/num_p_vert;
	
	local_Nx_hor = Nx/num_p_hor;
	local_Ny_hor = Ny/num_p_hor;
	local_Nz_hor = (Nz/2+1)/num_p_hor;
	
	local_Nx_start = 0;
	local_Ny_start = my_vert_pcoord*local_Ny_vert;
	local_Nz_start = my_hor_pcoord*local_Nz_hor;
	
	shape_x_array_3d = shape(local_Ny_vert,local_Nz_hor,Nx);
	shape_y_array_3d = shape(Ny,local_Nz_hor,local_Nx_vert);
	shape_z_array_3d = shape(local_Ny_hor,Nz/2+1,local_Nx_vert);

	X_3d.resize(shape_x_array_3d);	   // temp array
	A3d_interm.resize(shape_y_array_3d);	   // temp array

		
	//Initialize plans
	int Nx_plan[]={Nx};
	int Ny_plan[]={Ny};
	int Nz_plan[]={Nz};
	
	
	//Strided plans
	c2c_x_forward_plan = fftw_plan_many_dft(1, Nx_plan, local_Ny_vert*local_Nz_hor,
		reinterpret_cast<fftw_complex *>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<fftw_complex*>(X_3d.data()), NULL,
		1, Nx,
		FFTW_FORWARD, FFTW_MEASURE);

	c2c_x_inverse_plan = fftw_plan_many_dft(1, Nx_plan, local_Ny_vert*local_Nz_hor,
		reinterpret_cast<fftw_complex *>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<fftw_complex*>(X_3d.data()), NULL,
		1, Nx,
		FFTW_BACKWARD, FFTW_MEASURE);
	
	c2c_y_forward_plan = fftw_plan_many_dft(1, Ny_plan, local_Nx_vert*local_Nz_hor,
		reinterpret_cast<fftw_complex *>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<fftw_complex*>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		FFTW_FORWARD, FFTW_MEASURE);

	c2c_y_inverse_plan = fftw_plan_many_dft(1, Ny_plan, local_Nx_vert*local_Nz_hor,
		reinterpret_cast<fftw_complex *>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<fftw_complex*>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		FFTW_BACKWARD, FFTW_MEASURE);
	
	c2r_z_plan = fftw_plan_many_dft_c2r(1, Nz_plan, local_Nx_vert,
		reinterpret_cast<fftw_complex *>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		FFTW_MEASURE);
	
	r2c_z_plan = fftw_plan_many_dft_r2c(1, Nz_plan, local_Nx_vert,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		reinterpret_cast<fftw_complex *>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		FFTW_MEASURE);


	//Confugure Transpose
	Set_transpose_config<complx>(X, Y, MPI_COMM_VERT_SEGMENT, XY);
	YX=XY.conj();

	Set_transpose_config<complx>(Y, Z, MPI_COMM_HOR_SEGMENT, YZ);
	ZY=YZ.conj();
}
//*********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_FFF_PENCIL(Array<DP,3> Ar)
{
	Ar(Range::all(),Range(Nz,Nz+1),Range::all()) = 0.0;
}

//*********************************************************************************************

void SpectralTransform::Norm_FFF_PENCIL(Array<complx,3> A)
{
	A = A/(DP(Nx) * DP(Ny) * DP(Nz));
}

//*************************************************************************************


void SpectralTransform::Forward_transform_FFF_PENCIL(Array<DP,3> Ar, Array<complx,3> A)
{

	Zero_pad_last_plane_FFF_PENCIL(Ar);							 // Zero_pad the row=Nz

	FTr2c_z(Ar);

	Transpose_array(Ar,A3d_interm, ZY);
	
	FTc2c_y(A3d_interm);
	
	Transpose_array(A3d_interm, A, YX);
	
	FTc2c_x(A);

	Norm_FFF_PENCIL(A);
}
		
//*********************************************************************************************

void SpectralTransform::Inverse_transform_FFF_PENCIL(Array<complx,3> A, Array<DP,3> Ar)
{

	X_3d = A;

	IFTc2c_x(X_3d);
	
	Transpose_array(X_3d, A3d_interm, XY);
	
	IFTc2c_y(A3d_interm);
	
	Transpose_array(A3d_interm, Ar, YZ);
	
	FTc2r_z(Ar);

	Zero_pad_last_plane_FFF_PENCIL(Ar); 
}

//********************************	End of four_tr.cc *****************************************


