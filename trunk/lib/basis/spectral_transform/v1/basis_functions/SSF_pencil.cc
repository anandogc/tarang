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


	
//*********************************************************************************************
 void SpectralTransform::Init_SSF_PENCIL()
{

	local_Nx_vert = Nx/num_p_vert;
	local_Ny_vert = Ny/num_p_vert;
	local_Nz_vert = (Nz/2+1)/num_p_vert;
	
	local_Nx_hor = Nx/num_p_hor;
	local_Ny_hor = Ny/num_p_hor;
	local_Nz_hor = (Nz/2+1)/num_p_hor;
	
	local_Nx_start = my_vert_pcoord*local_Nx_vert;
	local_Ny_start = my_hor_pcoord*local_Ny_hor;
	local_Nz_start = 0; 
	
	shape_x_array_3d = shape(local_Ny_vert,local_Nz_hor,Nx);
	shape_y_array_3d = shape(Ny,local_Nz_hor,local_Nx_vert);
	shape_z_array_3d = shape(local_Ny_hor,Nz/2+1,local_Nx_vert);

	X_3d.resize(shape_z_array_3d);	   // temp array
	A3d_interm.resize(shape_y_array_3d);	   // intermediate array

		
	//Initialize plans
	int Nx_plan[]={Nx};
	int Ny_plan[]={Ny};
	int Nz_plan[]={Nz};
	
	
	//Strided plans
	fftw_r2r_kind kind[1];
	
	// Along X
	kind[0]=FFTW_RODFT10;
	sintr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, 2*local_Ny_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		kind, FFTW_MEASURE);
	
	kind[0]=FFTW_REDFT10;
	costr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, 2*local_Ny_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		kind, FFTW_MEASURE);
	
	kind[0]=FFTW_RODFT01;
	isintr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, 2*local_Ny_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		kind, FFTW_MEASURE);
	
	
	kind[0]=FFTW_REDFT01;
	icostr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, 2*local_Ny_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, Nx,
		kind, FFTW_MEASURE);
	
	
	// Along Y
	kind[0]=FFTW_RODFT10;
	sintr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, 2*local_Nx_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		kind, FFTW_MEASURE);
	
	kind[0]=FFTW_REDFT10;
	costr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, 2*local_Nx_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		kind, FFTW_MEASURE);
	
	kind[0]=FFTW_RODFT01;
	isintr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, 2*local_Nx_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		kind, FFTW_MEASURE);
	
	
	kind[0]=FFTW_REDFT01;
	icostr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, 2*local_Nx_vert*local_Nz_hor,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		2*local_Nx_vert*local_Nz_hor, 1,
		kind, FFTW_MEASURE);
	
	// Along-Z
	c2r_z_plan = FFTW_PLAN_MANY_DFT_C2R(1, Nz_plan, local_Nx_vert,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		FFTW_MEASURE);
	
	r2c_z_plan = FFTW_PLAN_MANY_DFT_R2C(1, Nz_plan, local_Nx_vert,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx_vert, 1,
		FFTW_MEASURE);

	//Confugure Transpose
	Set_transpose_config<DP>(X, Y, MPI_COMM_VERT_SEGMENT, XY);
	YX=XY.conj();

	Set_transpose_config<DP>(Y, Z, MPI_COMM_HOR_SEGMENT, YZ);
	ZY=YZ.conj();
}
//*********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_SSF_PENCIL(Array<DP,3> Ar)
{
	if (my_hor_pcoord == num_p_hor-1)
		Ar(Range::all(),Range(2*local_Nz_hor-2,toEnd),Range::all()) = 0.0;
}

//*********************************************************************************************

void SpectralTransform::Norm_SSF_PENCIL(Array<complx,3> A)
{
	A = A/(4*DP(Nx) * DP(Ny) * DP(Nz));
}

//*************************************************************************************


void SpectralTransform::Forward_transform_SSF_PENCIL(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A)
{

	Zero_pad_last_plane_SSF_PENCIL(Ar);							 // Zero_pad the row=Nz

	SinCostr_x(sincostr_switch[0], Ar);

	Transpose_array(Ar, A3d_interm, XY);
	
	SinCostr_y(sincostr_switch[1], A3d_interm);
	
	Transpose_array(A3d_interm, A, YZ);
	
	FTr2c_z(A);

	Norm_SSF_PENCIL(A);
}
		
//*********************************************************************************************

void SpectralTransform::Inverse_transform_SSF_PENCIL(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar)
{

	X_3d = A;
	
	FTc2r_z(X_3d);
	
	Transpose_array(X_3d, A3d_interm, ZY);
	
	ISinCostr_y(sincostr_switch[1], A3d_interm);
	
	Transpose_array(A3d_interm, Ar, YX);
	
	ISinCostr_x(sincostr_switch[0], Ar);

	Zero_pad_last_plane_SSF_PENCIL(Ar);
}

//********************************	End of four_tr.cc *****************************************


