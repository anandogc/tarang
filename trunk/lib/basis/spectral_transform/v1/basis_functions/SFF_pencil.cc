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
 void SpectralTransform::Init_SFF_PENCIL()
{

	local_Nx_vert = Nx/num_p_vert;
	local_Ny_vert = Ny/num_p_vert;
	local_Nz_vert = (Nz/2+1)/num_p_vert;
	
	local_Nx_hor = Nx/num_p_hor;
	local_Ny_hor = Ny/num_p_hor;
	local_Nz_hor = (Nz/2+1)/num_p_hor;
	
	local_Nx_start = my_vert_pcoord*local_Nx_vert;
	local_Ny_start = 0;
	local_Nz_start = my_hor_pcoord*local_Nz_hor;
	
	shape_x_array_3d = shape(local_Ny_hor,local_Nz_vert,Nx);
	shape_y_array_3d = shape(Ny,local_Nz_hor,local_Nx_vert);
	shape_z_array_3d = shape(local_Ny_hor,Nz/2+1,local_Nx_vert);

	X_3d.resize(shape_y_array_3d);	   // temp array
	A3d_interm.resize(shape_z_array_3d);	   // temp array

		
	//Initialize plans
	int Nx_plan[]={Nx};
	int Ny_plan[]={Ny};
	int Nz_plan[]={Nz};
	
	
	//Strided plans
	fftw_r2r_kind kind[1];
	
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
	
	
	c2c_y_forward_plan = FFTW_PLAN_MANY_DFT(1, Ny_plan, local_Nx_vert*local_Nz_hor,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		FFTW_FORWARD, FFTW_MEASURE);

	c2c_y_inverse_plan = FFTW_PLAN_MANY_DFT(1, Ny_plan, local_Nx_vert*local_Nz_hor,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		local_Nx_vert*local_Nz_hor, 1,
		FFTW_BACKWARD, FFTW_MEASURE);
	
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
	Set_transpose_config<complx>(Y, Z, MPI_COMM_HOR_SEGMENT, YZ);
	ZY=YZ.conj();

	Set_transpose_config<DP>(Z, X, MPI_COMM_VERT_SEGMENT, ZX);
	XZ=ZX.conj();
}
//*********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_SFF_PENCIL(Array<DP,3> Ar)
{
	if (my_vert_pcoord == num_p_vert-1)
		Ar(Range::all(),Range(2*local_Nz_vert-2,toEnd),Range::all()) = 0.0;
}

//*********************************************************************************************

void SpectralTransform::Norm_SFF_PENCIL(Array<complx,3> A)
{
	A = A/(2*DP(Nx) * DP(Ny) * DP(Nz));
}

//*************************************************************************************


void SpectralTransform::Forward_transform_SFF_PENCIL(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A)
{

	Zero_pad_last_plane_SFF_PENCIL(Ar);							 // Zero_pad the row=Nz

	SinCostr_x(sincostr_switch[0], Ar);

	Transpose_array(Ar, A3d_interm, XZ);
	
	FTr2c_z(A3d_interm);
	
	Transpose_array(A3d_interm, A, ZY);
	
	FTc2c_y(A);

	Norm_SFF_PENCIL(A);
}
		
//*********************************************************************************************

void SpectralTransform::Inverse_transform_SFF_PENCIL(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar)
{
	X_3d = A;

	IFTc2c_y(X_3d);
	
	Transpose_array(X_3d, A3d_interm, YZ);
	
	FTc2r_z(A3d_interm);
	
	Transpose_array(A3d_interm, Ar, ZX);
	
	ISinCostr_x(sincostr_switch[0], Ar);

	Zero_pad_last_plane_SFF_PENCIL(Ar);
}

//********************************	End of four_tr.cc *****************************************


