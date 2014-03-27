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


#include "sff_slab_alltoall_3d.h"


	
//*********************************************************************************************
SFF_slab_Alltoall_3D::SFF_slab_Alltoall_3D(int my_id, int numprocs, int Nx, int Ny, int Nz): SpectralPlan(my_id, numprocs, Nx, Ny, Nz)
{

	X_3d.resize(local_Nx, Ny, Nz/2+1);	   // temp array

	//Initialize plans
	int Nx_plan[]={Nx};
	int Nyz_plan[]={Ny,Nz};

	fftw_r2r_kind kind[1];

	kind[0]=FFTW_RODFT10;
	plan_sintr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT10;
	plan_costr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	plan_icostr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);


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

	X_3d.free();
}
//*********************************************************************************************

void SFF_slab_Alltoall_3D::Normalize(Array<complx,3> A)
{
	A /= (2.0 * DP(Nx) * DP(Ny) * DP(Nz));
}

//*************************************************************************************


void SFF_slab_Alltoall_3D::Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A)
{
	Zero_pad_last_plane(Ar);	// Zero_pad the plane Nz

	SinCostr_x(sincostr_option[0],Ar);

	Transpose(Ar,A);

	FT_r2c_yz(A);

	Normalize(A);
}


void SFF_slab_Alltoall_3D::Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar)
{
	IFT_c2r_yz(A);

	Transpose(A,Ar);

	ISinCostr_x(sincostr_option[0], Ar);

	Zero_pad_last_plane(Ar); 
}



//Transpose
void SFF_slab_Alltoall_3D::Transpose(Array<DP,3> Ar, Array<complx,3> A) {
	Alltoall(Ar,A,XZ);
}

void SFF_slab_Alltoall_3D::Transpose(Array<complx,3> A, Array<DP,3> Ar) {
	Alltoall(A,Ar,YZ);
}


void SFF_slab_Alltoall_3D::SinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
	if (sincostr_option == 'S') {
		FFTW_EXECUTE_R2R(plan_sintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'X');
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_costr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SFF_slab_Alltoall_3D::ISinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
    if (sincostr_option == 'S') {
		ArrayShiftLeft(Ar, 'X');
		FFTW_EXECUTE_R2R(plan_isintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_icostr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SFF_slab_Alltoall_3D::FT_r2c_yz(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_R2C(plan_ft_r2c_yz, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX*>(A.data()));
}


void SFF_slab_Alltoall_3D::IFT_c2r_yz(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_C2R(plan_ift_c2r_yz, reinterpret_cast<FFTW_COMPLEX*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

//********************************	End of four_tr.cc *****************************************


