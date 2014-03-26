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
 * along with Tarang-2; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */

/*! \file scft_tr.cc 
 * 
 * @sa scft_tr.h
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "spectral_transform.h"
#include "transpose.h"
 


//*********************************************************************************************
void SpectralTransform::Init_SSF_SLAB()
{
	if (Ny == 1) {
		if (my_id == 0) cerr << "ERROR: 2D Not implemented for SSF basis, Please use FFFW basis (uses original FFTW functions)" << endl;
		exit(1);
	}

	local_Nx = Nx/numprocs;
	local_Ny = Ny/numprocs;
	local_Nz = Nz/2+1;
	
	local_Nx_start = my_id * local_Nx;
	local_Ny_start = my_id * local_Ny;
	local_Nz_start = 0;

	shape_horizontal_array_3d = shape(Ny, Nz/2+1, local_Nx);
	shape_vertical_array_3d = shape(local_Ny, Nz/2+1, Nx);
	
	X_3d.resize(shape_horizontal_array_3d);  // temp array

	int Nx_plan[]={Nx};
	int Ny_plan[]={Ny};
	int Nz_plan[]={Nz};
	
	//Along Z
	r2c_z_plan = FFTW_PLAN_MANY_DFT_R2C(1, Nz_plan, local_Nx, reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx, 1, reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL, local_Nx, 1, FFTW_PLAN_FLAG);

	c2r_z_plan = FFTW_PLAN_MANY_DFT_C2R(1, Nz_plan, local_Nx, reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL, local_Nx, 1,
		reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx, 1, FFTW_PLAN_FLAG);
	

	fftw_r2r_kind kind[1];
	
	//Along X
	kind[0]=FFTW_RODFT10;
	sintr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, local_Ny*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT10;
	costr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, local_Ny*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	isintr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, local_Ny*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	icostr_x_plan = FFTW_PLAN_MANY_R2R(1, Nx_plan, local_Ny*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, kind, FFTW_PLAN_FLAG);


	//Along Y
	kind[0]=FFTW_RODFT10;
	sintr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, local_Nx*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT10;
	costr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, local_Nx*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	isintr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, local_Nx*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	icostr_y_plan = FFTW_PLAN_MANY_R2R(1, Ny_plan, local_Nx*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx*(Nz+2), 1, kind, FFTW_PLAN_FLAG);


	//Configure Transpose
    Set_transpose_config<DP>(ZX_PLANE, YZ_PLANE, MPI_COMM_WORLD, XY);
	YX=XY.conj();

}
//*********************************************************************************************



void SpectralTransform::Zero_pad_last_plane_SSF_SLAB(Array<DP,3> Ar)
{
   Ar(Range::all(),Range(Nz,Nz+1),Range::all()) = 0.0;
}

void SpectralTransform::Zero_pad_last_col_SSF_SLAB(Array<DP,2> Ar)
{
	if (my_id == 0) cerr << "ERROR: 2D Not implemented for SSF basis, Please use SSS basis" << endl;
}

//*********************************************************************************************

void SpectralTransform::Norm_SSF_SLAB(Array<complx,3> A)
{  
	A = A/(4*DP(Nx)*DP(Ny)*DP(Nz));
}

void SpectralTransform::Norm_SSF_SLAB(Array<complx,2> A)
{
	if (my_id == 0) cerr << "ERROR: 2D Not implemented for SSF basis, Please use SSS basis" << endl;
}


/**********************************************************************************************

	SFT(Ar) = A; Ar(Ny, Nx, Nz/2+1) is in transposed order

	Note: z= N[3]/2 plane does contain any real data (zero everywhere).
	
***********************************************************************************************/


void SpectralTransform::Forward_transform_SSF_SLAB(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A)
{
    
	Zero_pad_last_plane_SSF_SLAB(Ar);							 // Zero_pad the row=Nz/2
	
	SinCostr_x(sincostr_switch[0], Ar);
	
	Transpose_array(Ar, A, XY);
	
	SinCostr_y(sincostr_switch[1], A);

	FTr2c_z(A);
	
	Norm_SSF_SLAB(A);
}




/**********************************************************************************************

	Inverse SFT (A)  = Ar 
	IFT along perp dirn and SIN transform along x dirn of A 

***********************************************************************************************/

void SpectralTransform::Inverse_transform_SSF_SLAB(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar)
{
	X_3d = A;
	
	FTc2r_z(X_3d);

	ISinCostr_y(sincostr_switch[1], X_3d);
	
	Transpose_array(X_3d, Ar, YX);
	
	ISinCostr_x(sincostr_switch[0], Ar);

	Zero_pad_last_plane_SSF_SLAB(Ar);							 // Zero_pad the row=Nz/2

}

//***** 2D

//*** 2D
void SpectralTransform::Forward_transform_SSF_SLAB(string sincostr_switch, Array<DP,2> Ar, Array<complx,2> A)
{
	if (my_id == 0) cerr << "ERROR: 2D Not implemented for SSF basis, Please use SSS basis" << endl;
}


void SpectralTransform::Inverse_transform_SSF_SLAB(string sincostr_switch, Array<complx,2> A, Array<DP,2> Ar)
{
	if (my_id == 0) cerr << "ERROR: 2D Not implemented for SSF basis, Please use SSS basis" << endl;
}


//******************************** End of SSF_slab_tr.cc  *****************************************
