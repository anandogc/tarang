/* Tarang-2
 *
 * Copyright_PENCIL( (C) 2008, 2009  Mahendra K. Verma
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
 * SIZE of A(Nx, Ny, Nz/2+1) with Nx=2^n (FFTW convention);
 * FFTW plan uses (Nx+1, col_x, col_x).. Nx+1 points.
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
#include "set_transpose_config.h"
 

//*********************************************************************************************
void SpectralTransform::Init_ChFF_SLAB()
{
	if (Ny > 1) {
		local_Nx = Nx/numprocs;
		local_Ny = Ny/numprocs;
		local_Nz = Nz/2+1;
		
		local_Nx_start = my_id * local_Nx;
		local_Ny_start = my_id * local_Ny;
		local_Nz_start = 0;

		// for transpose only..
		shape_vertical_array_3d = shape(local_Ny, Nz/2+1, Nx);
		shape_horizontal_array_3d = shape(Ny, Nz/2+1, local_Nx);

		X_3d.resize(shape_vertical_array_3d*shape(1,2,1));  // temp array
		A3d_interm.resize(shape_horizontal_array_3d*shape(1,2,1));  // real array
	}
	
	else if (Ny == 1) {
		local_Nx = Nx/numprocs;
		local_Ny = 1;
		local_Nz = (Nz/2+1)/numprocs;
		
		local_Nx_start = my_id * local_Nx;
		local_Ny_start = 0;
		local_Nz_start = my_id * local_Nz;

		shape_vertical_array_2d = shape(local_Nz, Nx);
		shape_horizontal_array_2d = shape(Nz/2+1, local_Nx);
		
		X_2d.resize(shape_vertical_array_2d);		// temp array
		A2d_interm.resize(shape_horizontal_array_2d*shape(2,1));  // real array
	}
	
	cout << Nx << " " << Ny << " " << Nz << endl;
	cout << local_Nx << " " << local_Ny << " " << local_Nz << endl;

	if (Ny>1) {
		int Nyz_plan[]={Ny,Nz};
		int Nx_plan[]={Nx};
		
		r2c_yz_plan = fftw_plan_many_dft_r2c(2, Nyz_plan, local_Nx, reinterpret_cast<DP*>(A3d_interm.data()), NULL, local_Nx, 1, reinterpret_cast<fftw_complex *>(A3d_interm.data()), NULL, local_Nx, 1, FFTW_PLAN_FLAG);

		c2r_yz_plan = fftw_plan_many_dft_c2r(2, Nyz_plan, local_Nx, reinterpret_cast<fftw_complex *>(X_3d.data()), NULL, local_Nx, 1,
			reinterpret_cast<DP*>(X_3d.data()), NULL, local_Nx, 1, FFTW_PLAN_FLAG);
		

		fftw_r2r_kind kind[1];
		
		kind[0]=FFTW_REDFT00;
		Chebyshevtr_x_plan = fftw_plan_many_r2r(1, Nx_plan, local_Ny*(Nz+2), reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, reinterpret_cast<DP*>(X_3d.data()), NULL, 1, Nx, kind, FFTW_PLAN_FLAG);
	}
    
    
    else if (Ny==1) { // 2D
		int Nz_plan[]={Nz};
		int Nx_plan[]={Nx};
		
		r2c_z_plan = fftw_plan_many_dft_r2c(1, Nz_plan, local_Nx, reinterpret_cast<DP*>(X_2d.data()), NULL, local_Nx, 1, reinterpret_cast<fftw_complex *>(X_2d.data()), NULL, local_Nx, 1,
			FFTW_PLAN_FLAG);

		c2r_z_plan = fftw_plan_many_dft_c2r(1, Nz_plan, local_Nx, reinterpret_cast<fftw_complex *>(X_2d.data()), NULL, local_Nx, 1, reinterpret_cast<DP*>(X_2d.data()), NULL, local_Nx, 1,
			FFTW_PLAN_FLAG);
		

		fftw_r2r_kind kind[1];
		
		kind[0]=FFTW_REDFT00;
		Chebyshevtr_x_plan = fftw_plan_many_r2r(1, Nx_plan, 2*local_Nz, reinterpret_cast<DP*>(X_2d.data()), NULL, 1, Nx, reinterpret_cast<DP*>(X_2d.data()), NULL, 1, Nx, kind, FFTW_PLAN_FLAG);
	}


	if (Ny>1) {
		Set_transpose_config<DP>(YZ_PLANE, ZX_PLANE, MPI_COMM_WORLD, XY_real);
		YX_real = XY_real.conj();

		Set_transpose_config<complx>(YZ_PLANE, ZX_PLANE, MPI_COMM_WORLD, XY_complex);
		YX_complex=XY_complex.conj();
	}
	else {
		Set_transpose_config<DP>(X, Z, MPI_COMM_WORLD, XZ_real);
		ZX_real = XZ_real.conj();

		Set_transpose_config<complx>(X, Z, MPI_COMM_WORLD, XZ_complex);
		ZX_complex=XZ_complex.conj();
	}
}
//*************************************************************************************


void SpectralTransform::Zero_pad_last_plane_ChFF_SLAB(Array<DP,3> Ar)
{
    Ar(Range::all(),Range(Nz,Nz+1),Range::all()) = 0.0;
}


void SpectralTransform::Zero_pad_last_col_ChFF_SLAB(Array<DP,2> Ar)
{
    if (my_id == numprocs-1)
        Ar(Range(2*local_Nz-2,2*local_Nz-1), Range::all()) = 0.0;
    // Zero pad Ar; In real space.. 2*local_Nz cols
}



//*********************************************************************************************

void SpectralTransform::Norm_ChFF_SLAB(Array<complx,3> A) 
{  
	A = A/(2*DP(Nx-1)*DP(Ny)*DP(Nz));
}


void SpectralTransform::Norm_ChFF_SLAB(Array<complx,2> A)
{
	A = A/(2*DP(Nx-1)*DP(Nz));
}

/**********************************************************************************************

	SFT(Ar) = A; Ar(Ny, Nx, Nz/2+1) is in transposed order

	Note: z= Nz/2 plane does contain any real data (zero everywhere).
	
***********************************************************************************************/

void SpectralTransform::Forward_transform_ChFF_SLAB(Array<DP,3> Ar, Array<complx,3> A)
{
        
	Zero_pad_last_plane_ChFF_SLAB(Ar);							 // Zero_pad the row=Nz/2
  	
	Chebyshevtr_x(Ar);
  
	Transpose_array(Ar, A3d_interm, YX_real);
  
    FTr2c_yz(A3d_interm);

	Transpose_array(A3d_interm, A, XY_complex);
	
	Norm_ChFF_SLAB(A);
}




/**********************************************************************************************

	Inverse SFT (A)  = Ar 
	IFT along perp dirn and SIN transform along x dirn of A 

***********************************************************************************************/

void SpectralTransform::Inverse_transform_ChFF_SLAB(Array<complx,3> A, Array<DP,3> Ar)
{
	Transpose_array(A, A3d_interm, YX_complex);

	FTc2r_yz(A3d_interm);

	Transpose_array(A3d_interm, Ar, XY_real);

	IChebyshevtr_x(Ar);

	Zero_pad_last_plane_ChFF_SLAB(Ar);
}

// 2D ************
void SpectralTransform::Forward_transform_ChFF_SLAB(Array<DP,2> Ar, Array<complx,2> A)
{
	
	Zero_pad_last_col_ChFF_SLAB(Ar);			// Zero_pad the row=Nz/2
  	
	Chebyshevtr_x(Ar);
	
	Transpose_array(Ar, A2d_interm, XZ_real);
	
    FTr2c_z(A2d_interm);
	
	Transpose_array(A2d_interm, A, ZX_complex);
	
	Norm_ChFF_SLAB(A);
}

// 2D ************
void SpectralTransform::Inverse_transform_ChFF_SLAB(Array<complx,2> A, Array<DP,2> Ar)
{
	Transpose_array(A, A2d_interm, XZ_complex);
	
	FTc2r_z(A2d_interm);

	Transpose_array(A2d_interm, Ar, ZX_real);

	IChebyshevtr_x(Ar);

	Zero_pad_last_col_ChFF_SLAB(Ar);
}


//******************************** End of ChFF_SLAB_tr.cc  *****************************************
