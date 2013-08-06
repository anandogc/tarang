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

	
//*********************************************************************************************
 void SpectralTransform::Init_FFF_SLAB()
{
    
	if (Ny > 1) {
		local_Nx = Nx/numprocs;
		local_Ny = Ny/numprocs;
		local_Nz = Nz/2+1;
		
		local_Nx_start = my_id * local_Nx;
		local_Ny_start = my_id * local_Ny;
		local_Nz_start = 0;
	}
	
	else if (Ny == 1) {
		local_Nx = Nx/numprocs;
		local_Ny = 1;
		local_Nz = (Nz/2+1)/numprocs;
		
		local_Nx_start = my_id * local_Nx;
		local_Ny_start = 0;
		local_Nz_start = my_id * local_Nz;
	}
    
    X.resize(local_Nx, Ny, Nz/2+1);
    Xr.resize(local_Ny,Nx,Nz/2+1);

    
    if(fftw3D_switch)
		plane_xz.resize(local_Nx,Nz/2+1);
	else
		plane_xz.resize(Nx,Nz/2+1);
    
	col_y.resize(Ny);
	
	slab_transpose.AxyTR.resize(Ny,local_Nx);
	slab_transpose.AxyTR_recv.resize(Ny,local_Nx);
	
	slab_transpose.AyxTR.resize(Nx,local_Ny);
	slab_transpose.AyxTR_recv.resize(Nx,local_Ny);
	
	/*
	 // May not be required
	 global.temp_array..slab_transpose.AzxTR.resize(Nx,local_Nz);
	 global.temp_array..slab_transpose.AzxTR_recv.resize(Nx,local_Nz);
	 */
	
    //Initialize plans
	if (fftw3D_switch) {
		
		if (Ny > 1) {
			r2c_xyz_plan = FFTW_MPI_PLAN_DFT_R2C_3D_DP(Nx, Ny, Nz,reinterpret_cast<DP*>(X.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(X.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
			
			c2r_xyz_plan = FFTW_MPI_PLAN_DFT_C2R_3D_DP(Nx, Ny, Nz,reinterpret_cast<FFTW_COMPLEX_DP*>(X.data()), reinterpret_cast<DP*>(X.data()), MPI_COMM_WORLD, FFTW_PLAN_FLAG);
			
		}
		
		else if (Ny == 1) {
			r2c_xz_plan = FFTW_MPI_PLAN_DFT_R2C_2D_DP(Nx, Nz, reinterpret_cast<DP*>(plane_xz.data()),reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
			
			c2r_xz_plan = FFTW_MPI_PLAN_DFT_C2R_2D_DP(Nx, Nz,reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()), reinterpret_cast<DP*>(plane_xz.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
		}
	}
	
	else {
		r2c_xz_plan = FFTW_PLAN_DFT_R2C_2D_DP(Nx, Nz, reinterpret_cast<DP*>(plane_xz.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()), FFTW_PLAN_FLAG);
        
		c2r_xz_plan = FFTW_PLAN_DFT_C2R_2D_DP(Nx, Nz, reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()), reinterpret_cast<DP*>(plane_xz.data()),FFTW_PLAN_FLAG);
		
		c2c_y_forward_plan = FFTW_PLAN_DFT_1D_DP(Ny, reinterpret_cast<FFTW_COMPLEX_DP*>(col_y.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(col_y.data()), FFTW_FORWARD, FFTW_PLAN_FLAG);
		
		c2c_y_inverse_plan = FFTW_PLAN_DFT_1D_DP(Ny, reinterpret_cast<FFTW_COMPLEX_DP*>(col_y.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(col_y.data()), FFTW_BACKWARD, FFTW_PLAN_FLAG);
	}
}
//*********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_FFF_SLAB(Array<complx,3> Ar)
{
	Ar(Range::all(),Range::all(),Nz/2) = 0.0;
}

void SpectralTransform::Zero_pad_last_plane_FFF_SLAB(Array<complx,2> Ar)
{
	Ar(Range::all(),Nz/2) = 0.0;
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


//*********************************************************************************************

// Ar: yxz, A: xyz
void SpectralTransform::Forward_transform_array_transpose_order_FFF_SLAB(Array<complx,3> Ar, Array<complx,3> A)
{
		Zero_pad_last_plane_FFF_SLAB(Ar);							 // Zero_pad the row=Nz
		
		for (int ly=0; ly<local_Ny; ly++) {
			plane_xz = Ar(ly, Range::all(), Range::all());
           FTr2c_xz(plane_xz);
			Ar(ly, Range::all(), Range::all()) = plane_xz;
		}
		
		Transpose_array_SLAB(Ar, A, 'Y', 'X', 'Z', 'Z');
		
		for (int lx=0; lx<local_Nx; lx++)
			for (int lz=0; lz<=Nz/2; lz++) {
				col_y = A(lx, Range::all(), lz);
               FTc2c_y(col_y);
				A(lx, Range::all(), lz) = col_y;
			}
		
		Norm_FFF_SLAB(A);
}



// Ar: yxz, A: xyz
void SpectralTransform::Forward_transform_array_transpose_order_FFF_SLAB(Array<complx,2> Ar, Array<complx,2> A)
{

    cout << "ERROR:  ArrayFFT_FFF_transpose_order for 2D array not allowed in FOUR basis " << endl;
    exit(1);
}


//*********************************************************************************************

void SpectralTransform::Forward_transform_array_FFF_SLAB(Array<complx,3> A) 
{
	
	if (fftw3D_switch) {
		Zero_pad_last_plane_FFF_SLAB(A);
		
        FTr2c_xyz(A);
		
		Norm_FFF_SLAB(A);
	}
	
	// not original switch.. split in 2D and 1D ffts.
	else {
        Transpose_array_SLAB(A, Xr, 'X', 'Y', 'Z', 'Z');
        Forward_transform_array_transpose_order_FFF_SLAB(Xr, A);	// Norm done here
	}
	
}
        
        
        
void SpectralTransform::Forward_transform_array_FFF_SLAB(Array<complx,2> A)
{
    
    if (fftw3D_switch) {
        Zero_pad_last_plane_FFF_SLAB(A);
        FTr2c_xz(A);
        Norm_FFF_SLAB(A);
    }
    
    // not original switch.. split in 2D and 1D ffts.
    else {
        cout << "ERROR: USE ORIGINAL FFT FOR 2D. " << endl;
        exit(1);
    }
    
}

//*********************************************************************************************


void SpectralTransform::Inverse_transform_array_transpose_order_FFF_SLAB(Array<complx,3> A, Array<complx,3> Ar)
{
	
    for (int lx=0; lx<local_Nx; lx++) 
        for (int lz=0; lz<=Nz/2; lz++) {
            col_y = A(lx, Range::all(), lz);
            IFTc2c_y(col_y);
            X(lx, Range::all(), lz) = col_y;
        }
    
    Transpose_array_SLAB(X, Ar, 'X', 'Y', 'Z', 'Z');
    
    for (int ly=0; ly < local_Ny; ly++)	{
        plane_xz = Ar(ly, Range::all(), Range::all());
        FTc2r_xz(plane_xz);
        Ar(ly, Range::all(), Range::all()) = plane_xz;
    }
    
    Zero_pad_last_plane_FFF_SLAB(Ar);
}


void SpectralTransform::Inverse_transform_array_transpose_order_FFF_SLAB(Array<complx,2> A, Array<complx,2> Ar)
{
    cout << "ERROR:  ArrayIFFT_FFF_transpose_order for 2D array not allowed in FOUR basis " << endl;
    exit(1);
}


//*********************************************************************************************


void SpectralTransform::Inverse_transform_array_FFF_SLAB(Array<complx,3> A)
{
	
	if (fftw3D_switch) 
        FTc2r_xyz(A);
	
    // By splitting it up...
	else {
        Inverse_transform_array_transpose_order_FFF_SLAB(A, Xr);
        Transpose_array_SLAB(Xr, A, 'Y', 'X', 'Z', 'Z');
    }
}


void SpectralTransform::Inverse_transform_array_FFF_SLAB(Array<complx,2> A)
{
	
	if (fftw3D_switch) 
        FTc2r_xz(A);
	
    // By splitting it up...
	else {
			cout << "ERROR: USE ORIGINAL FFT FOR 2D. " << endl;
            exit(1);
	}
}




//********************************	End of four_tr.cc *****************************************


