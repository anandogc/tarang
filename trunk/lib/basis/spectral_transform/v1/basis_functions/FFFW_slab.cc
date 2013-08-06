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
 void SpectralTransform::Init_FFFW_SLAB()
{
    if (Ny > 1) {
        local_Nx = Nx/numprocs;
        local_Ny = Ny/numprocs;
        local_Nz = Nz/2+1;
        
        local_Nx_start = my_id * local_Nx;
        local_Ny_start = my_id * local_Ny;
        local_Nz_start = 0;

        shape_complex_array_3d = shape(local_Nx,Ny,Nz/2+1);
        shape_real_array_3d = shape(local_Ny,Nx,Nz/2+1);

        X.resize(shape_complex_array_3d);
        Xr.resize(shape_real_array_3d);
    }
    else if (Ny == 1) {
        local_Nx = Nx/numprocs;
        local_Ny = 1;
        local_Nz = (Nz/2+1)/numprocs;
        
        local_Nx_start = my_id * local_Nx;
        local_Ny_start = 0;
        local_Nz_start = my_id * local_Nz;

        X_2d.resize(shape_complex_array_2d);
        Xr_2d.resize(shape_real_array_2d);
    }

	
    //Initialize plans
		
    if (Ny > 1) {
        r2c_xyz_plan = FFTW_MPI_PLAN_DFT_R2C_3D_DP(Nx, Ny, Nz,reinterpret_cast<DP*>(Xr.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(X.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG | FFTW_MPI_TRANSPOSED_IN);
        
        c2r_xyz_plan = FFTW_MPI_PLAN_DFT_C2R_3D_DP(Nx, Ny, Nz,reinterpret_cast<FFTW_COMPLEX_DP*>(X.data()), reinterpret_cast<DP*>(Xr.data()), MPI_COMM_WORLD, FFTW_PLAN_FLAG | FFTW_MPI_TRANSPOSED_OUT);
        
    }
    
    else if (Ny == 1) {
    	r2c_xz_plan = FFTW_MPI_PLAN_DFT_R2C_2D_DP(Nx, Nz, reinterpret_cast<DP*>(plane_xz.data()),reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
        
        c2r_xz_plan = FFTW_MPI_PLAN_DFT_C2R_2D_DP(Nx, Nz,reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()), reinterpret_cast<DP*>(plane_xz.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
    }
}
//*********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_FFFW_SLAB(Array<complx,3> Ar)
{
	Ar(Range::all(),Range::all(),Nz/2) = 0.0;
}

void SpectralTransform::Zero_pad_last_col_FFFW_SLAB(Array<complx,2> Ar)
{
	Ar(Range::all(),Nz/2) = 0.0;
}

//*********************************************************************************************

void SpectralTransform::Norm_FFFW_SLAB(Array<complx,3> A)
{
	A = A/(DP(Nx) * DP(Ny) * DP(Nz));
}

void SpectralTransform::Norm_FFFW_SLAB(Array<complx,2> A)
{
	A = A/(DP(Nx) * DP(Nz));
}

//*************************************************************************************


void SpectralTransform::Forward_transform_FFFW_SLAB(Array<complx,3> Ar, Array<complx,3> A)
{
    FTr2c_xyz(Ar, A);
}
        

void SpectralTransform::Forward_transform_FFFW_SLAB(Array<complx,2> Ar, Array<complx,2> A)
{
    FTr2c_xz(Ar, A);
}

//*********************************************************************************************

void SpectralTransform::Inverse_transform_FFFW_SLAB(Array<complx,3> A, Array<complx,3> Ar)
{
    FTc2r_xyz(A, Ar);
}


//**************************************************************************************


void SpectralTransform::Inverse_transform_FFFW_SLAB(Array<complx,2> A, Array<complx,2> Ar)
{
	FTc2r_xz(A, Ar);
}




//********************************	End of four_tr.cc *****************************************


