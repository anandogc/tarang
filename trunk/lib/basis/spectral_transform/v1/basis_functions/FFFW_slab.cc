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
        local_Ny = Ny;
        local_Nz = (Nz/2+1);
        
        local_Nx_start = my_id * local_Nz;
        local_Ny_start = 0;
        local_Nz_start = 0;

        X_3d.resize(local_Nx,Ny,Nz/2+1);
        Xr_3d.resize(local_Nx,Ny,Nz+2);
    }
    else if (Ny == 1) {
        local_Nx = Nx/numprocs;
        local_Ny = 1;
        local_Nz = (Nz/2+1);
        
        local_Nx_start = my_id * local_Nx;
        local_Ny_start = 0;
        local_Nz_start = my_id * local_Nz;

        X_2d.resize(local_Nx,Nz/2+1);
        Xr_2d.resize(local_Nx,Nz);
    }

    //Initialize plans
		
    if (Ny > 1) {
        c2r_xyz_plan = FFTW_MPI_PLAN_DFT_C2R_3D(Nx, Ny, Nz,reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), reinterpret_cast<DP*>(Xr_3d.data()), MPI_COMM_WORLD, FFTW_PLAN_FLAG);

        r2c_xyz_plan = FFTW_MPI_PLAN_DFT_R2C_3D(Nx, Ny, Nz,reinterpret_cast<DP*>(Xr_3d.data()), reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
    }
    
    else if (Ny == 1) {
        c2r_xz_plan = FFTW_MPI_PLAN_DFT_C2R_2D(Nx, Nz,reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()), reinterpret_cast<DP*>(Xr_2d.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
        
        r2c_xz_plan = FFTW_MPI_PLAN_DFT_R2C_2D(Nx, Nz, reinterpret_cast<DP*>(Xr_2d.data()),reinterpret_cast<FFTW_COMPLEX*>(X_2d.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
    }
}
//*********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_FFFW_SLAB(Array<DP,3> Ar)
{
	Ar(Nz,Range::all(),Range::all()) = 0.0;
}

void SpectralTransform::Zero_pad_last_col_FFFW_SLAB(Array<DP,2> Ar)
{
	Ar(Nz,Range::all()) = 0.0;
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


void SpectralTransform::Forward_transform_FFFW_SLAB(Array<DP,3> Ar, Array<complx,3> A)
{
    FTr2c_xyz(Ar, A);
    Norm_FFFW_SLAB(A);
}
        
//*********************************************************************************************

void SpectralTransform::Inverse_transform_FFFW_SLAB(Array<complx,3> A, Array<DP,3> Ar)
{
    FTc2r_xyz(A, Ar);
}

//*********************************************************************************************

void SpectralTransform::Forward_transform_FFFW_SLAB(Array<DP,2> Ar, Array<complx,2> A)
{
    FTr2c_xz(Ar, A);
    Norm_FFFW_SLAB(A);
}

//**************************************************************************************


void SpectralTransform::Inverse_transform_FFFW_SLAB(Array<complx,2> A, Array<DP,2> Ar)
{
	FTc2r_xz(A, Ar);
}




//********************************	End of four_tr.cc *****************************************


