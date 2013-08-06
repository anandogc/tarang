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
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		ArrayIFFT
 */ 


#include "spectral_transform.h"


//*********************************************************************************************

// 3D: FFT original
void SpectralTransform::FTr2c_xyz(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_xyz_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::FTc2r_xyz(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_xyz_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

// 2D 
void SpectralTransform::FTr2c_xz(Array<complx,2> Plane)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_xz_plan, reinterpret_cast<DP*>(Plane.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Plane.data()));
}

void SpectralTransform::FTr2c_yz(Array<complx,2> Plane)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_yz_plan, reinterpret_cast<DP*>(Plane.data()),  reinterpret_cast<FFTW_COMPLEX_DP*>(Plane.data()));
}

void SpectralTransform::FTc2r_xz(Array<complx,2> Plane)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_xz_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(Plane.data()), reinterpret_cast<DP*>(Plane.data()));
}

void SpectralTransform::FTc2r_yz(Array<complx,2> Plane)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_yz_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(Plane.data()), reinterpret_cast<DP*>(Plane.data())); 
}

// 1D 
// c2c
void SpectralTransform::FTc2c_x(Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_DP(c2c_x_forward_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

void SpectralTransform::FTc2c_y(Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_DP(c2c_y_forward_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

void SpectralTransform::IFTc2c_x(Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_DP(c2c_x_inverse_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

void SpectralTransform::IFTc2c_y(Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_DP(c2c_y_inverse_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

// r2c, c2c FFF_PENCIL
void SpectralTransform::FTr2c_z(Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_z_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

void SpectralTransform::FTc2r_z(Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_z_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));
}

// sin/cos
void SpectralTransform::Sintr_x(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(sintr_x_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}

void SpectralTransform::Sintr_y(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(sintr_y_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}

void SpectralTransform::Sintr_z(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(sintr_z_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}

void SpectralTransform::Costr_x(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(costr_x_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}

void SpectralTransform::Costr_y(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(costr_y_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}

void SpectralTransform::Costr_z(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(costr_z_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}

// SpectralTransform::inverse
void SpectralTransform::ISintr_x(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(isintr_x_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}

void SpectralTransform::ISintr_y(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(isintr_y_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}	

void SpectralTransform::ISintr_z(Array<DP,1> Col)
{
	FFTW_EXECUTE_R2R_DP(isintr_z_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  
}	

void SpectralTransform::ICostr_x(Array<DP,1> Col)
{			 
	FFTW_EXECUTE_R2R_DP(icostr_x_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  		 
}

void SpectralTransform::ICostr_y(Array<DP,1> Col)
{			 
	FFTW_EXECUTE_R2R_DP(icostr_y_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  		 
}	

void SpectralTransform::ICostr_z(Array<DP,1> Col)
{			 
	FFTW_EXECUTE_R2R_DP(icostr_z_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));  		 
}	

// Chebyshev transform via fftw cosine transform
void SpectralTransform::Chebyshevtr_x(Array<DP,1> Col)
{
    FFTW_EXECUTE_R2R_DP(Chebyshevtr_x_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));
    // Postprocess (see definition)
    Col(Range(1,Nx-2)) *= 2.0;
    
}

// Inverse Chebyshev transform via fftw cosine transform
void SpectralTransform::IChebyshevtr_x(Array<DP,1> Col)
{
	// Preprocess (see definition)
    Col(Range(1,Nx-2)) /= 2.0;
    
	FFTW_EXECUTE_R2R_DP(Chebyshevtr_x_plan, reinterpret_cast<DP*>(Col.data()), reinterpret_cast<DP*>(Col.data()));
}


// GP
// 3D: FFT original
void SpectralTransform::FTc2c_xyz(Array<complx,3> A)
{
	FFTW_MPI_EXECUTE_DFT_DP(c2c_forward_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::IFTc2c_xyz(Array<complx,3> A)
{
	FFTW_MPI_EXECUTE_DFT_DP(c2c_inverse_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}



//****************************************************************************************
// After sin transform of A, shift right along the shift direction.
// A(Ny/Nx, Nx/Ny, Nz/2+1) except for SSS basis for which Ar(Ny, Nx, Nz/2)
// Shift direction != Xa
void SpectralTransform::ArrayShiftRight_SLAB(Array<complx,3> A, char Xa, char Xb, char Xc, char shift_dirn)
{
    int Nb;
    
    if (Xb == 'X')
        Nb = Nx;
    
    else if (Xb == 'Y')
        Nb = Ny;
    
    if (shift_dirn == Xb) {
        for (int ly=(Nb-2);  ly>=0; ly--) 
            A(Range::all(),ly+1,Range::all()) = A(Range::all(),ly,Range::all());
        
        A(Range::all(),0,Range::all()) = 0.0;
    }
    
    else if (shift_dirn == Xc) {  // for SSS basis: Xc = Z
        for  (int lz=(Nz/2-1); lz>=1; lz--) {
            imag(A(Range::all(),Range::all(),lz)) = real(A(Range::all(),Range::all(),lz));
            real(A(Range::all(),Range::all(),lz)) = imag(A(Range::all(),Range::all(),lz-1));
        }
        
        imag(A(Range::all(),Range::all(),0)) = real(A(Range::all(),Range::all(),0));
        real(A(Range::all(),Range::all(),0)) = 0.0;
    }
}
 
//****************************************************************************************
// Before sin transform, shift left
void SpectralTransform::ArrayShiftLeft_SLAB(Array<complx,3> A, char Xa, char Xb, char Xc, char shift_dirn)
{
    int Nb;
    
   if (Xb == 'X')
        Nb = Nx;
    
    else if(Xb == 'Y')
        Nb = Ny;
   
    if (shift_dirn == Xb) {
        for (int ly=1;  ly<=(Nb-1); ly++) 
            A(Range::all(),ly-1,Range::all()) = A(Range::all(),ly,Range::all());
        
        A(Range::all(),Nb-1,Range::all()) = 0.0;
    }
    
    else if (shift_dirn == Xc) {  // for SSS basis: Xc = Z
        for  (int lz=0; lz<=(Nz/2-2); lz++) {
            real(A(Range::all(),Range::all(),lz)) = imag(A(Range::all(),Range::all(),lz));
            imag(A(Range::all(),Range::all(),lz)) = real(A(Range::all(),Range::all(),lz+1));
        }
        
        real(A(Range::all(),Range::all(),Nz/2-1)) = imag(A(Range::all(),Range::all(),Nz/2-1));
        imag(A(Range::all(),Range::all(),Nz/2-1)) = 0.0;
    }
}


//****************************************************************************************
// 2D
// Shift direction = second direction
void SpectralTransform::ArrayShiftRight_SLAB(Array<complx,2> A, char shift_dirn)
{
    if (shift_dirn == 'X') { // A((Nz/2+1)/p, Nx)
        for (int lx=(Nx-2);  lx>=0; lx--) 
            A(Range::all(),lx+1) = A(Range::all(),lx);
        
        A(Range::all(),0) = 0.0;
    }
    
    else if (shift_dirn == 'Z') {  // for SSS basis: Xc = Z
        for  (int lz=(Nz/2-1); lz>=1; lz--) {
            imag(A(Range::all(),lz)) = real(A(Range::all(),lz));
            real(A(Range::all(),lz)) = imag(A(Range::all(),lz-1));
        }
        
        imag(A(Range::all(),0)) = real(A(Range::all(),0));
        real(A(Range::all(),0)) = 0.0;
    }
}

//****************************************************************************************
// 2D: shift_dirn = second dirn
void SpectralTransform::ArrayShiftLeft_SLAB(Array<complx,2> A, char shift_dirn)
{
    
    if (shift_dirn == 'X') {
        for (int lx=1;  lx<=(Nx-1); lx++) 
           A(Range::all(),lx-1) = A(Range::all(),lx);
        
        A(Range::all(),Nx-1) = 0.0;
    }
    
    else if (shift_dirn == 'Z') {  // for SSS basis: Xc = Z
        for  (int lz=0; lz<=(Nz/2-2); lz++) {
            real(A(Range::all(),lz)) = imag(A(Range::all(),lz));
            imag(A(Range::all(),lz)) = real(A(Range::all(),lz+1));
        }
        
        real(A(Range::all(),Nz/2-1)) = imag(A(Range::all(),Nz/2-1));
        imag(A(Range::all(),Nz/2-1)) = 0.0;
    }
}


//********************************	End of four_tr.cc *****************************************


