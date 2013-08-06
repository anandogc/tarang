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


#include "CFFF_slab.h"
#include "CFFF_slab_inline.h"

	
//*********************************************************************************************
/* void CFFF_SLAB::Init_fftw_plan(Array<complx,3> A)
{
		//	Transform::Init_fftw_plan_all(A);
}*/
//*********************************************************************************************

void CFFF_SLAB::Zero_pad_last_plane(Array<complx,3> Ar)
{
}

//*********************************************************************************************

void CFFF_SLAB::Norm(Array<complx,3> A)
{
	A = A/(DP(Nx) * DP(Ny) * DP(Nz));
}


//*********************************************************************************************

// Ar: yxz, A: xyz
void CFFF_SLAB::Forward_transform_array_transpose_order(Array<complx,3> Ar, Array<complx,3> A)
{
}


//*********************************************************************************************

void CFFF_SLAB::Forward_transform_array(Array<complx,3> A)
{
	
	if (global.fft.fftw3D_switch == true) {
        Transform::FTc2c_xyz(A);
        Norm(A);
	}
	
	// not original switch.. split in 2D and 1D ffts.
	else {
			cout << "ERROR: USE ORIGINAL FFTW " << endl;
            exit(1);
    }
	
}

//*********************************************************************************************


void CFFF_SLAB::Inverse_transform_array_transpose_order(Array<complx,3> A, Array<complx,3> Ar)
{
}


//*********************************************************************************************


void CFFF_SLAB::Inverse_transform_array(Array<complx,3> A)
{
	
	if (global.fft.fftw3D_switch == true) {
			Transform::IFTc2c_xyz(A);
	}
	
    // By splitting it up...
	else {
        cout << "ERROR: USE ORIGINAL FFT FOR 2D. " << endl;
        exit(1);
	}
	
}

//*********************************************************************************************

void CFFF_SLAB::Xderiv(Array<complx,3> A, Array<complx,3> B)
{
}	

void  CFFF_SLAB::Xderiv(Array<DP,3> A, Array<DP,3> B)
{
	cerr << "This is not defined for this basis. "<<endl;
}


//*********************************************************************************************



void CFFF_SLAB::Yderiv(Array<complx,3> A, Array<complx,3> B)
{
}

//*********************************************************************************************


void CFFF_SLAB::Zderiv(Array<complx,3> A, Array<complx,3> B)
{
}

//********************************	End of four_tr.cc *****************************************


