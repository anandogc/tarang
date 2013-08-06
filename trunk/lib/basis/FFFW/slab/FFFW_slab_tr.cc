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


#include "FFFW_slab.h"
#include "FFFW_slab_inline.h"

	
//*********************************************************************************************
/* void FFFW_SLAB::Init_fftw_plan(Array<complx,3> A)
{
		//	Transform::Init_fftw_plan_all(A);
}*/
//*********************************************************************************************

void FFFW_SLAB::Zero_pad_last_plane(Array<complx,3> Ar)
{
	Ar(Range::all(),Range::all(),Nz/2) = 0.0;
}

//*********************************************************************************************

void FFFW_SLAB::Norm(Array<complx,3> A) 
{
	A = A/(DP(Nx) * DP(Ny) * DP(Nz));
}


//*********************************************************************************************

// Ar: yxz, A: xyz
void FFFW_SLAB::Forward_transform_array_transpose_order(Array<complx,3> Ar, Array<complx,3> A)
{
	if (Ny > 1)  {
		Zero_pad_last_plane(Ar);							 // Zero_pad the row=Nz
		
		for (int ly=0; ly<local_Ny; ly++) {
			global.temp_array.plane_xz = Ar(ly, Range::all(), Range::all());
            Transform::FTr2c_xz(global.temp_array.plane_xz);
			Ar(ly, Range::all(), Range::all()) = global.temp_array.plane_xz;
		}
		
		Transform::Transpose_array_SLAB(Ar, A, 'Y', 'X', 'Z', 'Z');
		
		for (int lx=0; lx<local_Nx; lx++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				global.temp_array.col_y = A(lx, Range::all(), lz);
                Transform::FTc2c_y(global.temp_array.col_y);
				A(lx, Range::all(), lz) = global.temp_array.col_y;
			}
		
		Norm(A);
	}
	
	else if (Ny == 1) {
		cout << "ERROR:  ArrayFFT_FFFW_transpose_order for 2D array not" 
			 <<		"allowed in FOUR basis " << endl;
		exit(1);
	}
}


//*********************************************************************************************

void FFFW_SLAB::Forward_transform_array(Array<complx,3> A) 
{
	
	if (global.fft.fftw3D_switch == true) {
		Zero_pad_last_plane(A);
		
		if (Ny > 1){ 
            Transform::FTr2c_xyz(A);
        }
		
		else if (Ny == 1) {
			global.temp_array.plane_xz = A(Range::all(), 0, Range::all());			
            Transform::FTr2c_xz(global.temp_array.plane_xz); 
			A(Range::all(), 0, Range::all()) = global.temp_array.plane_xz;
		}
		
		Norm(A);
	}
	
	// not original switch.. split in 2D and 1D ffts.
	else {
		if (Ny > 1) {
			Transform::Transpose_array_SLAB(A, global.temp_array.Ar, 'X', 'Y', 'Z', 'Z');
			Forward_transform_array_transpose_order(global.temp_array.Ar, A);	// Norm done here
		}
		
		else if (Ny == 1) {
			cout << "ERROR: USE ORIGINAL FFT FOR 2D. " << endl;
            exit(1);
		}
	}
	
}

//*********************************************************************************************


void FFFW_SLAB::Inverse_transform_array_transpose_order(Array<complx,3> A, Array<complx,3> Ar)
{
	global.temp_array.X2 = A;
	
	if (Ny > 1) {
		for (int lx=0; lx<local_Nx; lx++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				global.temp_array.col_y = A(lx, Range::all(), lz);
				Transform::IFTc2c_y(global.temp_array.col_y);
				A(lx, Range::all(), lz) = global.temp_array.col_y;
			}
		
		Transform::Transpose_array_SLAB(A, Ar, 'X', 'Y', 'Z', 'Z');
		
		for (int ly=0; ly < local_Ny; ly++)	{
			global.temp_array.plane_xz = Ar(ly, Range::all(), Range::all());
			Transform::FTc2r_xz(global.temp_array.plane_xz);
			Ar(ly, Range::all(), Range::all()) = global.temp_array.plane_xz;
		}
		
		Zero_pad_last_plane(Ar);
	}
	
	else if (Ny == 1) {
		cout << "ERROR:  ArrayIFFT_FFFW_transpose_order for 2D array not allowed in FOUR basis " << endl;
		exit(1);
	}
	
	A = global.temp_array.X2;
}


//*********************************************************************************************


void FFFW_SLAB::Inverse_transform_array(Array<complx,3> A)
{
	
	if (global.fft.fftw3D_switch == true) {
		if (Ny > 1)
			Transform::FTc2r_xyz(A);
		
		else if (Ny == 1) {
			global.temp_array.plane_xz = A(Range::all(), 0, Range::all());
            Transform::FTc2r_xz(global.temp_array.plane_xz);
			A(Range::all(), 0, Range::all()) = global.temp_array.plane_xz;
		}
	}
	
    // By splitting it up...
	else {
		if (Ny > 1) {
			Inverse_transform_array_transpose_order(A, global.temp_array.Ar);
			Transform::Transpose_array_SLAB(global.temp_array.Ar, A, 'Y', 'X', 'Z', 'Z'); 
		}
		
		else if (Ny == 1) {
			cout << "ERROR: USE ORIGINAL FFT FOR 2D. " << endl;
            exit(1);
		}
	}
	
}

//*********************************************************************************************

void FFFW_SLAB::Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kx;
	
	for (int lx = 0; lx < local_Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		B(lx,Range::all(),Range::all()) = complex<DP>(0, Kx)* (A(lx,Range::all(),Range::all())); 	
	}
}	

void  FFFW_SLAB::Xderiv(Array<DP,3> A, Array<DP,3> B)
{
	cerr << "This is not defined for this basis. "<<endl;
}


//*********************************************************************************************



void FFFW_SLAB::Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Ky;
	
	if (Ny > 1)
		for (int ly=0; ly<Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			B(Range::all(),ly,Range::all()) = complex<DP>(0, Ky)* (A(Range::all(),ly,Range::all())); 
		}
	
	else 
		B = 0.0;
}

//*********************************************************************************************


void FFFW_SLAB::Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
    
    for (int lz=0; lz<=Nz/2; lz++) {
		Kz = lz*kfactor[3];
		B(Range::all(),Range::all(),lz) = complex<DP>(0, Kz)*(A(Range::all(),Range::all(),lz)); 	
	}
}

//********************************	End of four_tr.cc *****************************************


