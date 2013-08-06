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

#include "ChFF_slab.h"
#include "ChFF_slab_inline.h"



/**********************************************************************************************
 
 Forward SFT (A)  = Ar
 FT along perp dirn and SIN transform along x dirn of A
 
 ***********************************************************************************************/

void ChFF_SLAB::Forward_transform(Array<DP,3> Ar, Array<complx,3> A)
{
	if (Ny > 1)
        spectralTransform.Forward_transform_ChFF_SLAB(Ar, A);
    
    else if (Ny == 1)
        spectralTransform.Forward_transform_ChFF_SLAB(Ar(0,Range::all(),Range::all()), A(0,Range::all(),Range::all()));
	
}


/**********************************************************************************************
 
 Inverse SFT (A)  = Ar
 IFT along perp dirn and SIN transform along x dirn of A
 
 ***********************************************************************************************/


void ChFF_SLAB::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar)
{
    if (Ny > 1)
        spectralTransform.Inverse_transform_ChFF_SLAB(A, Ar);
    
    else if (Ny == 1)
		spectralTransform.Inverse_transform_ChFF_SLAB(A(0,Range::all(),Range::all()), Ar(0,Range::all(),Range::all()));
}



/**********************************************************************************************

	Derivative along x;  B(k) = Tr(derivative(A))

***********************************************************************************************/

void  ChFF_SLAB::Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	
    int kx;
	DP Kx;

    B(Range::all(), Range::all(), Nx-1) = 0;
    B(Range::all(), Range::all(), Nx-2) = (TWO*(Nx-1))*A(Range::all(), Range::all(), Nx-1);

	for (int lx=(Nx-3); lx>=1; lx--) 	{
		B(Range::all(),Range::all(), lx) = (TWO*(lx+1))*A(Range::all(),Range::all(),lx+1) + B(Range::all(),Range::all(),lx+2);
	}

    B(Range::all(),Range::all(),0) =  A(Range::all(),Range::all(),1) + B(Range::all(),Range::all(),2)/TWO;
    
    B = kfactor[1]*B;
}


void  ChFF_SLAB::Add_Xderiv(Array<complx,3> A, Array<complx,3> B)
{
 
	Xderiv(A, global.temp_array.X2);
	B += global.temp_array.X2;
}



void  ChFF_SLAB::Xderiv(Array<DP,3> A, Array<DP,3> B)
{
    int kx;
	DP Kx;
    
    B(Range::all(), Range::all(), Nx-1) = 0;
    B(Range::all(), Range::all(), Nx-2) = TWO*(Nx-1)*A(Range::all(), Range::all(), Nx-1);
    
	for (int lx=(Nx-3); lx>=1; lx--) 	{
		B(Range::all(),Range::all(),lx) = TWO*(lx+1)*A(Range::all(),Range::all(),lx+1) + B(Range::all(),Range::all(),lx+2);
	}
    
    B(Range::all(),Range::all(),0) = A(Range::all(),Range::all(),1) + B(Range::all(),Range::all(),2)/TWO;
    
    B = kfactor[1]*B;
}



/**********************************************************************************************

		Derivative along y;  B(k) = i*ky*A(k)

***********************************************************************************************/

// Note: In the first half- ky=i2;
// In the second half- i2=0:Ny/-1; fftw-index=(Ny/2 +1+i2); FT-index=fftw-index-N=(i2+1-Ny/2)
void ChFF_SLAB::Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Ky;
	
	if (Ny > 1) 
		for (int ly=0; ly<local_Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];

			B(ly,Range::all(),Range::all()) = complex<DP>(0, Ky)* (A(ly,Range::all(),Range::all()));
		}
	
	else // Ny = 1
		B = 0;
}

void ChFF_SLAB::Add_Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Ky;
	
	if (Ny > 1)
		for (int ly=0; ly<local_Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			
			B(ly,Range::all(),Range::all()) += complex<DP>(0, Ky)* (A(ly,Range::all(),Range::all()));
		}
}

/**********************************************************************************************

		Derivative along z; Derivative along z

***********************************************************************************************/


void ChFF_SLAB::Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
	
	for (int lz=0; lz<local_Nz; lz++) {
		Kz = Get_kz(lz)*kfactor[3];
		
		B(Range::all(),lz,Range::all()) = complex<DP>(0, Kz)*(A(Range::all(),lz,Range::all()));
	}   
}



void ChFF_SLAB::Add_Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
	
	for (int lz=0; lz<local_Nz; lz++) {
		Kz = Get_kz(lz)*kfactor[3];
		
		B(Range::all(),lz,Range::all()) += complex<DP>(0, Kz)*(A(Range::all(),lz,Range::all()));
	}
}


/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/


void ChFF_SLAB::Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
{
	Xderiv(A, global.temp_array.X2);
	Xderiv(global.temp_array.X2, B);
	
	if (abs(factor-1.0) > MYEPS)
		B *= factor;
	
	DP Kperp_sqr, Kperp_sqr_factor;
	
	for (int ly=0; ly<A.extent(0); ly++) {
		Kperp_sqr = my_pow(Get_ky(ly)*kfactor[2],2);
		
        for (int lz=0; lz<A.extent(1); lz++) {
			Kperp_sqr_factor = factor*(Kperp_sqr+my_pow(Get_lz(lz)*kfactor[3],2));
			B(ly,lz,Range::all()) -= Kperp_sqr_factor*A(ly,lz,Range::all());
		}
	}
}

/**********************************************************************************************
 
  B = B - factor*Laplacian(A) = B-factor*Dxx(A) - factor*K_perp^2 A
 
 ***********************************************************************************************/


void ChFF_SLAB::Subtract_Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
{
	Xderiv(A, global.temp_array.X);
	Xderiv(global.temp_array.X, global.temp_array.X2);  // could possibly use X for X2--- test it
	
	B -= factor*global.temp_array.X2;
	
	DP Kperp_sqr, Kperp_sqr_factor;
	
	for (int ly=0; ly<A.extent(0); ly++) {
		Kperp_sqr = my_pow(Get_ky(ly)*kfactor[2],2);
		
        for (int lz=0; lz<A.extent(1); lz++) {
			Kperp_sqr_factor = factor*(Kperp_sqr+my_pow(Get_lz(lz)*kfactor[3],2));
			B(ly,lz,Range::all()) += Kperp_sqr_factor*A(ly,lz,Range::all());
		}
	}

}


//******************************** End of SFF_slab_tr.cc  *****************************************
