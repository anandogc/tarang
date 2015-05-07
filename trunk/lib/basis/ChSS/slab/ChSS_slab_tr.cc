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

#include "ChSS_slab.h"
#include "ChSS_slab_inline.h"



/**********************************************************************************************
 
 Forward SFT (A)  = Ar
 FT along perp dirn and SIN transform along x dirn of A
 
 ***********************************************************************************************/

void ChSS_SLAB::Forward_transform(Array<Real,3> Ar, Array<Complex,3> A)
{
	if (Ny > 1)
        spectralTransform.Forward_transform_ChSS_SLAB(Ar, A);
    
    else if (Ny == 1)
        spectralTransform.Forward_transform_ChSS_SLAB(Ar(0,Range::all(),Range::all()), A(0,Range::all(),Range::all()));
	
}


/**********************************************************************************************
 
 Inverse SFT (A)  = Ar
 IFT along perp dirn and SIN transform along x dirn of A
 
 ***********************************************************************************************/


void ChSS_SLAB::Inverse_transform(Array<Complex,3> A, Array<Real,3> Ar)
{
    if (Ny > 1)
        spectralTransform.Inverse_transform_ChSS_SLAB(A, Ar);
    
    else if (Ny == 1)
		spectralTransform.Inverse_transform_ChSS_SLAB(A(0,Range::all(),Range::all()), Ar(0,Range::all(),Range::all()));
}



/**********************************************************************************************

	Derivative along x;  B(k) = Tr(derivative(A))

***********************************************************************************************/

void  ChSS_SLAB::Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{	
    int kx;
	Real Kx;

    B(Range::all(), Range::all(), Nx-1) = 0;
    B(Range::all(), Range::all(), Nx-2) = Complex(2.0*(Nx-1),0)*A(Range::all(), Range::all(), Nx-1);

	for (int lx=(Nx-3); lx>=1; lx--) 	{
		B(Range::all(),Range::all(), lx) = (Complex(2*(lx+1),0)*A(Range::all(),Range::all(),lx+1)) + B(Range::all(),Range::all(),lx+2);
	}

    B(Range::all(),Range::all(),0) =  A(Range::all(),Range::all(),1) + B(Range::all(),Range::all(),2)*Complex(0.5,0);
    
    B = kfactor[1]*B;
}


void  ChSS_SLAB::Add_Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{
 
	Xderiv(A, global.temp_array.X2);
	B += global.temp_array.X2;
}



void  ChSS_SLAB::Xderiv(Array<Real,3> A, Array<Real,3> B)
{
    int kx;
	Real Kx;
    
    B(Range::all(), Range::all(), Nx-1) = 0;
    B(Range::all(), Range::all(), Nx-2) = 2.0*(Nx-1)*A(Range::all(), Range::all(), Nx-1);
    
	for (int lx=(Nx-3); lx>=1; lx--) 	{
		B(Range::all(),Range::all(),lx) = TWO*(lx+1)*A(Range::all(),Range::all(),lx+1) + B(Range::all(),Range::all(),lx+2);
	}
    
    B(Range::all(),Range::all(),0) = A(Range::all(),Range::all(),1) + B(Range::all(),Range::all(),2)/2;
    
    B = kfactor[1]*B;
}


/**********************************************************************************************

		Derivative along y;  B(k) = i*ky*A(k)

***********************************************************************************************/

// Note: In the first half- ky=i2;
// In the second half- i2=0:Ny/-1; fftw-index=(Ny/2 +1+i2); FT-index=fftw-index-N=(i2+1-Ny/2)
void ChSS_SLAB::Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
	
	if (Ny > 1) 
		for (int ly=0; ly<Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			
			if (global.program.sincostr_switch[1] == 'S')
				B(ly,Range::all(),Range::all()) = complex<Real>(Ky, 0)* (A(ly,Range::all(),Range::all()));
			
			else if (global.program.sincostr_switch[1] == 'C')
				B(ly,Range::all(),Range::all()) = complex<Real>(-Ky, 0)* (A(ly,Range::all(),Range::all()));
		}
	
	else // Ny = 1
		B = 0;
}

void ChSS_SLAB::Add_Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
	
	if (Ny > 1)
		for (int ly=0; ly<Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			
			if (global.program.sincostr_switch[1] == 'S')
				B(ly,Range::all(),Range::all()) += complex<Real>(Ky, 0)* (A(ly,Range::all(),Range::all()));
			
			else if (global.program.sincostr_switch[1] == 'C')
				B(ly,Range::all(),Range::all()) += complex<Real>(-Ky, 0)* (A(ly,Range::all(),Range::all()));
		}
}

/**********************************************************************************************

		Derivative along z; Derivative along z

***********************************************************************************************/


void ChSS_SLAB::Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kz;
	
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_real_array, neverDeleteData);
	Array<Real,3> B_real(reinterpret_cast<Real*>(B.data()), shape_real_array, neverDeleteData);
	
	for (int lz=0; lz<Nz; lz++) {
		Kz = lz*kfactor[3];
		
		if (global.program.sincostr_switch[2] == 'S')
			B_real(Range::all(),lz,Range::all()) = Kz*A_real(Range::all(),lz,Range::all());
		
		else if (global.program.sincostr_switch[2] == 'C')
			B_real(Range::all(),lz,Range::all()) = -Kz*A_real(Range::all(),lz,Range::all());
		
	}
}



void ChSS_SLAB::Add_Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_real_array, neverDeleteData);
	Array<Real,3> B_real(reinterpret_cast<Real*>(B.data()), shape_real_array, neverDeleteData);
	
	Real Kz;
	
	for (int lz=0; lz<Nz; lz++) {
		Kz = lz*kfactor[3];
		
		if (global.program.sincostr_switch[2] == 'S')
			B_real(Range::all(),lz,Range::all()) += Kz*A_real(Range::all(),lz,Range::all());
		
		else if (global.program.sincostr_switch[2] == 'C')
			B_real(Range::all(),lz,Range::all()) += (-Kz)*A_real(Range::all(),lz,Range::all());
		
	}
}

/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/


void ChSS_SLAB::Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	Xderiv(A, global.temp_array.X);
	Xderiv(global.temp_array.X, B);
	
	B *= factor;
	
	Real Kperp_sqr, Kperp_sqr_factor;
	
	for (int ly=0; ly<A.extent(0); ly++) {
		Kperp_sqr = my_pow(Get_ky(ly)*kfactor[2],2);
		
        for (int lz=0; lz<A.extent(1); lz++) {
			Kperp_sqr_factor = factor*(Kperp_sqr+my_pow(Get_lz(lz)*kfactor[3],2));
			B(ly,lz,Range::all()) -= Complex(Kperp_sqr_factor,0)*A(ly,lz,Range::all());
		}
	}
}


/**********************************************************************************************
 
  B = B - factor*Laplacian(A) = B-factor*Dxx(A) - factor*K_perp^2 A
 
 ***********************************************************************************************/


void ChSS_SLAB::Subtract_Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	Xderiv(A, global.temp_array.X);
	Xderiv(global.temp_array.X, global.temp_array.X2);  // could possibly use X for X2--- test it
	
	B -= factor*global.temp_array.X2;
	
	Real Kperp_sqr, Kperp_sqr_factor;
	
	for (int ly=0; ly<A.extent(0); ly++) {
		Kperp_sqr = my_pow(Get_ky(ly)*kfactor[2],2);
		
        for (int lz=0; lz<A.extent(1); lz++) {
			Kperp_sqr_factor = factor*(Kperp_sqr +my_pow(Get_lz(lz)*kfactor[3],2));
			B(ly,lz,Range::all()) += Complex(Kperp_sqr_factor,0)*A(ly,lz,Range::all());
		}
	}

}


//******************************** End of SFF_slab_tr.cc  *****************************************
