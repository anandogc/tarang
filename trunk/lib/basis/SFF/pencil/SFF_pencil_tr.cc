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

#include "SFF_pencil.h"
#include "SFF_pencil_inline.h"

//*********************************************************************************************

// Ar: yxz, A: xyz
void SFF_PENCIL::Forward_transform(Array<Real,3> Ar, Array<Complex,3> A)
{
	if ( Ny>1 ) {
		fftk.Forward_transform("SFF", Ar, A);
	}
	else if (Ny == 1) {
		fftk.Forward_transform("SFF", Ar(Range::all(),0,Range::all()), A(Range::all(),0,Range::all()));
	}

}


//*********************************************************************************************


void SFF_PENCIL::Inverse_transform(Array<Complex,3> A, Array<Real,3> Ar)
{
	global.temp_array.X_transform = A;

    if ( Ny>1 ) {
		fftk.Inverse_transform("SFF", global.temp_array.X_transform, Ar);
	}
	else if (Ny == 1) {
		fftk.Inverse_transform("SFF", global.temp_array.X_transform(Range::all(),0,Range::all()), Ar(Range::all(),0,Range::all()));
	}
}


/**********************************************************************************************

	Derivative along x;  Bk_COS[i]=pi*i1*Ak_SIN[i]; Bk_SIN[]=-pi*i1*Ak_COS[i]

***********************************************************************************************/

void SFF_PENCIL::Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kx;
	
	for (int lx=0; lx<maxlx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
    
       if (global.program.sincostr_switch[0] == 'S')
            B(lx,Range::all(),Range::all()) = Kx * A(lx,Range::all(),Range::all());
		
         else if (global.program.sincostr_switch[0] == 'C')
            B(lx,Range::all(),Range::all()) = (-Kx) * A(lx,Range::all(),Range::all());
	}
}


void SFF_PENCIL::Add_Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kx;
	
	for (int lx=0; lx<maxlx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		
		if (global.program.sincostr_switch[0] == 'S')
            B(lx,Range::all(),Range::all()) += Kx * A(lx,Range::all(),Range::all());
		
		else if (global.program.sincostr_switch[0] == 'C')
            B(lx,Range::all(),Range::all()) += (-Kx) * A(lx,Range::all(),Range::all());
	}
}


void  SFF_PENCIL::Xderiv(Array<Real,3> A, Array<Real,3> B)
{
	cerr << "This is not defined for this basis. "<<endl;
}
/**********************************************************************************************

		Derivative along y;  B(k) = i*ky*A(k)

***********************************************************************************************/

// Note: In the first half- ky=i2;
// In the second half- i2=0:Ny/-1; fftw-index=(Ny/2 +1+i2); FT-index=fftw-index-N=(i2+1-Ny/2)
void SFF_PENCIL::Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
  
	for (int ly=0; ly<maxly; ly++) {
		Ky = Get_ky(ly)*kfactor[2];
		
		B(Range::all(),ly,Range::all()) = Complex(0, Ky)* (A(Range::all(),ly,Range::all()));
	}
}


void SFF_PENCIL::Add_Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
	
	for (int ly=0; ly<maxly; ly++) {
		Ky = Get_ky(ly)*kfactor[2];
		
		B(Range::all(),ly,Range::all()) += Complex(0, Ky)* (A(Range::all(),ly,Range::all()));
	}
}


/**********************************************************************************************

		Derivative along z; Derivative along z

***********************************************************************************************/


void SFF_PENCIL::Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kz;
	
	for (int lz=0; lz<maxlz; lz++) {
		Kz = Get_kz(lz)*kfactor[3];
		
		B(Range::all(),Range::all(),lz) = Complex(0, Kz)*(A(Range::all(),Range::all(),lz));
	}    
}



void SFF_PENCIL::Add_Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kz;
	
	for (int lz=0; lz<maxlz; lz++) {
		Kz = Get_kz(lz)*kfactor[3];
		
		B(Range::all(),Range::all(),lz) += Complex(0, Kz)*(A(Range::all(),Range::all(),lz));
	}
}


/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/

void SFF_PENCIL::Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kxsqr;    // Kx^2
	Real Kxysqr;
	Real Ksqr;
	
	for (int lx=0; lx<maxlx; lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
		
		for (int ly=0; ly<maxly; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);
			
			for (int lz=0; lz<maxlz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

				B(lx,ly,lz) = (-factor*Ksqr)*A(lx,ly,lz);
			}
		}
	}
}



/**********************************************************************************************
 
 B = B - factor*Laplacian(A) = B + factor*K^2 A
 
 ***********************************************************************************************/

void SFF_PENCIL::Subtract_Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kxsqr;    // Kx^2
	Real Kxysqr;
	Real Ksqr;
	
	for (int lx=0; lx<maxlx; lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
		
		for (int ly=0; ly<maxly; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);
			
			for (int lz=0; lz<maxlz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

				B(lx,ly,lz) += (factor*Ksqr)*A(lx,ly,lz);
			}
		}
	}
}



//******************************** End of SFF_pencil_tr.cc *****************************************
