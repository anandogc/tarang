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

// Ar: yxz, A: xyz
void FFFW_SLAB::Forward_transform(Array<Real,3> Ar, Array<Complex,3> A)
{
    if (Ny > 1) 
        spectralTransform.Forward_transform(Ar,A);
    
    else if (Ny == 1){
       spectralTransform.Forward_transform(Ar(Range::all(),0,Range::all()),A(Range::all(),0,Range::all()));
    }
}


//*********************************************************************************************


void FFFW_SLAB::Inverse_transform(Array<Complex,3> A, Array<Real,3> Ar)
{
	global.temp_array.X_transform = A;

    if (Ny > 1)
        spectralTransform.Inverse_transform(global.temp_array.X_transform, Ar);
    
    else if (Ny == 1){
    	spectralTransform.Inverse_transform(global.temp_array.X_transform(Range::all(),0,Range::all()),Ar(Range::all(),0,Range::all()));
   }
}



//*********************************************************************************************

void FFFW_SLAB::Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kx;
	
	for (int lx = 0; lx < local_Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		B(lx,Range::all(),Range::all()) = Complex(0, Kx)* (A(lx,Range::all(),Range::all())); 	
	}
}

void FFFW_SLAB::Add_Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kx;
	
	for (int lx = 0; lx < local_Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		B(lx,Range::all(),Range::all()) += Complex(0, Kx)* (A(lx,Range::all(),Range::all()));
	}
}

void  FFFW_SLAB::Xderiv(Array<Real,3> A, Array<Real,3> B)
{
	if (master) cerr <<  "Xderiv(real array)  is not defined for this basis. "<<endl;
}



//*********************************************************************************************



void FFFW_SLAB::Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
	
	if (Ny > 1)
		for (int ly=0; ly<Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			B(Range::all(),ly,Range::all()) = Complex(0, Ky)* (A(Range::all(),ly,Range::all())); 
		}
	
	else 
		B = 0.0;
}



void FFFW_SLAB::Add_Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
	
	if (Ny > 1)
		for (int ly=0; ly<Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			B(Range::all(),ly,Range::all()) += Complex(0, Ky)* (A(Range::all(),ly,Range::all()));
		}
}

//*********************************************************************************************


void FFFW_SLAB::Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kz;
    
    for (int lz=0; lz<=Nz/2; lz++) {
		Kz = lz*kfactor[3];
		B(Range::all(),Range::all(),lz) = Complex(0, Kz)*(A(Range::all(),Range::all(),lz)); 	
	}
}


void FFFW_SLAB::Add_Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kz;
    
    for (int lz=0; lz<=Nz/2; lz++) {
		Kz = lz*kfactor[3];
		B(Range::all(),Range::all(),lz) += Complex(0, Kz)*(A(Range::all(),Range::all(),lz));
	}
}

/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/

void FFFW_SLAB::Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	
	Real Ksqr;
	
	for (int lx=0; lx<A.extent(0); lx++) {
		Ksqr =  my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<A.extent(1); ly++) {
			Ksqr += my_pow(Get_ky(ly)*kfactor[2],2);
			
	        for (int lz=0; lz<A.extent(2); lz++) {
				Ksqr += my_pow(Get_lz(lz)*kfactor[3],2);
				
				B(lx,ly,lz) = (-factor*Ksqr)*A(lx,ly,lz);
			}
		}
	}
}

/**********************************************************************************************
 
 B = B - factor*Laplacian(A) = B + factor*K^2 A
 
 ***********************************************************************************************/

void FFFW_SLAB::Subtract_Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	
	Real Ksqr, Ksqr_factor;
	
    for (int lx=0; lx<A.extent(0); lx++) {
		Ksqr = my_pow(Get_lx(lx)*kfactor[1],2);

		for (int ly=0; ly<A.extent(1); ly++) {
			Ksqr += my_pow(Get_ky(ly)*kfactor[2],2);
			
	        for (int lz=0; lz<A.extent(2); lz++) {
				Ksqr_factor = factor * (Ksqr+my_pow(Get_kz(lz)*kfactor[3],2));
				
				B(lx,ly,lz) += Ksqr_factor*A(lx,ly,lz);
			}
		}
	}
}


//********************************	End of four_tr.cc *****************************************


