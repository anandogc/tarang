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
void FFFW_SLAB::Forward_transform(Array<DP,3> Ar, Array<complx,3> A)
{
    if (Ny > 1) 
        spectralTransform.Forward_transform_FFFW_SLAB(Ar,A);
    
    else if (Ny == 1){
     //   spectralTransform.Forward_transform_FFFW_SLAB(Ar(Range::all(),0,Range::all()),A(Range::all(),0,Range::all()));
    }
}


//*********************************************************************************************


void FFFW_SLAB::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar)
{
    if (Ny > 1)
        spectralTransform.Inverse_transform_FFFW_SLAB(A, Ar);
    
    else if (Ny == 1){
    //   spectralTransform.Inverse_transform_FFFW_SLAB(Ar(Range::all(),0,Range::all()),A(Range::all(),0,Range::all()));
   }
}



//*********************************************************************************************

void FFFW_SLAB::Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kx;
	
	for (int lx = 0; lx < Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		B(Range::all(),Range::all(),lx) = complx(0, Kx)* (A(Range::all(),Range::all(),lx)); 	
	}
}

void FFFW_SLAB::Add_Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kx;
	
	for (int lx = 0; lx < Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		B(Range::all(),Range::all(),lx) += complx(0, Kx)* (A(Range::all(),Range::all(),lx));
	}
}

void  FFFW_SLAB::Xderiv(Array<DP,3> A, Array<DP,3> B)
{
	if (master) cerr <<  "Xderiv(real array)  is not defined for this basis. "<<endl;
}



//*********************************************************************************************



void FFFW_SLAB::Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Ky;
	
	if (Ny > 1)
		for (int ly=0; ly<local_Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			B(ly,Range::all(),Range::all()) = complx(0, Ky)* (A(ly,Range::all(),Range::all())); 
		}
	
	else 
		B = 0.0;
}



void FFFW_SLAB::Add_Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Ky;
	
	if (Ny > 1)
		for (int ly=0; ly<local_Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			B(ly,Range::all(),Range::all()) += complx(0, Ky)* (A(ly,Range::all(),Range::all()));
		}
}

//*********************************************************************************************


void FFFW_SLAB::Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
    
    for (int lz=0; lz<=Nz/2; lz++) {
		Kz = lz*kfactor[3];
		B(Range::all(),lz,Range::all()) = complx(0, Kz)*(A(Range::all(),lz,Range::all())); 	
	}
}


void FFFW_SLAB::Add_Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
    
    for (int lz=0; lz<=Nz/2; lz++) {
		Kz = lz*kfactor[3];
		B(Range::all(),lz,Range::all()) += complx(0, Kz)*(A(Range::all(),lz,Range::all()));
	}
}

/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/

void FFFW_SLAB::Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
{
	
	DP Ksqr;
	
	for (int ly=0; ly<A.extent(0); ly++) {
		Ksqr = my_pow(Get_ky(ly)*kfactor[2],2);
		
        for (int lz=0; lz<A.extent(1); lz++) {
			Ksqr += my_pow(Get_lz(lz)*kfactor[3],2);
			
            for (int lx=0; lx<A.extent(2); lx++) {
				Ksqr +=  my_pow(Get_kx(lx)*kfactor[1],2);
				
				B(ly,lz,lx) = (-factor*Ksqr)*A(ly,lz,lx);
			}
		}
	}
}

/**********************************************************************************************
 
 B = B - factor*Laplacian(A) = B + factor*K^2 A
 
 ***********************************************************************************************/

void FFFW_SLAB::Subtract_Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
{
	
	DP Ksqr, Ksqr_factor;
	
	for (int ly=0; ly<A.extent(0); ly++) {
		Ksqr = my_pow(Get_ky(ly)*kfactor[2],2);
		
        for (int lz=0; lz<A.extent(1); lz++) {
			Ksqr += my_pow(Get_lz(lz)*kfactor[3],2);
			
            for (int lx=0; lx<A.extent(2); lx++) {
				Ksqr_factor = factor * (Ksqr+my_pow(Get_kx(lx)*kfactor[1],2));
				
				B(ly,lz,lx) += Ksqr_factor*A(ly,lz,lx);
			}
		}
	}
}


//********************************	End of four_tr.cc *****************************************


