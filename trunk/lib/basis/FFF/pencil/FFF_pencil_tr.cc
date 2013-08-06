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


#include "FFF_pencil.h"
#include "FFF_pencil_inline.h"

//*********************************************************************************************

// Ar: yxz, A: xyz
void FFF_PENCIL::Forward_transform(Array<DP,3> Ar, Array<complx,3> A)
{
	spectralTransform.Forward_transform_FFF_PENCIL(Ar,A);
}


//*********************************************************************************************


void FFF_PENCIL::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar)
{

    spectralTransform.Inverse_transform_FFF_PENCIL(A, Ar);
    
}

//*********************************************************************************************

void FFF_PENCIL::Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kx;
	
	for (int lx=0; lx<Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		
		B(Range::all(),Range::all(),lx) = complex<DP>(0, Kx)* (A(Range::all(),Range::all(),lx)); 	
	}
}


void FFF_PENCIL::Add_Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kx;
	
	for (int lx=0; lx<Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		
		B(Range::all(),Range::all(),lx) += complex<DP>(0, Kx)* (A(Range::all(),Range::all(),lx));
	}
}

void  FFF_PENCIL::Xderiv(Array<DP,3> A, Array<DP,3> B)
{
	cerr << "This is not defined for this basis. "<<endl;
}

//*********************************************************************************************


void FFF_PENCIL::Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Ky;
    
	for (int ly=0; ly<local_Ny; ly++) {
		Ky = Get_ky(ly)*kfactor[2];
		
		B(ly,Range::all(),Range::all()) = complx(0, Ky)* (A(ly,Range::all(),Range::all()));
	}
}


void FFF_PENCIL::Add_Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Ky;
    
	for (int ly=0; ly<local_Ny; ly++) {
		Ky = Get_ky(ly)*kfactor[2];
		
		B(ly,Range::all(),Range::all()) += complx(0, Ky)* (A(ly,Range::all(),Range::all()));
	}
}

//*********************************************************************************************


void FFF_PENCIL::Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
    
    for (int lz=0; lz<local_Nz; lz++) {
		Kz = Get_kz(lz)*kfactor[3];
		B(Range::all(),lz,Range::all()) = complx(0, Kz)*(A(Range::all(),lz,Range::all()));  	
	}
}

void FFF_PENCIL::Add_Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
    
    for (int lz=0; lz<local_Nz; lz++) {
		Kz = Get_kz(lz)*kfactor[3];
		B(Range::all(),lz,Range::all()) += complx(0, Kz)*(A(Range::all(),lz,Range::all()));
	}
}


/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/

void FFF_PENCIL::Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
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

void FFF_PENCIL::Subtract_Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
{
	
	DP Ksqr, Ksqr_factor;
	
	for (int ly=0; ly<A.extent(0); ly++) {
		Ksqr = my_pow(Get_ky(ly)*kfactor[2],2);
		
        for (int lz=0; lz<A.extent(1); lz++) {
			Ksqr +=  my_pow(Get_lz(lz)*kfactor[3],2);
			
            for (int lx=0; lx<A.extent(2); lx++) {
				Ksqr_factor = factor * (Ksqr+my_pow(Get_kx(lx)*kfactor[1],2));
				
				B(ly,lz,lx) += Ksqr_factor*A(ly,lz,lx);
			}
		}
	}
}

//********************************	End of four_tr.cc *****************************************


