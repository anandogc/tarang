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

#include "SSF_slab.h"
#include "SSF_slab_inline.h"


/**********************************************************************************************

Forward SSF (A)  = Ar
FT along perp dirn and SIN transform along x dirn of A

***********************************************************************************************/

void SSF_SLAB::Forward_transform(Array<Real,3> Ar, Array<Complex,3> A)
{
    spectralTransform.Forward_transform(global.program.sincostr_switch, Ar, A);
}


/**********************************************************************************************
 
 Inverse SFT (A)  = Ar
 IFT along perp dirn and SIN transform along x dirn of A
 
 ***********************************************************************************************/


void SSF_SLAB::Inverse_transform(Array<Complex,3> A, Array<Real,3> Ar)
{
	global.temp_array.X_transform = A;
    spectralTransform.Inverse_transform(global.program.sincostr_switch, global.temp_array.X_transform, Ar);
}




/**********************************************************************************************

	Derivative along x;  Bk_COS[i]=pi*i1*Ak_SIN[i]; Bk_SIN[]=-pi*i1*Ak_COS[i]

***********************************************************************************************/

void  SSF_SLAB::Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{	
	Real Kx;
	
	for (int lx = 0; lx < local_Nx; lx++) 	{
		Kx = Get_kx(lx)*kfactor[1];
        
		if (global.program.sincostr_switch[0] == 'S')
            B(lx,Range::all(),Range::all()) = Kx * A(lx,Range::all(),Range::all());
		
        else if (global.program.sincostr_switch[0] == 'C')
            B(lx,Range::all(),Range::all()) = (-Kx) * A(lx,Range::all(),Range::all());
	}
}

void  SSF_SLAB::Add_Xderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kx;
	
	for (int lx = 0; lx < local_Nx; lx++) 	{
		Kx = Get_kx(lx)*kfactor[1];
        
		if (global.program.sincostr_switch[0] == 'S')
            B(lx,Range::all(),Range::all()) += Kx * A(lx,Range::all(),Range::all());
		
        else if (global.program.sincostr_switch[0] == 'C')
            B(lx,Range::all(),Range::all()) += (-Kx) * A(lx,Range::all(),Range::all());
	}
}


void  SSF_SLAB::Xderiv(Array<Real,3> A, Array<Real,3> B)
{
	cerr << "This is not defined for this basis. "<<endl;
}

/**********************************************************************************************

		Derivative along y;  B(k) = i*ky*A(k)

***********************************************************************************************/

// Note: In the first half- ky=i2;
// In the second half- i2=0:Ny/-1; fftw-index=(Ny/2 +1+i2); FT-index=fftw-index-N=(i2+1-Ny/2)
void SSF_SLAB::Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
	
	for (int ly=0; ly<Ny; ly++) {
		Ky = Get_ky(ly)*kfactor[2];
				
		if (global.program.sincostr_switch[1] == 'S')
			B(Range::all(),ly,Range::all()) = Ky*(A(Range::all(),ly,Range::all()));
		
		else if (global.program.sincostr_switch[1] == 'C')
			B(Range::all(),ly,Range::all()) = (-Ky)*(A(Range::all(),ly,Range::all()));
	}
}

void SSF_SLAB::Add_Yderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Ky;
	
	for (int ly=0; ly<Ny; ly++) {
		Ky = Get_ky(ly)*kfactor[2];
		
		if (global.program.sincostr_switch[1] == 'S')
			B(Range::all(),ly,Range::all()) += Ky*A(Range::all(),ly,Range::all());
		
		else if (global.program.sincostr_switch[1] == 'C')
			B(Range::all(),ly,Range::all()) += (-Ky)*A(Range::all(),ly,Range::all());
	}
}

/**********************************************************************************************

		Derivative along z; Derivative along z

***********************************************************************************************/


void SSF_SLAB::Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kz;
	
	for (int lz=0; lz<=Nz/2; lz++) {
		Kz = lz*kfactor[3];
		
		B(Range::all(),Range::all(),lz) = Complex(0, Kz)*A(Range::all(),Range::all(),lz); 	
	}   
}



void SSF_SLAB::Add_Zderiv(Array<Complex,3> A, Array<Complex,3> B)
{
	Real Kz;
	
	for (int lz=0; lz<=Nz/2; lz++) {
		Kz = lz*kfactor[3];
		B(Range::all(),Range::all(),lz) += Complex(0, Kz)*A(Range::all(),Range::all(),lz);
	}
}

/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/

void SSF_SLAB::Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	
	Real Ksqr;
	
    for (int lx=0; lx<A.extent(0); lx++) {
		Ksqr =  my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<A.extent(1); ly++) {
			Ksqr += my_pow(Get_ky(ly)*kfactor[2],2);
			
	        for (int lz=0; lz<A.extent(2); lz++) {
				Ksqr += my_pow(Get_kz(lz)*kfactor[3],2);	
				
				B(lx,ly,lz) = (-factor*Ksqr)*A(lx,ly,lz);
			}
		}
	}
}


/**********************************************************************************************
 
 B = B - factor*Laplacian(A) = B + factor*K^2 A
 
 ***********************************************************************************************/

void SSF_SLAB::Subtract_Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B)
{
	
	Real Ksqr, Ksqr_factor;
			
    for (int lx=0; lx<A.extent(0); lx++) {
		Ksqr =  my_pow(Get_kx(lx)*kfactor[1],2);
	
		for (int ly=0; ly<A.extent(1); ly++) {
			Ksqr +=  my_pow(Get_ky(ly)*kfactor[2],2);
			
	        for (int lz=0; lz<A.extent(2); lz++) {
				Ksqr_factor = factor * (Ksqr+my_pow(Get_kz(lz)*kfactor[3],2));

				B(lx,ly,lz) += Ksqr_factor*A(lx,ly,lz);
			}
		}
	}
}


//******************************** End of SSF_slab_tr.cc  *****************************************
