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

#include "SSS_pencil.h"
#include "SSS_pencil_inline.h"

//*********************************************************************************************

// Ar: yxz, A: xyz
void SSS_PENCIL::Forward_transform(Array<DP,3> Ar, Array<complx,3> A)
{
	spectralTransform.Forward_transform(global.program.sincostr_switch,Ar,A);
}


//*********************************************************************************************


void SSS_PENCIL::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar)
{
	
	global.temp_array.X_transform = A;	
    spectralTransform.Inverse_transform(global.program.sincostr_switch,global.temp_array.X_transform, Ar);
    
}



/**********************************************************************************************

	Derivative along x;  Bk_COS[i]=pi*i1*Ak_SIN[i]; Bk_SIN[]=-pi*i1*Ak_COS[i]

***********************************************************************************************/

void  SSS_PENCIL::Xderiv(Array<complx,3> A, Array<complx,3> B)
{	
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Kx;


	for (int lx = 0; lx < local_Nx; lx++) 	{
		Kx = Get_kx(lx)*kfactor[1];
		
        if (global.program.sincostr_switch[0] == 'S')
            Br(lx,Range::all(),Range::all()) = Kx*Ar(lx,Range::all(),Range::all());
		
        else if (global.program.sincostr_switch[0] == 'C')
            Br(lx,Range::all(),Range::all()) = (-Kx)*Ar(lx,Range::all(),Range::all());
        
	}
}

void  SSS_PENCIL::Add_Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Kx;
	
	for (int lx = 0; lx < local_Nx; lx++) 	{
		Kx = Get_kx(lx)*kfactor[1];
		
        if (global.program.sincostr_switch[0] == 'S')
            Br(lx,Range::all(),Range::all()) += Kx*Ar(lx,Range::all(),Range::all());
		
        else if (global.program.sincostr_switch[0] == 'C')
            Br(lx,Range::all(),Range::all()) += (-Kx)*Ar(lx,Range::all(),Range::all());
        
	}
}

void  SSS_PENCIL::Xderiv(Array<DP,3> A, Array<DP,3> B)
{
	cerr << "SSS_PENCIL::Xderiv(Array<DP,3> A, Array<DP,3> B) is not defined for this basis. "<<endl;
}
/**********************************************************************************************

		Derivative along y;  B(k) = i*ky*A(k)

***********************************************************************************************/

// Note: In the first half- ky=i2;
// In the second half- i2=0:Ny/-1; fftw-index=(Ny/2 +1+i2); FT-index=fftw-index-N=(i2+1-Ny/2)
void SSS_PENCIL::Yderiv(Array<complx,3> A, Array<complx,3> B)
{
    Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Ky;
	
	for (int ly=0; ly<local_Ny; ly++)  	{
		Ky = Get_ky(ly)*kfactor[2];
        
        if (global.program.sincostr_switch[1] == 'S')
			Br(Range::all(),ly,Range::all()) = Ky * Ar(Range::all(),ly,Range::all());
		
		else if(global.program.sincostr_switch[1] == 'C')
			Br(Range::all(),ly,Range::all()) = (-Ky) * Ar(Range::all(),ly,Range::all());
	}
}

void SSS_PENCIL::Add_Yderiv(Array<complx,3> A, Array<complx,3> B)
{
    Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Ky;
	
	for (int ly=0; ly<local_Ny; ly++)  	{
		Ky = Get_ky(ly)*kfactor[2];
        
        if (global.program.sincostr_switch[1] == 'S')
			Br(Range::all(),ly,Range::all()) += Ky * Ar(Range::all(),ly,Range::all());
		
		else if(global.program.sincostr_switch[1] == 'C')
			Br(Range::all(),ly,Range::all()) += (-Ky) * Ar(Range::all(),ly,Range::all());
	}
}

/**********************************************************************************************

		Derivative along z; Derivative along z

***********************************************************************************************/


void SSS_PENCIL::Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
	
    Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	for (int lz=0; lz<local_Nz; lz++) {
		Kz = lz*kfactor[3];
		
		if (global.program.sincostr_switch[2] == 'S')
			Br(Range::all(),Range::all(),lz) = Kz*Ar(Range::all(),Range::all(),lz);
		
		else if (global.program.sincostr_switch[2] == 'C')
			Br(Range::all(),Range::all(),lz) = (-Kz)*Ar(Range::all(),Range::all(),lz);
		
	}
}



void SSS_PENCIL::Add_Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	DP Kz;
	
    Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	for (int lz=0; lz<local_Nz; lz++) {
		Kz = lz*kfactor[3];
		
		if (global.program.sincostr_switch[2] == 'S')
			Br(Range::all(),Range::all(),lz) += Kz*Ar(Range::all(),Range::all(),lz);
		
		else if (global.program.sincostr_switch[2] == 'C')
			Br(Range::all(),Range::all(),lz) += (-Kz)*Ar(Range::all(),Range::all(),lz);
		
	}
}


/**********************************************************************************************
 
 B =  factor*Laplacian(A)
 
 ***********************************************************************************************/

void SSS_PENCIL::Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
{
	
	Array<DP,3> Ar(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Ksqr;
	
	for (int lx=0; lx<local_Nx; lx++) {
		Ksqr = my_pow(Get_kx(lx)*kfactor[1],2);
	
		for (int ly=0; ly<local_Ny; ly++) {
			Ksqr += my_pow(Get_ky(ly)*kfactor[2],2);
			
	        for (int lz=0; lz<local_Nz; lz++) {
				Ksqr += my_pow(lz*kfactor[3],2);
				Br(lx,ly,lz) = (-factor*Ksqr)*Ar(lx,ly,lz);
			}
		}
	}
}


/**********************************************************************************************
 
 B = B - factor*Laplacian(A) = B + factor*K^2 A
 
 ***********************************************************************************************/

void SSS_PENCIL::Subtract_Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B)
{
	
	Array<DP,3> Ar(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Ksqr, Ksqr_factor;
	
	for (int lx=0; lx<local_Nx; lx++) {
		Ksqr = my_pow(Get_kx(lx)*kfactor[1],2);
	
		for (int ly=0; ly<local_Ny; ly++) {
			Ksqr += my_pow(Get_ky(ly)*kfactor[2],2);
			
	        for (int lz=0; lz<local_Nz; lz++) {
				Ksqr_factor = factor*(Ksqr+my_pow(lz*kfactor[3],2));

				Br(lx,ly,lz) += Ksqr_factor*Ar(lx,ly,lz);
			}
		}
	}
}

//******************************** End of scft_tr.cc  *****************************************






