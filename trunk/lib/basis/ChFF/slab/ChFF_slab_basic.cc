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

/*! \file scft_basic.cc
 *
 * @sa scft_basic.h
 *
 * @author  M. K. Verma
 * @version 4.0 
 * @date 30 August 2008
 * @bug  No known bug
 */
 

#include "ChFF_slab.h"
#include "ChFF_slab_inline.h"



/**********************************************************************************************
comment
 
 WORK ON IT
 
***********************************************************************************************/

void ChFF_SLAB::Last_component(int kx, int ky, int kz, Real& dvxdx, Real& Vy, Real& Vz)
{
}

void ChFF_SLAB::Last_component(int kx, int ky, int kz, Complex& dvxdx, Complex& Vy, Complex& Vz)
{
	Real Kx = kx*kfactor[1];
	Real Ky = ky*kfactor[2];
	Real Kz = kz*kfactor[3];
	
	if (kz != 0)
		Vz = (dvxdx + Complex(0,Ky)*Vy)/Complex(0,-Kz);
	// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) { // 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/Complex(0,-Ky);
		}
		else if (abs(dvxdx) > MYEPS) {	// k = (kx,0,0); input fields are (Vy, Vz)
			cerr << "MyERROR: Initcond with modes: u^(1)(kx,0,0) is nonzero!! " << endl;
		}
	}
}

/**********************************************************************************************
 
Dealias
 
 ***********************************************************************************************/

void ChFF_SLAB::Dealias(Array<Complex,3> A)
{
	Array<int,1> Ay_filter(Ny);
	
	Ay_filter = 0;
	
	Ay_filter(Range(Ny/3+1,2*Ny/3-1)) = 1;
	
	int first_y = first(Ay_filter(Range(my_id*local_Ny,(my_id+1)*local_Ny-1)) == 1);
	int last_y = last(Ay_filter(Range(my_id*local_Ny,(my_id+1)*local_Ny-1)) == 1);
	
	A(Range(first_y,last_y), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0.0;
}

// Data resides till outer_radius in k-space
bool ChFF_SLAB::Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > 2*Nx/3) || (ky_max > Ny/3) || (kz_max > Nz/3))
		return true;
	else
		return false;
}


/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
 
f(m, -ky, 0) = conj(f(m, ky, 0))- do for kz=N[3]/2
 
 ***********************************************************************************************/

void ChFF_SLAB::Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A)
{
	int array_index_minus_ky;
	
	// for kz=0
	ArrayOps::Get_XY_plane(A, global.temp_array.plane_xy, 0);
	
    // For a given (minusky, minuskx), locate (ky,ky) and then subst.
    // A(minusky, 0, minuskx) = conj(A(ky,0,kx))
	for (int ly=0; ly<local_Ny; ly++) {
		if (Get_ky(ly) < 0) {  // Do not apply for ky=Ny/2
			array_index_minus_ky =  Get_iy(-Get_ky(ly));  // minusky = -Get_ky(ly) > 0;
			
			A(ly, 0, Range::all()) = conj(global.temp_array.plane_xy(array_index_minus_ky,Range::all()));	
		}
	}
	
	// for kz=Nz/2;
	A(Range::all(),Nz/2,Range::all()) = 0.0;
}

void ChFF_SLAB::Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A)
{
	// for kz=Nz/2;
	A(Range::all(),Nz/2,Range::all()) = 0.0;
}


void ChFF_SLAB::Test_reality_condition_in_Array(Array<Complex,3> A)
{
	int array_index_minus_ky;
	
	// for kz=0
	ArrayOps::Get_XY_plane(A, global.temp_array.plane_xy, 0);
	
    // For a given (minusky, minuskx), locate (ky,ky) and then subst.
    // A(minusky, 0, minuskx) = conj(A(ky,0,kx))
	for (int ly=0; ly<local_Ny; ly++) {
		if (Get_ky(ly) < 0) {  // Do not apply for ky=Ny/2
			array_index_minus_ky =  Get_iy(-Get_ky(ly));  // minusky = -Get_ky(ly) > 0;
			
			for (int lx=0; lx<Nx; lx++) {
				if (abs(A(ly,0,lx)-conj(global.temp_array.plane_xy(array_index_minus_ky,lx))) > MYEPS2)
					cout << "Reality condition voilated for (kx,ky,kz)=(" << lx <<  "," << Get_ky(ly) << "," << 0 << ")" << endl;
			}
		}
	}
	
	// for kz=Nz/2;
	A(Range::all(),Nz/2,Range::all()) = 0.0;
}

//*********************************  End of scft_basic.cc *************************************


