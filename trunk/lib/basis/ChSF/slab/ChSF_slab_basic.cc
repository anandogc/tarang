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
 

#include "ChSF_slab.h"
#include "ChSF_slab_inline.h"



/**********************************************************************************************
comment
 
 WORK ON IT
 
***********************************************************************************************/

void ChSF_SLAB::Last_component(int kx, int ky, int kz, DP& dvxdx, DP& Vy, DP& Vz)
{
}

void ChSF_SLAB::Last_component(int kx, int ky, int kz,  complx &dvxdx, complx &Vy, complx &Vz)
{
	DP Kx = kx*kfactor[1];
	DP Ky = ky*kfactor[2];
	DP Kz = kz*kfactor[3];
	
	complx dvydy;
	int kysign;
	
	global.program.sincostr_switch = sincostr_switch_Vy;
	if (global.program.sincostr_switch[1] == 'S') {
		dvydy = complx(Ky,0)*Vy;
		kysign = 1;
	}
	else if (global.program.sincostr_switch[1] == 'C') {
		dvydy = -complx(Ky,0)*Vy;
		kysign = -1;
	}

	
	if (kz != 0)
		Vz = (dvxdx+dvydy)/complx(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) { // 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/complx(0,-Ky*kysign); 
		}
		else if (abs(dvxdx) > MYEPS) {	// k = (kx,0,0); input fields are (Vy, Vz)
			cerr << "MyERROR: Initcond with modes: u^(1)(kx,0,0) is nonzero!! " << endl;
		}
	}
}

/**********************************************************************************************
 
Dealias
 
 ***********************************************************************************************/

void ChSF_SLAB::Dealias(Array<complx,3> A)
{
	Array<int,1> Ay_filter(Ny);
	
	Ay_filter = 0;
	
	Ay_filter(Range(Ny/3+1,2*Ny/3-1)) = 1;
	
	int first_y = first(Ay_filter(Range(my_id*local_Ny,(my_id+1)*local_Ny-1)) == 1);
	int last_y = last(Ay_filter(Range(my_id*local_Ny,(my_id+1)*local_Ny-1)) == 1);
	
	A(Range(first_y,last_y), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0.0;
}

// Data resides till outer_radius in k-space
bool ChSF_SLAB::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
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

void ChSF_SLAB::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
	return;
}

void ChSF_SLAB::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{
	return;
}


void ChSF_SLAB::Test_reality_condition_in_Array(Array<complx,3> A)
{
	return;
}


//*********************************  End of scft_basic.cc *************************************


