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
 

#include "SSF_slab.h"
#include "SSF_slab_inline.h"

/**********************************************************************************************
comment
 
***********************************************************************************************/

void SSF_SLAB::Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz)
{}

void SSF_SLAB::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
{
	DP Kx = kx*kfactor[1];
	DP Ky = ky*kfactor[2];
	DP Kz = kz*kfactor[3];
	
	complx dvxdx, dvydy;
	int kysign;
	
    global.program.sincostr_switch = sincostr_switch_Vx;
	if (global.program.sincostr_switch[0] == 'S')
		dvxdx = complx(Kx,0)*Vx;
	else if (global.program.sincostr_switch[0] == 'C')
		dvxdx = -complx(Kx,0)*Vx;
	
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
		if (ky != 0) {// 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/complx(0,-Ky*kysign); 
		}
		else { // k = (kx,0,0); input fields are (Vy, Vz) field is purely real
			if ( (abs(imag(Vx)) < MYEPS) || (abs(imag(Vx)) < MYEPS))
				cout << "MYERROR: SSF_SLAB::Last_component(); For (kx,0,0), the modes are purely real; Setting imag part to zero. " << endl;
			Vz = complx(real(Vy), 0);
			Vy = complx(real(Vx), 0);
			Vx = complx(0,0);
		}
	}
}


/*****************************************************************************************
 
Dealias
 step 1: A(:,:,Range(Nx/3+1,toEnd)) = 0;
 
 step 2: A(:,Range(Nz/3+1,Nz/2),:) = 0;
 
 step 3: A(Range(Ny/3+1,toEnd),Range(0,Nz/3),:) = 0;
 
 *****************************************************************************************/

void SSF_SLAB::Dealias(Array<complx,3> A)
{
	
	Array<int,1> Ax_filter(Nx);
	
	Ax_filter = 0;
	
	Ax_filter(Range(2*Nx/3+1,toEnd)) = 1;
	
	int first_x = first(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
	int last_x = last(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
	
	A(Range(2*Ny/3+1,toEnd), Range(Nz/3,toEnd), Range(first_x,last_x)) = 0.0;
}

// Data resides till outer_radius in k-space
bool SSF_SLAB::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > 2*Nx/3) || (ky_max > 2*Ny/3) || (kz_max > Nz/3))
		return true;
	else
		return false;
}
/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 A is real aleardy.. do nothing.
 
 ***********************************************************************************************/

void SSF_SLAB::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
    return;  // Do nothing
}

void SSF_SLAB::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{
    return;  // Do nothing
}


void SSF_SLAB::Test_reality_condition_in_Array(Array<complx,3> A)
{
     return;  // Do nothing
}



//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SSF_SLAB::Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
    // lx = 0 reside in master node
    
    global.program.sincostr_switch = sincostr_switch_Vx;
    
    if (global.program.sincostr_switch == "SCF") {
        if (master)
            Ax(Range::all(),Range::all(),0) = 0.0;
        
        Ay(0,Range::all(),Range::all()) = 0.0;
    }
    
    else if (global.program.sincostr_switch == "CSF") {
        Ax(0,Range::all(),Range::all()) = 0.0;
        
        if (master)
            Ay(Range::all(),Range::all(),0) = 0.0;
        
        if (master)
            Az(0,Range::all(),0) = 0.0;
    }
    
    else if (global.program.sincostr_switch == "SSF") {
        if (master)
            Ax(0,Range::all(),0) = 0.0;
        
        Az(0,Range::all(),Range::all()) = 0.0;
    }
    
    else if (global.program.sincostr_switch == "CCF") {
        if (master)
            Ay(0,Range::all(),0) = 0.0;
        
        if (master)
            Az(Range::all(),Range::all(),0) = 0.0;
    }
    
}
    
/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SSF_SLAB::Zero_modes(Array<complx,3> F)
{
    // lx = 0 reside in master node
    global.program.sincostr_switch = sincostr_switch_F;
    
	if (global.program.sincostr_switch == "SCF") {
        if (master)
            F(Range::all(),Range::all(),0) = 0.0;
    }
    
    else if (global.program.sincostr_switch == "CSF")
        F(0,Range::all(),Range::all()) = 0.0;
    
    else if (global.program.sincostr_switch == "SSF") {
        if (master)
            F(0,Range::all(),0) = 0.0;
    }
}


//*********************************  End of scft_basic.cc *************************************


