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
 

#include "SFF_slab.h"
#include "SFF_slab_inline.h"



/**********************************************************************************************
comment
 
***********************************************************************************************/

void SFF_SLAB::Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz)
{}

void SFF_SLAB::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
{
	DP Kx = kx*kfactor[1];
	DP Ky = ky*kfactor[2];
	DP Kz = kz*kfactor[3];
	
	complx dvxdx;
	
	if (global.program.sincostr_switch[0] == 'S')
		dvxdx = Kx*Vx;
	
	else if (global.program.sincostr_switch[0] == 'C') 
		dvxdx = (-Kx)*Vx;
	
	
	if (kz != 0) 
		Vz = (dvxdx + complx(0,Ky)*Vy)/complx(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) {// 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/complx(0,-Ky); 
		}
		else { // k = (kx,0,0); input fields are (Vy, Vz) field is purely real
			if ( (abs(imag(Vx)) > MYEPS) || (abs(imag(Vx)) > MYEPS))
				cout << "MYERROR: SFF_SLAB::Last_component(); For (kx,0,0), the modes are purely real; Setting imag part to zero. " << endl;
			Vz = complx(real(Vy), 0);
			Vy = complx(real(Vx), 0);
			Vx = complx(0,0);
		}
	}
	
}

/**********************************************************************************************
 
Dealias
 
 A(Range(2*Ny/3+1,toEnd), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0;

 
 ***********************************************************************************************/

void SFF_SLAB::Dealias(Array<complx,3> A)
{
	Array<int,1> Ax_filter(Nx);
	
	Ax_filter = 0;
	
	Ax_filter(Range(2*Nx/3+1,toEnd)) = 1;
	
	int first_x = first(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
	int last_x = last(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
	
	A(Range(2*Ny/3+1,toEnd), Range(Nz/3,toEnd), Range(first_x,last_x)) = 0.0;
	
}

// Data resides till outer_radius in k-space
bool SFF_SLAB::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
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
 
// ky=Ny/2: set A(kx,ky,Nz/2)=0; quite efficient without losing many data points
 
 For 2d: f(m, 0, 0) = conj(f(m, 0, 0)) => set Im() = 0.
 
 ***********************************************************************************************/

void SFF_SLAB::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
	
	int array_index_minus_ky;
	
    // For a given (., minusky), locate (.,ky) and then subst.
    // A(minuskx, minusky,0) = conj(A(kx,ky,0))
	for (int ly=Ny/2+1; ly<Ny; ly++)
		for (int lx=0; lx<local_Nx; lx++) {
            array_index_minus_ky = -Get_ky(ly);        // minusky = Get_ky(ly)
            
            A(ly,0,lx) = conj( A(array_index_minus_ky,0,lx) );
        }
    
	// for (ky=0; kz=0) line
	imag(A(0,0,Range::all())) = 0.0;
    
    // for kz=Nz/2
	A(Range::all(),Nz/2,Range::all()) = 0.0;
}


void SFF_SLAB::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{

	A(Range::all(),Nz/2,Range::all()) = 0.0;
}


void SFF_SLAB::Test_reality_condition_in_Array(Array<complx,3> A)
{
	int array_index_ky;
	
	for (int ly=Ny/2+1; ly<Ny; ly++) {
		for (int lx=0; lx<local_Nx; lx++) {
            array_index_ky = -Get_ky(ly);        // minusky = Get_ky(ly)
            
            if (abs(A(ly,0,lx)-conj(A(array_index_ky,0,lx))) > MYEPS2)
                cout << "Reality condition voilated for (kx,ky,kz)=(" << lx <<  "," << ly << "," << 0 << ")" << endl; 

        }
    }
	
	// for (ky=0; kz=0) line
	if (sum(abs(imag(A(0,0,Range::all())))) > MYEPS2)
		cout << "Reality condition voilated for (kx,ky,kz)=(:" <<  "," << 0 << "," << 0 << ")" << endl;
    
    // for kz=Nz/2 plane
	for (int ly=0; ly<Ny; ly++)
		for (int lx=0; lx<local_Nx; lx++)
            if (abs(A(ly,Nz/2,lx)) > MYEPS)
                cout << "Reality condition voilated for (kx,ky,kz)=(" << Get_kx(lx) <<  "," << ly << "," << Nz/2 << ")" << endl;
}



//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SFF_SLAB::Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
    // lx = 0 reside in master node
    global.program.sincostr_switch = sincostr_switch_Vx;
    
	if (master) {
        if (global.program.sincostr_switch == "SFF")
            Ax(Range::all(),Range::all(),0) = 0.0;
        
        else if (global.program.sincostr_switch == "CFF") {
            Ay(Range::all(),Range::all(),0) = 0.0;
            Az(Range::all(),Range::all(),0) = 0.0;
        }
    } 
}

/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SFF_SLAB::Zero_modes(Array<complx,3> F)
{
    // lx = 0 reside in master node
    global.program.sincostr_switch = sincostr_switch_F;
    
    if ((master) && (global.program.sincostr_switch == "SFF"))
        F(Range::all(),Range::all(),0) = 0.0;
}

//*********************************  End of scft_basic.cc *************************************


