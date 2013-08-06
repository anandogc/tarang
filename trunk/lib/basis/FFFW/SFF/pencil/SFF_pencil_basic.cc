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
 

#include "SFF_pencil.h"
#include "SFF_pencil_inline.h"


 
/**********************************************************************************************
 
 COmment
 
***********************************************************************************************/
void SFF_PENCIL::Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz)
{}

void SFF_PENCIL::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
{
	DP Kx = kx*kfactor[1];
	DP Ky = ky*kfactor[2];
	DP Kz = kz*kfactor[3];
	
	complx dvxdx;
	
	if (global.program.sincostr_switch[0] == 'S')
		dvxdx = complx(Kx,0)*Vx;
	
	else if (global.program.sincostr_switch[0] == 'C') 
		dvxdx = complx(-Kx,0)*Vx;
	
	if (kz != 0) 
		Vz = (dvxdx + complx(0,Ky)*Vy)/complx(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) {// 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/complx(0,-Ky); 
		}
		else { // k = (kx,0,0)
			if ( (abs(imag(Vx)) > MYEPS) || (abs(imag(Vx)) > MYEPS))
				cout << "MYERROR: SCFT_SLAB::Last_component(); For (kx,0,0), the modes are purely real; Setting imag part to zero. " << endl;
			Vz = complx(real(Vy), 0);
			Vy = complx(real(Vx), 0);
			Vx = complx(0,0);
		}
	}
	
}


/*******************************************************************************
 
 Dealias

	A(Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0;
 
 *******************************************************************************/

void SFF_PENCIL::Dealias(Array<complx,3> A)
{
	Assign_sub_array(Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd), A, complx(0,0));
}


// Data resides till outer_radius in k-space
bool SFF_PENCIL::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > 2*Nx/3) || (ky_max > Ny/3) || (kz_max > Nz/3))
		return true;

	return false;
}


/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
 No need for data transfer; local..
 
f(m, -ky, 0) = conj(f(m, ky, 0))- do for kz=N[3]/2
 
 ***********************************************************************************************/

void SFF_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
	int array_index_minus_ky;
    
    // For a given (., minusky), locate (.,ky) and then subst.
    // A(., minusky,0) = conj(A(.,ky,0))
	if (my_z_pcoord == 0) {
        for (int lx=0; lx<local_Nx; lx++)
            for (int ly=Ny/2+1; ly<Ny; ly++) {
                array_index_minus_ky = -Get_ky(ly);        // minusky = Get_ky(ly)
                
                A(ly,0,lx) = conj(A(array_index_minus_ky,0,lx));
            }
            // for (ky=0,kz=0) line
            imag(A(0,0,Range::all())) = 0.0;
	}
	
	if (my_z_pcoord == num_z_procs-1)
         A(local_Nz-1,Range::all(),Range::all()) = 0.0;
}

void SFF_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{
    if (my_z_pcoord == num_z_procs-1)
        A(Range::all(),local_Nz,Range::all()) = 0.0;
}



void SFF_PENCIL::Test_reality_condition_in_Array(Array<complx,3> A)
{
    int array_index_minus_ky;
	
	if (my_z_pcoord == 0) {
        for (int lx=0; lx<local_Nx; lx++) {
           for (int ly=Ny/2+1; ly<Ny; ly++) {
                array_index_minus_ky = -Get_ky(ly);
               
                if (abs(A(ly,0,lx)-conj(A(array_index_minus_ky,0,lx))) > MYEPS2)
                    cout << "Reality condition voilated for (kx,ky,kz)=(" << lx <<  "," << ly << "," << 0 << ")" << endl; 
            }
            // for (ky=0,kz=0) line
            if (abs(imag(A(0,0,lx))) > MYEPS2)
                cout << "Reality condition voilated for (kx,ky,kz)=(" << lx <<  "," << 0 << "," << 0 << ")" << endl;
        }
	}
	
    // for kz=Nz/2
    int last_index=local_Nz-1;
	if (my_z_pcoord == num_z_procs-1) {
        for (int lx=0; lx<local_Nx; lx++)
			for (int ly=0; ly<Ny; ly++)
                if (abs(A(last_index,lx,ly)) > MYEPS)
                    cout << "Reality condition voilated for (kx,ky,kz)=(" << Get_kx(lx) <<  "," << Get_ky(ly) << "," << Nz/2 << ")" << endl;
    }

}


//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SFF_PENCIL::Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
    // lx = 0 reside in master node
    
    Range zero(0,0);
    global.program.sincostr_switch = sincostr_switch_Vx;

	if (my_vert_pcoord == 0) {
        if (global.program.sincostr_switch == "SFF")
        	Assign_sub_array(Range::all(),Range::all(),zero,Ax,complx(0,0));
        
        else if (global.program.sincostr_switch == "CFF") {
        	Assign_sub_array(Range::all(),Range::all(),zero,Ay,complx(0,0));
        	Assign_sub_array(Range::all(),Range::all(),zero,Az,complx(0,0));
        }
    }
}

/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SFF_PENCIL::Zero_modes(Array<complx,3> F)
{
    // lx = 0 reside in master node
    Range zero(0,0);
    global.program.sincostr_switch = sincostr_switch_F;
    if ((my_vert_pcoord == 0) && (global.program.sincostr_switch == "SFF"))
    	Assign_sub_array(Range::all(),Range::all(),zero,F,complx(0,0));
}

void SFF_PENCIL::Assign_sub_array(Range y_range, Range z_range, Range x_range, Array<complx,3> A, complx value)
{
	
	static Array<int,1> z_filter(Nz/2+1);
	static Array<int,1> x_filter(Nx);
	
	z_filter = 0;
	x_filter = 0;
	
	z_filter(z_range)=1;
	x_filter(x_range)=1;
	
	
	static Range z_apply, x_apply;
	
		
	z_apply = Range(first(z_filter(Range(my_z_pcoord*local_Nz,(my_z_pcoord+1)*local_Nz-1)) == 1),
					 last(z_filter(Range(my_z_pcoord*local_Nz,(my_z_pcoord+1)*local_Nz-1)) == 1));
	
	x_apply = Range(first(x_filter(Range(my_x_pcoord*local_Nx,(my_x_pcoord+1)*local_Nx-1)) == 1),
					 last(x_filter(Range(my_x_pcoord*local_Nx,(my_x_pcoord+1)*local_Nx-1)) == 1));

	
	
	if ( (z_apply(0)>=0) && (x_apply(0)>=0))
		A(y_range, z_apply, x_apply) = value;
}

//*********************************  End of scft_basic.cc *************************************


