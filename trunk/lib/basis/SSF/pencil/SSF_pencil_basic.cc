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
 

#include "SSF_pencil.h"
#include "SSF_pencil_inline.h"

/**********************************************************************************************
comment
 
***********************************************************************************************/

void SSF_PENCIL::Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz)
{}

void SSF_PENCIL::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
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
	if (global.program.sincostr_switch[1] == 'S')
		dvydy = complx(Ky,0)*Vy;
	else if (global.program.sincostr_switch[1] == 'C')
		dvydy = complx(-Ky,0)*Vy;
	
	if (global.program.sincostr_switch[1] == 'S')
		kysign = 1;
	else if (global.program.sincostr_switch[1] == 'C')
		kysign = -1;

	
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
				cout << "MYERROR: SSF_PENCIL::Last_component(); For (kx,0,0), the modes are purely real; Setting imag part to zero. " << endl;
			Vx = complx(0,0);
			Vy = complx(real(Vx), 0);
			Vz = complx(real(Vy), 0);
		}
	}
	
}

/**********************************************************************************************
 
Dealias
	A(Range(2*Ny/3+1,toEnd), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0;
 
 ***********************************************************************************************/

void SSF_PENCIL::Dealias(Array<complx,3> A)
{
	Assign_sub_array(Range(2*Nx/3+1,toEnd), Range(2*Ny/3+1,toEnd), Range(Nz/3+1,toEnd), A, complx(0,0));
}

// Data resides till outer_radius in k-space
bool SSF_PENCIL::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > 2*Nx/3) || (ky_max > 2*Ny/3) || (kz_max > Nz/3))
		return true;
	
	return false;
}
/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 A is real aleardy.. do nothing.
 
 ***********************************************************************************************/

void SSF_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
	return; // Do nothing
}

void SSF_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{
	return; // Do nothing
}

void SSF_PENCIL::Test_reality_condition_in_Array(Array<complx,3> A)
{
	return; // Do nothing
}


//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SSF_PENCIL::Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
	// lx = 0 reside in master node
	
	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_Vx;
	
	if (global.program.sincostr_switch == "SCF") {
		Assign_sub_array(zero,Range::all(),Range::all(),Ax,complx(0,0));
		Assign_sub_array(Range::all(),zero,Range::all(),Ay,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "CSF") {
		if (Ny>1)
			Ax(Range::all(),0,Range::all()) = 0.0;

		Assign_sub_array(Range::all(),Range::all(),zero,Ax,complx(0,0));


		Assign_sub_array(zero,Range::all(),Range::all(),Ay,complx(0,0));
		Assign_sub_array(Range::all(),Range::all(),zero,Ay,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "SSF") {
		Assign_sub_array(zero,Range::all(),zero,Ax,complx(0,0));
		if (Ny>1)
			Ax(Range::all(),0,Range::all()) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CCF") {
		Assign_sub_array(Range::all(),Range::all(),zero,Ax,complx(0,0));

		Assign_sub_array(zero,Range::all(),Range::all(),Ay,complx(0,0));
		if (Ny>1)
			Ay(Range::all(),0,Range::all()) = 0.0;	
	}
	
}
	
/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SSF_PENCIL::Zero_modes(Array<complx,3> F)
{
	// lx = 0 reside in master node
	
	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_F;
	
	if (global.program.sincostr_switch == "SCF")
		Assign_sub_array(zero,Range::all(),Range::all(),F,complx(0,0));

	else if (global.program.sincostr_switch == "CSF")
		Assign_sub_array(zero,Range::all(),Range::all(),F,complx(0,0));

	else if (global.program.sincostr_switch == "SSF")
		Assign_sub_array(zero,Range::all(),zero,F,complx(0,0));
}

void SSF_PENCIL::Assign_sub_array(Range x_range, Range y_range, Range z_range, Array<complx,3> A, complx value)
{
	static Array<int,1> x_filter(Nx);
	static Array<int,1> y_filter(Ny);
	
	x_filter = 0;
	y_filter = 0;
	
	x_filter(x_range)=1;
	y_filter(y_range)=1;
	
	
	static Range x_apply, y_apply;
	
	
	x_apply = Range(first(x_filter(Range(my_x_pcoord*local_Nx,(my_x_pcoord+1)*local_Nx-1)) == 1),
	                 last(x_filter(Range(my_x_pcoord*local_Nx,(my_x_pcoord+1)*local_Nx-1)) == 1));
	
	y_apply = Range(first(y_filter(Range(my_y_pcoord*local_Ny,(my_y_pcoord+1)*local_Ny-1)) == 1),
	                 last(y_filter(Range(my_y_pcoord*local_Ny,(my_y_pcoord+1)*local_Ny-1)) == 1));
	
	
	if ( (x_apply(0)>=0) && (y_apply(0)>=0))
		A(x_apply, y_apply, z_range) = value;
}

int SSF_PENCIL::Read(Array<complx,3> A, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	int err = BasicIO::Read(global.temp_array.Xr.data(), plan, file_name, dataset_name);
	spectralTransform.Transpose(global.temp_array.Xr, A);
	return err;
}

int SSF_PENCIL::Read(Array<DP,3> Ar, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	return BasicIO::Read(Ar.data(), plan, file_name, dataset_name);
}


int SSF_PENCIL::Write(Array<complx,3> A, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	spectralTransform.Transpose(A, global.temp_array.Xr);
	return BasicIO::Write(global.temp_array.Xr.data(), plan, folder_name, file_name, dataset_name);
}

int SSF_PENCIL::Write(Array<DP,3> Ar, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	return BasicIO::Write(Ar.data(), plan, folder_name, file_name, dataset_name);  
}

//*********************************  End of scft_basic.cc *************************************


