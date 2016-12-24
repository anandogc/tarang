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

void SSF_PENCIL::Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz)
{}

void SSF_PENCIL::Last_component(int kx, int ky, int kz, Complex &Vx, Complex &Vy, Complex &Vz)
{
	Real Kx = kx*kfactor[1];
	Real Ky = ky*kfactor[2];
	Real Kz = kz*kfactor[3];
	
	Complex dvxdx, dvydy;
	int kysign;
	
	global.program.sincostr_switch = sincostr_switch_Vx;
	if (global.program.sincostr_switch[0] == 'S')
		dvxdx = Complex(Kx,0)*Vx;
	else if (global.program.sincostr_switch[0] == 'C')
		dvxdx = -Complex(Kx,0)*Vx;
	
	global.program.sincostr_switch = sincostr_switch_Vy;
	if (global.program.sincostr_switch[1] == 'S')
		dvydy = Complex(Ky,0)*Vy;
	else if (global.program.sincostr_switch[1] == 'C')
		dvydy = Complex(-Ky,0)*Vy;
	
	if (global.program.sincostr_switch[1] == 'S')
		kysign = 1;
	else if (global.program.sincostr_switch[1] == 'C')
		kysign = -1;

	
	if (kz != 0) 
		Vz = (dvxdx+dvydy)/Complex(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) {// 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/Complex(0,-Ky*kysign); 
		}
		else { // k = (kx,0,0); input fields are (Vy, Vz) field is purely real
			if ( (abs(imag(Vx)) < MYEPS) || (abs(imag(Vx)) < MYEPS))
				cout << "MYERROR: SSF_PENCIL::Last_component(); For (kx,0,0), the modes are purely real; Setting imag part to zero. " << endl;
			Vx = Complex(0,0);
			Vy = Complex(real(Vx), 0);
			Vz = Complex(real(Vy), 0);
		}
	}
	
}

/**********************************************************************************************
 
Dealias
	A(Range(2*Ny/3+1,toEnd), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0;
 
 ***********************************************************************************************/

void SSF_PENCIL::Dealias(Array<Complex,3> A)
{
	Assign_sub_array(Range(2*Nx/3+1,toEnd), Range(2*Ny/3+1,toEnd), Range(Nz/3+1,toEnd), A, Complex(0,0));
}

// Data resides till outer_radius in k-space
bool SSF_PENCIL::Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius)
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

void SSF_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}

void SSF_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}

void SSF_PENCIL::Test_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}


//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SSF_PENCIL::Zero_modes(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	// lx = 0 reside in master node
	
	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_Vx;
	
	if (global.program.sincostr_switch == "SCF") {
		Assign_sub_array(zero,Range::all(),Range::all(),Ax,Complex(0,0));
		Assign_sub_array(Range::all(),zero,Range::all(),Ay,Complex(0,0));
	}
	
	else if (global.program.sincostr_switch == "CSF") {
		if (Ny>1)
			Ax(Range::all(),0,Range::all()) = 0.0;

		Assign_sub_array(Range::all(),Range::all(),zero,Ax,Complex(0,0));


		Assign_sub_array(zero,Range::all(),Range::all(),Ay,Complex(0,0));
		Assign_sub_array(Range::all(),Range::all(),zero,Ay,Complex(0,0));
	}
	
	else if (global.program.sincostr_switch == "SSF") {
		Assign_sub_array(zero,Range::all(),zero,Ax,Complex(0,0));
		if (Ny>1)
			Ax(Range::all(),0,Range::all()) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CCF") {
		Assign_sub_array(Range::all(),Range::all(),zero,Ax,Complex(0,0));

		Assign_sub_array(zero,Range::all(),Range::all(),Ay,Complex(0,0));
		if (Ny>1)
			Ay(Range::all(),0,Range::all()) = 0.0;	
	}
	
}
	
/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SSF_PENCIL::Zero_modes(Array<Complex,3> F)
{
	// lx = 0 reside in master node
	
	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_F;
	
	if (global.program.sincostr_switch == "SCF")
		Assign_sub_array(zero,Range::all(),Range::all(),F,Complex(0,0));

	else if (global.program.sincostr_switch == "CSF")
		Assign_sub_array(zero,Range::all(),Range::all(),F,Complex(0,0));

	else if (global.program.sincostr_switch == "SSF")
		Assign_sub_array(zero,Range::all(),zero,F,Complex(0,0));
}

void SSF_PENCIL::Assign_sub_array(Range x_range, Range y_range, Range z_range, Array<Complex,3> A, Complex value)
{
	static Array<int,1> y_filter(Ny);
	static Array<int,1> z_filter(Nz/2+1);
	
	//Sanitize ranges. if last index is less than first index, modify the range (this happens for 2D)
	if (y_range.last() < y_range.first())
		y_range = Range(y_range.first(), y_range.first());

	if (z_range.last() < z_range.first())
		z_range = Range(z_range.first(), z_range.first());

	y_filter = 0;
	z_filter = 0;

	y_filter(y_range)=1;
	z_filter(z_range)=1;

	
	static Range y_apply, z_apply;
	
	
	y_apply = Range(first(y_filter(Range(my_y_pcoord*maxly,(my_y_pcoord+1)*maxly-1)) == 1),
					 last(y_filter(Range(my_y_pcoord*maxly,(my_y_pcoord+1)*maxly-1)) == 1));
	
	z_apply = Range(first(z_filter(Range(my_z_pcoord*maxlz,(my_z_pcoord+1)*maxlz-1)) == 1),
					 last(z_filter(Range(my_z_pcoord*maxlz,(my_z_pcoord+1)*maxlz-1)) == 1));
	
	
	if ( (y_apply(0)>=0) && (z_apply(0)>=0))
			A(x_range, y_apply, z_apply) = value;
}

int SSF_PENCIL::Read(Array<Complex,3> A, h5::Plan plan, string file_name, string dataset_name)
{
	if (Ny > 1) {
		if (my_z_pcoord == 0)  
			BasicIO::Read(global.temp_array.Xr_slab.data(), plan, file_name, dataset_name);

		fftk.To_pencil(global.temp_array.Xr_slab, global.temp_array.Xr);
		fftk.Transpose(global.temp_array.Xr, A);

	}
	if (Ny == 1) {
		BasicIO::Read(global.temp_array.Xr.data(), plan, file_name, dataset_name);
		fftk.Transpose(global.temp_array.Xr(Range::all(),0,Range::all()), A(Range::all(),0,Range::all()));
	}
	return 0;
}

int SSF_PENCIL::Read(Array<Real,3> Ar, h5::Plan plan, string file_name, string dataset_name)
{
	if (Ny > 1) {
		if (my_z_pcoord == 0)  
			BasicIO::Read(global.temp_array.Xr_slab.data(), plan, file_name, dataset_name);

		fftk.To_pencil(global.temp_array.Xr_slab, Ar);
	}
	if (Ny == 1) {
		BasicIO::Read(Ar.data(), plan, file_name, dataset_name);
	}
	return 0;
}


int SSF_PENCIL::Write(Array<Complex,3> A, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
{
	if (Ny > 1) {
		fftk.Transpose(A, global.temp_array.Xr);
		fftk.To_slab(global.temp_array.Xr, global.temp_array.Xr_slab);
		if (my_z_pcoord == 0)  
			BasicIO::Write(global.temp_array.Xr_slab.data(), plan, access_mode, folder_name, file_name, dataset_name);

	}
	if (Ny == 1) {
		fftk.Transpose(A(Range::all(),0,Range::all()), global.temp_array.Xr(Range::all(),0,Range::all()));
		BasicIO::Write(global.temp_array.Xr.data(), plan, access_mode, folder_name, file_name, dataset_name);

	}

	return 0;
}

int SSF_PENCIL::Write(Array<Real,3> Ar, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
{
	if (Ny > 1) {
		fftk.To_slab(Ar, global.temp_array.Xr_slab);
		if (my_z_pcoord == 0)  
			BasicIO::Write(global.temp_array.Xr_slab.data(), plan, access_mode, folder_name, file_name, dataset_name);

	}
	if (Ny == 1) {
		BasicIO::Write(Ar.data(), plan, access_mode, folder_name, file_name, dataset_name);
	}
	return 0;
}

//*********************************  End of scft_basic.cc *************************************


