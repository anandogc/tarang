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
 

#include "SSS_pencil.h"
#include "SSS_pencil_inline.h"

/**********************************************************************************************
 
 COmment
 
***********************************************************************************************/


void SSS_PENCIL::Print_large_Fourier_elements(Array<Complex,3> A, string array_name)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	for (int lx=0; lx<Ar.extent(0); lx++)
		for (int ly=0; ly<Ar.extent(1); ly++)
			for (int lz=0; lz<Ar.extent(2); lz++) {
				if (abs(Ar(lx,ly,lz)) > MYEPS2)
					cout << "my_id = " << my_id <<  " " << array_name <<"(" << Get_kx(lx) << "," << Get_ky(ly) << "," << Get_kz(lz) <<") = " << Ar(lx, ly, lz) << '\n';
			}
	
	cout << endl;
}


void SSS_PENCIL::Last_component(int kx, int ky, int kz, Complex &Vx, Complex &Vy, Complex &Vz)
{}

void SSS_PENCIL::Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz)
{
	Real Kx = kx*kfactor[1];
	Real Ky = ky*kfactor[2];
	Real Kz = kz*kfactor[3];
	
	Real dvxdx, dvydy;
	int kysign, kzsign;
	
	global.program.sincostr_switch = sincostr_switch_Vx;
	if (global.program.sincostr_switch[0] == 'S')
		dvxdx = Kx*Vx;
	else if (global.program.sincostr_switch[0] == 'C')
		dvxdx = -Kx*Vx;
	
	global.program.sincostr_switch = sincostr_switch_Vy;
	if (global.program.sincostr_switch[1] == 'S')
		dvydy = Ky*Vy;
	else if (global.program.sincostr_switch[1] == 'C')
		dvydy = -Ky*Vy;
	else if (global.program.sincostr_switch[1] == '0')
		dvydy = 0;
	
	if (global.program.sincostr_switch[1] == 'S')
		kysign = 1;
	else if (global.program.sincostr_switch[1] == 'C')
		kysign = -1;
	
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	if (global.program.sincostr_switch[2] == 'S')
		kzsign = 1;
	else if (global.program.sincostr_switch[2] == 'C')
		kzsign = -1;
	
	if (kz != 0)
		Vz = (dvxdx+dvydy)/(-Kz*kzsign);
	// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) {// 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/(-Ky*kysign);
		}
		else { // 2D: Input fields are Vy, Vz
			Vz = Vy;
			Vy = Vx;
			Vx = 0.0;
		}
	}
}

/**********************************************************************************************
 
Dealias
 
 A(Range(2*Ny/3+1,toEnd), Range(2*Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0;

 In complex array A: z-range = Nz/3

 
 ***********************************************************************************************/

void SSS_PENCIL::Dealias(Array<Complex,3> A)
{
	Array<Real,3> Ar(reinterpret_cast<Real*>(A.data()), A.shape()*shape(1,1,2), neverDeleteData);

	Assign_sub_array(Range(2*Nx/3+1,toEnd), Range(2*Ny/3+1,toEnd), Range(Nz/3,toEnd), Ar, 0);
}


// Data resides till outer_radius in k-space
bool SSS_PENCIL::Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > 2*Nx/3) || (ky_max > 2*Ny/3) || (kz_max > 2*Nz/3))
		return true;

	return false;
}


/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 A is real aleardy.. do nothing.
 
 ***********************************************************************************************/

void SSS_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}

void SSS_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}

void SSS_PENCIL::Test_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}


//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SSS_PENCIL::Zero_modes(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	Array<Real,3> Axr(reinterpret_cast<Real*>(Ax.data()), Ax.shape()*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ayr(reinterpret_cast<Real*>(Ay.data()), Ay.shape()*shape(1,1,2), neverDeleteData);
	Array<Real,3> Azr(reinterpret_cast<Real*>(Az.data()), Az.shape()*shape(1,1,2), neverDeleteData);

	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_Vx;
	
	if (global.program.sincostr_switch == "SCC") {
		Assign_sub_array(zero,Range::all(),Range::all(),Axr,0);
		Assign_sub_array(Range::all(),zero,Range::all(),Ayr,0);
		Assign_sub_array(Range::all(),Range::all(),zero,Azr,0);
	}
	
	else if (global.program.sincostr_switch == "CSS") {
		Assign_sub_array(Range::all(),zero,zero,Axr,0);
		Assign_sub_array(zero,zero,Range::all(),Ayr,0);
		Assign_sub_array(zero,zero,Range::all(),Azr,0);
	}
	
	else if (global.program.sincostr_switch == "CCS") {
		Assign_sub_array(Range::all(),Range::all(),zero,Axr,0);
		Assign_sub_array(zero,zero,zero,Ayr,0);
		Assign_sub_array(zero,Range::all(),Range::all(),Azr,0);
	}
	
	else if (global.program.sincostr_switch == "SSC") {
		Assign_sub_array(zero,zero,Range::all(),Axr,0);
		Assign_sub_array(Range::all(),zero,zero,Azr,0);
	}
	
	else if (global.program.sincostr_switch == "CSC") {
		Assign_sub_array(Range::all(),zero,Range::all(),Axr,0);
		Assign_sub_array(zero,Range::all(),Range::all(),Ayr,0);
		Assign_sub_array(zero,zero,zero,Azr,0);
	}
	
	else if (global.program.sincostr_switch == "SCS") {
		Assign_sub_array(zero,zero,Range::all(),Axr,0);
		Assign_sub_array(Range::all(),zero,zero,Ayr,0);
	}
	
	else if (global.program.sincostr_switch == "CCC") {
		Assign_sub_array(zero,zero,Range::all(),Ayr,0);
		Assign_sub_array(zero,Range::all(),zero,Azr,0);
	}
	
	else if (global.program.sincostr_switch == "SSS") {
		Assign_sub_array(zero,zero,zero,Axr,0);
		Assign_sub_array(Range::all(),Range::all(),zero,Ayr,0);
		Assign_sub_array(Range::all(),zero,Range::all(),Azr,0);
	}
}
	
/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SSS_PENCIL::Zero_modes(Array<Complex,3> F)
{
	Array<Real,3> Fr(reinterpret_cast<Real*>(F.data()), F.shape()*shape(1,1,2), neverDeleteData);

	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_F;
	
	if (global.program.sincostr_switch == "SCC")
		Assign_sub_array(zero,Range::all(),Range::all(),Fr,0);
	
	else if (global.program.sincostr_switch == "CSS")
		Assign_sub_array(Range::all(),zero,zero,Fr,0);
	
	else if (global.program.sincostr_switch == "CCS")
		Assign_sub_array(Range::all(),Range::all(),zero,Fr,0);

	else if (global.program.sincostr_switch == "SSC")
		Assign_sub_array(zero,zero,Range::all(),Fr,0);
	
	else if (global.program.sincostr_switch == "CSC")
		Assign_sub_array(Range::all(),zero,Range::all(),Fr,0);
	
	else if (global.program.sincostr_switch == "SCS")
		Assign_sub_array(zero,zero,Range::all(),Fr,0);
	
	else if (global.program.sincostr_switch == "CCC")
		;
	
	else if (global.program.sincostr_switch == "SSS")
		Assign_sub_array(zero,zero,zero,Fr,0);
}

void SSS_PENCIL::Assign_sub_array(Range x_range, Range y_range, Range z_range, Array<Real,3> Ar, Real value)
{

	static Array<int,1> y_filter(Ny);
	static Array<int,1> z_filter(Nz);

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
	
	z_apply = Range(first(z_filter(Range(my_z_pcoord*2*maxlz,(my_z_pcoord+1)*2*maxlz-1)) == 1),
					 last(z_filter(Range(my_z_pcoord*2*maxlz,(my_z_pcoord+1)*2*maxlz-1)) == 1));
	
	
	if ( (y_apply(0)>=0) && (z_apply(0)>=0))
		Ar(x_range, y_apply, z_apply) = value;
}


/**********************************************************************************************

			Replaces A(k) by A(k)*K^2.

***********************************************************************************************/


void SSS_PENCIL::Array_mult_ksqr(Array<Complex,3> A)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kxsqr;    // Kx^2
	Real Kxysqr;
	Real Ksqr;
	
	for (int lx=0; lx<maxlx; lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<maxly; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);	
			
			for (int lz=0; lz<2*maxlz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

				Ar(lx, ly, lz) *= Ksqr;
		  }
		}
	}
}


/**********************************************************************************************

	Replaces A(k) by A(k)/K^2 with K^2 = sum_i [ki*kfactor(i)]^2 ; A(0)=0

***********************************************************************************************/



void SSS_PENCIL::Array_divide_ksqr(Array<Complex,3> A)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kxsqr;    // Kx^2
	Real Kxysqr;
	Real Ksqr;
	
	for (int lx=0; lx<maxlx; lx++) {
		Kxsqr= my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<maxly; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);    
			
			for (int lz=0; lz<2*maxlz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);
			
				Ar(lx, ly, lz) /= Ksqr;
			}
		}
	}  
	   
	// To avoid division by zero
	if (master)
		Ar(0,0,0) = 0.0;
}


/**********************************************************************************************

	Replaces A(k) by A(k)*exp(factor*K^2).

***********************************************************************************************/


void SSS_PENCIL::Array_exp_ksqr(Array<Complex,3> A, Real factor)
{	
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kxsqr;    // Kx^2
	Real Kxysqr;
	Real Ksqr;
	
	for (int lx=0; lx<maxlx; lx++) {
		Kxsqr= my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<maxly; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);    
		
			for (int lz=0; lz<2*maxlz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

			  Ar(lx, ly, lz) *= exp(factor*Ksqr);
		  }
		}
	}  
}



/**********************************************************************************************

	Replaces A(k) by A(k)*[ exp( factor*K^2+hyper_factor*K^4 )]

***********************************************************************************************/

void SSS_PENCIL::Array_exp_ksqr(Array<Complex,3> A, Real factor, Real hyper_factor, int hyper_exponent)
{
	
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);

	Real Kxsqr;    // Kx^2
	Real Kxysqr;
	Real Ksqr;
	Real Kpownm2;	// K^{q-2} where q = hyper_exponent
	
	for (int lx=0; lx<maxlx; lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
	
		for (int ly=0; ly<maxly; ly++) {
			Kxysqr = my_pow(Get_ky(ly)*kfactor[2],2);
			
			for (int lz=0; lz<2*maxlz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);
				
				if (hyper_exponent == 4)
					Kpownm2 = Ksqr;
				else
					Kpownm2 = my_pow(Ksqr,(hyper_exponent-2)/2);
				
				Ar(lx, ly, lz) *= exp((factor+hyper_factor*Kpownm2)* Ksqr);
			}
		}
	}
}
 


/**********************************************************************************************

	Replaces A(k) by A(k)*(V0.K)^2 / K^2.

***********************************************************************************************/


void SSS_PENCIL::Array_mult_V0_khat_sqr(Array<Complex,3> A, TinyVector<Real,3> V0)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kx, Ky, Kz;
	Real Kxsqr;    // Ky^2
	Real Kxysqr;
	Real Ksqr;
  
	Real V0x = V0(0);
	Real V0y = V0(1);
	Real V0z = V0(2);
	
	for (int lx=0; lx<maxlx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		Kxsqr = my_pow(Kx,2);
		
		for (int ly=0; ly<maxly; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			Kxysqr = Kxsqr + my_pow(Ky,2);
			
			for (int lz=0; lz<2*maxlz; lz++) {
				Kz = Get_kz(lz)*kfactor[3];
				Ksqr = Kxysqr + my_pow(Kz,2);

				Ar(lx, ly, lz) *= my_pow(V0x*Kx+V0y*Ky+V0z*Kz, 2)/Ksqr;
			}
		}
	}
	
	// To avoid division by zero
	if (my_id == master_id)
		Ar(0,0,0) = 0.0;
}

//*********************************************************************************************


// kz=0 is already filled.

void SSS_PENCIL::Fill_Vz(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	if (global.io.input_vx_vy_switch && global.field.incompressible) {
		
		Array<Real,3> Axr=Array<Real,3>(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		Array<Real,3> Ayr=Array<Real,3>(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		Array<Real,3> Azr=Array<Real,3>(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		
		int kx, ky, kz;

		int lz_min=0;
		if (my_z_pcoord==0)
			lz_min=1;
		
		for (int lx=0; lx<maxlx; lx++)
			for (int ly=0; ly<maxly; ly++)
				for (int lz=lz_min; lz<2*maxlz; lz++) {
					kx = Get_kx(lx);
					ky = Get_ky(ly);
					kz = Get_kz(lz);
					
					Last_component(kx, ky, kz, Axr(lx,ly,lz), Ayr(lx,ly,lz), Azr(lx,ly,lz));
				}
	}
}

int SSS_PENCIL::Read(Array<Complex,3> A, h5::Plan plan, string file_name, string dataset_name)
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

int SSS_PENCIL::Read(Array<Real,3> Ar, h5::Plan plan, string file_name, string dataset_name)
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


int SSS_PENCIL::Write(Array<Complex,3> A, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
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

int SSS_PENCIL::Write(Array<Real,3> Ar, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
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


