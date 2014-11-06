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


void SSS_PENCIL::Print_large_Fourier_elements(Array<complx,3> A, string array_name)
{
	Array<DP,3> B=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				if (abs(B(lx,ly,lz)) > MYEPS2)
					cout << "my_id = " << my_id <<  " vect(k) = (" << Get_kx(lx) << "," << ly  << "," << lz <<");  " << array_name << "(k) = " <<  B(lx, ly, lz) << '\n';
			}
	
	cout << endl;
}


void SSS_PENCIL::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
{}

void SSS_PENCIL::Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz)
{
	DP Kx = kx*kfactor[1];
	DP Ky = ky*kfactor[2];
	DP Kz = kz*kfactor[3];
	
	DP dvxdx, dvydy;
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

void SSS_PENCIL::Dealias(Array<complx,3> A)
{
	Assign_sub_array(Range(2*Nx/3+1,toEnd), Range(2*Ny/3+1,toEnd), Range(Nz/3,toEnd), A, complx(0,0));
}


// Data resides till outer_radius in k-space
bool SSS_PENCIL::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
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

void SSS_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
	return; // Do nothing
}

void SSS_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{
	return; // Do nothing
}

void SSS_PENCIL::Test_reality_condition_in_Array(Array<complx,3> A)
{
	return; // Do nothing
}


//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SSS_PENCIL::Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
	// lx = 0 reside in master node
	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_Vx;
	
	if (global.program.sincostr_switch == "SCC") {
		Assign_sub_array(zero,Range::all(),Range::all(),Ax,complx(0,0));
		Assign_sub_array(Range::all(),zero,Range::all(),Ay,complx(0,0));
		Assign_sub_array(Range::all(),Range::all(),zero,Az,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "CSS") {
		Assign_sub_array(Range::all(),zero,zero,Ax,complx(0,0));
		Assign_sub_array(zero,zero,Range::all(),Ay,complx(0,0));
		Assign_sub_array(zero,zero,Range::all(),Az,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "CCS") {
		Assign_sub_array(Range::all(),Range::all(),zero,Ax,complx(0,0));
		Assign_sub_array(zero,zero,zero,Ay,complx(0,0));
		Assign_sub_array(zero,Range::all(),Range::all(),Az,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "SSC") {
		Assign_sub_array(zero,zero,Range::all(),Ax,complx(0,0));
		Assign_sub_array(Range::all(),zero,zero,Az,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "CSC") {
		Assign_sub_array(Range::all(),zero,Range::all(),Ax,complx(0,0));
		Assign_sub_array(zero,Range::all(),Range::all(),Ay,complx(0,0));
		Assign_sub_array(zero,zero,zero,Az,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "SCS") {
		Assign_sub_array(zero,zero,Range::all(),Ax,complx(0,0));
		Assign_sub_array(Range::all(),zero,zero,Ay,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "CCC") {
		Assign_sub_array(zero,zero,Range::all(),Ay,complx(0,0));
		Assign_sub_array(zero,Range::all(),zero,Az,complx(0,0));
	}
	
	else if (global.program.sincostr_switch == "SSS") {
		Assign_sub_array(zero,zero,zero,Ax,complx(0,0));
		Assign_sub_array(Range::all(),Range::all(),zero,Ay,complx(0,0));
		Assign_sub_array(Range::all(),zero,Range::all(),Az,complx(0,0));
	}
}
	
/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SSS_PENCIL::Zero_modes(Array<complx,3> F)
{
	// lx = 0 reside in master node
	Range zero(0,0);
	global.program.sincostr_switch = sincostr_switch_F;
	
	if (global.program.sincostr_switch == "SCC")
		Assign_sub_array(zero,Range::all(),Range::all(),F,complx(0,0));
	
	else if (global.program.sincostr_switch == "CSS")
		Assign_sub_array(Range::all(),zero,zero,F,complx(0,0));
	
	else if (global.program.sincostr_switch == "CCS")
		Assign_sub_array(Range::all(),Range::all(),zero,F,complx(0,0));

	else if (global.program.sincostr_switch == "SSC")
		Assign_sub_array(zero,zero,Range::all(),F,complx(0,0));
	
	else if (global.program.sincostr_switch == "CSC")
		Assign_sub_array(Range::all(),zero,Range::all(),F,complx(0,0));
	
	else if (global.program.sincostr_switch == "SCS")
		Assign_sub_array(zero,zero,Range::all(),F,complx(0,0));
	
	else if (global.program.sincostr_switch == "CCC")
		;
	
	else if (global.program.sincostr_switch == "SSS")
		Assign_sub_array(zero,zero,zero,F,complx(0,0));
}

void SSS_PENCIL::Assign_sub_array(Range x_range, Range y_range, Range z_range, Array<complx,3> A, complx value)
{
    Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);


	static Array<int,1> x_filter(Nx);
	static Array<int,1> y_filter(Ny);
	
	x_filter = 0;
	y_filter = 0;
	
	x_filter(x_range)=1;
	y_filter(y_range)=1;
	
	
	static Range y_apply, x_apply;
	
	
	x_apply = Range(first(x_filter(Range(my_x_pcoord*local_Nx,(my_x_pcoord+1)*local_Nx-1)) == 1),
					 last(x_filter(Range(my_x_pcoord*local_Nx,(my_x_pcoord+1)*local_Nx-1)) == 1));
	
	y_apply = Range(first(y_filter(Range(my_y_pcoord*local_Ny,(my_y_pcoord+1)*local_Ny-1)) == 1),
					 last(y_filter(Range(my_y_pcoord*local_Ny,(my_y_pcoord+1)*local_Ny-1)) == 1));
	
	// if (my_x_pcoord==0) {
	// 	cout << "apply shape = " << x_apply << " " << y_apply << " " <<z_range << endl;
	// }
	if ( (x_apply(0)>=0) && (y_apply(0)>=0) ) {
		// Print_large_Fourier_elements(A, "Before");
		Ar(x_apply, y_apply, z_range) = real(value);
		// Print_large_Fourier_elements(A, "After");
	}
}


/**********************************************************************************************

			Replaces A(k) by A(k)*K^2.

***********************************************************************************************/


void SSS_PENCIL::Array_mult_ksqr(Array<complx,3> A)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Kxsqr;    // Kx^2
	DP Kxysqr;
	DP Ksqr;
	
	// #pragma omp parallel for private(Kysqr,Kyzsqr,Ksqr) 
	for (int lx=0; lx<local_Nx; lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<local_Ny; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);	
			
			for (int lz=0; lz<local_Nz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

				Ar(lx, ly, lz) *= Ksqr;
		  }
		}
	}
}


/**********************************************************************************************

	Replaces A(k) by A(k)/K^2 with K^2 = sum_i [ki*kfactor(i)]^2 ; A(0)=0

***********************************************************************************************/



void SSS_PENCIL::Array_divide_ksqr(Array<complx,3> A)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Kxsqr;    // Kx^2
	DP Kxysqr;
	DP Ksqr;
	
	// #pragma omp parallel for private(Kysqr,Kyzsqr,Ksqr) 
	for (int lx=0; lx<local_Nx; lx++) {
		Kxsqr= my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<local_Ny; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);    
			
			for (int lz=0; lz<local_Nz; lz++) {
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


void SSS_PENCIL::Array_exp_ksqr(Array<complx,3> A, DP factor)
{	
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Kxsqr;    // Kx^2
	DP Kxysqr;
	DP Ksqr;
	
	// #pragma omp parallel for private(Kysqr,Kyzsqr,Ksqr) 
	for (int lx=0; lx<local_Nx; lx++) {
		Kxsqr= my_pow(Get_kx(lx)*kfactor[1],2);

		for (int ly=0; ly<local_Ny; ly++) {
			Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);    
		
			for (int lz=0; lz<local_Nz; lz++) {
				Ksqr = Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

			  Ar(lx, ly, lz) *= exp(factor*Ksqr);
		  }
		}
	}  
}



/**********************************************************************************************

	Replaces A(k) by A(k)*[ exp( factor*K^2+hyper_factor*K^4 )]

***********************************************************************************************/

void SSS_PENCIL::Array_exp_ksqr(Array<complx,3> A, DP factor, DP hyper_factor, int hyper_exponent)
{
	
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);

	DP Kxsqr;    // Kx^2
	DP Kxysqr;
	DP Ksqr;
	DP Kpownm2;	// K^{q-2} where q = hyper_exponent
	
	for (int lx=0; lx<local_Nx; lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
	
		for (int ly=0; ly<local_Ny; ly++) {
			Kxysqr = my_pow(Get_ky(ly)*kfactor[2],2);
			
			for (int lz=0; lz<local_Nz; lz++) {
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


void SSS_PENCIL::Array_mult_V0_khat_sqr(Array<complx,3> A, TinyVector<DP,3> V0)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Kx, Ky, Kz;
	DP Kxsqr;    // Ky^2
	DP Kxysqr;
	DP Ksqr;
  
	DP V0x = V0(0);
	DP V0y = V0(1);
	DP V0z = V0(2);
	
	for (int lx=0; lx<local_Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		Kxsqr = my_pow(Kx,2);
		
		for (int ly=0; ly<local_Ny; ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			Kxysqr = Kxsqr + my_pow(Ky,2);
			
			for (int lz=0; lz<local_Nz; lz++) {
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

void SSS_PENCIL::Fill_Vz(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
	if (global.io.input_vx_vy_switch && global.field.incompressible) {
		
		Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		
		int kx, ky, kz;
		DP vz;
		
		// #pragma omp parallel for
		for (int lx=0; lx<local_Nx; lx++)
			for (int ly=0; ly<local_Ny; ly++)
				for (int lz=1; lz<local_Nz; lz++) {
					kx = Get_kx(lx);
					ky = Get_ky(ly);
					kz = Get_kz(lz);
					
					Last_component(kx, ky, kz, Axr(lx,ly,lz), Ayr(lx,ly,lz), vz);
					
					Assign_spectral_field(kx, ky, kz, Az, vz);
				}
	}
}

int SSS_PENCIL::Read(Array<complx,3> A, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	int err = BasicIO::Read(global.temp_array.Xr.data(), plan, file_name, dataset_name);
	spectralTransform.Transpose(global.temp_array.Xr, A);
	return err;
}

int SSS_PENCIL::Read(Array<DP,3> Ar, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	return BasicIO::Read(Ar.data(), plan, file_name, dataset_name);
}


int SSS_PENCIL::Write(Array<complx,3> A, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	spectralTransform.Transpose(A, global.temp_array.Xr);
	return BasicIO::Write(global.temp_array.Xr.data(), plan, folder_name, file_name, dataset_name);
}

int SSS_PENCIL::Write(Array<DP,3> Ar, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	return BasicIO::Write(Ar.data(), plan, folder_name, file_name, dataset_name);  
}

//*********************************  End of scft_basic.cc *************************************


