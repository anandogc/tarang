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
 

#include "SSS_slab.h"
#include "SSS_slab_inline.h"


/**********************************************************************************************
 
 COmment
 
***********************************************************************************************/



void SSS_SLAB::Print_large_Fourier_elements(Array<Complex,3> A, string array_name)
{
	Array<Real,3> B=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	for (int lx=0; lx<B.extent(0); lx++)
		for (int ly=0; ly<B.extent(1); ly++)
			for (int lz=0; lz<B.extent(2); lz++) {
				if (abs(B(lx,ly,lz)) > MYEPS2)
					cout << "my_id = " << my_id <<  " vect(k) = (" << Get_kx(lx) << "," << ly  << "," << lz <<");  "<< array_name << "(k) = " <<  B(lx, ly, lz) << '\n';
			}

	cout << endl;
}

//**************************************************************************************

void SSS_SLAB::Last_component(int kx, int ky, int kz, Complex &Vx, Complex &Vy, Complex &Vz)
{}

void SSS_SLAB::Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz)
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


/*****************************************************************************************
 
Dealias
 
 A(Range(2*Ny/3+1,toEnd), Range(2*Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0;
 
 In complex array A: z-range = Nz/3

 
 *****************************************************************************************/

void SSS_SLAB::Dealias(Array<Complex,3> A)
{
	Array<int,1> Ax_filter(Nx);
	
	Ax_filter = 0;
	
	Ax_filter(Range(2*Nx/3+1,toEnd)) = 1;
	
	int first_x = first(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
	int last_x = last(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
	
	A(Range(first_x,last_x), Range(2*Ny/3+1,toEnd), Range(Nz/3,toEnd)) = 0.0;
}

// Data resides till outer_radius in k-space
bool SSS_SLAB::Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > 2*Nx/3) || (ky_max > 2*Ny/3) || (kz_max > 2*Nz/3))
		return true;
	else
		return false;
}


/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 A is real aleardy.. do nothing.
 
 ***********************************************************************************************/

void SSS_SLAB::Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}

void SSS_SLAB::Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}

void SSS_SLAB::Test_reality_condition_in_Array(Array<Complex,3> A)
{
	return; // Do nothing
}


//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SSS_SLAB::Zero_modes(Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz)
{
	Array<Real,3> Ax(reinterpret_cast<Real*>(Bx.data()), Bx.shape()*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ay(reinterpret_cast<Real*>(By.data()), By.shape()*shape(1,1,2), neverDeleteData);
	Array<Real,3> Az(reinterpret_cast<Real*>(Bz.data()), Bz.shape()*shape(1,1,2), neverDeleteData);

	// lx = 0 reside in master node
	global.program.sincostr_switch = sincostr_switch_Vx;
	
	if (global.program.sincostr_switch == "SCC") {
		if (master)
			Ax(0,Range::all(),Range::all()) = 0.0;

		if (Ny>1)
			Ay(Range::all(),0,Range::all()) = 0.0;

		Az(Range::all(),Range::all(),0) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CSS") {
		if (Ny>1)
			Ax(Range::all(),0,Range::all()) = 0.0;

		Ax(Range::all(),Range::all(),0) = 0.0;
		
		if (master)
			Ay(0,Range::all(),Range::all()) = 0.0;
		Ay(Range::all(),Range::all(),0) = 0.0;
		
		if (master)
			Az(0,Range::all(),Range::all()) = 0.0;
		
		if (Ny>1)
			Az(Range::all(),0,Range::all()) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CCS") {
		Ax(Range::all(),Range::all(),0)  = 0.0;
		
		if (master)
			Ay(0,Range::all(),Range::all()) = 0.0;
		if (Ny>1)
			Ay(Range::all(),0,Range::all()) = 0.0;
		Ay(Range::all(),Range::all(),0) = 0.0;
	   
		if (master)
			Az(0,Range::all(),Range::all()) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "SSC") {
		if (master)
			Ax(0,Range::all(),Range::all()) = 0.0;
		if (Ny>1)
			Ax(Range::all(),0,Range::all()) = 0.0;
		
		if (Ny>1)
			Az(Range::all(),0,Range::all()) = 0.0;
		Az(Range::all(),Range::all(),0) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CSC") {
		if (Ny>1)
			Ax(Range::all(),0,Range::all())  = 0.0;
		
		if (master)
			Ay(0,Range::all(),Range::all()) = 0.0;
		
		if (master)
			Az(0,Range::all(),Range::all()) = 0.0;
		if (Ny>1)
			Az(Range::all(),0,Range::all()) = 0.0;
		Az(Range::all(),Range::all(),0) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "SCS") {
		if (master)
			Ax(0,Range::all(),Range::all()) = 0.0;
		Ax(Range::all(),Range::all(),0) = 0.0;
		
		if (Ny>1)
			Ay(Range::all(),0,Range::all()) = 0.0;
		Ay(Range::all(),Range::all(),0) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CCC") {
		if (master)
			Ay(0,Range::all(),Range::all()) = 0.0;
		if (Ny>1)
			Ay(Range::all(),0,Range::all()) = 0.0;
		
		if (master)
			Az(0,Range::all(),Range::all()) = 0.0;
		Az(Range::all(),Range::all(),0) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "SSS") {
		if (master)
			Ax(0,Range::all(),Range::all()) = 0.0;
		if (Ny>1)
			Ax(Range::all(),0,Range::all()) = 0.0;
		Ax(Range::all(),Range::all(),0) = 0.0;
		
		if (master)
			Ay(Range::all(),Range::all(),0) = 0.0;
		
		if (Ny>1)
			Az(Range::all(),0,Range::all()) = 0.0;
	}
}
	
/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SSS_SLAB::Zero_modes(Array<Complex,3> F)
{
	Array<Real,3> Fr(reinterpret_cast<Real*>(F.data()), F.shape()*shape(1,1,2), neverDeleteData);

	// lx = 0 reside in master node
	
	global.program.sincostr_switch = sincostr_switch_F;
	
	if (global.program.sincostr_switch == "SCC") {
		if (master)
			Fr(0,Range::all(),Range::all()) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CSS"){
		if (Ny>1)
			Fr(Range::all(),0,Range::all()) = 0.0;
		Fr(Range::all(),Range::all(),0) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CCS")
		F(Range::all(),Range::all(),0)  = 0.0;
		
	else if (global.program.sincostr_switch == "SSC") {
		if (master)
			Fr(0,Range::all(),Range::all()) = 0.0;
		if (Ny>1)
			Fr(Range::all(),0,Range::all()) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CSC") {
		if (Ny>1)
			Fr(Range::all(),0,Range::all())  = 0.0;
	}
			  
	else if (global.program.sincostr_switch == "SCS") {
		if (master)
			Fr(0,Range::all(),Range::all()) = 0.0;
		Fr(Range::all(),Range::all(),0) = 0.0;
	}
	
	else if (global.program.sincostr_switch == "CCC")
		;
	
	else if (global.program.sincostr_switch == "SSS") {
		if (master)
			Fr(0,Range::all(),Range::all()) = 0.0;
		if (Ny>1)
			Fr(Range::all(),0,Range::all()) = 0.0;
		Fr(Range::all(),Range::all(),0) = 0.0;
	} 

}


/**********************************************************************************************

			Replaces A(k) by A(k)*K^2.

***********************************************************************************************/


void SSS_SLAB::Array_mult_ksqr(Array<Complex,3> A)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kxsqr;
	Real Kxysqr;
	Real Ksqr;
	
	//#pragma omp parallel for private(Kysqr,Kyzsqr,Ksqr) 
	for (int lx=0; lx<Ar.extent(0); lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
		
        for (int ly=0; ly<Ar.extent(1); ly++) {
            Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);

            for (int lz=0; lz<Ar.extent(2); lz++) {
                Ksqr= Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

				Ar(lx, ly, lz) *= Ksqr;
        	}
        }
    }
}


/**********************************************************************************************

	Replaces A(k) by A(k)/K^2 with K^2 = sum_i [ki*kfactor(i)]^2 ; A(0)=0

***********************************************************************************************/



void SSS_SLAB::Array_divide_ksqr(Array<Complex,3> A)
{
	
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kxsqr;
	Real Kxysqr;
	Real Ksqr;
	
	//#pragma omp parallel for private(Kysqr,Kyzsqr,Ksqr) 
	for (int lx=0; lx<Ar.extent(0); lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
		
        for (int ly=0; ly<Ar.extent(1); ly++) {
            Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);

            for (int lz=0; lz<Ar.extent(2); lz++) {
                Ksqr= Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

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


void SSS_SLAB::Array_exp_ksqr(Array<Complex,3> A, Real factor)
{
	
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);

	Real Kxsqr;
	Real Kxysqr;
	Real Ksqr;
	
	//#pragma omp parallel for private(Kysqr,Kyzsqr,Ksqr) 
	for (int lx=0; lx<Ar.extent(0); lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
		
        for (int ly=0; ly<Ar.extent(1); ly++) {
            Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);

            for (int lz=0; lz<Ar.extent(2); lz++) {
                Ksqr= Kxysqr + my_pow(Get_kz(lz)*kfactor[3],2);

				Ar(lx, ly, lz) *= exp(factor*Ksqr);
        	}
        }
    }
}



/**********************************************************************************************

	Replaces A(k) by A(k)*[ exp( factor*K^2+hyper_factor*K^4 )]

***********************************************************************************************/

void SSS_SLAB::Array_exp_ksqr(Array<Complex,3> A, Real factor, Real hyper_factor, int hyper_exponent)
{
	
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);

	Real Kxsqr;
	Real Kxysqr;
	Real Ksqr;
	Real Kpownm2;	// K^{q-2} where q = hyper_exponent
	
	//#pragma omp parallel for private(Kysqr,Kyzsqr,Ksqr) 
	for (int lx=0; lx<Ar.extent(0); lx++) {
		Kxsqr = my_pow(Get_kx(lx)*kfactor[1],2);
		
        for (int ly=0; ly<Ar.extent(1); ly++) {
            Kxysqr = Kxsqr + my_pow(Get_ky(ly)*kfactor[2],2);

            for (int lz=0; lz<Ar.extent(2); lz++) {
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


void SSS_SLAB::Array_mult_V0_khat_sqr(Array<Complex,3> A, TinyVector<Real,3> V0)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kx, Ky, Kz;
	Real Kxsqr;    // Ky^2
	Real Kxysqr;
	Real Ksqr;
  
	Real V0x = V0(0);
	Real V0y = V0(1);
	Real V0z = V0(2);

	for (int lx=0; lx<Ar.extent(0); lx++) {
		Kx = Get_kx(lx)*kfactor[1];
		Kxsqr = my_pow(Kx,2);
		
		for (int ly=0; ly<Ar.extent(1); ly++) {
			Ky = Get_ky(ly)*kfactor[2];
			Kxysqr = Kxsqr + my_pow(Ky,2);
			
			for (int lz=0; lz<Ar.extent(2); lz++) {
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

void SSS_SLAB::Fill_Vz(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	if (global.io.input_vx_vy_switch && global.field.incompressible) {
		
		Array<Real,3> Axr=Array<Real,3>(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		Array<Real,3> Ayr=Array<Real,3>(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		Array<Real,3> Azr=Array<Real,3>(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
		
		int kx, ky, kz;
		Real vz;
		
		//#pragma omp parallel for
		for (int lx=0; lx<Axr.extent(0); lx++)
			for (int ly=0; ly<Axr.extent(1); ly++)
				for (int lz=1; lz<Axr.extent(2); lz++) {
					kx = Get_kx(lx);
					ky = Get_ky(ly);
					kz = Get_kz(lz);
					
					Last_component(kx, ky, kz, Axr(lx,ly,lz), Ayr(lx,ly,lz), vz);
					
					Assign_spectral_field(kx, ky, kz, Az, vz);
				}
	}
}

int SSS_SLAB::Read(Array<Complex,3> A, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	return BasicIO::Read(A.data(), plan, file_name, dataset_name);
}

int SSS_SLAB::Read(Array<Real,3> Ar, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	int err = BasicIO::Read(global.temp_array.X.data(), plan, file_name, dataset_name);
	if (Ny>1)
		spectralTransform.Transpose(global.temp_array.X, Ar);
	else
		spectralTransform.Transpose(global.temp_array.X(Range::all(),0,Range::all()), Ar(Range::all(),0,Range::all()));

	return err;
}


int SSS_SLAB::Write(Array<Complex,3> A, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	return BasicIO::Write(A.data(), plan, folder_name, file_name, dataset_name);
}

int SSS_SLAB::Write(Array<Real,3> Ar, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	if (Ny>1)
		spectralTransform.Transpose(Ar, global.temp_array.X);
	else
		spectralTransform.Transpose(Ar(Range::all(),0,Range::all()), global.temp_array.X(Range::all(),0,Range::all()));

	return BasicIO::Write(global.temp_array.X.data(), plan, folder_name, file_name, dataset_name);  
}

//*********************************  End of scft_basic.cc *************************************


