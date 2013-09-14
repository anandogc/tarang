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

/*! \file  Cvf.cc
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @sa Cvf.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bugs No known bugs
 */

#include "FluidVF.h"

					
/**********************************************************************************************

   Class constructor: Allocates for Vi, Ek etc.; Initializes the arrays.
   
**********************************************************************************************/

FluidVF::FluidVF
(
	DP dissipation_coefficient, 
	DP hyper_dissipation_coefficient, 
	int hyper_dissipation_exponent,
	bool force_switch,
	string field_name
 ) : PlainFluidVF(field_name)
{			

	this->dissipation_coefficient = dissipation_coefficient;
	this->hyper_dissipation_coefficient = hyper_dissipation_coefficient;
	this->hyper_dissipation_exponent = hyper_dissipation_exponent;
	
	if (hyper_dissipation_exponent == 0)
		hyper_dissipation_switch = false;
	
	this->force_switch = force_switch;
    
    nlin1.resize(shape_complex_array);
    nlin2.resize(shape_complex_array);
	nlin3.resize(shape_complex_array);
	
	
	// Memory allocation if force_switch is on
	if (force_switch) {
        Force1.resize(shape_complex_array);
        Force2.resize(shape_complex_array);
        Force3.resize(shape_complex_array);
	}
}
	

	//********************************************************************************************* 


void FluidVF::Compute_divergence_nlin(Array<complx,3> div)
{
    DP total_abs_div;  // not reqd for this function.
    universal->Compute_divergence(nlin1, nlin2, nlin3, div, "nlin", total_abs_div, false);
    // PS: global.temp_array.X is used in Compute_divergence_field.. so
    // be careful.
    // Used in pressure computation.
}

//********************************************************************************************* 


/** @brief Compute divergence of VF \f$ \vec{V} \f$.
 *
 *	@note	Divergence is stored in CVF of class NLIN.  We use the same class to store
 *			pressure also.  Hence divergence function should be called in the beginning 
 *			or end of loop when *F is free.
 *
 *  @return  \f$ *F = \mathcal{F}(D_i V_i) \f$. 
 */
void FluidVF::Compute_divergence_field(Array<complx,3> div, DP &total_abs_div, bool print_switch)
{
    universal->Compute_divergence(cvf.V1, cvf.V2, cvf.V3, div, "field", total_abs_div, print_switch);
    // PS: global.temp_array.X is used in Compute_divergence_field.. so
    // be careful.  We call div using global.temp_array.X2
}
	
//*********************************************************************************************	

DP	FluidVF::Get_mag_V0()
{
	DP modV0;
	
	if (my_id == master_id)
		modV0 = sqrt( pow2(real((cvf.V1)(0,0,0)))  +  pow2(real((cvf.V2)(0,0,0)))  +  pow2(real((cvf.V3)(0,0,0))) ) ;	

	MPI_Bcast( &modV0, 1, MPI_DP, master_id, MPI_COMM_WORLD);

	return modV0;														
}	
	
	
//*********************************************************************************************	

void FluidVF::Satisfy_strong_reality_condition_field()
{
	universal->Satisfy_strong_reality_condition_in_Array(cvf.V1);
	universal->Satisfy_strong_reality_condition_in_Array(cvf.V2);
	universal->Satisfy_strong_reality_condition_in_Array(cvf.V3);
}


//*********************************************************************************************	

void FluidVF::Satisfy_strong_reality_condition_force_field()
{

	universal->Satisfy_strong_reality_condition_in_Array(Force1);
	universal->Satisfy_strong_reality_condition_in_Array(Force2);
	universal->Satisfy_strong_reality_condition_in_Array(Force3);

}

//*********************************************************************************************

void FluidVF::Satisfy_weak_reality_condition_field()
{
	universal->Satisfy_weak_reality_condition_in_Array(cvf.V1);
	universal->Satisfy_weak_reality_condition_in_Array(cvf.V2);
	universal->Satisfy_weak_reality_condition_in_Array(cvf.V3);
}


//*********************************************************************************************

void FluidVF::Satisfy_weak_reality_condition_force_field()
{
    
	universal->Satisfy_weak_reality_condition_in_Array(Force1);
	universal->Satisfy_weak_reality_condition_in_Array(Force2);
	universal->Satisfy_weak_reality_condition_in_Array(Force3);
    
}
	

//******************************************************************************
void FluidVF::Test_reality_condition_field()
{
	universal->Test_reality_condition_in_Array(cvf.V1);
	universal->Test_reality_condition_in_Array(cvf.V2);
	universal->Test_reality_condition_in_Array(cvf.V3);
}
	

void FluidVF::Test_reality_condition_force_field()
{
    
	universal->Test_reality_condition_in_Array(Force1);
	universal->Test_reality_condition_in_Array(Force2);
	universal->Test_reality_condition_in_Array(Force3);
    
}

//******************************************************************************
void FluidVF::Dealias_force_field()
{
	if (global.program.alias_option == "DEALIAS") {
		universal->Dealias(Force1);
		universal->Dealias(Force2);
		universal->Dealias(Force3);
	}
}

//*********************************************************************************************	

void FluidVF::Zero_Prandtl_number_compute_temperature(FluidSF& T)
{
	T.csf.F = cvf.V1;
	T.csf.Divide_ksqr();
}


void FluidVF::Infinite_Prandtl_number_compute_velocity(FluidSF& T)
{
	//global.temp_array.X = T.csf.F;
	global.program.sincostr_switch = sincostr_switch_Vx; // cvf.V1 and T.csf.F has same basis
	universal->Xderiv(T.csf.F, global.temp_array.X);
	universal->Array_divide_ksqr(global.temp_array.X);  //X contains -Pressure
	
	// Compute Vy = -Dy(P)/(K^2); Vz = -Dz(P)/(K^2); Vx = (-Dx(P)+theta)/(K^2)
	// Adjusting the sign..
	if (N[2] > 1) {
		global.program.sincostr_switch = global.program.sincostr_switch_divergence;
		universal->Yderiv(global.temp_array.X, cvf.V2); 
		universal->Array_divide_ksqr(cvf.V2);
	}
	
	global.program.sincostr_switch = global.program.sincostr_switch_divergence;
	universal->Zderiv(global.temp_array.X, cvf.V3);
	universal->Array_divide_ksqr(cvf.V3);
	
	global.program.sincostr_switch = global.program.sincostr_switch_divergence;
	universal->Xderiv(global.temp_array.X, cvf.V1);   // V1 = Dx(P)
	cvf.V1 = cvf.V1 + T.csf.F;
	universal->Array_divide_ksqr(cvf.V1);
	
	if (global.PHYSICS.Uscaling == "ULARGE") {
		cvf.V1 = cvf.V1*sqrt(global.PHYSICS.Rayleigh);
		cvf.V2 = cvf.V2*sqrt(global.PHYSICS.Rayleigh);
		cvf.V3 = cvf.V3*sqrt(global.PHYSICS.Rayleigh);
	}
    
    else if (global.PHYSICS.Uscaling == "USMALL") {
		cvf.V1 = cvf.V1*global.PHYSICS.Rayleigh;
		cvf.V2 = cvf.V2*global.PHYSICS.Rayleigh;
		cvf.V3 = cvf.V3*global.PHYSICS.Rayleigh;
	}
    
}


/**********************************************************************************************
 Copy_field_to(W):	  W.V <- V
 Copy_field_from(W):   V <-W. V
***********************************************************************************************/

void FluidVF::Copy_field_to(PlainCVF& W)
{	
	W.V1 = cvf.V1;  W.V2 = cvf.V2;	W.V3 = cvf.V3;
} 

void FluidVF::Copy_field_to(CVF& W)
{	
	W.V1 = cvf.V1;  W.V2 = cvf.V2;	W.V3 = cvf.V3;
}

void FluidVF::Copy_field_from(PlainCVF& W)
{	
	cvf.V1 = W.V1;  cvf.V2 = W.V2;	cvf.V3 = W.V3;
}
  

void FluidVF::Copy_field_from(CVF& W)
{	
	cvf.V1 = W.V1;  cvf.V2 = W.V2;	cvf.V3 = W.V3;
}


/**********************************************************************************************
 
Function 1:  V = V + nlin*dt;  
Function 2: for CVF: Y = Y+ factor*dt*U.nlin; 
 
***********************************************************************************************/

void FluidVF::Add_nlin_factor_dt(DP factor)
{	
	if (abs(factor) > MYEPS2) {
		cvf.V1 += (factor*global.time.dt)*(nlin1);
		cvf.V2 += (factor*global.time.dt)*(nlin2);
		cvf.V3 += (factor*global.time.dt)*(nlin3);
	}
	
}


void FluidVF::Add_field_nlin_factor_dt(PlainCVF& Y, DP factor)
{	
	if (abs(factor) > MYEPS2) {
		Y.V1 += (factor*global.time.dt)* (nlin1);         
		Y.V2 += (factor*global.time.dt)* (nlin2);
		Y.V3 += (factor*global.time.dt)* (nlin3);
	}
}

//*********************************************************************************************

/** @brief Multiply vector field by \f$ \exp(-\nu K^2 dt) \f$ or 
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ V(k) = V(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ V(k) = V(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */
void FluidVF::Mult_field_exp_ksqr_dt(DP a)
{
	if (abs(a) > MYEPS) {
		if (!hyper_dissipation_switch) {
			universal->Array_exp_ksqr(cvf.V1, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(cvf.V2, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(cvf.V3, -dissipation_coefficient*a*global.time.dt);
		}
	
		else {
			universal->Array_exp_ksqr(cvf.V1, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(cvf.V2, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(cvf.V3, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
		}
	}
	
}



//*********************************************************************************************
/** @brief Multiply nonlinear field by \f$ \exp(-\nu K^2 dt) \f$ or 
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ N(k) = N(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ N(k) = N(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */

void FluidVF::Mult_nlin_exp_ksqr_dt(DP a)
{
	if (abs(a) > MYEPS) {
		if (!hyper_dissipation_switch) {
			universal->Array_exp_ksqr(nlin1, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(nlin2, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(nlin3, -dissipation_coefficient*a*global.time.dt);
		}
	
		else {
			universal->Array_exp_ksqr(nlin1, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(nlin2, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(nlin3, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
		}
	}
}



//******************************************************************************

/** @brief D=3: Add complex conjugate at \f$ \vec{k} \f$ 
 *				once \f$ \vec{V}(kx, ky, kz) \f$ is given.
 *  
 *	@bug  Add_complex_conj does not work across the processors.
 * @note  For ky = 0 and ky = N(2)/2 planes only.
 * 
 * @return  FFF: \f$ \vec{V}(-kx, -ky, kz) = \vec{V}^*(kx, ky, kz)\f$
 * @return  SCFT: For (ky != 0) -> \f$ \vec{V}(kx, -ky, kz) = \vec{V}^*(kx, ky, kz)\f$.
 * @return  SCFT: For (ky = 0) -> \f$ imag(V_y(kx, 0, kz)) = imag(V_z(kx, 0, kz)) = 0\f$,
 *					and \f$ V_x(kx, 0, kz) = 0 \f$.
 *				
 */
void FluidVF::Add_complex_conj(int kx, int ky, int kz, complx Vx, complx Vy, complx Vz)			
{	
	TinyVector<complx,3> localV;
	
		// On kz=0 or kz=N[3]/2 plane
	if ((basis_type == "FFF") && (kz == 0)) {
		localV = conj(Vx), conj(Vy), conj(Vz);
		universal->Assign_spectral_field(-kx, -ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);
		cout << "Complex-conj(V) added for k = (" << -kx << "," << -ky << "," << kz << ")" << endl;	
	}
	
	else if ((basis_type == "SFF") && (kz == 0)) {          
		localV = conj(Vx), conj(Vy), conj(Vz);
		universal->Assign_spectral_field(kx, -ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);
		cout << "Complex-conj(V) added for k = (" << kx << "," << -ky << "," << kz << ")" << endl;	
	}
}

void FluidVF::Assign_field_and_comp_conj(int kx, int ky, int kz, complx Vx, complx Vy, complx Vz)
{	
    
    TinyVector<complx,3> localV(Vx, Vy, Vz);
	
	universal->Assign_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);	// in appropriate proc.
	
	Add_complex_conj(kx, ky, kz, Vx, Vy, Vz);	
}

/// 3D: Generate random vector at (kx,ky,kz) with range = rand_range

void FluidVF::Assign_random_complex_vector(int kx, int ky, int kz, DP rand_range)
{
	complx  Vx, Vy, Vz;
	
	Vx = complx( 2*rand_range*(SPECrand.random()-0.5), 2*rand_range*(SPECrand.random()-0.5) );
	
	Vy = complx( 2*rand_range*(SPECrand.random()-0.5), 2*rand_range*(SPECrand.random()-0.5) );
    
   // universal->Last_component(kx, ky, kz, Vx, Vy, Vz);
	
	Assign_field_and_comp_conj(kx, ky, kz, Vx, Vy, Vz);
}

void FluidVF::Assign_random_real_vector(int kx, int ky, int kz, DP rand_range)
{
	DP Vx, Vy, Vz;
	
	Vx = 2*rand_range*(SPECrand.random()-0.5);
	
	Vy = 2*rand_range*(SPECrand.random()-0.5);
	
	universal->Last_component(kx, ky, kz, Vx, Vy, Vz);
	
	//Assign_field(kx, ky, kz, Vx, Vy, Vz);
}


void FluidVF::Assign_field(int kx, int ky, int kz, DP Vx, DP Vy, DP Vz)
{	
    
    TinyVector<DP,3> localV(Vx, Vy, Vz);
	
	universal->Assign_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);	// in appropriate proc.
}

//*********************************************************************************************
/// 3D: Return Tk = Real(-nlin(k). conj(V(k))  for vector V
DP FluidVF::Get_Tk(int kx, int ky, int kz)
{	
	TinyVector<complx,3> localV_complex, localnlin_complex;
	TinyVector<DP,3> localV_real, localnlin_real;
	
	if (universal->Probe_in_me(kx,ky,kz)) {
		if (basis_type == "SSS") {
			localV_real = real(universal->Get_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3));
			localnlin_real = real(universal->Get_spectral_field(kx, ky, kz, nlin1, nlin2, nlin3));
			return dot(localV_real, localnlin_real);
		}
		
		else {
			localV_complex = universal->Get_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3);
			localnlin_complex = universal->Get_spectral_field(kx, ky, kz, nlin1, nlin2, nlin3);
			return mydot(localV_complex, localnlin_complex);
		}
	}

	return 0;
}


//*********************************************************************************************

void FluidVF::Get_local_max_real_space(DP local_ur[])
{	
	local_ur[0] = max(abs(rvf.V1r));
	local_ur[2] = max(abs(rvf.V3r));
		
	if (Ny > 1) 
		local_ur[1] = max(abs(rvf.V2r));
	else
		local_ur[1] = 0.0;
}

//*********************************************************************************************
/*
DP FluidVF::Get_dt()
{
	DP local_ux,local_uy,local_uz;
	DP ux, uy, uz;
	
	if (global.program.dt_option == 0) {
		return global.time.dt_fixed;
	}
	
	else if (global.program.dt_option == 1) {
		Get_local_max_real_space(local_ux,local_uy,local_uz);
		
		MPI_Reduce(&local_ux, &ux, 1, MPI_DP, MPI_MAX, master_id, MPI_COMM_WORLD);
		MPI_Reduce(&local_uy, &uy, 1, MPI_DP, MPI_MAX, master_id, MPI_COMM_WORLD);
		MPI_Reduce(&local_uz, &uz, 1, MPI_DP, MPI_MAX, master_id, MPI_COMM_WORLD);
		
		DP sum_ubydx = ux/global.field.Delta_x[1] +  uz/global.field.Delta_x[3];
		
		if (Ny > 1)
			sum_ubydx += uy/global.field.Delta_x[2];
		
		DP dt = global.time.Courant_no/sum_ubydx;
		
		MPI_Bcast(&dt, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
		
		return min(dt, global.time.dt_fixed);
	}
}	
*/



DP FluidVF::Get_dt()
{
	DP local_ur[3];
	DP ur[3];
	
	if (global.program.dt_option == 0) {
		return global.time.dt_fixed;
	}
	
	else if (global.program.dt_option == 1) {
		
		Get_local_max_real_space(local_ur);
		
		MPI_Allreduce(local_ur, ur, 3, MPI_DP, MPI_MAX, MPI_COMM_WORLD);
		
		DP sum_ubydx = ur[0]/global.field.Delta_x[1] +  ur[2]/global.field.Delta_x[3];
		
		if (Ny > 1)
			sum_ubydx += ur[1]/global.field.Delta_x[2];
		
		DP dt = global.time.Courant_no/sum_ubydx;
		
		return min(dt, global.time.dt_fixed);
	}

	return -1;
}	


DP FluidVF::Get_dt(FluidSF& T)
{
	return Get_dt();
}


// Alfven waves move with the speed of local mean magnetic field.. So we look for Bmax
// in the real space.
DP FluidVF::Get_dt(FluidVF& W)
{
	
	
	if (global.program.dt_option == 0) {
		return global.time.dt_fixed;
	}
	
	else if (global.program.dt_option == 1) {
		
		DP local_ur[3], local_wr[3];
		
		DP local_uwr[6];
		DP uwr[6];
		
		Get_local_max_real_space(local_ur);
		W.Get_local_max_real_space(local_wr);
		
		local_uwr[0] = local_ur[0];
		local_uwr[1] = local_ur[1];
		local_uwr[2] = local_ur[2];
		local_uwr[3] = local_wr[0]; // w fields
		local_uwr[4] = local_wr[1];
		local_uwr[5] = local_wr[2];
		
		MPI_Allreduce(local_uwr, uwr, 6, MPI_DP, MPI_MAX, MPI_COMM_WORLD);
		
		// for fluid
		DP sum_ubydx = (uwr[0]+uwr[3])/global.field.Delta_x[1] +  (uwr[2]+uwr[5])/global.field.Delta_x[3];
		
		if (Ny > 1)
			sum_ubydx += (uwr[1]+uwr[4])/global.field.Delta_x[2];
		
		DP dt = global.time.Courant_no/sum_ubydx;
		
		return min(dt, global.time.dt_fixed);
	}

	return -1;
}


DP FluidVF::Get_dt(FluidVF& W, FluidSF& T)
{
	return Get_dt(W);
}


//*********************************************************************************************



//************************ END of CVF class Definitions ***************************************
