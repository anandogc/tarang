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
/*! \file  compute_force_RB.cc
 * 
 * @brief Force for RB convection
 *
 *		For Pr_large, u_small:		 F1 = Ra*Pr*T(k);  
 *		For Pr_large, u_large:		 F1 = Ra*T(k);    	
 *		For Pr_small or 0, u_small:  F1 = Ra*T(k);	   
 *		For Pr_small or 0, u_large:  F1 = Pr*T(k)
 *
 *		For Pr_large				Ftheta = V1;    	
 *		For Pr_small				Ftheta = V1/Pr;	   
 *		For Pr = 0					Ftheta = 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 
#include "FORCE.h"


//*********************************************************************************************

void FORCE::Compute_force_RBC_basic_assign(FluidVF& U, FluidSF& T)
{
	
	if (global.program.basis_type == "ChFF") { // box size (2,Ly,Lz)
        U.Force1 =  0.125*(global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)*(T.csf.F);
        U.Force2 = 0.0;
        U.Force3 = 0.0;
        T.Force = 0.5*(U.cvf.V1);
    }
    
    else {
	// For the velocity field
		if (global.PHYSICS.Pr_option == "PRLARGE") {
			
			if (global.PHYSICS.Uscaling == "USMALL") 		
				U.Force1 =  (global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)*(T.csf.F);	// (u.grad)u-Ra*Pr*theta
			
			else if (global.PHYSICS.Uscaling  == "ULARGE") 		
				U.Force1 =  T.csf.F;					// (u.grad)u-theta			
		}
		
		else if ((global.PHYSICS.Pr_option == "PRSMALL") || (global.PHYSICS.Pr_option == "PRZERO")) {
			
			if (global.PHYSICS.Uscaling == "USMALL") 
				U.Force1 = (global.PHYSICS.Rayleigh)* (T.csf.F);				// (u.grad)u-theta
			
			else if (global.PHYSICS.Uscaling == "ULARGE") 		
				U.Force1 = (global.PHYSICS.Prandtl)*(T.csf.F);				// (u.grad)u-theta
		}
		
		else if (global.PHYSICS.Pr_option == "PRINFTY") {
			;  // Force only Temperature; do nothing here
		}

		
		// For the temperature field

		if (global.PHYSICS.Pr_option == "PRZERO") 
			T.Force = 0.0;  
		// Do nothing; Nonlinear term (u.grad)T does not exist- see single-time step
		
		else  {
			if (global.PHYSICS.Pr_option == "PRLARGE")
			T.Force = (U.cvf.V1);
			
			// F(T) = globalvar_temperature_grad * ux(k)	
			// globalvar_temperature_grad = +1 for RB, and -1 for stratified flows
			
			else if (global.PHYSICS.Pr_option == "PRSMALL") {
				T.Force =   (1/global.PHYSICS.Prandtl)*(U.cvf.V1);
			}
			//F(T) =  globalvar_temperature_grad * ux(k)/Pr
			
			else if (global.PHYSICS.Pr_option == "PRINFTY")
					T.Force =  (U.cvf.V1);
			
			if (global.PHYSICS.temperature_grad == -1)  // for 
				T.Force = -T.Force;
		}
	}
}

void FORCE::Compute_force_RBC_basic_add(FluidVF& U, FluidSF& T)
{
	
	if (global.program.basis_type == "ChFF") { // box size (2,Ly,Lz)
        U.Force1 +=  0.125*(global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)*(T.csf.F);
        T.Force += 0.5*(U.cvf.V1);
    }
    
    else {
		// For the velocity field
		if (global.PHYSICS.Pr_option == "PRLARGE") {
			
			if (global.PHYSICS.Uscaling == "USMALL") 		
				U.Force1 +=  (global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)*(T.csf.F);	// (u.grad)u-Ra*Pr*theta
			
			else if (global.PHYSICS.Uscaling  == "ULARGE") 		
				U.Force1 +=  T.csf.F;					// (u.grad)u-theta			
		}
		
		else if ((global.PHYSICS.Pr_option == "PRSMALL") || (global.PHYSICS.Pr_option == "PRZERO")) {
			
			if (global.PHYSICS.Uscaling == "USMALL") 
				U.Force1 += (global.PHYSICS.Rayleigh)* (T.csf.F);				// (u.grad)u-theta
			
			else if (global.PHYSICS.Uscaling == "ULARGE") 		
				U.Force1 += (global.PHYSICS.Prandtl)*(T.csf.F);				// (u.grad)u-theta
		}
		
		else if (global.PHYSICS.Pr_option == "PRINFTY") {
			;  // Force only Temperature; do nothing here
		}

		
		// For the temperature field

		//if (global.PHYSICS.Pr_option == "PRZERO") 
		//	T.Force = 0.0;  // No point in adding Zero
		// Do nothing; Nonlinear term (u.grad)T does not exist- see single-time step
		
		if (global.PHYSICS.Pr_option != "PRZERO")   {
			if (global.PHYSICS.Pr_option == "PRLARGE")
			T.Force += (U.cvf.V1);
			
			// F(T) = globalvar_temperature_grad * ux(k)	
			// globalvar_temperature_grad = +1 for RB, and -1 for stratified flows
			
			else if (global.PHYSICS.Pr_option == "PRSMALL") {
				T.Force += (1/global.PHYSICS.Prandtl)*(U.cvf.V1);
			}
			//F(T) =  globalvar_temperature_grad * ux(k)/Pr
			
			else if (global.PHYSICS.Pr_option == "PRINFTY")
					T.Force +=  (U.cvf.V1);
			
			if (global.PHYSICS.temperature_grad == -1)  // for 
				T.Force += -T.Force;
		}
	}
			
}

void FORCE::Compute_force_RBC(FluidVF& U, FluidSF& T)
{
	Compute_force_RBC_basic_assign(U, T);
}


//*********************************************************************************************


void FORCE::Compute_force_RBC(FluidVF& U, FluidVF& W, FluidSF& T)
{
	Compute_force_RBC_basic_assign(U, T);
	
	if (W.force_switch) {
		W.Force1 = 0.0; 
		W.Force2 = 0.0; 
		W.Force3 = 0.0;	
	}
}


//*********************************************************************************************
//*********************************************************************************************


void FORCE::Compute_force_RBC_rotation(FluidVF& U, FluidSF& T)
{
	
	Compute_force_RBC_basic_assign(U, T);
	
	int omega_components = global.force.int_para(0);
	int rotation_direction = global.force.int_para(1);
	
	if (omega_components == 1) {
		int rotation_direction = global.force.int_para(1);
		DP two_omega =  2*global.force.double_para(0);
		
		Compute_force_Coriolis_basic_add(U, rotation_direction, two_omega);
	}
	
	else if (omega_components == 3) {
		DP two_omega1 = 2*global.force.double_para(0);
		DP two_omega2 = 2*global.force.double_para(1);
		DP two_omega3 = 2*global.force.double_para(2);
		
		Compute_force_Coriolis_basic_add(U, two_omega1, two_omega2, two_omega3);
	}
	
}

//*********************************************************************************************
void FORCE::Compute_force_RBC_rotation(FluidVF& U, FluidVF& W, FluidSF& T)
{
	Compute_force_RBC_rotation(U, T);
}


//*********************************************************************************************
void FORCE::Compute_force_stratified_random(FluidVF& U, FluidSF& T)
{
	DP inner_radius =  global.force.double_para(0);
	DP outer_radius = global.force.double_para(1);
	DP alpha_factor = global.force.double_para(3);
	
	Compute_force_RBC(U, T);
	
	U.Force1  = -U.Force1;
	universal->Subtract_Laplacian(U.dissipation_coefficient, U.cvf.V1, U.Force1); // Force1 = -F_buoyancy- nu*laplacian(V1)
	universal->Subtract_Laplacian(U.dissipation_coefficient, U.cvf.V2, U.Force2); // Force2 = - nu*laplacian(V2)
	universal->Subtract_Laplacian(U.dissipation_coefficient, U.cvf.V3, U.Force3); // Force3 = - nu*laplacian(V3)
	
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	DP modal_energy;

	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF") || (basis_type == "SFF"))
		ky_min = -ky_max;
	
	int lx, ly, lz;
	DP Kmag, FdotV, alpha_k, beta_k, sk;
	bool add_flag;
	
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = ky_min; ky <= ky_max; ky++)
			for (int kz = 0; kz <= kz_max; kz++) {
				
				if (universal->Probe_in_me(kx,ky,kz))  {
					lx = universal->Get_lx(kx);
					ly = universal->Get_ly(ky);
					lz = universal->Get_kz(kz);
					
					Kmag = universal->Kmagnitude(lx, ly, lz);
					if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
						
						modal_energy = U.cvf.Modal_energy(lx, ly, lz);
						if (modal_energy > MYEPS) {
							FdotV = mydot(U.Force1(ly,lz,lx), U.Force2(ly,lz,lx), U.Force3(ly,lz,lx), U.cvf.V1(ly,lz,lx), U.cvf.V2(ly,lz,lx), U.cvf.V3(ly,lz,lx));
							alpha_k = alpha_factor*FdotV/modal_energy;
							beta_k = 0.0;

							add_flag= false;
							Const_energy_supply_alpha_beta(U, lx, ly, lz, alpha_k, beta_k, add_flag); 
						}
					}
				}	//  of if (Probe_in_me())
			}		// of for 

	T.Force = 0.0;
	Compute_force_RBC_basic_add(U, T);
	
	
	//Force_energy_helicity_supply_or_level_basic_assign(U, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply, epsh_by_k_epse);
//	Compute_force_const_energy_helicity(U, T);
//	Compute_force_RBC_basic_add(U, T);
	/*U.Force1 = 0;
	U.Force2 = 0;
	U.Force3 = 0;*/
}

//************************ End of compute_force_RB.cc *****************************************

