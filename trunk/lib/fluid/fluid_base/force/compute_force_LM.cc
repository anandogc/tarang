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

/*! \file  compute_force_LM.cc
 * 
 * @brief  Set up force using Liquid metal for under strong mean magnetic field.
 *
 * @note F =  (hat(B0).hat(K))^2 V
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "FORCE.h"


//*********************************************************************************************


void FORCE::Compute_force_Liquid_metal_basic_assign(FluidVF& U, Real B0x, Real B0y, Real B0z)
{	
	
	TinyVector <Real,3> B0;
	B0 = B0x, B0y, B0z;
	
	U.Force1 = U.cvf.V1;
	U.Force2 = U.cvf.V2; 
	U.Force3 = U.cvf.V3;
	
	universal->Array_mult_V0_khat_sqr(U.Force1, B0);
	universal->Array_mult_V0_khat_sqr(U.Force2, B0);
	universal->Array_mult_V0_khat_sqr(U.Force3, B0);
	
	U.Force1 = Complex(-1,0) * (U.Force1);
	U.Force2 = Complex(-1,0) * (U.Force2);
	U.Force3 = Complex(-1,0) * (U.Force3); 
}


void FORCE::Compute_force_Liquid_metal_basic_add(FluidVF& U, Real B0x, Real B0y, Real B0z)
{	
	
	TinyVector <Real,3> B0;
	B0 = B0x, B0y, B0z;
		
	global.temp_array.X = U.cvf.V1;
	universal->Array_mult_V0_khat_sqr(global.temp_array.X, B0);
	U.Force1 += Complex(-1,0) * (global.temp_array.X);
	
	global.temp_array.X = U.cvf.V2;
	universal->Array_mult_V0_khat_sqr(global.temp_array.X, B0);
	U.Force2 += Complex(-1,0) * (global.temp_array.X);
	
	global.temp_array.X = U.cvf.V3;
	universal->Array_mult_V0_khat_sqr(global.temp_array.X, B0);
	U.Force3 += Complex(-1,0) * (global.temp_array.X);
}
void FORCE::Compute_force_Liquid_metal(FluidVF& U)
{
	
	Real B0x = global.force.double_para(0);
	Real B0y = global.force.double_para(1);
	Real B0z = global.force.double_para(2);
	
	Compute_force_Liquid_metal_basic_assign(U, B0x, B0y, B0z);
}


void FORCE::Compute_force_Liquid_metal_const_energy_supply(FluidVF& U)
{	
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_supply = global.force.double_para(2);
	Real epsh_by_k_epse = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	
	Real B0x = global.force.double_para(4);
	Real B0y = global.force.double_para(5);
	Real B0z = global.force.double_para(6);
	
	// first feed const eps force;
	Force_energy_helicity_supply_or_level_basic_assign(U, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply, epsh_by_k_epse); 
	
	// LM force
	
	Compute_force_Liquid_metal_basic_add(U, B0x, B0y, B0z);
}	

//*********************************************************************************************

void FORCE::Compute_force_Liquid_metal(FluidVF& U, FluidSF& T)
{
	Real B0x = global.force.double_para(0);
	Real B0y = global.force.double_para(1);
	Real B0z = global.force.double_para(2);
	
	Compute_force_RBC_basic_assign(U, T);
	
	Compute_force_Liquid_metal_basic_add(U, B0x, B0y, B0z);
}


//*********************************************************************************************

// Eqn 2D NS with additional forcing: -u/Rh

/** @brief Kolmogorov flow forcing
 * 
 * @note Force_x = amp*sin(k0 x) cos(k0 y) cos(k0 z) - V1/Rh
 *		Force_y = 0
 *		Force_z = -amp*cos(k0 x)sin(k0 y) cos(k0 z) -V3/Rh
 *
 *  @param k0
 *	@param amp  
 *
 */
void FORCE::Compute_force_Kolmogorov_flow_basic_assign(FluidVF& U, Real k0, Real force_amp, Real Rh)
{	
		
	// -u/Rh
	U.Force1 = Complex(-1/Rh,0) * (U.cvf.V1);
	U.Force2 = 0.0; 
	U.Force3 = Complex(-1/Rh,0) * (U.cvf.V3);
	
	// CHOOSE SCC type
	if (basis_type == "SSS") {
		universal->Assign_spectral_field(k0, 0, k0, U.Force1, force_amp/4);
		universal->Assign_spectral_field(k0, 0, k0, U.Force3, -force_amp/4);
	}
	
	// CHOOSE SCF type
	else if (basis_type == "SSF") {
		cout << "MYERROR:  In Compute_force_Kolmogorov_flow(), SSF basis not allowed " << endl;
		exit(1);
	}
	
	// CHOOSE SFF type
	else if (basis_type == "SFF") {
		
		universal->Assign_spectral_field(k0, 0, k0, U.Force1, Complex(force_amp/2, 0));
		universal->Assign_spectral_field(k0, 0, k0, U.Force3, Complex(-force_amp/2, 0));
	}
	
	else if (basis_type == "FFF" || basis_type == "FFFW") {
		universal->Assign_spectral_field(k0, 0, k0, U.Force1, Complex(0, -force_amp/4));
		universal->Assign_spectral_field(-k0, 0, k0, U.Force1, Complex(0, force_amp/4));
		
		universal->Assign_spectral_field(k0, 0, k0, U.Force3, Complex(0, force_amp/4));
		universal->Assign_spectral_field(-k0, 0, k0, U.Force3, Complex(0, force_amp/4));
	}
}

void FORCE::Compute_force_Kolmogorov_flow_basic_add(FluidVF& U, Real k0, Real force_amp, Real Rh)
{	
		
	// -u/Rh
	U.Force1 += Complex(-1/Rh,0) * (U.cvf.V1);
	//U.Force2 += 0.0; 
	U.Force3 += Complex(-1/Rh,0) * (U.cvf.V3);
	
	// CHOOSE SCC type
	if (basis_type == "SSS") {
		universal->Add_spectral_field(k0, 0, k0, U.Force1, force_amp/4);
		universal->Add_spectral_field(k0, 0, k0, U.Force3, -force_amp/4);
	}
	
	// CHOOSE SCF type
	else if (basis_type == "SSF") {
		cout << "MYERROR:  In Compute_force_Kolmogorov_flow(), SSF basis not allowed " << endl;
		exit(1);
	}
	
	// CHOOSE SFF type
	else if (basis_type == "SFF") {
		
		universal->Add_spectral_field(k0, 0, k0, U.Force1, Complex(force_amp/2, 0));
		universal->Add_spectral_field(k0, 0, k0, U.Force3, Complex(-force_amp/2, 0));
	}
	
	else if (basis_type == "FFF" || basis_type == "FFFW") {
		universal->Add_spectral_field(k0, 0, k0, U.Force1, Complex(0, -force_amp/4));
		universal->Add_spectral_field(-k0, 0, k0, U.Force1, Complex(0, force_amp/4));
		
		universal->Add_spectral_field(k0, 0, k0, U.Force3, Complex(0, force_amp/4));
		universal->Add_spectral_field(-k0, 0, k0, U.Force3, Complex(0, force_amp/4));
	}
}


void FORCE::Compute_force_Kolmogorov_flow(FluidVF& U)
{
	int k0 = ((int) global.force.double_para(0));
	Real force_amp = global.force.double_para(1);
	Real Rh = global.force.double_para(2);
	
	Compute_force_Kolmogorov_flow_basic_assign(U, k0, force_amp, Rh);
	
}



//****************************** End of compute_force_LM.cc ***********************************


