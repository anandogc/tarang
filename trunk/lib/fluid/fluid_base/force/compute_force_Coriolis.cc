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


/*! \file  compute_force_TG.cc
 * 
 * @brief  Set up force using Taylor Green theory.
 *
 * @note Fx = F*sin(k0 x) cos(k0 y) cos(k0 z) <BR>
 *		Fy = -F*cos(k0 x)sin(k0 y) cos(k0 z) <BR>
 *		Fz = 0
 *
 *  @note is_force_field_para_read is used to read only once.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "FORCE.h"

//*********************************************************************************************

// omega = 2 omega(orig) L/nu
void FORCE::Compute_force_Coriolis_basic_assign(FluidVF& U, int rotation_direction, Real two_omega)
{
	if (rotation_direction == 1) {	// omega along x
		U.Force2 = two_omega*(U.cvf.V3);
		U.Force3 = -two_omega*(U.cvf.V2);
	}
	else if (rotation_direction == 2) {// omega along y
		U.Force1 = -two_omega*(U.cvf.V3);
		U.Force3 = two_omega*(U.cvf.V1);
	}
	else if (rotation_direction == 3) {// omega along z
	   U.Force1 = two_omega*(U.cvf.V2);
		U.Force2 = -two_omega*(U.cvf.V1);
}
}
void FORCE::Compute_force_Coriolis_basic_add(FluidVF& U, int rotation_direction, Real two_omega)
{
	if (rotation_direction == 1) {	// omega along x
		U.Force2 += two_omega*(U.cvf.V3);
		U.Force3 += -two_omega*(U.cvf.V2);
	}
	else if (rotation_direction == 2) {// omega along y
		U.Force1 += -two_omega*(U.cvf.V3);
		U.Force3 += two_omega*(U.cvf.V1);
	}
	else if (rotation_direction == 3) {// omega along z
		U.Force1 += two_omega*(U.cvf.V2);
		U.Force2 += -two_omega*(U.cvf.V1);
	}	
}


void FORCE::Compute_force_Coriolis_basic_assign(FluidVF& U, Real two_omega1, Real two_omega2, Real two_omega3)
{
	U.Force1 = -two_omega2*(U.cvf.V3) + two_omega3*(U.cvf.V2);
	U.Force2 = -two_omega3*(U.cvf.V1) + two_omega1*(U.cvf.V3);
	U.Force3 = -two_omega1*(U.cvf.V2) + two_omega2*(U.cvf.V1);

}

void FORCE::Compute_force_Coriolis_basic_add(FluidVF& U, Real two_omega1, Real two_omega2, Real two_omega3)
{
	U.Force1 += -two_omega2*(U.cvf.V3) + two_omega3*(U.cvf.V2);
	U.Force2 += -two_omega3*(U.cvf.V1) + two_omega1*(U.cvf.V3);
	U.Force3 += -two_omega1*(U.cvf.V2) + two_omega2*(U.cvf.V1);
}


// omega = 2 omega(orig) L/nu
// derived fn
void FORCE::Compute_force_Coriolis(FluidVF& U)
{
	U.Force1 = 0.0;
	U.Force2 = 0.0;
	U.Force3 = 0.0;
	
	if (U.force_switch) {
		int omega_components = global.force.int_para(0);
		
		if (omega_components == 1) {
			int rotation_direction = global.force.int_para(1);
			Real two_omega =  2*global.force.double_para(0);
			
			Compute_force_Coriolis_basic_assign(U, rotation_direction, two_omega);
		}
		
		else if (omega_components == 3) {
			Real two_omega1 = 2*global.force.double_para(0);
			Real two_omega2 = 2*global.force.double_para(1);
			Real two_omega3 = 2*global.force.double_para(2);
			
			Compute_force_Coriolis_basic_assign(U, two_omega1, two_omega2, two_omega3);
		}
	}
}

// derived fn
void FORCE::Compute_force_Coriolis(FluidVF& U, FluidSF& T) 
{
	Compute_force_Coriolis(U);
}

// derived fn
void FORCE::Compute_force_Coriolis(FluidVF& U, FluidVF& W) 
{
	Compute_force_Coriolis(U);
}

// derived fn
void FORCE::Compute_force_Coriolis(FluidVF& U, FluidVF& W, FluidSF& T) 
{
	Compute_force_Coriolis(U);
}


//******************************************************************************


void FORCE::Compute_force_Ekman_friction_basic_assign(FluidVF& U, Real alpha)
{		
	U.Force1 = Complex(-alpha,0) * (U.cvf.V1);
	
	if (N[2] > 1)	
		U.Force2 = Complex(-alpha,0) * (U.cvf.V2);
	
	U.Force3 = Complex(-alpha,0) * (U.cvf.V3);
}

void FORCE::Compute_force_Ekman_friction_basic_add(FluidVF& U, Real alpha)
{
	U.Force1 += Complex(-alpha,0) * (U.cvf.V1);
	
	if (N[2] > 1)	
		U.Force2 += Complex(-alpha,0) * (U.cvf.V2);
	
	U.Force3 += Complex(-alpha,0) * (U.cvf.V3);
}

// derived fn
void FORCE::Compute_force_Ekman_friction(FluidVF& U)
{
	Real alpha = global.force.double_para(0);
	
	U.Force1 = 0.0;
	U.Force2 = 0.0;
	U.Force3 = 0.0;
	
	Compute_force_Ekman_friction_basic_assign(U, alpha);
}


// derived fn
void FORCE::Compute_force_Ekman_friction_const_energy_supply(FluidVF& U)
{	
	
	Real inner_radius =  global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_supply = global.force.double_para(2);
	Real epsh_by_k_epse = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	
	Real alpha = global.force.double_para(4);
	
	U.Force1 = 0.0;
	U.Force2 = 0.0;
	U.Force3 = 0.0;
	
		// first feed const eps force;
	Force_energy_helicity_supply_or_level_basic_assign(U, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply, epsh_by_k_epse);
	
	Compute_force_Ekman_friction_basic_add(U, alpha);
}


//******************************************************************************


void FORCE::Compute_force_Keplerian_basic_assign(FluidVF& U, Real omega_keplerian, Real q_keplerian)
{
	
	if (U.force_switch == 1) {
		U.Force1 = 2*omega_keplerian*(U.cvf.V2);
		U.Force2 = (q_keplerian-2)*omega_keplerian*(U.cvf.V1);
		U.Force3 = 0.0;
		
		universal->Yderiv(U.cvf.V1, global.temp_array.X);
		U.Force1 += global.temp_array.X;
		
		universal->Yderiv(U.cvf.V2, global.temp_array.X);
		U.Force2 += global.temp_array.X;
		
		universal->Yderiv(U.cvf.V3, global.temp_array.X);
		U.Force3 += global.temp_array.X;
	}
}

void FORCE::Compute_force_Keplerian_basic_add(FluidVF& U, Real omega_keplerian, Real q_keplerian)
{
	
	if (U.force_switch == 1) {
		U.Force1 += 2*omega_keplerian*(U.cvf.V2);
		U.Force2 += (q_keplerian-2)*omega_keplerian*(U.cvf.V1);
		//	U.Force3 = 0.0;
		
		universal->Yderiv(U.cvf.V1, global.temp_array.X);
		U.Force1 += global.temp_array.X;
		
		universal->Yderiv(U.cvf.V2, global.temp_array.X);
		U.Force2 += global.temp_array.X;
		
		universal->Yderiv(U.cvf.V3, global.temp_array.X);
		U.Force3 += global.temp_array.X;
	}
}

// derived fn
void FORCE::Compute_force_Keplerian(FluidVF& U)
{
	
	Real omega_keplerian = global.force.double_para(0);
	Real q_keplerian = global.force.double_para(1);
	
	Compute_force_Keplerian_basic_assign(U, omega_keplerian, q_keplerian);
	
}

void FORCE::Compute_force_Keplerian_basic_assign(FluidVF& U, FluidVF& W, Real omega_keplerian, Real q_keplerian)
{
	
	if (U.force_switch == 1) {
		U.Force1 = 2*omega_keplerian*(U.cvf.V2);
		U.Force2 = (q_keplerian-2)*omega_keplerian*(U.cvf.V1);
		U.Force3 = 0.0;
	}
	
	if (W.force_switch == 1) {
		W.Force1 = 0.0;
		W.Force2 = -q_keplerian*omega_keplerian*(W.cvf.V1);
		W.Force3 = 0.0;
	};
}

void FORCE::Compute_force_Keplerian_basic_add(FluidVF& U, FluidVF& W, Real omega_keplerian, Real q_keplerian)
{
	
	if (U.force_switch == 1) {
		U.Force1 += 2*omega_keplerian*(U.cvf.V2);
		U.Force2 += (q_keplerian-2)*omega_keplerian*(U.cvf.V1);
	//	U.Force3 = 0.0;
	}
	
	if (W.force_switch == 1) {
	//	W.Force1 = 0.0;
		W.Force2 += -q_keplerian*omega_keplerian*(W.cvf.V1);
	//	W.Force3 = 0.0;
	};
}

// derived fn
void FORCE::Compute_force_Keplerian(FluidVF& U, FluidVF& W)
{
	
	Real omega_keplerian = global.force.double_para(0);
	Real q_keplerian = global.force.double_para(1);
	
	Compute_force_Keplerian_basic_assign(U, W, omega_keplerian, q_keplerian);
	
}

// Sujit_Banibrata -- add random noise
void FORCE::Compute_force_Keplerian_SB(FluidVF& U)
{
	
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real force_spectrum_amplitude = global.force.double_para(2);
	Real force_spectrum_exponent = global.force.double_para(3);
	Real hk_by_kek = global.force.double_para(4);
	// Real omega_keplerian = global.force.double_para(5);
	Real q_keplerian = global.force.double_para(6);
	Real x0 = global.force.double_para(7);
	
	Real omega_keplerian = 1/q_keplerian;
	
	
	Compute_force_using_random_energy_helicity_spectrum_basic_assign(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek);
	
	Compute_force_Keplerian_basic_assign(U, omega_keplerian, q_keplerian);
	
	universal->Yderiv(U.cvf.V1, global.temp_array.X);
	U.Force1 += x0*global.temp_array.X;
	
	universal->Yderiv(U.cvf.V2, global.temp_array.X);
	U.Force2 += x0*global.temp_array.X;
	
	universal->Yderiv(U.cvf.V3, global.temp_array.X);
	U.Force3 += x0*global.temp_array.X;
	
}

// Sujit_Banibrata -- add random noise
void FORCE::Compute_force_Keplerian_SB(FluidVF& U, FluidVF& W)
{
	
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real force_spectrum_amplitude = global.force.double_para(2);
	Real force_spectrum_exponent = global.force.double_para(3);
	Real hk_by_kek = global.force.double_para(4);
	Real Wforce_spectrum_amplitude = global.force.double_para(5);
	Real Wforce_spectrum_exponent = global.force.double_para(6);
	Real Whk_by_kek = global.force.double_para(7);
	//Real omega_keplerian = global.force.double_para(8);
	Real q_keplerian = global.force.double_para(9);
	Real omega_keplerian = 1/q_keplerian;
	
	Compute_force_using_random_energy_helicity_spectrum_basic_assign(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek);
	
	Compute_force_using_random_energy_helicity_spectrum_basic_assign(W, inner_radius, outer_radius, Wforce_spectrum_amplitude, Wforce_spectrum_exponent, Whk_by_kek);
	
	Compute_force_Keplerian_basic_add(U, W, omega_keplerian, q_keplerian);
	
	universal->Yderiv(U.cvf.V1, global.temp_array.X);
	U.Force1 += global.temp_array.X;
	
	universal->Yderiv(U.cvf.V2, global.temp_array.X);
	U.Force2 += global.temp_array.X;
	
	universal->Yderiv(U.cvf.V3, global.temp_array.X);
	U.Force3 += global.temp_array.X;
	
	universal->Yderiv(W.cvf.V1, global.temp_array.X);
	W.Force1 += global.temp_array.X;
	
	universal->Yderiv(W.cvf.V2, global.temp_array.X);
	W.Force2 += global.temp_array.X;
	
	universal->Yderiv(W.cvf.V3, global.temp_array.X);
	W.Force3 += global.temp_array.X;

}


//****************************** End of compute_force_TG.cc ***********************************


