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

/*! \file  force_model.cc
 * 
 * @brief Forcing when low k modes are forced: TG, ABC flows etc.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */



#ifndef _H_Force
#define _H_Force

#include "fluid_base.h"

//*********************************************************************************************

class FORCE
{
public: 
	ifstream force_field_in_file;
	
	Array<Real,3> temp1;
	Array<Complex,3> temp2;
	
public:
	void Compute_force(FluidVF& U);
	void Compute_force(FluidVF& U, FluidSF& T);
	void Compute_force(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Compute_force(FluidVF& U, FluidVF& W);
	void Compute_force(FluidVF& U, FluidVF& W, FluidSF& T);
	void Compute_force(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	
	void Compute_force_decay(FluidVF& U);
	void Compute_force_decay(FluidVF& U, FluidSF& T);
	void Compute_force_decay(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Compute_force_decay(FluidVF& U, FluidVF& W);
	void Compute_force_decay(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force(FluidSF& T);
	
	void Compute_force_pressure_grad(FluidVF& U);
	
	void Const_energy_supply_alpha_beta(FluidVF& U, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag);
	void Const_energy_alpha_beta(FluidVF& U, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag);
	
	void Const_energy_supply_alpha(FluidSF& T, int lx, int ly, int lz, Real alpha, bool add_flag);
	void Const_energy_alpha(FluidSF& T, int lx, int ly, int lz, Real alpha, bool add_flag);
	
	void Force_energy_helicity_supply_or_level_basic(FluidVF& U, string force_type, Real inner_radius, Real outer_radius, Real para1, Real para2, bool add_flag);
	void Force_energy_helicity_supply_or_level_basic(FluidSF& T, string force_type, Real inner_radius, Real outer_radius, Real para, bool add_flag);
	void Force_energy_helicity_supply_or_level_basic_assign(FluidVF& U, string force_type, Real inner_radius, Real outer_radius, Real para1, Real para2);
	void Force_energy_helicity_supply_or_level_basic_add(FluidVF& U, string force_type, Real inner_radius, Real outer_radius, Real para1, Real para2);
	void Force_energy_helicity_supply_or_level_basic_assign(FluidSF& T, string force_type, Real inner_radius, Real outer_radius, Real para);
	void Force_energy_helicity_supply_or_level_basic_add(FluidSF& T, string force_type, Real inner_radius, Real outer_radius, Real para);

	void Force_energy_helicity_supply_or_level(FluidVF& U, string force_type, Real inner_radius, Real outer_radius, Real para1, Real para2);
	void Force_energy_helicity_supply_or_level(FluidSF& T, string force_type, Real inner_radius, Real outer_radius, Real para);
	
	void Compute_force_const_energy_helicity_supply(FluidVF& U);
	void Compute_force_const_energy_helicity_supply(FluidVF& U, FluidSF& T);
	void Compute_force_const_energy_helicity_supply(FluidVF& U, FluidVF& W);
	void Compute_force_const_energy_helicity_supply(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_const_energy_helicity(FluidVF& U);			
	void Compute_force_const_energy_helicity(FluidVF& U, FluidSF& T);
	void Compute_force_const_energy_helicity(FluidVF& U, FluidVF& W);						
	void Compute_force_const_energy_helicity(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Add_complex_conj_force(FluidVF& U, int kx, int ky, int kz, Complex Fx, Complex Fy, Complex Fz);
	void Add_complex_conj_force(FluidSF& T, int kx, int ky, int kz, Complex localG);
	void Assign_force_and_comp_conj(FluidVF& U, int kx, int ky, int kz, Complex Fx, Complex Fy, Complex Fz);
	void Assign_force_and_comp_conj(FluidSF& T, int kx, int ky, int kz, Complex localG);
	void Assign_force(FluidVF& U, int kx, int ky, int kz, Real Fx, Real Fy, Real Fz);
	void Assign_force(FluidSF& T, int kx, int ky, int kz, Real localG);
	
	void Compute_force_given_modes(FluidVF& U);
	void Compute_force_given_modes(FluidVF& U, FluidSF& T);
	void Compute_force_given_modes(FluidVF& U, FluidVF& W);
	void Compute_force_given_modes(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Setup_Taylor_Green_force_field(FluidVF& U, int k0, Real amp);								
	void Setup_ABC_force_field(FluidVF& U, int k0, Real amp, Real A, Real B, Real C);
	void Setup_SIX_MODE_force_field(FluidVF& U, int k0, Real amp101, Real amp011, Real amp112, Real h);
	
	
	void Compute_force_Taylor_Green(FluidVF& U);
	void Compute_force_Taylor_Green(FluidVF& U, FluidSF& T);
	void Compute_force_Taylor_Green(FluidVF& U, FluidVF& W);
	void Compute_force_Taylor_Green(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_ABC(FluidVF& U);
	void Compute_force_ABC(FluidVF& U, FluidSF& T);
	void Compute_force_ABC(FluidVF& U, FluidVF& W);
	void Compute_force_ABC(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_RBC_basic_assign(FluidVF& U, FluidSF& T);
	void Compute_force_RBC_basic_add(FluidVF& U, FluidSF& T);
	void Compute_force_RBC(FluidVF& U, FluidSF& T);
	void Compute_force_RBC(FluidVF& U, FluidVF& W, FluidSF& T);
	void Compute_force_RBC_rotation(FluidVF& U, FluidSF& T);
	void Compute_force_RBC_rotation(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_astro(FluidVF& U, FluidVF& W, FluidSF& T);
	void Compute_force_astro(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	
	void Compute_force_MRBC(FluidVF& U, FluidSF& T1, FluidSF& T2);
	
	void Compute_force_DYNAMO_SIX_MODE(FluidVF& U, FluidVF& W);
	
	void Compute_force_Coriolis_basic_assign(FluidVF& U, int rotation_direction, Real two_omega);
	void Compute_force_Coriolis_basic_assign(FluidVF& U, Real two_omega1, Real two_omega2, Real two_omega3);

	void Compute_force_Coriolis_basic_add(FluidVF& U, int rotation_direction, Real two_omega);
	void Compute_force_Coriolis_basic_add(FluidVF& U, Real two_omega1, Real two_omega2, Real two_omega3);

	
	void Compute_force_Coriolis(FluidVF& U);
	void Compute_force_Coriolis(FluidVF& U, FluidSF& T);
	void Compute_force_Coriolis(FluidVF& U, FluidVF& W);
	void Compute_force_Coriolis(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_Ekman_friction(FluidVF& U);
	void Compute_force_Ekman_friction_basic_assign(FluidVF& U, Real alpha);
	void Compute_force_Ekman_friction_basic_add(FluidVF& U, Real alpha);
	void Compute_force_Ekman_friction_const_energy_supply(FluidVF& U);
	
	void Compute_force_Keplerian_basic_assign(FluidVF& U, Real omega_keplerian, Real q_keplerian);
	void Compute_force_Keplerian_basic_add(FluidVF& U, Real omega_keplerian, Real q_keplerian);
	void Compute_force_Keplerian(FluidVF& U);
	void Compute_force_Keplerian_SB(FluidVF& U);
	
	void Compute_force_Keplerian_basic_assign(FluidVF& U, FluidVF& W, Real omega_keplerian, Real q_keplerian);
	void Compute_force_Keplerian_basic_add(FluidVF& U, FluidVF& W, Real omega_keplerian, Real q_keplerian);
	void Compute_force_Keplerian(FluidVF& U, FluidVF& W);
	void Compute_force_Keplerian_SB(FluidVF& U, FluidVF& W);
	
	
	void Compute_force_Liquid_metal_basic_assign(FluidVF& U, Real B0x, Real B0y, Real B0z);
	void Compute_force_Liquid_metal_basic_add(FluidVF& U, Real B0x, Real B0y, Real B0z);
	void Compute_force_Liquid_metal(FluidVF& U);
	void Compute_force_Liquid_metal(FluidVF& U, FluidSF& T);
	void Compute_force_Liquid_metal_const_energy_supply(FluidVF& U);
	
	void Compute_force_Kolmogorov_flow_basic_assign(FluidVF& U, Real k0, Real force_amp, Real Rh);
	void Compute_force_Kolmogorov_flow_basic_add(FluidVF& U, Real k0, Real force_amp, Real Rh);
	void Compute_force_Kolmogorov_flow(FluidVF& U);

	void Compute_force_stratified_random(FluidVF& U, FluidSF& T);

	void Model_force_spectrum(Real force_spectrum_amplitude, Real force_spectrum_exponent, Array<Real,1> Sk);
	void Put_force_amp_phase_comp_conj(FluidVF U, int lx, int ly, int lz,  Real amp, Real phase1, Real phase2, Real phase3, bool add_flag);
	void Put_force_amp_phase_comp_conj(FluidSF& T, int lx, int ly, int lz, Real amp, Real phase, bool add_flag);
	
	void Compute_force_using_random_energy_helicity_spectrum_basic(FluidVF& U, Real inner_radius, Real outer_radius, Real force_spectrum_amplitude, Real force_spectrum_exponent, Real hk_by_kek, bool add_flag);
	void Compute_force_using_random_energy_helicity_spectrum_basic_assign(FluidVF& U, Real inner_radius, Real outer_radius, Real force_spectrum_amplitude, Real force_spectrum_exponent, Real hk_by_kek);
	void Compute_force_using_random_energy_helicity_spectrum_basic_add(FluidVF& U, Real inner_radius, Real outer_radius, Real force_spectrum_amplitude, Real force_spectrum_exponent, Real hk_by_kek);
	void Compute_force_using_random_energy_spectrum_basic(FluidSF& T, Real inner_radius, Real outer_radius, Real force_spectrum_amplitude, Real force_spectrum_exponent, bool add_flag);
	void Compute_force_using_random_energy_spectrum_basic_assign(FluidSF& T, Real inner_radius, Real outer_radius, Real force_spectrum_amplitude, Real force_spectrum_exponent);
	void Compute_force_using_random_energy_spectrum_basic_add(FluidSF& T, Real inner_radius, Real outer_radius, Real force_spectrum_amplitude, Real force_spectrum_exponent);
	
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U);
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidSF& T);
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidVF& W);
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
							
};

#endif

//*******************************  End of force_model.cc **************************************

	
