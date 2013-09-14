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
	
	Array<DP,3> temp1;
	Array<complx,3> temp2;
	
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
	
	void Const_energy_supply_alpha_beta(FluidVF& U, int lx, int ly, int lz, DP alpha, DP beta, bool add_flag);
	void Const_energy_alpha_beta(FluidVF& U, int lx, int ly, int lz, DP alpha, DP beta, bool add_flag);
	
	void Const_energy_supply_alpha(FluidSF& T, int lx, int ly, int lz, DP alpha, bool add_flag);
	void Const_energy_alpha(FluidSF& T, int lx, int ly, int lz, DP alpha, bool add_flag);
	
	void Force_energy_helicity_supply_or_level_basic(FluidVF& U, string force_type, DP inner_radius, DP outer_radius, DP para1, DP para2, bool add_flag);
	void Force_energy_helicity_supply_or_level_basic(FluidSF& T, string force_type, DP inner_radius, DP outer_radius, DP para, bool add_flag);
	void Force_energy_helicity_supply_or_level_basic_assign(FluidVF& U, string force_type, DP inner_radius, DP outer_radius, DP para1, DP para2);
	void Force_energy_helicity_supply_or_level_basic_add(FluidVF& U, string force_type, DP inner_radius, DP outer_radius, DP para1, DP para2);
	void Force_energy_helicity_supply_or_level_basic_assign(FluidSF& T, string force_type, DP inner_radius, DP outer_radius, DP para);
	void Force_energy_helicity_supply_or_level_basic_add(FluidSF& T, string force_type, DP inner_radius, DP outer_radius, DP para);

	void Force_energy_helicity_supply_or_level(FluidVF& U, string force_type, DP inner_radius, DP outer_radius, DP para1, DP para2);
	void Force_energy_helicity_supply_or_level(FluidSF& T, string force_type, DP inner_radius, DP outer_radius, DP para);
	
	void Compute_force_const_energy_helicity_supply(FluidVF& U);
	void Compute_force_const_energy_helicity_supply(FluidVF& U, FluidSF& T);
	void Compute_force_const_energy_helicity_supply(FluidVF& U, FluidVF& W);
	void Compute_force_const_energy_helicity_supply(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_const_energy_helicity(FluidVF& U);			
	void Compute_force_const_energy_helicity(FluidVF& U, FluidSF& T);
	void Compute_force_const_energy_helicity(FluidVF& U, FluidVF& W);						
	void Compute_force_const_energy_helicity(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Add_complex_conj_force(FluidVF& U, int kx, int ky, int kz, complx Fx, complx Fy, complx Fz);
	void Add_complex_conj_force(FluidSF& T, int kx, int ky, int kz, complx localG);
	void Assign_force_and_comp_conj(FluidVF& U, int kx, int ky, int kz, complx Fx, complx Fy, complx Fz);
	void Assign_force_and_comp_conj(FluidSF& T, int kx, int ky, int kz, complx localG);
	void Assign_force(FluidVF& U, int kx, int ky, int kz, DP Fx, DP Fy, DP Fz);
	void Assign_force(FluidSF& T, int kx, int ky, int kz, DP localG);
	
	void Compute_force_given_modes(FluidVF& U);
	void Compute_force_given_modes(FluidVF& U, FluidSF& T);
	void Compute_force_given_modes(FluidVF& U, FluidVF& W);
	void Compute_force_given_modes(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Setup_Taylor_Green_force_field(FluidVF& U, int k0, DP amp);								
	void Setup_ABC_force_field(FluidVF& U, int k0, DP amp, DP A, DP B, DP C);
	void Setup_SIX_MODE_force_field(FluidVF& U, int k0, DP amp101, DP amp011, DP amp112, DP h);
	
	
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
	
	void Compute_force_Coriolis_basic_assign(FluidVF& U, int rotation_direction, DP two_omega);
	void Compute_force_Coriolis_basic_assign(FluidVF& U, DP two_omega1, DP two_omega2, DP two_omega3);

	void Compute_force_Coriolis_basic_add(FluidVF& U, int rotation_direction, DP two_omega);
	void Compute_force_Coriolis_basic_add(FluidVF& U, DP two_omega1, DP two_omega2, DP two_omega3);

	
	void Compute_force_Coriolis(FluidVF& U);
	void Compute_force_Coriolis(FluidVF& U, FluidSF& T);
	void Compute_force_Coriolis(FluidVF& U, FluidVF& W);
	void Compute_force_Coriolis(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_Ekman_friction(FluidVF& U);
	void Compute_force_Ekman_friction_basic_assign(FluidVF& U, DP alpha);
	void Compute_force_Ekman_friction_basic_add(FluidVF& U, DP alpha);
	void Compute_force_Ekman_friction_const_energy_supply(FluidVF& U);
	
	void Compute_force_Keplerian_basic_assign(FluidVF& U, DP omega_keplerian, DP q_keplerian);
	void Compute_force_Keplerian_basic_add(FluidVF& U, DP omega_keplerian, DP q_keplerian);
	void Compute_force_Keplerian(FluidVF& U);
	void Compute_force_Keplerian_SB(FluidVF& U);
	
	void Compute_force_Keplerian_basic_assign(FluidVF& U, FluidVF& W, DP omega_keplerian, DP q_keplerian);
	void Compute_force_Keplerian_basic_add(FluidVF& U, FluidVF& W, DP omega_keplerian, DP q_keplerian);
	void Compute_force_Keplerian(FluidVF& U, FluidVF& W);
	void Compute_force_Keplerian_SB(FluidVF& U, FluidVF& W);
	
	
	void Compute_force_Liquid_metal_basic_assign(FluidVF& U, DP B0x, DP B0y, DP B0z);
	void Compute_force_Liquid_metal_basic_add(FluidVF& U, DP B0x, DP B0y, DP B0z);
	void Compute_force_Liquid_metal(FluidVF& U);
	void Compute_force_Liquid_metal(FluidVF& U, FluidSF& T);
	void Compute_force_Liquid_metal_const_energy_supply(FluidVF& U);
	
	void Compute_force_Kolmogorov_flow_basic_assign(FluidVF& U, DP k0, DP force_amp, DP Rh);
	void Compute_force_Kolmogorov_flow_basic_add(FluidVF& U, DP k0, DP force_amp, DP Rh);
	void Compute_force_Kolmogorov_flow(FluidVF& U);

	void Compute_force_stratified_random(FluidVF& U, FluidSF& T);

	void Model_force_spectrum(DP force_spectrum_amplitude, DP force_spectrum_exponent, Array<DP,1> Sk);
	void Put_force_amp_phase_comp_conj(FluidVF U, int lx, int ly, int lz,  DP amp, DP phase1, DP phase2, DP phase3, bool add_flag);
	void Put_force_amp_phase_comp_conj(FluidSF& T, int lx, int ly, int lz, DP amp, DP phase, bool add_flag);
	
	void Compute_force_using_random_energy_helicity_spectrum_basic(FluidVF& U, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, DP hk_by_kek, bool add_flag);
	void Compute_force_using_random_energy_helicity_spectrum_basic_assign(FluidVF& U, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, DP hk_by_kek);
	void Compute_force_using_random_energy_helicity_spectrum_basic_add(FluidVF& U, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, DP hk_by_kek);
	void Compute_force_using_random_energy_spectrum_basic(FluidSF& T, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, bool add_flag);
	void Compute_force_using_random_energy_spectrum_basic_assign(FluidSF& T, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent);
	void Compute_force_using_random_energy_spectrum_basic_add(FluidSF& T, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent);
	
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U);
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidSF& T);
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidVF& W);
	void  Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
							
};

#endif

//*******************************  End of force_model.cc **************************************

	
