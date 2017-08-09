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
    
    // For the shells of energy_supply rate, helicity_supply rate etc.
    Real inner_radius, outer_radius;
    
    // arrays temp1, temp2 used in compute_force_astro.cc
	Array<Real,3> temp1;
	Array<Complex,3> temp2;
	
public:
 //   void FORCE();
    void Initialize();
	void Compute_force(FluidVF& U);
	void Compute_force(FluidSF& T);
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
	
	
	
	void Compute_force_pressure_grad(FluidVF& U);
	
	void Const_energy_supply_alpha_beta(FluidVF& U, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag);
	void Const_energy_supply_alpha(FluidSF& T, int lx, int ly, int lz, Real alpha, bool add_flag);
    
    void Const_energy_alpha_beta(FluidVF& U, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag);
	void Const_energy_alpha(FluidSF& T, int lx, int ly, int lz, Real alpha, bool add_flag);
	
    
    void Compute_force_Carati_scheme_basic(FluidVF& U, string force_type, bool global_alpha_beta, bool add_flag);
    
    void Compute_force_Carati_scheme_energy_supply(FluidVF& U, bool global_alpha_beta, bool add_flag);
    
    void Compute_force_Carati_scheme_const_energy(FluidVF& U, bool global_alpha_beta, bool add_flag);
    
    

    void Compute_force_Carati_scheme_basic(FluidSF& T, string force_type, bool global_alpha_beta, bool add_flag);
    
    void Compute_force_Carati_scheme(FluidVF& U, string force_type,  bool global_alpha_beta);

    void Compute_force_Carati_scheme(FluidVF& U, FluidSF& T, string force_type, bool global_alpha_beta);
    
    void Compute_force_Carati_scheme(FluidVF& U, FluidVF& W,  string force_type, bool global_alpha_beta);

    
    void Compute_force_Carati_scheme(FluidVF& U);
    
    void Compute_force_Carati_scheme(FluidVF& U, FluidSF& T);
    
    void Compute_force_Carati_scheme(FluidVF& U, FluidVF& W);
    
    void Compute_force_Carati_scheme(FluidVF& U, FluidVF& W, FluidSF& T);

    
    void Compute_force_Carati_scheme_crosshelicity_basic(FluidVF& U, FluidVF& W,  bool global_alpha_beta, bool add_flag);
	
    void Compute_force_Carati_scheme_crosshelicity(FluidVF& U, FluidVF& W);

	
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
	
	
	void Compute_force_RBC_basic(FluidVF& U, FluidSF& T);
	void Compute_force_RBC(FluidVF& U, FluidSF& T);
	void Compute_force_RBC(FluidVF& U, FluidVF& W, FluidSF& T);
	void Compute_force_RBC_rotation(FluidVF& U, FluidSF& T);
	void Compute_force_RBC_rotation(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Compute_force_astro(FluidVF& U, FluidVF& W, FluidSF& T);
	void Compute_force_astro(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	
	void Compute_force_MRBC(FluidVF& U, FluidSF& T1, FluidSF& T2);
	
	void Compute_force_DYNAMO_SIX_MODE(FluidVF& U, FluidVF& W);
	
	

	
	void Compute_force_Coriolis(FluidVF& U);
	void Compute_force_Coriolis(FluidVF& U, Real two_omega1, Real two_omega2, Real two_omega3);
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
    void Model_dissipation_spectrum(Real dissipation_spectrum_amplitude, Real dissipation_spectrum_exponent, Array<Real,1> Sk);
	
    void Put_or_add_force_vector(FluidVF U, int lx, int lz,  Real amp, Real phase, bool add_flag); // for 2D

    void Put_or_add_force_vector(FluidVF U, int lx, int ly, int lz,  Real amp_plus, Real amp_minus, Real phase_plus, Real phase_minus, bool add_flag);
    
    void Put_or_add_force_scalar(FluidSF& T, int lx, int ly, int lz, Real amp, Real phase, bool add_flag);
    
    
    
    void Compute_force_using_random_noise_basic(FluidVF& U,  Real random_interval, Uniform<Real> rand_struct, int rand_seed, bool add_flag);
    
    void Compute_force_using_random_noise_assign(FluidVF& U, Real random_interval);
    
    void Compute_force_using_random_noise_add(FluidVF& U,  Real random_interval);
    
    void Compute_force_using_random_noise_basic(FluidSF& T,  Real random_interval, Uniform<Real> rand_struct, int rand_seed, bool add_flag);
    
    void Compute_force_using_random_noise_assign(FluidVF& U, FluidSF& T, Real random_interval);
    
    void Compute_force_using_random_noise_add(FluidVF& U, FluidSF& T,  Real random_interval);

    
    void Compute_force_using_random_noise_assign(FluidVF& U, FluidVF& W,  Real random_interval);
    
    void  Compute_force_using_random_noise_add(FluidVF& U, FluidVF& W,  Real random_interval);
	
    void Compute_force_using_random_noise_assign(FluidVF& U, FluidVF& W, FluidSF& T, Real random_interval);
    
    void Compute_force_using_random_noise_add(FluidVF& U, FluidVF& W, FluidSF& T, Real random_interval);
    
	void  Compute_force_using_random_noise(FluidVF& U);
	void  Compute_force_using_random_noise(FluidVF& U, FluidSF& T);
	void  Compute_force_using_random_noise(FluidVF& U, FluidVF& W);
    void  Compute_force_using_random_noise(FluidVF& U, FluidVF& W, FluidSF& T);
    
    void Compute_force_user_defined1(FluidVF& U);
    void Compute_force_user_defined2(FluidVF& U);
    void Compute_force_user_defined1(FluidVF& U, FluidSF& T);
    void Compute_force_user_defined2(FluidVF& U, FluidSF& T);
    void Compute_force_user_defined1(FluidVF& U, FluidVF& W);
    void Compute_force_user_defined2(FluidVF& U, FluidVF& W);
							
};

#endif

//*******************************  End of force_model.cc **************************************

	
