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

/*! \file  init_cond_energy.ccinit_cond_double_para
 * 
 * @brief Construct energy spectrum given parameter values.  Given spectrum construct
 *			random initial condition (phases).
 *
 * @note 	Given  Model energy spectrum for initial condition
 *		Sk(k) = a k^4 exp(-b k^1.1) / (k^4 + q^4)^(2.8/12)
 *		with q = 1.5, b = 0.02
 *
 * @note:   Satisfy reality condition is critical here for kz=0 and N[3]/2 planes.   
 *			Do not remove this function.
 *
 *	Notation:  (Ki =) kki = ki * kfactor[i]
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "FORCE.h"


extern Uniform<DP> SPECrand;

//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters eps_u = IC(1); sk = Hk/(k*ek) = IC(2).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */

// force spectrum F(k) = A k^{-e} saved in U.cvf.shell_ek;  hk_by_kek=h_F(k)/(k F(k))  
void  FORCE::Compute_force_using_random_energy_helicity_spectrum_basic(FluidVF& U, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, DP hk_by_kek, bool add_flag)
{
	DP Kmag, ek, amp, phase1, phase2, phase3;
	int index;
    
	Model_force_spectrum(force_spectrum_amplitude, force_spectrum_exponent, Correlation::shell_ek);
    
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
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
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = ky_min; ky <= ky_max; ky++)
			for (int kz = 0; kz <= kz_max; kz++) {
				
				if (universal->Probe_in_me(kx,ky,kz))  {
					lx = universal->Get_lx(kx);
					ly = universal->Get_ly(ky);
					lz = universal->Get_kz(kz);

					Kmag = universal->Kmagnitude(lx, ly, lz);
					
					if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
						index = (int) ceil(Kmag);
					
						ek = Correlation::shell_ek(index)/ universal->Approx_number_modes_in_shell(index);
						amp = sqrt(2*ek);
						
						phase1 = 2*M_PI * SPECrand.random();
						
						if (abs(hk_by_kek) > MYEPS) {		// helical
							phase2 = phase1 + M_PI/2.0;
							phase3 = asin(hk_by_kek)/2.0;			// zeta
						}
						
						else { // no helicity
							phase2 = 2*M_PI * SPECrand.random();
							phase3 = 2*M_PI * SPECrand.random();
						}
						
						Put_force_amp_phase_comp_conj(U, lx, ly, lz, amp, phase1, phase2, phase3, add_flag);
                        
                    //    cout << "li, amp = " << lx << " " << ly << " " << lz << " " << Kmag << " " << amp << endl;
						
					} // of if(Kmag > MYEPS)
				}
			}
}

void  FORCE::Compute_force_using_random_energy_helicity_spectrum_basic_assign(FluidVF& U, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, DP hk_by_kek)
{
	Compute_force_using_random_energy_helicity_spectrum_basic(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek, false);
}


void  FORCE::Compute_force_using_random_energy_helicity_spectrum_basic_add(FluidVF& U, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, DP hk_by_kek)
{
	Compute_force_using_random_energy_helicity_spectrum_basic(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek, true);
}

//*********************************************************************************************

void  FORCE::Compute_force_using_random_energy_spectrum_basic(FluidSF& T, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent, bool add_flag)
{
	DP Kmag, ek, amp, phase;
	int index;
	
	Model_force_spectrum(force_spectrum_amplitude, force_spectrum_exponent, Correlation::shell_ek);
	
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
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
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = ky_min; ky <= ky_max; ky++)
			for (int kz = 0; kz <= kz_max; kz++) {
				if (universal->Probe_in_me(kx,ky,kz))  {
					lx = universal->Get_lx(kx);
					ly = universal->Get_ly(ky);
					lz = universal->Get_kz(kz);
					
					Kmag = universal->Kmagnitude(lx, ly, lz);
					
					if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
						index = (int) ceil(Kmag);
						
						ek = Correlation::shell_ek(index)/ universal->Approx_number_modes_in_shell(index);
						amp = sqrt(2*ek);
						phase = 2*M_PI * SPECrand.random();
						
						Put_force_amp_phase_comp_conj(T, lx, ly, lz, amp, phase, add_flag);
					} // of if(Kmag > MYEPS)
				}
			}
	
	if (my_id == master_id)
		T.csf.F(0,0,0) = 0.0;
}


void  FORCE::Compute_force_using_random_energy_spectrum_basic_assign(FluidSF& T, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent)
{
	Compute_force_using_random_energy_spectrum_basic(T, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, false);
}

void  FORCE::Compute_force_using_random_energy_spectrum_basic_add(FluidSF& T, DP inner_radius, DP outer_radius, DP force_spectrum_amplitude, DP force_spectrum_exponent)
{
	Compute_force_using_random_energy_spectrum_basic(T, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, true);
}


//*********************************************************************************************

void  FORCE::Compute_force_using_random_energy_helicity_spectrum(FluidVF& U) 
{
	DP inner_radius = global.force.double_para(0);
	DP outer_radius = global.force.double_para(1);
	DP force_spectrum_amplitude = global.force.double_para(2);
	DP force_spectrum_exponent = global.force.double_para(3);
	DP hk_by_kek = global.force.double_para(4);
    
    Compute_force_using_random_energy_helicity_spectrum_basic_assign(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek);
}

//*********************************************************************************************
/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters eps_u = IC(1); sk = Hk/(k*ek) = IC(2); eps_T = IC(3),
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */

void  FORCE::Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidSF& T)
{    
	DP inner_radius = global.force.double_para(0);
	DP outer_radius = global.force.double_para(1);
    DP force_spectrum_amplitude = global.force.double_para(2);
	DP force_spectrum_exponent = global.force.double_para(3);
	DP hk_by_kek = global.force.double_para(4);
	DP Tforce_spectrum_amplitude = global.force.double_para(5);
	DP Tforce_spectrum_exponent = global.force.double_para(6);
	
	Compute_force_using_random_energy_helicity_spectrum_basic_assign(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek);
	Compute_force_using_random_energy_spectrum_basic_assign(T, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent);
}


//*********************************************************************************************
/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters eps_U = IC(1); sk = IC(2); 
 *			for W: eps_W = IC(3), skW =  HM(k) * k/ EW(k) = IC(4), 
 *					h = 2*Hc / (amp, ampW) = IC(5).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */
void  FORCE::Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidVF& W)
{
	DP inner_radius = global.force.double_para(0);
	DP outer_radius = global.force.double_para(1);
	DP force_spectrum_amplitude = global.force.double_para(2);
	DP force_spectrum_exponent = global.force.double_para(3);
	DP hk_by_kek = global.force.double_para(4);
	DP Wforce_spectrum_amplitude = global.force.double_para(5);
	DP Wforce_spectrum_exponent = global.force.double_para(6);
	DP Whk_by_kek = global.force.double_para(7);
    
    Compute_force_using_random_energy_helicity_spectrum_basic_assign(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek);
	
	Compute_force_using_random_energy_helicity_spectrum_basic_assign(W, inner_radius, outer_radius, Wforce_spectrum_amplitude, Wforce_spectrum_exponent, Whk_by_kek);
}
											  


//*********************************************************************************************
/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters a = IC(1); sk = IC(2); 
 *			for W: a = IC(3), skW =  HM(k) * k/ EW(k) = IC(4), 
 *					h = 2*Hc / (amp, ampW) = IC(5).
 *			for T: a = IC(6).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */
void  FORCE::Compute_force_using_random_energy_helicity_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{
	DP inner_radius = global.force.double_para(0);
	DP outer_radius = global.force.double_para(1);
	DP force_spectrum_amplitude = global.force.double_para(2);
	DP force_spectrum_exponent = global.force.double_para(3);
	DP hk_by_kek = global.force.double_para(4);
	
	DP Wforce_spectrum_amplitude = global.force.double_para(5);
	DP Wforce_spectrum_exponent = global.force.double_para(6);
	DP Whk_by_kek = global.force.double_para(7);
	
	DP Tforce_spectrum_amplitude = global.force.double_para(8);
	DP Tforce_spectrum_exponent = global.force.double_para(9);
    
    Compute_force_using_random_energy_helicity_spectrum_basic_assign(U, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent, hk_by_kek);
	
	Compute_force_using_random_energy_helicity_spectrum_basic_assign(W, inner_radius, outer_radius, Wforce_spectrum_amplitude, Wforce_spectrum_exponent, Whk_by_kek);
	
	Compute_force_using_random_energy_spectrum_basic_assign(T, inner_radius, outer_radius, force_spectrum_amplitude, force_spectrum_exponent);
	
	
}

//********************************** init_cond_energy.cc **************************************



  
