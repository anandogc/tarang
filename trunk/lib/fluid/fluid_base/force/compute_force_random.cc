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


//*********************************************************************************************

/** @brief Create initial condition based on given energy supply and kinetic helicity supply.
 *			The force is updated every random_interval = 0.01, between which the same force is
 *          applied at every time step.  Note that the random seed remains the same betwween the interval.
 *
 * @note  The supply rates are divided equally among all the modes in the shell.
 * @note  helicity = U.curl U
 * 
 * @note	The mean mode has zero energy.
 */


void  FORCE::Compute_force_using_random_noise_basic(FluidVF& U, Real random_interval, Uniform<Real> rand_struct, int rand_seed, bool add_flag)
{
    
    // No forcing....
    Real total_energy_supply = sum(U.energy_supply_spectrum);
    Real total_helicity_supply = sum(U.helicity_supply_spectrum);
    if (abs(total_energy_supply)+abs(total_helicity_supply) < MYEPS)
        return;
    
    // force ..
    Real Kmag, energy_supply_k, helicity_supply_k;
    Real amp_ch1, phase_ch1;  // Amplitudes in Craya-Herring basis.
    Real amp_plus, amp_minus, phase_plus, phase_minus;
    int index;
    
    rand_struct.seed(rand_seed);
    
    int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
    kx_max = (int) ceil(outer_radius/kfactor[1]);
    
    if (Ny > 1)
        ky_max = (int) ceil(outer_radius/kfactor[2]);
    else
        ky_max = 0;
    
    kz_max = (int) ceil(outer_radius/kfactor[3]);
    
    kx_min = ky_min = kz_min = 0;
    
    if (basis_type == "FFF" || basis_type == "FFFW")
        kx_min = -kx_max;
    
    if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
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
                        energy_supply_k = U.energy_supply_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                        
                        if ((Ny == 1) && global.program.two_dimension) {  // 2D 2C
                            amp_ch1 = sqrt(2*energy_supply_k/random_interval);
                            phase_ch1 = 2*M_PI * rand_struct.random();
                            
                            Put_or_add_force_vector(U, lx, lz, amp_ch1, phase_ch1, add_flag);
                        }
                        else { //  for 3D & 2D3C
                            helicity_supply_k = U.helicity_supply_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                            
                            if (abs(Kmag*energy_supply_k - abs(helicity_supply_k)) > MYEPS) {
                                amp_plus = sqrt((energy_supply_k+helicity_supply_k/Kmag)/random_interval);
                                amp_minus = sqrt((energy_supply_k-helicity_supply_k/Kmag)/random_interval);
                                phase_plus = 2*M_PI * rand_struct.random();
                                phase_minus = 2*M_PI * rand_struct.random();
                                
                                Put_or_add_force_vector(U, lx, ly, lz, amp_plus, amp_minus, phase_plus, phase_minus, add_flag);
                            }
                            else {
                                if (master) cout << "EXITING; ERROR: You have |helicity_supply_k| <= k* energy_supply_k|";
                                exit(1);
                            }
                        }
                    }
                }
            }
}

void  FORCE::Compute_force_using_random_noise_assign(FluidVF& U, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed;
    rand_seed = (unsigned int) my_id*time(0);
    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed = (unsigned int) my_id*time(0);
        random_next += random_interval;
        //   if (master) cout <<" time random_next  = "<<global.time.now<<" "<< random_next << " " << random_interval << endl;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed, false);
}


void  FORCE::Compute_force_using_random_noise_add(FluidVF& U, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed;
    rand_seed = (unsigned int) my_id*time(0);
    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed = (unsigned int) my_id*time(0);
        random_next += random_interval;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed, true);
}

//*********************************************************************************************

/** @brief Create initial condition based on given energy supply and kinetic helicity supply.
 *			The force is updated every random_interval = 0.01, between which the same force is
 *          applied at every time step.  Note that the random seed remains the same betwween the interval.
 *
 * @note  The supply rates are divided equally among all the modes in the shell.
 * @note  helicity = U.curl U
 *
 * @note	The mean mode has zero energy.
 */

void  FORCE::Compute_force_using_random_noise_basic(FluidSF& T, Real random_interval, Uniform<Real> rand_struct, int rand_seed,  bool add_flag)
{

    // No random forcing if supplies are zero.
    Real total_energy_supply = sum(T.energy_supply_spectrum);
    if (abs(total_energy_supply) < MYEPS)
        return;

    // Force now
    Real Kmag, energy_supply_k;
    Real amp, phase; // For scalar
    int index;
    static Real random_next = 0;
    
    rand_struct.seed(rand_seed);
    
    int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
    kx_max = (int) ceil(outer_radius/kfactor[1]);
    
    if (Ny > 1)
        ky_max = (int) ceil(outer_radius/kfactor[2]);
    else
        ky_max = 0;
    
    kz_max = (int) ceil(outer_radius/kfactor[3]);
    
    kx_min = ky_min = kz_min = 0;
    
    if (basis_type == "FFF" || basis_type == "FFFW")
        kx_min = -kx_max;
    
    if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
        ky_min = -ky_max;
    
    int lx, ly, lz;
    for (int kx = kx_min; kx <= kx_max; kx++)
        for (int ky = ky_min; ky <= ky_max; ky++)
            for (int kz = 0; kz <= kz_max; kz++)
                
                if (universal->Probe_in_me(kx,ky,kz))  {
                    lx = universal->Get_lx(kx);
                    ly = universal->Get_ly(ky);
                    lz = universal->Get_kz(kz);
                    
                    Kmag = universal->Kmagnitude(lx, ly, lz);
                    
                    if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
                        index = (int) ceil(Kmag);
                        
                        if (T.force_switch) {
                            energy_supply_k = T.energy_supply_spectrum(index)/  global.spectrum.shell.modes_in_shell(index);
                            amp = sqrt(2*energy_supply_k/random_interval);
                            phase = 2*M_PI * rand_struct.random();
                            
                            Put_or_add_force_scalar(T, lx, ly, lz, amp, phase, add_flag);
                        }
                    }
                }
}




void  FORCE::Compute_force_using_random_noise_assign(FluidVF& U, FluidSF& T, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed1, rand_seed2;
    rand_seed1 = ((unsigned int) my_id*time(0));
    rand_seed2 = ((unsigned int) my_id*time(0)+1);
    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed1 = ((unsigned int) my_id*time(0));
        rand_seed2 = ((unsigned int) my_id*time(0)+1);
        random_next += random_interval;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed1, false);
    
    if (T.force_switch)
        Compute_force_using_random_noise_basic(T, random_interval, force_rand,  rand_seed2, false);
}

void  FORCE::Compute_force_using_random_noise_add(FluidVF& U, FluidSF& T, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed1, rand_seed2;
    rand_seed1 = ((unsigned int) my_id*time(0));
    rand_seed2 = ((unsigned int) my_id*time(0)+1);

    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed1 = ((unsigned int) my_id*time(0));
        rand_seed2 = ((unsigned int) my_id*time(0)+1);
        random_next += random_interval;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed1, true);
    
    if (T.force_switch)
        Compute_force_using_random_noise_basic(T, random_interval, force_rand, rand_seed2, true);
}

//*********************************************************************************************


void  FORCE::Compute_force_using_random_noise_assign(FluidVF& U, FluidVF& W, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed1, rand_seed2;
    rand_seed1 = ((unsigned int) my_id*time(0));
    rand_seed2 = ((unsigned int) my_id*time(0)+1);
    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed1 = ((unsigned int) my_id*time(0));
        rand_seed2 = ((unsigned int) my_id*time(0)+1);
        random_next += random_interval;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed1, false);
    Compute_force_using_random_noise_basic(W, random_interval, force_rand, rand_seed2, false);
}


void  FORCE::Compute_force_using_random_noise_add(FluidVF& U, FluidVF& W, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed1, rand_seed2;
    rand_seed1 = ((unsigned int) my_id*time(0));
    rand_seed2 = ((unsigned int) my_id*time(0)+1);
    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed1 = ((unsigned int) my_id*time(0));
        rand_seed2 = ((unsigned int) my_id*time(0)+1);
        random_next += random_interval;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed1, true);
    Compute_force_using_random_noise_basic(W, random_interval, force_rand, rand_seed2, true);
}

//*********************************************************************************************


void  FORCE::Compute_force_using_random_noise_assign(FluidVF& U, FluidVF& W, FluidSF& T, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed1, rand_seed2, rand_seed3;
    rand_seed1 = ((unsigned int) my_id*time(0));
    rand_seed2 = ((unsigned int) my_id*time(0)+1);
    rand_seed3 = ((unsigned int) my_id*time(0)+2);
    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed1 = ((unsigned int) my_id*time(0));
        rand_seed2 = ((unsigned int) my_id*time(0)+1);
        rand_seed3 = ((unsigned int) my_id*time(0)+2);
        random_next += random_interval;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed1, false);
    Compute_force_using_random_noise_basic(W, random_interval, force_rand, rand_seed2, false);
    Compute_force_using_random_noise_basic(T, random_interval, force_rand, rand_seed3, false);
}


void  FORCE::Compute_force_using_random_noise_add(FluidVF& U, FluidVF& W, FluidSF& T, Real random_interval)
{
    static Real random_next = 0;
    
    static Uniform<Real> force_rand;
    static int rand_seed1, rand_seed2, rand_seed3;
    rand_seed1 = ((unsigned int) my_id*time(0));
    rand_seed2 = ((unsigned int) my_id*time(0)+1);
    rand_seed3 = ((unsigned int) my_id*time(0)+2);
    
    if (global.time.now >=random_next) { // Get new rand_seed
        rand_seed1 = ((unsigned int) my_id*time(0));
        rand_seed2 = ((unsigned int) my_id*time(0)+1);
        rand_seed3 = ((unsigned int) my_id*time(0)+2);
        random_next += random_interval;
    } // else rand_seed is unchanged.
    
    Compute_force_using_random_noise_basic(U, random_interval, force_rand, rand_seed1, true);
    Compute_force_using_random_noise_basic(W, random_interval, force_rand, rand_seed2, true);
    Compute_force_using_random_noise_basic(T, random_interval, force_rand, rand_seed3, true);
}

//*********************************************************************************************

void  FORCE::Compute_force_using_random_noise(FluidVF& U)
{
    inner_radius = global.force.double_para(0);
    outer_radius = global.force.double_para(1);
    
    Real random_interval = global.force.double_para(2);
    
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;
    
    Compute_force_using_random_noise_assign(U, random_interval);
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

void  FORCE::Compute_force_using_random_noise(FluidVF& U, FluidSF& T)
{    
    inner_radius = global.force.double_para(0);
    outer_radius = global.force.double_para(1);
    
    Real random_interval = global.force.double_para(2);
    
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    T.energy_supply_spectrum = global.force.T_energy_supply_spectrum;
    
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;
	
	Compute_force_using_random_noise_assign(U, T, random_interval);
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
void  FORCE::Compute_force_using_random_noise(FluidVF& U, FluidVF& W)
{
	inner_radius = global.force.double_para(0);
	outer_radius = global.force.double_para(1);
    
    Real random_interval = global.force.double_para(2);
    
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    W.energy_supply_spectrum = global.force.W_energy_supply_spectrum;
    
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;
    W.helicity_supply_spectrum = global.force.W_helicity_supply_spectrum;
    W.crosshelicity_supply_spectrum = global.force.W_crosshelicity_supply_spectrum;
    
    Compute_force_using_random_noise_assign(U, W, random_interval);
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
void  FORCE::Compute_force_using_random_noise(FluidVF& U, FluidVF& W, FluidSF& T)
{
    inner_radius = global.force.double_para(0);
    outer_radius = global.force.double_para(1);
    
    Real random_interval = global.force.double_para(2);
    
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    W.energy_supply_spectrum = global.force.W_energy_supply_spectrum;
    T.energy_supply_spectrum = global.force.T_energy_supply_spectrum;
    
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;
    W.helicity_supply_spectrum = global.force.W_helicity_supply_spectrum;
    W.crosshelicity_supply_spectrum = global.force.W_crosshelicity_supply_spectrum;
    
    Compute_force_using_random_noise_assign(U, W, T, random_interval);
}

//********************************** init_cond_energy.cc **************************************



  
