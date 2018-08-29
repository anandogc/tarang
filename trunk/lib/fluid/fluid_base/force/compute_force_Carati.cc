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

/*! \file  compute_force_ek_hk_supply.cc
 * 
 * @brief Compute force when ek and hk supply rate is given
 *
 * @note 2D:   F(k) = alpha * V(k)
 * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *				alpha and beta determined from the supply rates
 *
 * @note:   Satisfy reality condition is critical here for kz=0 and N[3]/2 planes.   
 *			Do not remove this function.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug sk=1, -1 needs to be handled separately.
 */

#include "FORCE.h"


//*********************************************************************************************
// Real energy_supply, Real epsh_by_k_epse
// Real energy_level, h_by_k_E

void FORCE::Compute_force_Carati_scheme_basic(FluidVF& U, string force_type, bool global_alpha_beta, bool add_flag)
{
    
    if (force_type == "ENERGY_SUPPLY")
        Compute_force_Carati_scheme_energy_supply(U, global_alpha_beta, add_flag);
    else if (force_type == "CONSTANT_ENERGY")
        Compute_force_Carati_scheme_const_energy(U, global_alpha_beta, add_flag);
}

void FORCE::Compute_force_Carati_scheme_energy_supply(FluidVF& U, bool global_alpha_beta, bool add_flag)
{

    // No forcing....
    Real total_energy_supply = sum(U.energy_supply_spectrum);
    Real total_helicity_supply = sum(U.helicity_supply_spectrum);
    if ((abs(total_energy_supply) + abs(total_helicity_supply)) < MYEPS)
        return;
    
    // Force now
    Real energy_supply_k, helicity_supply_k;  // for the mode {\bf k}
    Real denr;
    Real epsh_by_k_epse;
    
    int lx, ly, lz;
    int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;

    Real Kmag, alpha_k, beta_k, sk;
    Real modal_energy, modal_helicity;
	Real temp, temp1, temp2, temp3;
    int index;
    
	
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
	
    if (global_alpha_beta) {
        Real total_Ek_in_force_shell = 0;
        Real total_ksqrEk_in_force_shell = 0;
        Real total_Hk_in_force_shell = 0;
        
        Array<Real,1> alpha_k_shell(global.spectrum.shell.no_shells);
        Array<Real,1> beta_k_shell(global.spectrum.shell.no_shells);
        alpha_k_shell=0;
        beta_k_shell=0;
        
        for (int i=int(inner_radius)+1; i<= int(outer_radius); i++) {
            total_Ek_in_force_shell =  total_ksqrEk_in_force_shell =  total_Hk_in_force_shell = 0;
            
            Correlation::Compute_shell_spectrum(U);
            total_Ek_in_force_shell += Correlation::shell_ek1_force(i) + Correlation::shell_ek2_force(i) + Correlation::shell_ek3_force(i);
            total_ksqrEk_in_force_shell += Correlation::shell_dissk1(i) + Correlation::shell_dissk2(i) + Correlation::shell_dissk3(i);
            
            total_ksqrEk_in_force_shell /= TWO;   // shell_dissk1 contains 2*k^2*E(k), so div by 2
            
            Correlation::Compute_shell_spectrum_helicity(U);
            total_Hk_in_force_shell += Correlation::shell_ek1_force(i) + Correlation::shell_ek2_force(i) + Correlation::shell_ek3_force(i);
            
            denr = total_Ek_in_force_shell*total_ksqrEk_in_force_shell - my_pow(total_Hk_in_force_shell,2);
            
            if (denr > MYEPS) {
                alpha_k_shell(i) = (U.energy_supply_spectrum(i)*total_ksqrEk_in_force_shell - U.helicity_supply_spectrum(i)*total_Hk_in_force_shell)/(2*denr);
                beta_k_shell(i) = (U.helicity_supply_spectrum(i)*total_Ek_in_force_shell - U.energy_supply_spectrum(i)*total_Hk_in_force_shell)/(2*denr);

            }

            else if ((total_Ek_in_force_shell > MYEPS2) && (total_ksqrEk_in_force_shell > MYEPS2)) { // helical
                alpha_k_shell(i) = U.energy_supply_spectrum(i)/(4*total_Ek_in_force_shell);
                beta_k_shell(i) = U.energy_supply_spectrum(i)/(4*sqrt(total_Ek_in_force_shell*total_ksqrEk_in_force_shell));
            }
            
        }
        for (int kx = kx_min; kx <= kx_max; kx++)
            for (int ky = ky_min; ky <= ky_max; ky++)
				for (int kz = 0; kz <= kz_max; kz++) {
					if (universal->Probe_in_me(kx,ky,kz))  {
						lx = universal->Get_lx(kx);
                        ly = universal->Get_ly(ky);
                        lz = universal->Get_lz(kz);
                            
                        Kmag = universal->Kmagnitude(lx, ly, lz);
						if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
							index = (int) ceil(Kmag);
                            
							Const_energy_supply_alpha_beta(U, lx, ly, lz, alpha_k_shell(index), beta_k_shell(index), add_flag);
						}
					}
				}
        
       // universal->Print_large_Fourier_elements(U.cvf.V1);
        return;  // Done... global_alpha_beta
    } // The above scheme works for 2D2C, 2D3C, and 3D.
    
    // Now alpha_beta computations for !global_alpha_beta
    
    for (int kx = kx_min; kx <= kx_max; kx++)
        for (int ky = ky_min; ky <= ky_max; ky++)
            for (int kz = 0; kz <= kz_max; kz++)
                if (universal->Probe_in_me(kx,ky,kz))  {  // for the modes inside the proc
                    lx = universal->Get_lx(kx);
                    ly = universal->Get_ly(ky);
                    lz = universal->Get_lz(kz);
                    
                    Kmag = universal->Kmagnitude(lx, ly, lz);
                    if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
                        
                        index = (int) ceil(Kmag);
                        modal_energy = U.cvf.Modal_energy(lx, ly, lz);
                        if (modal_energy > MYEPS) {
                            
                            energy_supply_k = U.energy_supply_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                            
                            if ((Ny == 1) && global.program.two_dimension) {  // 2D 2C
                                alpha_k = energy_supply_k/(2*modal_energy);
                                beta_k = 0;
                            }
                            else { // 3D or 2D3C
                                helicity_supply_k = U.helicity_supply_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                                modal_helicity = U.cvf.Modal_helicity(lx,ly,lz);
                                
                                denr = TWO*(my_pow(Kmag*modal_energy,2) - my_pow(modal_helicity, 2));
                                epsh_by_k_epse = helicity_supply_k/(Kmag*energy_supply_k);
                                
                                if (abs(denr) > MYEPS) {
                                    alpha_k = (my_pow(Kmag,2)*energy_supply_k*modal_energy - helicity_supply_k*modal_helicity)/denr;
                                    
                                    beta_k = (helicity_supply_k*modal_energy - energy_supply_k*modal_helicity)/denr;
                                }
                                else if (abs(epsh_by_k_epse-1) < MYEPS) {
                                    alpha_k = energy_supply_k/(4*modal_energy);
                                    beta_k = alpha_k/Kmag;
                                }
                                else {
                                    if (master)
                                        cout << "ERROR: Max helicity supply case: Hk approx k*Ek" << endl;
                                    exit(1);
                                } // of error
                            }   // 3D or 2D3C
                            
                            Const_energy_supply_alpha_beta(U, lx, ly, lz, alpha_k, beta_k, add_flag);
                            
                        } // (modal_energy > MYEPS)  F computation for both 2D and 3D
                    } // Kmag cond
                } // for the modes inside the proc
                
}

void FORCE::Compute_force_Carati_scheme_const_energy(FluidVF& U, bool global_alpha_beta, bool add_flag)
{
    // TO DO
  /*  for (int kx = kx_min; kx <= kx_max; kx++)
        for (int ky = ky_min; ky <= ky_max; ky++)
            for (int kz = 0; kz <= kz_max; kz++) {
                
                if (universal->Probe_in_me(kx,ky,kz))  {
                    lx = universal->Get_lx(kx);
                    ly = universal->Get_ly(ky);
                    lz = universal->Get_lz(kz);
                    
                    Kmag = universal->Kmagnitude(lx, ly, lz);
                    if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
                        index = (int) ceil(Kmag);
                        modal_energy = U.cvf.Modal_energy(lx, ly, lz);
                        if (modal_energy > MYEPS) {
                            
                            if (force_type == "ENERGY_SUPPLY") {
                                // if (!global_alpha_beta) {
                                energy_supply_k = U.energy_supply_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                                helicity_supply_k = U.helicity_supply_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                                
                                if ((Ny == 1) && global.program.two_dimension) {  // 2D 2C
                                    alpha_k = energy_supply_k/(2*modal_energy);
                                    beta_k = 0;
                                }
                                else { // 3D or 2D3C
                                    sk = U.cvf.Modal_helicity(lx,ly,lz)/ (Kmag*modal_energy);
                                    
                                    if (energy_supply_k > MYEPS2) { // energy_supply > 0
                                        temp = energy_supply_k/ (2*modal_energy);
                                        if (abs(sk*sk-1) > MYEPS2) {
                                            epsh_by_k_epse = helicity_supply_k/(Kmag*energy_supply_k);
                                            alpha_k = temp * (1-sk*epsh_by_k_epse) / (1 - sk*sk);
                                            beta_k = temp * (epsh_by_k_epse - sk) / ((1 - sk*sk)*Kmag);
                                        }
                                        
                                        else { // max helicity
                                            alpha_k = temp/2;
                                            beta_k = alpha_k/(sk*Kmag);
                                        }
                                    }
                                    else { // energy_supply = 0, force using helicity_supply
                                        if (abs(sk*sk-1) > MYEPS2) {
                                            temp = helicity_supply_k/ (2*Kmag*modal_energy);
                                            alpha_k = -temp * sk / (1 - sk*sk);
                                            beta_k = temp/ ((1 - sk*sk)*Kmag);
                                        }
                                        
                                        else {
                                            if (master)
                                                cout << "ERROR: Max heliicty supply case: Hk approx k*Ek" << endl;
                                            exit(1);
                                        }
                                    }
                                }
                                //   }
                                
                                Const_energy_supply_alpha_beta(U, lx, ly, lz, alpha_k, beta_k, add_flag);
                                
                            } // end of (force_type == "ENERGY_SUPPLY")
                            
                            else if (force_type == "CONSTANT_ENERGY") {
                                /*	temp1 = sqrt(energy_per_mode/modal_energy);
                                 
                                 if (abs(h_by_k_E) < MYEPS) { // No helical forcing
                                 alpha_k = temp1;
                                 beta_k = 0.0;
                                 }
                                 
                                 else {
                                 sk = U.cvf.Modal_helicity(lx,ly,lz)/ (Kmag*modal_energy);
                                 
                                 if (abs(sk*sk-1) > MYEPS2) {
                                 temp2 = sqrt((1+h_by_k_E)/ (1+sk));
                                 temp3 = sqrt((1-h_by_k_E)/ (1-sk));
                                 
                                 alpha_k = (temp1/2) * (temp2 + temp3);
                                 beta_k =  (temp1/(2*Kmag)) * (temp2 - temp3);
                                 }
                                 
                                 else {
                                 alpha_k = temp1/2;
                                 beta_k = alpha_k/(sk*Kmag);
                                 }	
                                 }
                                 Const_energy_alpha_beta(U, lx, ly, lz, alpha_k, beta_k, add_flag);
                                 */
    
}


// Shubhadeep & A. G. Chatterjee	
void FORCE::Compute_force_Carati_scheme(FluidVF& U, string force_type,  bool global_alpha_beta)
{

	Compute_force_Carati_scheme_basic(U, force_type, global_alpha_beta, true);
}
		


// Shubhadeep & A. G. Chatterjee
void FORCE::Compute_force_Carati_scheme(FluidVF& U, FluidVF& W, string force_type, bool global_alpha_beta)
{

    if (U.force_switch)
        Compute_force_Carati_scheme_basic(U, force_type, global_alpha_beta, true);
    
    if (W.force_switch)
        Compute_force_Carati_scheme_basic(W, force_type, global_alpha_beta, true);
}

//*********************************************************************************************
//
//	Scalar
//

void FORCE::Compute_force_Carati_scheme_basic(FluidSF& T, string force_type, bool global_alpha_beta, bool add_flag)
{

    // No forcing....
    Real total_energy_supply_level = sum(T.energy_supply_spectrum);
    if (abs(total_energy_supply_level) < MYEPS)
        return;
	
    if (T.force_switch) {
        int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
        int lx, ly, lz;

        Real modal_energy;
        Real alpha_k, Kmag, energy_supply_k;
        int index;
        
        T.Force = 0.0;
        
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
        
        for (int kx = kx_min; kx <= kx_max; kx++)
            for (int ky = ky_min; ky <= ky_max; ky++)
                for (int kz = 0; kz <= kz_max; kz++) {
                    
                    if (universal->Probe_in_me(kx, ky, kz)) {
                        lx = universal->Get_lx(kx);
                        ly = universal->Get_ly(ky);
                        lz = universal->Get_lz(kz);
                        
                        Kmag = universal->Kmagnitude(lx, ly, lz);
                        if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
                            index = (int) ceil(Kmag);
                            energy_supply_k = T.energy_supply_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                            
                            modal_energy = T.csf.Modal_energy(lx, ly, lz);
                            
                            if (modal_energy > MYEPS) {
                                if (force_type == "ENERGY_SUPPLY") {
                                    alpha_k = energy_supply_k/ (2*modal_energy);
                                    Const_energy_supply_alpha(T, lx, ly, lz, alpha_k, add_flag);
                                }
                                
                                else if (force_type == "CONSTANT_ENERGY") {
                                    //	alpha_k = sqrt(energy_per_mode/modal_energy);
                                    //	Const_energy_alpha(T, lx, ly, lz, alpha_k, add_flag);
                                }
                            }
                        }
                    }   // of (Probe_in_me())						
                }		// of for
    }
}	


// Shubhadeep & A. G. Chatterjee
void FORCE::Compute_force_Carati_scheme(FluidVF& U, FluidSF& T, string force_type,   bool  global_alpha_beta)
{
    if (U.force_switch)
        Compute_force_Carati_scheme_basic(U, force_type, global_alpha_beta, true);
    
    if (T.force_switch)
        Compute_force_Carati_scheme_basic(T, force_type, global_alpha_beta, true);
}
//*********************************************************************************************
// derived fn
/*void FORCE::Compute_force_Carati_scheme(FluidVF& U)
{
    char force_type[80];
    bool global_alpha_beta;
    global.force.Get_para("%s, %d, %f, %f",force_type, &global_alpha_beta, &inner_radius, &outer_radius);

    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;   
	if (U.force_switch)
        Compute_force_Carati_scheme_assign(U, force_type, global_alpha_beta);
}*/

void FORCE::Compute_force_Carati_scheme(FluidVF& U)
{
	// Shubhadeep & A. G. Chatterjee
    char force_type[80];
    bool global_alpha_beta;
    global.force.Get_para("%s %d %f %f",force_type, &global_alpha_beta, &inner_radius, &outer_radius);
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;   
	if (U.force_switch)
        Compute_force_Carati_scheme(U, force_type, global_alpha_beta);
}


//*********************************************************************************************
// derived fn
void FORCE::Compute_force_Carati_scheme(FluidVF& U, FluidSF& T)
{
	char force_type[80];
    bool global_alpha_beta;
    global.force.Get_para("%s %d %f %f",force_type, &global_alpha_beta, &inner_radius, &outer_radius);   
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    T.energy_supply_spectrum = global.force.T_energy_supply_spectrum;
    
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;
	
    Compute_force_Carati_scheme(U, T, force_type, global_alpha_beta);
}

//*********************************************************************************************
// derived fn
void FORCE::Compute_force_Carati_scheme(FluidVF& U, FluidVF& W)
{
		
    /*
    inner_radius = global.force.double_para(0);
	outer_radius = global.force.double_para(1);
    string force_type = global.force.string_para(0);
    
    bool global_alpha_beta = !!(global.force.int_para(0));
    */
    // Shubhadeep & A. G. Chatterjee
    char force_type[80];
    bool global_alpha_beta;
    global.force.Get_para("%s %d %f %f",force_type, &global_alpha_beta, &inner_radius, &outer_radius);
    
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    W.energy_supply_spectrum = global.force.W_energy_supply_spectrum;
    
    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;
    W.helicity_supply_spectrum = global.force.W_helicity_supply_spectrum;
    W.crosshelicity_supply_spectrum = global.force.W_crosshelicity_supply_spectrum;
    
    Compute_force_Carati_scheme(U, W, force_type,  global_alpha_beta);
}




//*********************************************************************************************
// derived fn: assign
void FORCE::Compute_force_Carati_scheme(FluidVF& U, FluidVF& W, FluidSF& T)
{

    inner_radius = global.force.double_para(0);
    outer_radius = global.force.double_para(1);
    string force_type = global.force.string_para(0);
    
    bool global_alpha_beta = !!(global.force.int_para(0));
 
    // TO WORK IN GLOBAL READ>>>>
    U.energy_supply_spectrum = global.force.U_energy_supply_spectrum;
    W.energy_supply_spectrum = global.force.W_energy_supply_spectrum;
    T.energy_supply_spectrum = global.force.T_energy_supply_spectrum;

    U.helicity_supply_spectrum = global.force.U_helicity_supply_spectrum;
    W.helicity_supply_spectrum = global.force.W_helicity_supply_spectrum;
    W.crosshelicity_supply_spectrum = global.force.W_crosshelicity_supply_spectrum;

    Compute_force_Carati_scheme(U, W, force_type,  global_alpha_beta);
    
    if (T.force_switch)
        Compute_force_Carati_scheme_basic(T, force_type, global_alpha_beta, true);
}


//*********************************************************************************************
// derived fn
void FORCE::Compute_force_Carati_scheme_crosshelicity(FluidVF& U, FluidVF& W)
{
    
/*  inner_radius = global.force.double_para(0);
    outer_radius = global.force.double_para(1);
    Real eps = global.force.double_para(2);
    Real eps_hc = global.force.double_para(3);				// epsh(k)/(k*eps(k))
    
    if (U.force_switch)
        Compute_force_Carati_scheme_crosshelicity_basic(U, W);
 */
}

//New crosshelical forsing (15/03/2017) which wa written with Stepanov R
void FORCE::Compute_force_Carati_scheme_crosshelicity_basic(FluidVF& U, FluidVF& W,  bool global_alpha_beta, bool add_flag)
{
    /*	Real eps = para1;
     Real eps_hc = para2;
     
     //cout << "para1,2 = " << para1 << " " << para2 << endl;
     
     int nf;
     int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
     Real modal_energy;
     Real temp, temp1, temp2, temp3;
     
     nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
     
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
     Real Kmag, alpha_k, beta_k, sk;
     Real eu, eb, hc, delta_c = 0.000000001; // With small delta_c forcing will work.
     
     for (int kx = kx_min; kx <= kx_max; kx++)
     for (int ky = ky_min; ky <= ky_max; ky++)
     for (int kz = 0; kz <= kz_max; kz++) {
     
     if (universal->Probe_in_me(kx,ky,kz))  { //if this mode exists for this proc
					lx = universal->Get_lx(kx);
					ly = universal->Get_ly(ky);
					lz = universal->Get_kz(kz);
     
					Kmag = universal->Kmagnitude(lx, ly, lz);
					if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
     eu = U.cvf.Modal_energy(lx, ly, lz);
     eb = W.cvf.Modal_energy(lx, ly, lz);
     hc = 0.5*mydot(U.cvf.V1(lx, ly, lz), U.cvf.V2(lx, ly, lz), U.cvf.V3(lx, ly, lz), W.cvf.V1(lx, ly, lz), W.cvf.V2(lx, ly, lz), W.cvf.V3(lx, ly, lz));
     
     if (fabs(eu*eb - hc*hc) > MYEPS) {
     alpha_k = (0.5*eps*eb - eps_hc*hc)/(eu*eb - hc*hc);
     beta_k = (0.5*eps*hc - eps_hc*eu)/(hc*hc - eu*eb);
     }
     else {
     alpha_k = (0.5*eps*eb - eps_hc*hc)/(delta_c);
     beta_k = (0.5*eps*hc - eps_hc*eu)/(delta_c);
     }
     
     
     
     Const_energy_supply_alpha_beta(U, W, lx, ly, lz, alpha_k, beta_k, add_flag);
					}
     }	//  of if (Probe_in_me())						
     }		// of for	
     */
    
}				



//***********************  End of compute_force_ek_hk_supply.cc *******************************


