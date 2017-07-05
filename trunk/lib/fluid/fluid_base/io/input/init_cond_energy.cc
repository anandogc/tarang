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


#include "FluidIO.h"


// extern Uniform<Real> SPECrand;

//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters eps_u = IC(1); sk = Hk/(k*ek) = IC(2).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */


//New helicity spectrum function --> Satyajit
void  FluidIO::Init_cond_energy_helicity_spectrum_k_ne0(FluidVF& U)
{
    
    Real Kmag, ek, hk;
    Real amp_u_ch1, phase_u_ch1;  // Amplitude along Craya-herring basis vector e1.
    Real amp_u_plus, amp_u_minus, phase_u_plus, phase_u_minus; // Amplitudes & phases along the helical basis
    int index;
    
    Uniform<Real> IC_rand;
    IC_rand.seed((unsigned int) my_id*time(0));
    
    for (int lx=0; lx<global.field.maxlx; lx++)
        for (int ly=0; ly<global.field.maxly; ly++)
            for (int lz=0; lz<global.field.maxlz; lz++) {
                Kmag = universal->Kmagnitude(lx, ly, lz);
                
                if (Kmag > MYEPS) {
                    index = (int) ceil(Kmag);
                    ek = U.IC_energy_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                    
                    if (ek < MYEPS) {
                        U.cvf.V1(lx,ly,lz) = 0;
                        U.cvf.V2(lx,ly,lz) = 0;
                        U.cvf.V3(lx,ly,lz) = 0;
                    }
                    else {
                        if ((Ny == 1) && global.program.two_dimension) {  // 2D 2C
                            amp_u_ch1 = sqrt(2*ek);
                            phase_u_ch1 = 2*M_PI * IC_rand.random();
                            U.cvf.Put_or_add_vector(lx, lz, amp_u_ch1, amp_u_ch1, 0);
                            // The last arg=0: do not add to exisiting field
                        }
                        else { // 3D or 2D3C
                            hk = U.IC_helicity_spectrum(index)/ global.spectrum.shell.modes_in_shell(index);
                            
                            if ( (Kmag*ek-abs(hk)) > MYEPS ) {
                                amp_u_plus = sqrt((ek + hk/Kmag)/TWO);
                                amp_u_minus = sqrt((ek - hk/Kmag)/TWO);
                                phase_u_plus = 2*M_PI * IC_rand.random();
                                phase_u_minus = 2*M_PI * IC_rand.random();
                                
                                U.cvf.Put_or_add_vector(lx, ly, lz, amp_u_plus, amp_u_minus, phase_u_plus, phase_u_minus, 0);
                                
                                // The last arg=0: do not add to exisiting field
                            }
                            else {
                                if (master) cout << "EXITING; ERROR: You have |Hk| > k* U_energy_supply_k|" << endl;
                                exit(1);
                            }
                        }
                    }
                }  // of if(Kmag > MYEPS)
            }
}

//*********************************************************************************************

void  FluidIO::Init_cond_energy_helicity_spectrum_k_ne0(FluidSF& T)
{
	Real Kmag, ek, amp, phase;
	int index;
	
    Uniform<Real> IC_rand;
    IC_rand.seed((unsigned int) my_id*time(0));
	
    for (int lx=0; lx<global.field.maxlx; lx++)
        for (int ly=0; ly<global.field.maxly; ly++)
            for (int lz=0; lz<global.field.maxlz; lz++) {
				Kmag = universal->Kmagnitude(lx, ly, lz);
				
				if (Kmag > MYEPS) {
					index = (int) ceil(Kmag);
					
					ek = T.IC_energy_spectrum(index)/ universal->Approx_number_modes_in_shell(index);
					amp = sqrt(2*ek);
					phase = 2*M_PI * IC_rand.random();
					
					T.csf.Put_or_add_vector(lx, ly, lz, amp, phase, 0);
                    // The last arg=0: do not add to exisiting field
				} // of if(Kmag > MYEPS)			
			}
}

//*********************************************************************************************

void FluidIO::Init_cond_energy_helicity_spectrum(FluidVF& U) 
{
    Real Ux000 = 0.0;
    Real Uy000 = 0.0;
    Real Uz000 = 0.0;

    if (global.io.double_para.size() == 3) {
        Ux000 = global.io.double_para(0);
        Uy000 = global.io.double_para(1);
        Uz000 = global.io.double_para(2);
    }
    
    U.IC_energy_spectrum = global.io.U_IC_energy_spectrum;
    U.IC_helicity_spectrum = global.io.U_IC_helicity_spectrum;
    
    cout << "here" << sum(U.IC_energy_spectrum) << U.IC_energy_spectrum << endl;
    
    cout << "here" << sum(U.IC_helicity_spectrum) << U.IC_helicity_spectrum << endl;
    
    Init_cond_energy_helicity_spectrum_k_ne0(U);
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
          U.cvf.V1(0,0,0).real(Ux000);
          U.cvf.V2(0,0,0).real(Uy000);
          U.cvf.V3(0,0,0).real(Uz000);
        }
        else {
            U.cvf.V1(0,0,0) = Complex(Ux000,0.0);
            U.cvf.V2(0,0,0) = Complex(Uy000,0.0);
            U.cvf.V3(0,0,0) = Complex(Uz000,0.0);
        }
    }
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

void  FluidIO::Init_cond_energy_helicity_spectrum(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR") 
		Init_cond_energy_helicity_spectrum_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Init_cond_energy_helicity_spectrum_RBC(U, T);
}


void  FluidIO::Init_cond_energy_helicity_spectrum_scalar(FluidVF& U, FluidSF& T)
{
    Real Ux000 = 0.0;
    Real Uy000 = 0.0;
    Real Uz000 = 0.0;
    Real T000 = 0.0;
    
   if (global.io.double_para.size() == 4) {
        Ux000 = global.io.double_para(0);
        Uy000 = global.io.double_para(1);
        Uz000 = global.io.double_para(2);
        T000 = global.io.double_para(3);
    }
	
    U.IC_energy_spectrum = global.io.U_IC_energy_spectrum;
    T.IC_energy_spectrum = global.io.T_IC_energy_spectrum;

    U.IC_helicity_spectrum = global.io.U_IC_helicity_spectrum;
    
	Init_cond_energy_helicity_spectrum_k_ne0(U);
	Init_cond_energy_helicity_spectrum_k_ne0(T);
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
          U.cvf.V1(0,0,0).real(Ux000);
          U.cvf.V2(0,0,0).real(Uy000);
          U.cvf.V3(0,0,0).real(Uz000);
          T.csf.F(0,0,0).real(T000);
        }
        else {
            U.cvf.V1(0,0,0) = Complex(Ux000,0.0);
            U.cvf.V2(0,0,0) = Complex(Uy000,0.0);
            U.cvf.V3(0,0,0) = Complex(Uz000,0.0);
            T.csf.F(0,0,0)  = Complex(T000,0.0);
        }
    }
}


void  FluidIO::Init_cond_energy_helicity_spectrum_RBC(FluidVF& U, FluidSF& T)
{
	
	if (global.PHYSICS.Pr_option == "PRZERO") {
		Init_cond_energy_helicity_spectrum_k_ne0(U);
		U.Zero_Prandtl_number_compute_temperature(T);		
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
        Init_cond_energy_helicity_spectrum_k_ne0(T);
        
        U.Infinite_Prandtl_number_compute_velocity(T);
	}
	
    else {
        Init_cond_energy_helicity_spectrum_scalar(U, T);
    }
		
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

void  FluidIO::Init_cond_energy_helicity_spectrum_k_ne0(FluidVF& U, FluidVF& W)
{
    if (abs(abs(sum(W.IC_crosshelicity_spectrum))) < MYEPS) {  // zero cross helicity
        Init_cond_energy_helicity_spectrum_k_ne0(U);
        Init_cond_energy_helicity_spectrum_k_ne0(W);
    }
    
    else {  // nonzero cross helicity  // TO REDO
        
        /*		Real Kmag, ek, hk_by_k, u_plus_mod, u_minus_mod;
         Real eWk, hWk_times_k, W_plus_mod, W_minus_mod;
         Real Hck, phi_u_plus, phi_u_minus, phi_W_plus, phi_W_minus;
         Real alpha, delta_plus, delta_minus, corrected_phi_u_plus, corrected_phi_u_minus, corrected_phi_W_plus, corrected_phi_W_minus;
         
         for (int lx=0; lx<global.field.maxlx; lx++)
         for (int ly=0; ly<global.field.maxly; ly++)
         for (int lz=0; lz<global.field.maxlz; lz++) {
         
         Kmag = universal->Kmagnitude(lx, ly, lz);
         
         if (Kmag > MYEPS) {
         index = int(ceil(Kmag));
         
         ek = Correlation::shell_ek1(index)/ universal->Approx_number_modes_in_shell(index);
         hk_by_k = sk*ek;
         
         eWk = Correlation::shell_ek2(index)/ universal->Approx_number_modes_in_shell(index);
         hWk_times_k = Wsk*eWk;
         
         u_plus_mod = sqrt((ek + hk_by_k)/2);
         u_minus_mod = sqrt(fabs((ek - hk_by_k)/2));
         
         W_plus_mod = sqrt((eWk + hWk_times_k)/2);
         W_minus_mod = sqrt(fabs((eWk - hWk_times_k)/2));
         phi_u_plus = 2*M_PI * SPECrand.random();
         phi_u_minus = 2*M_PI * SPECrand.random();
         phi_W_plus = 2*M_PI * SPECrand.random();
         phi_W_minus = 2*M_PI * SPECrand.random();
         if (hc==1){
         alpha=0;
         delta_plus = 0.5*(alpha + phi_W_plus - phi_u_plus);
         delta_minus = 0.5*(alpha + phi_W_minus - phi_u_minus);
         
         }else if(hc==-1){
         alpha=M_PI;
         delta_plus = 0.5*(alpha + phi_W_plus - phi_u_plus);
         delta_minus = 0.5*(alpha + phi_W_minus - phi_u_minus);
         }
         else{
         Hck = hc*(ek + eWk)/2.0;
         alpha = acos(Hck/(u_plus_mod * W_plus_mod + u_minus_mod * W_minus_mod));
         delta_plus = 0.5*(alpha + phi_W_plus - phi_u_plus);
         delta_minus = 0.5*(alpha + phi_W_minus - phi_u_minus);
         }
         corrected_phi_u_plus = phi_u_plus + delta_plus;
         corrected_phi_u_minus = phi_u_minus + delta_minus;
         corrected_phi_W_plus = phi_W_plus - delta_plus;
         corrected_phi_W_minus = phi_W_minus - delta_minus;
         
         Put_vector_amp_phase_comp_conj_cross_helicity(U, lx, ly, lz, u_plus_mod, u_minus_mod, corrected_phi_u_plus, corrected_phi_u_minus);
         Put_vector_amp_phase_comp_conj_cross_helicity(W, lx, ly, lz, W_plus_mod, W_minus_mod, corrected_phi_W_plus, corrected_phi_W_minus);
         }	// of if (kkmax > MYEPS)
         }  // of for loop
         */
    }
}

void  FluidIO::Init_cond_energy_helicity_spectrum(FluidVF& U, FluidVF& W)
{
    
    Real Ux000 = 0.0;
    Real Uy000 = 0.0;
    Real Uz000 = 0.0;
    Real Wx000 = 0.0;
    Real Wy000 = 0.0;
    Real Wz000 = 0.0;
    
  if (global.io.double_para.size() == 6) {
        Ux000 = global.io.double_para(0);
        Uy000 = global.io.double_para(1);
        Uz000 = global.io.double_para(2);
        
        Wx000 = global.io.double_para(3);
        Wy000 = global.io.double_para(4);
        Wz000 = global.io.double_para(5);
    }
    
    U.IC_energy_spectrum = global.io.U_IC_energy_spectrum;
    W.IC_energy_spectrum = global.io.W_IC_energy_spectrum;
	
    U.IC_helicity_spectrum = global.io.U_IC_helicity_spectrum;
    W.IC_helicity_spectrum = global.io.W_IC_helicity_spectrum;
    W.IC_crosshelicity_spectrum = global.io.W_IC_crosshelicity_spectrum;
    
    Init_cond_energy_helicity_spectrum_k_ne0(U,W);
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
          U.cvf.V1(0,0,0).real(Ux000);
          U.cvf.V2(0,0,0).real(Uy000);
          U.cvf.V3(0,0,0).real(Uz000);
          
          W.cvf.V1(0,0,0).real(Wx000);
          W.cvf.V2(0,0,0).real(Wy000);
          W.cvf.V3(0,0,0).real(Wz000);
        }
        else {
            U.cvf.V1(0,0,0) = Complex(Ux000,0.0);
            U.cvf.V2(0,0,0) = Complex(Uy000,0.0);
            U.cvf.V3(0,0,0) = Complex(Uz000,0.0);
            
            W.cvf.V1(0,0,0) = Complex(Wx000,0.0);
            W.cvf.V2(0,0,0) = Complex(Wy000,0.0);
            W.cvf.V3(0,0,0) = Complex(Wz000,0.0);
        }
    }
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
void  FluidIO::Init_cond_energy_helicity_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{
    
    Real Ux000 = 0.0;
    Real Uy000 = 0.0;
    Real Uz000 = 0.0;
    Real Wx000 = 0.0;
    Real Wy000 = 0.0;
    Real Wz000 = 0.0;
    Real T000 = 0.0;
    
    
    if (global.io.double_para.size() == 7) {
        Ux000 = global.io.double_para(0);
        Uy000 = global.io.double_para(1);
        Uz000 = global.io.double_para(2);
        
        Wx000 = global.io.double_para(3);
        Wy000 = global.io.double_para(4);
        Wz000 = global.io.double_para(5);
        T000 = global.io.double_para(6);
    }
    
    // TO WORK in global.cc
    U.IC_energy_spectrum = global.io.U_IC_energy_spectrum;
    W.IC_energy_spectrum = global.io.W_IC_energy_spectrum;
    T.IC_energy_spectrum = global.io.T_IC_energy_spectrum;
    
    U.IC_helicity_spectrum = global.io.U_IC_helicity_spectrum;
    W.IC_helicity_spectrum = global.io.W_IC_helicity_spectrum;
    W.IC_crosshelicity_spectrum = global.io.W_IC_crosshelicity_spectrum;

    
    Init_cond_energy_helicity_spectrum_k_ne0(U, W);
	
    Init_cond_energy_helicity_spectrum_k_ne0(T);
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
          
          U.cvf.V1(0,0,0).real(Ux000);
          U.cvf.V2(0,0,0).real(Uy000);
          U.cvf.V3(0,0,0).real(Uz000);
          
          W.cvf.V1(0,0,0).real(Wx000);
          W.cvf.V2(0,0,0).real(Wy000);
          W.cvf.V3(0,0,0).real(Wz000);
          T.csf.F(0,0,0).real(T000);
        }
        else {
            U.cvf.V1(0,0,0) = Complex(Ux000,0.0);
            U.cvf.V2(0,0,0) = Complex(Uy000,0.0);
            U.cvf.V3(0,0,0) = Complex(Uz000,0.0);
            
            W.cvf.V1(0,0,0) = Complex(Wx000,0.0);
            W.cvf.V2(0,0,0) = Complex(Wy000,0.0);
            W.cvf.V3(0,0,0) = Complex(Wz000,0.0);
            T.csf.F(0,0,0)  = Complex(T000,0.0);
        }
    }
}

//********************************** init_cond_energy.cc **************************************



  
