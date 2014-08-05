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

void  FluidIO::Initialize_using_energy_helicity_spectrum(FluidVF& U, DP epsilon, DP hk_by_kek)
{
	DP Kmag, ek, amp, phase1, phase2, phase3;
	int index;
	
	
    
	Model_initial_using_shell_spectrum_Pope(U.dissipation_coefficient, epsilon, Correlation::shell_ek1);
    
	for (int lx=0; lx<=global.field.maxlx; lx++)
		for (int ly=0; ly<=global.field.maxly; ly++)
			for (int lz=0; lz<=global.field.maxlz; lz++) {

				Kmag = universal->Kmagnitude(lx, ly, lz);
				
				if (Kmag > MYEPS) {
					index = (int) ceil(Kmag);
				
					ek = Correlation::shell_ek1(index)/ universal->Approx_number_modes_in_shell(index);
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
					
					Put_vector_amp_phase_comp_conj(U, lx, ly, lz, amp, phase1, phase2, phase3);
                    
				} // of if(Kmag > MYEPS)			
			}
}

//*********************************************************************************************

void  FluidIO::Initialize_using_energy_helicity_spectrum(FluidSF& T, DP epsilon)
{
	DP Kmag, ek, amp, phase;
	int index;
	
	Model_initial_using_shell_spectrum_Pope(T.diffusion_coefficient, epsilon, Correlation::shell_ek);
	
    for (int lx=0; lx<=global.field.maxlx; lx++)
        for (int ly=0; ly<=global.field.maxly; ly++)
            for (int lz=0; lz<=global.field.maxlz; lz++) {
				Kmag = universal->Kmagnitude(lx, ly, lz);
				
				if (Kmag > MYEPS) {
					index = (int) ceil(Kmag);
					
					ek = Correlation::shell_ek(index)/ universal->Approx_number_modes_in_shell(index);
					amp = sqrt(2*ek);
					phase = 2*M_PI * SPECrand.random();
					
					Put_scalar_amp_phase_comp_conj(T, lx, ly, lz, amp, phase);
				} // of if(Kmag > MYEPS)			
			}	
	
	if (my_id == master_id)
		T.csf.F(0,0,0) = 0.0;
}

//*********************************************************************************************

void  FluidIO::Init_cond_energy_helicity_spectrum(FluidVF& U) 
{
    DP epsilon, sk;
    DP Ux000 = 0.0;
    DP Uy000 = 0.0;
    DP Uz000 = 0.0;
    
    if (global.io.double_para.size() == 2) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);  // sk = Hk/(k*ek) = IC(2)
    }
    
    else if (global.io.double_para.size() == 5) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);  
        Ux000 = global.io.double_para(2);
        Uy000 = global.io.double_para(3);
        Uz000 = global.io.double_para(4);
    }
    
    Initialize_using_energy_helicity_spectrum(U, epsilon, sk);
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
            real(U.cvf.V1(0,0,0)) = Ux000;
            real(U.cvf.V2(0,0,0)) = Uy000;
            real(U.cvf.V3(0,0,0)) = Uz000;
        }
        else {
            U.cvf.V1(0,0,0) = complx(Ux000,0.0);
            U.cvf.V2(0,0,0) = complx(Uy000,0.0);
            U.cvf.V3(0,0,0) = complx(Uz000,0.0);
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

    DP epsilon, sk, Tepsilon;
    DP Ux000 = 0.0;
    DP Uy000 = 0.0;
    DP Uz000 = 0.0;
    DP T000 = 0.0;
    
    if (global.io.double_para.size() == 3) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);  // sk = Hk/(k*ek) = IC(2)
        Tepsilon = global.io.double_para(2);
    }
    
    else if (global.io.double_para.size() == 7) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);  
        Tepsilon = global.io.double_para(2);
        Ux000 = global.io.double_para(3);
        Uy000 = global.io.double_para(4);
        Uz000 = global.io.double_para(5);
        T000 = global.io.double_para(6);
    }
	
	Initialize_using_energy_helicity_spectrum(U, epsilon, sk);
	Initialize_using_energy_helicity_spectrum(T, Tepsilon);
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
            real(U.cvf.V1(0,0,0)) = Ux000;
            real(U.cvf.V2(0,0,0)) = Uy000;
            real(U.cvf.V3(0,0,0)) = Uz000;
            real(T.csf.F(0,0,0))  = T000;
        }
        else {
            U.cvf.V1(0,0,0) = complx(Ux000,0.0);
            U.cvf.V2(0,0,0) = complx(Uy000,0.0);
            U.cvf.V3(0,0,0) = complx(Uz000,0.0);
            T.csf.F(0,0,0)  = complx(T000,0.0);
        }
    }
}


void  FluidIO::Init_cond_energy_helicity_spectrum_RBC(FluidVF& U, FluidSF& T)
{
	
	if (global.PHYSICS.Pr_option == "PRZERO") {
		Init_cond_energy_helicity_spectrum(U);
		U.Zero_Prandtl_number_compute_temperature(T);		
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
		DP Tepsilon = global.io.double_para(0);
        Initialize_using_energy_helicity_spectrum(T, Tepsilon);
        
        U.Infinite_Prandtl_number_compute_velocity(T);
	}
	
	else
		Init_cond_energy_helicity_spectrum_scalar(U, T);
		
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
void  FluidIO::Init_cond_energy_helicity_spectrum(FluidVF& U, FluidVF& W)
{
    DP epsilon, sk, Wepsilon, Wsk, h;
    
    DP Ux000 = 0.0;
    DP Uy000 = 0.0;
    DP Uz000 = 0.0;
    DP Wx000 = 0.0;
    DP Wy000 = 0.0;
    DP Wz000 = 0.0;
    
    
    if (global.io.double_para.size() == 5) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);   // sk = Hk/(k*ek) = IC(2)
        
        Wepsilon = global.io.double_para(2);
        Wsk = global.io.double_para(3);  // sk = HkW/(k*ekW) = IC(2)
        h = global.io.double_para(4);	 // = 2*Hc/(amp*ampW)
    }
    
    else if (global.io.double_para.size() == 11) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);
        Wepsilon = global.io.double_para(2);
        Wsk = global.io.double_para(3);  
        h = global.io.double_para(4);
        
        Ux000 = global.io.double_para(5);
        Uy000 = global.io.double_para(6);
        Uz000 = global.io.double_para(7);
        
        Wx000 = global.io.double_para(8);
        Wy000 = global.io.double_para(9);
        Wz000 = global.io.double_para(10);
    }
    
	int index;
	DP temp;
	
	if (abs(h) < MYEPS) {  // zero cross helicity
		Initialize_using_energy_helicity_spectrum(U, epsilon, sk);
		Initialize_using_energy_helicity_spectrum(W, Wepsilon, Wsk);
	}
	
	else {  // nonzero cross helicity
		Model_initial_using_shell_spectrum_Pope(U.dissipation_coefficient, epsilon, Correlation::shell_ek1);
		Model_initial_using_shell_spectrum_Pope(W.dissipation_coefficient, Wepsilon, Correlation::shell_ek2);
		
		DP Kmag, ek, amp, phase1, phase2, phase3;
		DP ekW, ampW, phase1W, phase2W, phase3W;
		
        for (int lx=0; lx<=global.field.maxlx; lx++)
            for (int ly=0; ly<=global.field.maxly; ly++)
                for (int lz=0; lz<=global.field.maxlz; lz++) {
					
					Kmag = universal->Kmagnitude(lx, ly, lz);
					
					if (Kmag > MYEPS) {
						index = int(ceil(Kmag));
						
						ek = Correlation::shell_ek1(index)/ universal->Approx_number_modes_in_shell(index);
						amp = sqrt(2*ek);
						
						phase1 = 2*M_PI * SPECrand.random();
						phase2 = phase1 + M_PI/2.0;
						phase3 = asin(sk)/2.0;			// zeta
						
						Put_vector_amp_phase_comp_conj(U, lx, ly, lz, amp, phase1, phase2, phase3);
					
						// W field
						ekW = Correlation::shell_ek2(index)/ universal->Approx_number_modes_in_shell(index);
						ampW = sqrt(2*ekW);
						
						phase3W = asin(Wsk)/2.0;			// zeta_b				
						temp = h / cos(phase3-phase3W);
						// cos(phase1 - phase1W)
											 
						phase1W = phase1 - acos(temp);
						phase2W = phase1W + M_PI/2.0;
						
						Put_vector_amp_phase_comp_conj(W, lx, ly, lz, ampW, phase1W, phase2W, phase3W);
					}	// of if (kkmax > MYEPS)		
				}  // of for loop
	}
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
            real(U.cvf.V1(0,0,0)) = Ux000;
            real(U.cvf.V2(0,0,0)) = Uy000;
            real(U.cvf.V3(0,0,0)) = Uz000;
            
            real(W.cvf.V1(0,0,0)) = Wx000;
            real(W.cvf.V2(0,0,0)) = Wy000;
            real(W.cvf.V3(0,0,0)) = Wz000;
        }
        else {
            U.cvf.V1(0,0,0) = complx(Ux000,0.0);
            U.cvf.V2(0,0,0) = complx(Uy000,0.0);
            U.cvf.V3(0,0,0) = complx(Uz000,0.0);
            
            W.cvf.V1(0,0,0) = complx(Wx000,0.0);
            W.cvf.V2(0,0,0) = complx(Wy000,0.0);
            W.cvf.V3(0,0,0) = complx(Wz000,0.0);
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
	DP epsilon, sk, Wepsilon, Wsk, h, Tepsilon;
    
    DP Ux000 = 0.0;
    DP Uy000 = 0.0;
    DP Uz000 = 0.0;
    DP Wx000 = 0.0;
    DP Wy000 = 0.0;
    DP Wz000 = 0.0;
    DP T000 = 0.0;
    
    
    if (global.io.double_para.size() == 6) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);   // sk = Hk/(k*ek) = IC(2)
        
        Wepsilon = global.io.double_para(2);
        Wsk = global.io.double_para(3);  // sk = HkW/(k*ekW) = IC(2)
        h = global.io.double_para(4);	 // = 2*Hc/(amp*ampW)
        Tepsilon = global.io.double_para(5);
    }
    
    else if (global.io.double_para.size() == 13) {
        epsilon = global.io.double_para(0);
        sk = global.io.double_para(1);
        Wepsilon = global.io.double_para(2);
        Wsk = global.io.double_para(3);
        h = global.io.double_para(4);
        Tepsilon = global.io.double_para(5);
        
        Ux000 = global.io.double_para(6);
        Uy000 = global.io.double_para(7);
        Uz000 = global.io.double_para(8);
        
        Wx000 = global.io.double_para(9);
        Wy000 = global.io.double_para(10);
        Wz000 = global.io.double_para(11);
        T000 = global.io.double_para(12);
    }
    
	int index;
	DP temp;
	
	if (abs(h) < MYEPS) {  // zero cross helicity
		Initialize_using_energy_helicity_spectrum(U, epsilon, sk);
		Initialize_using_energy_helicity_spectrum(W, Wepsilon, Wsk);
	}
	
	else {  // nonzero cross helicity
		Model_initial_using_shell_spectrum_Pope(U.dissipation_coefficient, epsilon, Correlation::shell_ek1);
		Model_initial_using_shell_spectrum_Pope(W.dissipation_coefficient, Wepsilon, Correlation::shell_ek2);
		
		DP Kmag, ek, amp, phase1, phase2, phase3;
		DP ekW, ampW, phase1W, phase2W, phase3W;
		
        for (int lx=0; lx<=global.field.maxlx; lx++)
            for (int ly=0; ly<=global.field.maxly; ly++)
                for (int lz=0; lz<=global.field.maxlz; lz++) {
					
					Kmag = universal->Kmagnitude(lx, ly, lz);
					
					if (Kmag > MYEPS) {
						index = int(ceil(Kmag));
						
						ek = Correlation::shell_ek1(index)/ universal->Approx_number_modes_in_shell(index);
						amp = sqrt(2*ek);
						
						phase1 = 2*M_PI * SPECrand.random();
						phase2 = phase1 + M_PI/2.0;
						phase3 = asin(sk)/2.0;			// zeta
						
						Put_vector_amp_phase_comp_conj(U, lx, ly, lz, amp, phase1, phase2, phase3);
                        
						// W field
						ekW = Correlation::shell_ek2(index)/ universal->Approx_number_modes_in_shell(index);
						ampW = sqrt(2*ekW);
						
						phase3W = asin(Wsk)/2.0;			// zeta_b
						temp = h / cos(phase3-phase3W);
						// cos(phase1 - phase1W)
                        
						phase1W = phase1 - acos(temp);
						phase2W = phase1W + M_PI/2.0;
						
						Put_vector_amp_phase_comp_conj(W, lx, ly, lz, ampW, phase1W, phase2W, phase3W);
					}	// of if (kkmax > MYEPS)
				}  // of for loop
	}
    
    Initialize_using_energy_helicity_spectrum(T, Tepsilon);
    
    if (my_id == master_id) {
        if (basis_type == "SSS") {
            real(U.cvf.V1(0,0,0)) = Ux000;
            real(U.cvf.V2(0,0,0)) = Uy000;
            real(U.cvf.V3(0,0,0)) = Uz000;
            
            real(W.cvf.V1(0,0,0)) = Wx000;
            real(W.cvf.V2(0,0,0)) = Wy000;
            real(W.cvf.V3(0,0,0)) = Wz000;
            real(T.csf.F(0,0,0)) = T000;
        }
        else {
            U.cvf.V1(0,0,0) = complx(Ux000,0.0);
            U.cvf.V2(0,0,0) = complx(Uy000,0.0);
            U.cvf.V3(0,0,0) = complx(Uz000,0.0);
            
            W.cvf.V1(0,0,0) = complx(Wx000,0.0);
            W.cvf.V2(0,0,0) = complx(Wy000,0.0);
            W.cvf.V3(0,0,0) = complx(Wz000,0.0);
            T.csf.F(0,0,0)  = complx(T000,0.0);
        }
    }
}

//********************************** init_cond_energy.cc **************************************



  
