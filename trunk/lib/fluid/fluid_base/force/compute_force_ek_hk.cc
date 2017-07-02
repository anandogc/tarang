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
// for global alpha, beta, set global_alpha_beta= true
// global alpha, beta: Based on suggestions by Rodion Stepanov

void FORCE::Force_energy_helicity_supply_or_level_basic(FluidVF& U, string force_type, Real inner_radius, Real outer_radius, Real para1, Real para2, bool add_flag)
{ 
	
	Real energy_supply = 0.0;
	Real epsh_by_k_epse = 0.0;
	Real energy_level = 0.0;
	Real h_by_k_E = 0.0;
    Real kinetic_helicity_supply = 0.0;
    
    bool global_alpha_beta = false;
	
	if (force_type == "ENERGY_SUPPLY") {
		energy_supply = para1;
		if (global_alpha_beta)
            kinetic_helicity_supply = para2;
        else
            epsh_by_k_epse = para2;
	}
	
	else if (force_type == "CONSTANT_ENERGY") {
		energy_level = para1;
		h_by_k_E = para2;
	}
	
	// No forcing....
	if ( (force_type == "ENERGY_SUPPLY") && (abs(energy_supply) < MYEPS) )
		return;
	
	if ( (force_type == "CONSTANT_ENERGY") && (abs(energy_level) < MYEPS) )
		return;
	
	int nf; 
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	Real modal_energy;
	Real temp, temp1, temp2, temp3;
	
	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
	Real energy_supply_per_mode, energy_per_mode;
	
	if (force_type == "ENERGY_SUPPLY")
		energy_supply_per_mode = energy_supply / nf;
	
	else if (force_type == "CONSTANT_ENERGY") 
		energy_per_mode = energy_level / nf;
	
	
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
    
    if (global_alpha_beta) {
        Real total_Ek_in_force_shell = 0;
        Real total_ksqrEk_in_force_shell = 0;

        Correlation::Compute_shell_spectrum(U);
        for (int i=int(inner_radius)+1; i<= int(outer_radius); i++) {
            total_Ek_in_force_shell += Correlation::shell_ek1(i) + Correlation::shell_ek2(i) + Correlation::shell_ek3(i);
            total_ksqrEk_in_force_shell += Correlation::shell_dissk1(i) + Correlation::shell_dissk3(i) + Correlation::shell_dissk3(i);
        }
        total_ksqrEk_in_force_shell /= TWO;   // shell_dissk1 contains 2*k^2*E(k), so div by 2
        
        Real total_Hk_in_force_shell = 0;
        Correlation::Compute_shell_spectrum_helicity(U);
        for (int i=int(inner_radius)+1; i<= int(outer_radius); i++)
            total_Hk_in_force_shell += Correlation::shell_ek1(i) + Correlation::shell_ek2(i) + Correlation::shell_ek3(i);
        
        Real denr = total_Ek_in_force_shell*total_ksqrEk_in_force_shell - my_pow(total_Hk_in_force_shell,2);
        
        alpha_k = (energy_supply*total_ksqrEk_in_force_shell - kinetic_helicity_supply*total_Hk_in_force_shell)/(2*denr);
        beta_k = (kinetic_helicity_supply*total_Ek_in_force_shell - energy_supply*total_Hk_in_force_shell)/(2*denr);
    }


	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = ky_min; ky <= ky_max; ky++)  
			for (int kz = 0; kz <= kz_max; kz++) {

				if (universal->Probe_in_me(kx,ky,kz))  {
					lx = universal->Get_lx(kx);
					ly = universal->Get_ly(ky);
					lz = universal->Get_kz(kz);
				
					Kmag = universal->Kmagnitude(lx, ly, lz);
					if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
						modal_energy = U.cvf.Modal_energy(lx, ly, lz);
						
						if (modal_energy > MYEPS) {
							
							if (force_type == "ENERGY_SUPPLY" && !global_alpha_beta) { // for local alpha-beta
								temp = energy_supply_per_mode/ (2*modal_energy);
                                sk = U.cvf.Modal_helicity(lx,ly,lz)/ (Kmag*modal_energy);
								
								if (abs(epsh_by_k_epse) < MYEPS) {  // No helicity forcing, Only Ek forcing
                                    alpha_k = temp * (1-sk*epsh_by_k_epse) / (1 - sk*sk);
                                    beta_k = temp * (epsh_by_k_epse - sk) / ((1 - sk*sk)*Kmag);
								}
								
								else { // helical									
									if (abs(sk*sk-1) > MYEPS2) {
										alpha_k = temp * (1-sk*epsh_by_k_epse) / (1 - sk*sk);
										beta_k = temp * (epsh_by_k_epse - sk) / ((1 - sk*sk)*Kmag);
									}
									
									else { // max helicity
										alpha_k = temp/2;
										beta_k = alpha_k/(sk*Kmag);
									}	
								}
							} // end of (force_type == "ENERGY_SUPPLY")
							
							else if (force_type == "CONSTANT_ENERGY") {
								temp1 = sqrt(energy_per_mode/modal_energy);
								
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
							} // end of (force_type == "CONSTANT_ENERGY")	
							
							if (force_type == "ENERGY_SUPPLY")
								Const_energy_supply_alpha_beta(U, lx, ly, lz, alpha_k, beta_k, add_flag);

							else if (force_type == "CONSTANT_ENERGY") 
								Const_energy_alpha_beta(U, lx, ly, lz, alpha_k, beta_k, add_flag);
						}		
					}
				}	//  of if (Probe_in_me())						
			}		// of for
						
}


void FORCE::Force_energy_helicity_supply_or_level_basic_assign(FluidVF& U, string force_type, Real inner_radius, Real outer_radius, Real para1, Real para2)
{
	Force_energy_helicity_supply_or_level_basic(U, force_type, inner_radius, outer_radius, para1, para2, false);
}

void FORCE::Force_energy_helicity_supply_or_level_basic_add(FluidVF& U, string force_type, Real inner_radius, Real outer_radius, Real para1, Real para2)
{
	Force_energy_helicity_supply_or_level_basic(U, force_type, inner_radius, outer_radius, para1, para2, true);
}
			
//New crosshelical forsing (15/03/2017) which wa written with Stepanov R
void FORCE::Force_energy_crosshelicity_supply_or_level_basic(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real para1, Real para2, bool add_flag)
{ 
	Real eps = para1;
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

}				


void FORCE::Force_energy_crosshelicity_supply_or_level_basic_assign(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real para1, Real para2)
{
	Force_energy_crosshelicity_supply_or_level_basic(U, W, inner_radius, outer_radius, para1, para2, false);
}

void FORCE::Force_energy_crosshelicity_supply_or_level_basic_add(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real para1, Real para2)
{
	Force_energy_crosshelicity_supply_or_level_basic(U, W, inner_radius, outer_radius, para1, para2, true);
}

//*********************************************************************************************
//
//	Scalar
//

void FORCE::Force_energy_helicity_supply_or_level_basic(FluidSF& T, string force_type, Real inner_radius, Real outer_radius, Real para, bool add_flag)
{
	
	Real energy_supply = 0.0;
	Real energy_level = 0.0;
	
	if (force_type == "ENERGY_SUPPLY")
		energy_supply = para;
	
	else if (force_type == "CONSTANT_ENERGY")
		energy_level = para;
	
	T.Force = 0.0;
	
	// No forcing....
	if ( (force_type == "ENERGY_SUPPLY") && (abs(energy_supply) < MYEPS) )
		return;
	
	if ( (force_type == "CONSTANT_ENERGY") && (abs(energy_level) < MYEPS) )
		return;
	
	int nf; 
	Real energy_supply_per_mode;
    Real energy_per_mode;
	Real modal_energy;
	
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	
	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
	
	
	if (force_type == "ENERGY_SUPPLY")
		energy_supply_per_mode = energy_supply/nf;
	
	else if (force_type == "CONSTANT_ENERGY")
		energy_per_mode = energy_level / nf;
    
	
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
	Real Kmag, alpha_k;
		
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = ky_min; ky <= ky_max; ky++)  
			for (int kz = 0; kz <= kz_max; kz++) {
				
				if (universal->Probe_in_me(kx, ky, kz)) {
					lx = universal->Get_lx(kx);
					ly = universal->Get_ly(ky);
					lz = universal->Get_kz(kz);

					Kmag = universal->Kmagnitude(lx, ly, lz);
					if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
						modal_energy = T.csf.Modal_energy(lx, ly, lz);
						
						if (modal_energy > MYEPS) {
							if (force_type == "ENERGY_SUPPLY") {
								alpha_k = energy_supply_per_mode/ (2*modal_energy);
								Const_energy_supply_alpha(T, lx, ly, lz, alpha_k, add_flag);
							}
							
							else if (force_type == "CONSTANT_ENERGY") {
								alpha_k = sqrt(energy_per_mode/modal_energy);
								Const_energy_alpha(T, lx, ly, lz, alpha_k, add_flag);
							}
							
							
							
						}
					}
				}   // of (Probe_in_me())						
			}		// of for
}	

void FORCE::Force_energy_helicity_supply_or_level_basic_assign(FluidSF& T, string force_type, Real inner_radius, Real outer_radius, Real para)
{
	Force_energy_helicity_supply_or_level_basic(T, force_type, inner_radius, outer_radius, para, false);
}

void FORCE::Force_energy_helicity_supply_or_level_basic_add(FluidSF& T, string force_type, Real inner_radius, Real outer_radius, Real para)
{
	Force_energy_helicity_supply_or_level_basic(T, force_type, inner_radius, outer_radius, para, true);
}
//*********************************************************************************************
// derived fn
void FORCE::Compute_force_const_energy_helicity_supply(FluidVF& U)
{
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_supply = global.force.double_para(2);
	Real epsh_by_k_epse = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply, epsh_by_k_epse);
	
}


//*********************************************************************************************
// derived fn
void FORCE::Compute_force_const_energy_helicity_supply(FluidVF& U, FluidSF& T)
{
	
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_supply = global.force.double_para(2);
	Real epsh_by_k_epse = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	Real energy_supply_scalar = global.force.double_para(4);
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply, epsh_by_k_epse);
	
	if (T.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(T, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply_scalar);
}

//*********************************************************************************************
// derived fn
void FORCE::Compute_force_const_energy_helicity_supply(FluidVF& U, FluidVF& W)
{
		
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_supply = global.force.double_para(2);
	Real epsh_by_k_epse = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	Real energy_supply_W = global.force.double_para(4);
	Real epsh_by_k_epse_W = global.force.double_para(5);			// W.epsh(k)/(k*W.eps(k))
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply, epsh_by_k_epse);
	
	if (W.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(W, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply_W, epsh_by_k_epse_W);
}


//*********************************************************************************************
// derived fn
void FORCE::Compute_force_crosshelicity_supply(FluidVF& U, FluidVF& W)
{
		
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real eps = global.force.double_para(2);
	Real eps_hc = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	
	if (U.force_switch) 
		Force_energy_crosshelicity_supply_or_level_basic_assign(U, W, inner_radius, outer_radius, eps, eps_hc);
}

//*********************************************************************************************
// derived fn
void FORCE::Compute_force_const_energy_helicity_supply(FluidVF& U, FluidVF& W, FluidSF& T)
{

	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_supply = global.force.double_para(2);
	Real epsh_by_k_epse = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	Real energy_supply_W = global.force.double_para(4);
	Real epsh_by_k_epse_W = global.force.double_para(5);			// W.epsh(k)/(k*W.eps(k))
	Real energy_supply_scalar = global.force.double_para(6);
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply, epsh_by_k_epse);
		
	
	if (W.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(W, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply_W, epsh_by_k_epse_W);
		
	if (T.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(T, "ENERGY_SUPPLY", inner_radius, outer_radius, energy_supply_scalar);
}

//*********************************************************************************************
void FORCE::Compute_force_const_energy_helicity(FluidVF& U)
{
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_level = global.force.double_para(2);
	Real h_by_k_E = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level, h_by_k_E);
}

//*********************************************************************************************
// derived fn
void FORCE::Compute_force_const_energy_helicity(FluidVF& U, FluidSF& T)
{
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_level = global.force.double_para(2);
	Real h_by_k_E = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	Real energy_level_scalar = global.force.double_para(4);
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level, h_by_k_E);

	
	if (T.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(T, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level_scalar);
		
}


//*********************************************************************************************
// derived fn
void FORCE::Compute_force_const_energy_helicity(FluidVF& U, FluidVF& W)
{
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_level = global.force.double_para(2);
	Real h_by_k_E = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	Real energy_level_W = global.force.double_para(4);
	Real h_by_k_E_W = global.force.double_para(5);
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level, h_by_k_E);
	
	if (W.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(W, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level_W, h_by_k_E_W);
}

//*********************************************************************************************
// derived fn
void FORCE::Compute_force_const_energy_helicity(FluidVF& U, FluidVF& W, FluidSF& T)
{
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real energy_level = global.force.double_para(2);
	Real h_by_k_E = global.force.double_para(3);				// epsh(k)/(k*eps(k))
	Real energy_level_W = global.force.double_para(4);
	Real h_by_k_E_W = global.force.double_para(5);
	Real energy_level_scalar = global.force.double_para(6);
	
	if (U.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(U, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level, h_by_k_E);
		
	
	if (W.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(W, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level_W, h_by_k_E_W);
		
	if (T.force_switch) 
		Force_energy_helicity_supply_or_level_basic_assign(T, "CONSTANT_ENERGY", inner_radius, outer_radius, energy_level_scalar);
}




//***********************  End of compute_force_ek_hk_supply.cc *******************************


