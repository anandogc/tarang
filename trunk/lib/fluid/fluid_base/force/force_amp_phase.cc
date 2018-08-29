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

/*! \file  ic_amp_phase.cc
 * 
 * @brief Initial conditions where fields are chosen randomly.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "FORCE.h"


//******************************************************************************
/** @brief Put a Force vector at (lx,0,lz) given the amp and phase of F_c1 (Craya-Herring)
 *  For 2D
 *
 *  @param lx, lz  (ly=0): location where the F is to be assigned.
 *	@param amp		Amplitude of the vector F_c1 (Craya-Herring).
 *	@param phase Phase of the vector F_c1 (Craya-Herring)
 *  @param add_flag: if yes,it adds the vector to existing F, otherwise assigns F at (lx,ly,lz
 *
 *  @detail If k=(0,0,0), set Force=(0,0,0)
 *  @detail For 2D, (xz) plane of the code is to be treaed as (x,y) of Craya-Herring decomposition.  For 2D, the anisotropy direction is ignored, and Fx = F1 sin(phi), Fz= -F1 cos(phi).
 */
void FORCE::Put_or_add_force_vector(FluidVF U, int lx, int lz,  Real amp, Real phase, bool add_flag)
{
    Complex F_ch1; // Amplitude along Craya-herring basis vector e1.
    Real kx, kz, Kmag, phi;
    
    TinyVector<Complex,3> F_FOUR, Flocal_complex;
    TinyVector<Real,3> Flocal_real;
    
    // k=(0,0,0) case first
    Kmag  = universal->Kmagnitude(lx, 0, lz);
    if (Kmag < MYEPS)
        if (my_id == master_id) { // origin lies in the master node
            U.Force1(0,0,0) = 0;
            U.Force3(0,0,0) = 0;
            return;
        }
    
    // k \ne (0,0,0)
    F_ch1 = amp * exp(I*phase);
    kx = universal->Get_kx(lx) * kfactor[1];
    kz = universal->Get_kz(lz) * kfactor[3];
    phi = Get_azimuthal_angle(kx, kz);
    
    F_FOUR = F_ch1*sin(phi), 0, -F_ch1*cos(phi);
    
    if (basis_type == "FFF" || basis_type == "FFFW") {
        Flocal_complex = F_FOUR;
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, 0, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
        else
            universal->Add_local_spectral_field(lx, 0, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
    }
    
    else if (basis_type == "SSS") {
        Convert_from_Fourier_space(F_FOUR, Flocal_real);
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, 0, lz, U.Force1, U.Force2, U.Force3, Flocal_real);
        else
            universal->Add_local_spectral_field(lx, 0, lz, U.Force1, U.Force2, U.Force3, Flocal_real);
    }
    
    else {
        Convert_from_Fourier_space(F_FOUR, Flocal_complex);
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, 0, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
        else
            universal->Add_local_spectral_field(lx, 0, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
    }
}


//******************************************************************************
/** @brief Put a Force vector at (lx,ly,lz) given the amps and phases of F+ and F- (helical basis)
*
*   For 3D and 2D3C
*
*   @param lx, ly, lz: : location where the F is to be assigned.
*	@param amp_plus, amp_minus:		Amplitudes of the vectors F+, F- (helical basis).
*	@param phase_plus, phase_minus:  Phases of the vectors F_c1 F+, F- (helical basis).
*   @param add_flag: if yes,it adds the vector to existing F, otherwise assigns F at (lx,ly,lz)
*/

void FORCE::Put_or_add_force_vector(FluidVF U, int lx, int ly, int lz,  Real amp_plus, Real amp_minus, Real phase_plus, Real phase_minus, bool add_flag)
{
    Complex F_plus, F_minus;  // Amplitudes along helical basis
    Complex F_ch1, F_ch2;     // Amplitudes along Craya-Herring basis
	Complex fpll, fh1, fh2;
	Real kx,kz,Kmag, phi;
	
	TinyVector<Complex,3> F_FOUR, Flocal_complex;
	TinyVector<Real,3> Flocal_real;
    
    // k=(0,0,0) case first
    Kmag  = universal->Kmagnitude(lx, ly, lz);
    if (Kmag < MYEPS)
        if (my_id == master_id) { // origin lies in the master node
            U.Force1(0,0,0) = 0;
            U.Force2(0,0,0) = 0;
            U.Force3(0,0,0) = 0;
            return;
        }

    
    // k \ne (0,0,0)
    F_plus = amp_plus * exp(I*phase_plus);
    F_minus = amp_minus * exp(I*phase_minus);
    universal->Helical_to_Craya(F_plus, F_minus, F_ch1, F_ch2);
    // Gets components F_ch1, F_ch2 along the Craya-Herring basis.
    
    if (Ny == 1) {  // 2D3C
        kx = universal->Get_kx(lx) * kfactor[1];
        kz = universal->Get_kz(lz) * kfactor[3];
        phi = Get_azimuthal_angle(kx, kz);
        
        F_FOUR = F_ch1*sin(phi), F_ch2, -F_ch1*cos(phi);
    }
    else { // 3D
        universal->Craya_to_cartesian(lx, ly, lz, F_ch1, F_ch2, fpll, fh1, fh2);
        // Gets components pll, fh1, fh2 in Cartesian basis.
        
        if (global.field.anisotropy_dirn == 1)
            F_FOUR = fpll, fh1, fh2;
        
        else if (global.field.anisotropy_dirn == 2)
            F_FOUR = fh2, fpll, fh1;
        
        else if (global.field.anisotropy_dirn == 3)
            F_FOUR = fh1, fh2, fpll;
    }
		
        
	if (basis_type == "FFF" || basis_type == "FFFW") {
		Flocal_complex = F_FOUR;
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
        else
            universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
	}
	
	else if (basis_type == "SSS") {
        Convert_from_Fourier_space(F_FOUR, Flocal_real);
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, Flocal_real);
        else
            universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, Flocal_real);
	}
	
	else {
		Convert_from_Fourier_space(F_FOUR, Flocal_complex);
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
        else
            universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, Flocal_complex);
	}
}



//*********************************************************************************************

/** @brief Put a scalar field at (lx,ly,lz) given its amplitude and phase, 
 *			and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param lx, ly, lz
 *	@param amp  Amplitude of the vector.
 *
 */
void FORCE::Put_or_add_force_scalar(FluidSF& T, int lx, int ly, int lz, Real amp, Real phase, bool add_flag)
{
	
	Real Glocal_real;
	Complex Glocal_complex;
    Real Kmag;
	
    // k=(0,0,0) case first
    Kmag  = universal->Kmagnitude(lx, ly, lz);
    if (Kmag < MYEPS)
        if (my_id == master_id) { // origin lies in the master node
            T.Force(0,0,0) = 0;
            return;
        }
    
    
    // k \ne (0,0,0)
    Complex G_Four = amp * exp(I * phase);
	
	if (basis_type == "FFF" || basis_type == "FFFW") {
		Glocal_complex = G_Four;
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, ly, lz, T.Force, Glocal_complex);
        else
            universal->Add_local_spectral_field(lx, ly, lz, T.Force, Glocal_complex);
    }
	
	else if (basis_type == "SSS") {
		Convert_from_Fourier_space(G_Four, Glocal_real);
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, ly, lz, T.Force, Glocal_real);
        else
            universal->Add_local_spectral_field(lx, ly, lz, T.Force, Glocal_real);
            
	}
	
	else {
		Convert_from_Fourier_space(G_Four, Glocal_complex);
        
        if (!add_flag)
            universal->Assign_local_spectral_field(lx, ly, lz, T.Force, Glocal_complex);
        else
            universal->Add_local_spectral_field(lx, ly, lz, T.Force, Glocal_complex);
	}
}


//********************************** ic_random.cc *********************************************

	
