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

extern Uniform<Real> SPECrand;


//******************************************************************************
/** @brief Put a velocity vector at (lx,ly,lz) given its amp and phases
 *			and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param lx, ly, lz
 *	@param amp		Amplitude of the vector.
 *	@param phase1, phase2, phase3	Phases
 *
 */
void FORCE::Put_force_amp_phase_comp_conj(FluidVF U, int lx, int ly, int lz,  Real amp, Real phase1, Real phase2, Real phase3, bool add_flag)
{
	Complex fperp1, fperp2;
	Complex fpll, fh1, fh2;
	Real theta = 0;
	Real phi = 0;
	Real kkmag, kkperp;
	
	TinyVector<Complex,3> F_FOUR, Flocal_complex;
	TinyVector<Real,3> Flocal_real;

	if (Ny >1) {
		fperp1 = amp * exp(I * phase1) * cos(phase3);
		fperp2 = amp * exp(I * phase2) * sin(phase3);
		
		theta = universal->AnisKvect_polar_angle(lx, ly, lz);
		phi = universal->AnisKvect_azimuthal_angle(lx, ly, lz);
	}
	
	else {
		fperp1 = 0.0;
		fperp2 = amp * exp(I * phase2);
		
		theta = universal->AnisKvect_polar_angle(lx, ly, lz);
		int kx = universal->Get_kx(lx);
		phi	  = (kx >=0) ? 0 : M_PI;  
        // phi = 0 or pi depending i1>0 or i1 < 0.
	}
	
	kkmag  = universal->Kmagnitude(lx, ly, lz);
	kkperp = universal->AnisKperp(lx, ly, lz);
	
	if (kkmag > MYEPS) {
		if ( kkperp > MYEPS) {	
			fpll = -fperp2 * sin(theta);
			fh1 =  fperp2 * cos(theta)*cos(phi) + fperp1 * sin(phi); 
			fh2 =  fperp2 * cos(theta)*sin(phi) - fperp1 * cos(phi); 	
		}
		
		else {	// k along pll axis.  V on the perp plane.	
				fpll = 0.0;
				fh1 = fperp2;
				fh2 = -fperp1;
		}
	}
	
	else {
        if (my_id == master_id) { // origin lies in the master node
            U.Force1(ly,lz,lx) = 0;
            U.Force2(ly,lz,lx) = 0;
            U.Force3(ly,lz,lx) = 0;
            return;
        }
	}
	
	if (global.field.anisotropy_dirn == 1)
		F_FOUR = fpll, fh1, fh2;
		//F_FOUR = fpll, fh2, fh1; //Patch for 2D
	
	else if (global.field.anisotropy_dirn == 2)
		F_FOUR = fh2, fpll, fh1;
	
	else if (global.field.anisotropy_dirn == 3)
		F_FOUR = fh1, fh2, fpll;
        
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
void FORCE::Put_force_amp_phase_comp_conj(FluidSF& T, int lx, int ly, int lz, Real amp, Real phase, bool add_flag)
{
	
	Real Glocal_real;
	Complex Glocal_complex;
	
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

	
