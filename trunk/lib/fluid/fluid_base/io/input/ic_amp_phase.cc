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


#include "FluidIO.h"

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

void FluidIO::Put_vector_amp_phase_comp_conj(FluidVF U, int lx, int ly, int lz,  Real amp, Real phase1, Real phase2, Real phase3)
{
	Complex uperp1, uperp2;
	Complex vpll, vh1, vh2;
	Real theta = 0;
	Real phi = 0;
	Real kkmag, kkperp;
	
	TinyVector<Complex,3> VFour, Vlocal_complex;
	TinyVector<Real,3> Vlocal_real;

	if (Ny >1) {
		uperp1 = amp * exp(I * phase1) * cos(phase3);
		uperp2 = amp * exp(I * phase2) * sin(phase3);
		
		theta = universal->AnisKvect_polar_angle(lx, ly, lz);
		phi = universal->AnisKvect_azimuthal_angle(lx, ly, lz);
	}
	
	else {
		uperp1 = 0.0;
		uperp2 = amp * exp(I * phase2);
		
		theta = universal->AnisKvect_polar_angle(lx, ly, lz);
		int kx = universal->Get_kx(lx);
		phi	  = (kx >=0) ? 0 : M_PI;  
        // phi = 0 or pi depending i1>0 or i1 < 0.
	}
	
	kkmag  = universal->Kmagnitude(lx, ly, lz);
	kkperp = universal->AnisKperp(lx, ly, lz);
	
	if (kkmag > MYEPS) {
		if ( kkperp > MYEPS) {	
			vpll = -uperp2 * sin(theta);
			vh1 =  uperp2 * cos(theta)*cos(phi) + uperp1 * sin(phi); 
			vh2 =  uperp2 * cos(theta)*sin(phi) - uperp1 * cos(phi); 	
		}
		
		else {	// k along pll axis.  V on the perp plane.	
				vpll = 0.0;
				vh1 = uperp2;
				vh2 = -uperp1;
		}
	}
	
	else {
        if (my_id == master_id) { // origin lies in the master node
            U.cvf.V1(ly,lz,lx) = 0;
            U.cvf.V2(ly,lz,lx) = 0;
            U.cvf.V3(ly,lz,lx) = 0;
            return;
        }
	}
	
	if (global.field.anisotropy_dirn == 1)
		VFour = vpll, vh1, vh2;
		//VFour = vpll, vh2, vh1; //Patch for 2D

	else if (global.field.anisotropy_dirn == 2)
		VFour = vh2, vpll, vh1;
	
	else if (global.field.anisotropy_dirn == 3)
		VFour = vh1, vh2, vpll;
        
	if (basis_type == "FFF" || basis_type == "FFFW") {
		Vlocal_complex = VFour;
		universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_complex);
	}
	
	else if (basis_type == "SSS") {
        Convert_from_Fourier_space(VFour, Vlocal_real);
        universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_real);
	}
	
	else {
		Convert_from_Fourier_space(VFour, Vlocal_complex);
		universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_complex); 
	}
}
//New helicity spectrum function --> Satyajit
void FluidIO::Put_vector_amp_phase_comp_conj_helicity(FluidVF U, int lx, int ly, int lz,  Real ek, Real hk_by_k, Real phase_plus, Real phase_minus)
{
  Complex u_plus, u_minus;
  Complex u1, u2;
  Complex vpll, vh1, vh2;
  Real theta = 0;
  Real phi = 0;
  Real kkmag, kkperp;
  
  TinyVector<Complex,3> VFour, Vlocal_complex;
  TinyVector<Real,3> Vlocal_real;
  
  if (Ny >1) {
    u_plus = (sqrt(abs((ek + hk_by_k) * 0.5))) * exp(I*phase_plus);
    u_minus = (sqrt(abs((ek - hk_by_k) * 0.5))) * exp(I*phase_minus);
    universal->Helical_to_Craya(u_plus,u_minus,u1,u2);
  }
  
  else {
    cout<<"Helicity 2D yet to be implemented"<<endl;
    /*uperp1 = 0.0;
    uperp2 = amp * exp(I * phase2);
    
    theta = universal->AnisKvect_polar_angle(lx, ly, lz);
    int kx = universal->Get_kx(lx);
    phi	  = (kx >=0) ? 0 : M_PI;
    // phi = 0 or pi depending i1>0 or i1 < 0.*/
  }
  
  kkmag  = universal->Kmagnitude(lx, ly, lz);
  kkperp = universal->AnisKperp(lx, ly, lz);
  
  if (kkmag > MYEPS) {
    if ( kkperp > MYEPS) {
      universal->Craya_to_cartesian(lx,ly,lz,u1,u2,vpll,vh1,vh2);
    }
    
    else {	// k along pll axis.  V on the perp plane.
      vpll = 0.0;
      vh1 = u2;
      vh2 = -u1;
    }
  }
  
  else {
    if (my_id == master_id) { // origin lies in the master node
      U.cvf.V1(ly,lz,lx) = 0;
      U.cvf.V2(ly,lz,lx) = 0;
      U.cvf.V3(ly,lz,lx) = 0;
      return;
    }
  }
  
  if (global.field.anisotropy_dirn == 1)
    VFour = vpll, vh1, vh2;
		//VFour = vpll, vh2, vh1; //Patch for 2D
  
  else if (global.field.anisotropy_dirn == 2)
    VFour = vh2, vpll, vh1;
  
  else if (global.field.anisotropy_dirn == 3)
    VFour = vh1, vh2, vpll;
  
  if (basis_type == "FFF" || basis_type == "FFFW") {
    Vlocal_complex = VFour;
    universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_complex);
  }
  
  else if (basis_type == "SSS") {
    Convert_from_Fourier_space(VFour, Vlocal_real);
    universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_real);
  }
  
  else {
    Convert_from_Fourier_space(VFour, Vlocal_complex);
    universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_complex);
  }
}


//New cross helicity spectrum function --> Satyajit
void FluidIO::Put_vector_amp_phase_comp_conj_cross_helicity(FluidVF U, int lx, int ly, int lz,  Real u_plus_mod, Real u_minus_mod, Real phase_plus, Real phase_minus)
{
  Complex u_plus, u_minus;
  Complex u1, u2;
  Complex vpll, vh1, vh2;
  Real theta = 0;
  Real phi = 0;
  Real kkmag, kkperp;
  
  TinyVector<Complex,3> VFour, Vlocal_complex;
  TinyVector<Real,3> Vlocal_real;
  
  if (Ny >1) {
      u_plus = u_plus_mod * exp(I*phase_plus);
      u_minus = u_minus_mod * exp(I*phase_minus);
      universal->Helical_to_Craya(u_plus,u_minus,u1,u2);
  }
  
  else {
    cout<<"Helicity 2D yet to be implemented"<<endl;
    /*uperp1 = 0.0;
    uperp2 = amp * exp(I * phase2);
    
    theta = universal->AnisKvect_polar_angle(lx, ly, lz);
    int kx = universal->Get_kx(lx);
    phi	  = (kx >=0) ? 0 : M_PI;
    // phi = 0 or pi depending i1>0 or i1 < 0.*/
  }
  
  kkmag  = universal->Kmagnitude(lx, ly, lz);
  kkperp = universal->AnisKperp(lx, ly, lz);
  
  if (kkmag > MYEPS) {
    if ( kkperp > MYEPS) {
      universal->Craya_to_cartesian(lx,ly,lz,u1,u2,vpll,vh1,vh2);
      
    }
    
    else {	// k along pll axis.  V on the perp plane.
      vpll = 0.0;
      vh1 = u2;
      vh2 = -u1;
    }
  }
  
  else {
    if (my_id == master_id) { // origin lies in the master node
      U.cvf.V1(ly,lz,lx) = 0;
      U.cvf.V2(ly,lz,lx) = 0;
      U.cvf.V3(ly,lz,lx) = 0;
      return;
    }
  }
  
  if (global.field.anisotropy_dirn == 1)
    VFour = vpll, vh1, vh2;
		//VFour = vpll, vh2, vh1; //Patch for 2D
  
  else if (global.field.anisotropy_dirn == 2)
    VFour = vh2, vpll, vh1;
  
  else if (global.field.anisotropy_dirn == 3)
    VFour = vh1, vh2, vpll;
  
  if (basis_type == "FFF" || basis_type == "FFFW") {
    Vlocal_complex = VFour;
    universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_complex);
  }
  
  else if (basis_type == "SSS") {
    Convert_from_Fourier_space(VFour, Vlocal_real);
    universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_real);
  }
  
  else {
    Convert_from_Fourier_space(VFour, Vlocal_complex);
    universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, Vlocal_complex);
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
void FluidIO::Put_scalar_amp_phase_comp_conj(FluidSF& T, int lx, int ly, int lz, Real amp, Real phase)
{
	
	Real Glocal_real;
	Complex Glocal_complex;
	
	Complex G_Four = amp * exp(I * phase);
	
	if (basis_type == "FFF" || basis_type == "FFFW") {
		Glocal_complex = G_Four;
		universal->Assign_local_spectral_field(lx, ly, lz, T.csf.F, Glocal_complex);
    }
	
	else if (basis_type == "SSS") {
		Convert_from_Fourier_space(G_Four, Glocal_real);
		universal->Assign_local_spectral_field(lx, ly, lz, T.csf.F, Glocal_real);
	}
	
	else {
		Convert_from_Fourier_space(G_Four, Glocal_complex);	
		universal->Assign_local_spectral_field(lx, ly, lz, T.csf.F, Glocal_complex);
	}
}


//********************************** ic_random.cc *********************************************

	
