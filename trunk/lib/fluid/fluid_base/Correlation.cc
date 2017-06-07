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
/*! \file  nusselt.cc
 * 
 * @brief Computes nusselt number
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "Correlation.h"

Array<Real,1> Correlation::shell_ek;
Array<Real,1> Correlation::shell_dissk;

Array<Real,1> Correlation::shell_ek1;
Array<Real,1> Correlation::shell_ek2;
Array<Real,1> Correlation::shell_ek3;

Array<Real,1> Correlation::shell_dissk1;
Array<Real,1> Correlation::shell_dissk2;
Array<Real,1> Correlation::shell_dissk3;

Array<Real,2> Correlation::ring_ek;
Array<Real,2> Correlation::ring_dissk;

Array<Real,2> Correlation::ring_ek1;
Array<Real,2> Correlation::ring_ek2;
Array<Real,2> Correlation::ring_ek3;

Array<Real,2> Correlation::ring_hk1;
Array<Real,2> Correlation::ring_hk2;
Array<Real,2> Correlation::ring_hk3;

Array<Real,2> Correlation::ring_dissk1;
Array<Real,2> Correlation::ring_dissk2;
Array<Real,2> Correlation::ring_dissk3;

Array<Real,2> Correlation::cylindrical_ring_ek;
Array<Real,2> Correlation::cylindrical_ring_dissk;

Array<Real,2> Correlation::cylindrical_ring_ek1;
Array<Real,2> Correlation::cylindrical_ring_ek2;

Array<Real,2> Correlation::cylindrical_ring_dissk1;
Array<Real,2> Correlation::cylindrical_ring_dissk2;


void Correlation::Initialize()
{
	shell_ek1.resize(global.spectrum.shell.no_shells);
	shell_ek2.resize(global.spectrum.shell.no_shells);
	shell_ek3.resize(global.spectrum.shell.no_shells);

	if (global.spectrum.shell.turnon){
		shell_ek.resize(global.spectrum.shell.no_shells);
		shell_dissk.resize(global.spectrum.shell.no_shells);

		shell_dissk1.resize(global.spectrum.shell.no_shells);
		shell_dissk2.resize(global.spectrum.shell.no_shells);
		shell_dissk3.resize(global.spectrum.shell.no_shells);
	}

	if (global.spectrum.ring.turnon){
		ring_ek.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
		ring_dissk.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);

		ring_ek1.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
		ring_ek2.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
		ring_ek3.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);

		ring_hk1.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
		ring_hk2.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
		ring_hk3.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);

		ring_dissk1.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
		ring_dissk2.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
		ring_dissk3.resize(global.spectrum.ring.no_shells, global.spectrum.ring.no_sectors+1);
	}

	if (global.spectrum.cylindrical_ring.turnon) {
			cylindrical_ring_ek.resize(global.spectrum.cylindrical_ring.no_shells, global.spectrum.cylindrical_ring.no_slabs+1);
			cylindrical_ring_dissk.resize(global.spectrum.cylindrical_ring.no_shells, global.spectrum.cylindrical_ring.no_slabs+1);

			cylindrical_ring_ek1.resize(global.spectrum.cylindrical_ring.no_shells, global.spectrum.cylindrical_ring.no_slabs+1);
			cylindrical_ring_ek2.resize(global.spectrum.cylindrical_ring.no_shells, global.spectrum.cylindrical_ring.no_slabs+1);

			cylindrical_ring_dissk1.resize(global.spectrum.cylindrical_ring.no_shells, global.spectrum.cylindrical_ring.no_slabs+1);
			cylindrical_ring_dissk2.resize(global.spectrum.cylindrical_ring.no_shells, global.spectrum.cylindrical_ring.no_slabs+1);
	}
}

//******************************************************************************
/** @brief Computes Nusselt number
 *
 *  @param IncSF T
 *  @param Ra  Rayleigh number
 *  @param Pr  Prandtl number
 *  @param  string Pr_option, string RB_Uscaling
 *
 * @return \f$ Nu = 1 + \sum u_z T \f$
 * @return For Pr=0, Nu = 1, but the function returns \f$ Nu = \sum u_z T \f$
 */
Real Correlation::Get_Nusselt_no(FluidVF& U, FluidSF& T)
{
    if (!(basis_type.find("Ch") != string::npos)) {
        // Actually Nu = 1 for Pr=0.  Here we just report the product W*theta.
        if (global.PHYSICS.Pr_option == "PRZERO")
            return ( 2*universal->Get_total_energy(U.cvf.V1, T.csf.F) );
        
        
        else if (global.PHYSICS.Pr_option == "PRLARGE") {
            if (global.PHYSICS.Uscaling == "USMALL")
                return ( 1 + 2*universal->Get_total_energy(U.cvf.V1, T.csf.F) );
            
            else if (global.PHYSICS.Uscaling == "ULARGE")
                return ( 1 + 2*sqrt(global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)
                        * universal->Get_total_energy(U.cvf.V1, T.csf.F) );
        }
        
        else if (global.PHYSICS.Pr_option == "PRSMALL") {
            if (global.PHYSICS.Uscaling == "USMALL")
                return ( 1 + 2*pow2(global.PHYSICS.Prandtl)
                        * universal->Get_total_energy(U.cvf.V1, T.csf.F) );
            
            else if (global.PHYSICS.Uscaling == "ULARGE")
                return ( 1 + 2*global.PHYSICS.Prandtl*sqrt(global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)
                        * universal->Get_total_energy(U.cvf.V1, T.csf.F) );
        }
        
        else if (global.PHYSICS.Pr_option == "PRINFTY") {
            if (global.PHYSICS.Uscaling == "USMALL")
                return ( 1 + 2*universal->Get_total_energy(U.cvf.V1, T.csf.F) );
            
            else if (global.PHYSICS.Uscaling == "ULARGE")
                return ( 1 + 2* sqrt(global.PHYSICS.Rayleigh) *universal->Get_total_energy(U.cvf.V1, T.csf.F) );
        }
    }
    
    // Non-chebyshev
    
    else {
		// Actually Nu = 1 for Pr=0.  Here we just report the product W*theta.
        if (global.PHYSICS.Pr_option == "PRZERO")
            return ( 4* universal->Get_total_energy_real_space(U.rvf.V1r, T.rsf.Fr) );
        
        
        else if (global.PHYSICS.Pr_option == "PRLARGE") {
            if (global.PHYSICS.Uscaling == "USMALL")
                return ( 1 + 4* universal->Get_total_energy_real_space(U.rvf.V1r, T.rsf.Fr) );
            
            else if (global.PHYSICS.Uscaling == "ULARGE")
                return ( 1 + 2*sqrt(global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)
                        * universal->Get_total_energy_real_space(U.rvf.V1r, T.rsf.Fr) );
        }
        
        else if (global.PHYSICS.Pr_option == "PRSMALL") {
            if (global.PHYSICS.Uscaling == "USMALL")
                return ( 1 + 2*pow2(global.PHYSICS.Prandtl)
                        * universal->Get_total_energy_real_space(U.rvf.V1r, T.rsf.Fr) );
            
            else if (global.PHYSICS.Uscaling == "ULARGE")
                return ( 1 + 2*global.PHYSICS.Prandtl*sqrt(global.PHYSICS.Rayleigh*global.PHYSICS.Prandtl)
                        * universal->Get_total_energy_real_space(U.rvf.V1r, T.rsf.Fr) );
        }
        
        else if (global.PHYSICS.Pr_option == "PRINFTY") {
            if (global.PHYSICS.Uscaling == "USMALL")
                return ( 1 + 2*universal->Get_total_energy_real_space(U.rvf.V1r, T.rsf.Fr) );
            
            else if (global.PHYSICS.Uscaling == "ULARGE")
                return ( 1 + 2* sqrt(global.PHYSICS.Rayleigh) *universal->Get_total_energy_real_space(U.rvf.V1r, T.rsf.Fr) );
        }
    }
 	
 	return 0; //To avoid compiler complain   
}

/** @brief Computes cross helicity V.W/2.
 *
 *  @param IncVF W
 *
 * @return \f$ Hc = \sum V.W/2 \f$ <BR>
 */
Real Correlation::Get_cross_helicity(FluidVF& U, FluidVF& W)
{
	return (universal->Get_total_energy(U.cvf.V1, W.cvf.V1)
			+ universal->Get_total_energy(U.cvf.V2, W.cvf.V2)
			+ universal->Get_total_energy(U.cvf.V3, W.cvf.V3));
}

//******************************************************************************



//******************************************************************************
//******** Shell spectrum
void Correlation::Compute_shell_spectrum(FluidVF& U)
{
	
	universal->Compute_shell_spectrum(U.cvf.V1, 0, shell_ek1);
	
	if (!global.program.two_dimension)
		universal->Compute_shell_spectrum(U.cvf.V2, 0, shell_ek2);
	
	universal->Compute_shell_spectrum(U.cvf.V3, 0, shell_ek3);
	
	
	universal->Compute_shell_spectrum(U.cvf.V1, 2, shell_dissk1);
	
	if (!global.program.two_dimension)
		universal->Compute_shell_spectrum(U.cvf.V2, 2, shell_dissk2);
	
	universal->Compute_shell_spectrum(U.cvf.V3, 2, shell_dissk3);
	
	shell_dissk1 *= TWO;
	shell_dissk2 *= TWO;
	shell_dissk3 *= TWO;
}


void Correlation::Compute_shell_spectrum(FluidVF& U, FluidVF& W)
{
	universal->Compute_shell_spectrum(U.cvf.V1, W.cvf.V1, 0, shell_ek1);
	
	if (!global.program.two_dimension)
		universal->Compute_shell_spectrum(U.cvf.V2, W.cvf.V2, 0, shell_ek2);
	
	universal->Compute_shell_spectrum(U.cvf.V3, W.cvf.V3, 0, shell_ek3);
}


void Correlation::Compute_shell_spectrum(FluidVF& U, FluidSF& T)
{
	
	universal->Compute_shell_spectrum(U.cvf.V1, T.csf.F, 0, shell_ek1);
    
	if (!global.program.two_dimension)
		universal->Compute_shell_spectrum(U.cvf.V2, T.csf.F, 0, shell_ek2);
    
	universal->Compute_shell_spectrum(U.cvf.V3, T.csf.F, 0, shell_ek3);
}


//**** ring
void Correlation::Compute_ring_spectrum(FluidVF& U)
{
		
	universal->Compute_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, 0, ring_ek1, ring_ek2, ring_ek3);
	
	universal->Compute_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, 2, ring_dissk1, ring_dissk2, ring_dissk3);
	
	ring_dissk1 *= TWO;
	ring_dissk2 *= TWO;
	ring_dissk3 *= TWO;
}

void Correlation::Compute_ring_spectrum(FluidVF& U, FluidVF& W)
{		
	universal->Compute_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, W.cvf.V1, W.cvf.V2, W.cvf.V3, 0, ring_ek1, ring_ek2, ring_ek3);
		
}


void Correlation::Compute_ring_spectrum(FluidVF&U, FluidSF& T)
{
	
	universal->Compute_ring_spectrum(U.cvf.V1, T.csf.F, 0, ring_ek1);
	
	universal->Compute_ring_spectrum(U.cvf.V2, T.csf.F, 0, ring_ek2);
	
	universal->Compute_ring_spectrum(U.cvf.V3, T.csf.F, 0, ring_ek3);
}


//Helical_ring

void Correlation::Compute_helical_ring_spectrum(FluidVF& U, FluidVF& helicalU)
{
		
	universal->Compute_helical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, ring_hk1, ring_hk2, ring_hk3);
	
	universal->Compute_helical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 2, ring_dissk1, ring_dissk2, ring_dissk3);
	
	ring_dissk1 *= TWO;
	ring_dissk2 *= TWO;
	ring_dissk3 *= TWO;
}



//***cylindrical_ring
void Correlation::Compute_cylindrical_ring_spectrum(FluidVF& U)
{
		
	universal->Compute_cylindrical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, 0, cylindrical_ring_ek1, cylindrical_ring_ek2);
	
	universal->Compute_cylindrical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, 2, cylindrical_ring_dissk1, cylindrical_ring_dissk2);
}


void Correlation::Compute_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W)
{
		
	universal->Compute_cylindrical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, W.cvf.V1, W.cvf.V2, W.cvf.V3, 0, cylindrical_ring_ek1, cylindrical_ring_ek2);
		
}

void Correlation::Compute_cylindrical_ring_spectrum(FluidVF& U, FluidSF& T)
{
	universal->Compute_cylindrical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, T.csf.F, 0, cylindrical_ring_ek1, cylindrical_ring_ek2);
}



//*** helicity


void Correlation::Compute_shell_spectrum_helicity(FluidVF& U)
{
	
	universal->Compute_shell_spectrum_helicity(U.cvf.V1, U.cvf.V2, U.cvf.V3, shell_ek1, shell_ek2, shell_ek3);
}

void Correlation::Compute_shell_spectrum_helicity2(FluidVF& U)
{
  
  universal->Compute_shell_spectrum_helicity2(U.cvf.V1, U.cvf.V2, U.cvf.V3, shell_ek1, shell_ek2, shell_ek3);
}

void Correlation::Compute_ring_spectrum_helicity(FluidVF& U)
{
	universal->Compute_ring_spectrum_helicity(U.cvf.V1, U.cvf.V2, U.cvf.V3, ring_ek1, ring_ek2, ring_ek3);
}


void Correlation::Compute_cylindrical_ring_spectrum_helicity(FluidVF& U)
{
	
	universal->Compute_cylindrical_ring_spectrum_helicity(U.cvf.V1, U.cvf.V2, U.cvf.V3, cylindrical_ring_ek1, cylindrical_ring_ek2);
}


// For scalar fields
void Correlation::Compute_shell_spectrum(FluidSF& T)
{
	universal->Compute_shell_spectrum(T.csf.F, 0, shell_ek);
	
	universal->Compute_shell_spectrum(T.csf.F, 2, shell_dissk);
	
	shell_dissk *= TWO;
}


void Correlation::Compute_shell_spectrum(FluidSF& T1, FluidSF& T2)
{
	universal->Compute_shell_spectrum(T1.csf.F, T2.csf.F, 0, shell_ek);
}


void Correlation::Compute_ring_spectrum(FluidSF& T)
{
	
	universal->Compute_ring_spectrum(T.csf.F, 0, ring_ek);
	
	universal->Compute_ring_spectrum(T.csf.F, 2, ring_dissk);
	
	ring_dissk *= TWO;
}

void Correlation::Compute_ring_spectrum(FluidSF& T1, FluidSF& T2)
{
	
	universal->Compute_ring_spectrum(T1.csf.F, T2.csf.F, 0, ring_ek);
}


void Correlation::Compute_cylindrical_ring_spectrum(FluidSF& T)
{
	
	universal->Compute_cylindrical_ring_spectrum(T.csf.F, 0, cylindrical_ring_ek);
	
	universal->Compute_cylindrical_ring_spectrum(T.csf.F, 2, cylindrical_ring_dissk);
}


void Correlation::Compute_cylindrical_ring_spectrum(FluidSF& T1, FluidSF& T2)
{

	universal->Compute_cylindrical_ring_spectrum(T1.csf.F, T2.csf.F, 0, cylindrical_ring_ek);
}

// Force spectrum

//******************************************************************************

/** @brief Computes Force spectrum (Force(K).V(K)^*), no factor of 1/2
 *
 * @return \f$ (*shell_spectrum_force_Vk) = \Re(F(K) \cdot V(K)^*)\f$ <BR>
 */
void Correlation::Compute_force_shell_spectrum(FluidVF& U)
{
	universal->Compute_shell_spectrum(U.cvf.V1, U.Force1, 0, shell_ek1);
	
	if (!global.program.two_dimension)
		universal->Compute_shell_spectrum(U.cvf.V2, U.Force2, 0, shell_ek2);
	
	universal->Compute_shell_spectrum(U.cvf.V3, U.Force3, 0, shell_ek3);
	
	shell_ek1 *= TWO;
	shell_ek2 *= TWO;
	shell_ek3 *= TWO;
}


void Correlation::Compute_force_ring_spectrum(FluidVF& U)
{
	universal->Compute_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, U.Force1, U.Force2, U.Force3, 0, ring_ek1, ring_ek2, ring_ek3);
	
	ring_ek1 *= TWO;
	ring_ek2 *= TWO;
	ring_ek3 *= TWO;
}

void Correlation::Compute_force_cylindrical_ring_spectrum(FluidVF& U)
{
	
	universal->Compute_cylindrical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, U.Force1, U.Force2, U.Force3, 0, cylindrical_ring_ek1, cylindrical_ring_ek2);
	
	cylindrical_ring_ek1 *= TWO;
	cylindrical_ring_ek2 *= TWO;
}



void Correlation::Compute_force_shell_spectrum(FluidSF& T)
{
	universal->Compute_shell_spectrum(T.csf.F, T.Force, 0, shell_ek);
	
	shell_ek *= TWO;
}

void Correlation::Compute_force_ring_spectrum(FluidSF& T)
{
	
	universal->Compute_ring_spectrum(T.csf.F, T.Force, 0, ring_ek);
	
	ring_ek *= TWO;
}

void Correlation::Compute_force_cylindrical_ring_spectrum(FluidSF& T)
{
	
	universal->Compute_cylindrical_ring_spectrum(T.csf.F, T.Force, 0, cylindrical_ring_ek);
	
	cylindrical_ring_ek *= TWO;
	
}


//******************************************************************************

/** @brief Computes Tk spectrum (nlin(K).V(K)^*), no factor of 1/2
 *
 * @return \f$ (*shell_spectrum_force_Vk) = \Re(F(K) \cdot V(K)^*)\f$ <BR>
 */
void Correlation::Compute_Tk_shell_spectrum(FluidVF& U)
{
	
	universal->Compute_shell_spectrum(U.cvf.V1, U.nlin1, 0, shell_ek1);
    
	if (!global.program.two_dimension)
		universal->Compute_shell_spectrum(U.cvf.V2, U.nlin2, 0, shell_ek2);
    
	universal->Compute_shell_spectrum(U.cvf.V3, U.nlin3, 0, shell_ek3);
    
	shell_ek1 *= TWO;
	shell_ek2 *= TWO;
	shell_ek3 *= TWO;
}


void Correlation::Compute_Tk_ring_spectrum(FluidVF& U)
{
	universal->Compute_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, U.nlin1, U.nlin2, U.nlin3, 0, ring_ek1, ring_ek2, ring_ek3);
	
	ring_ek1 *= TWO;
	ring_ek2 *= TWO;
	ring_ek3 *= TWO;
	
}

void Correlation::Compute_Tk_cylindrical_ring_spectrum(FluidVF& U)
{
	universal->Compute_cylindrical_ring_spectrum(U.cvf.V1, U.cvf.V2, U.cvf.V3, U.nlin1, U.nlin2, U.nlin3, 0, cylindrical_ring_ek1, cylindrical_ring_ek2);
	
	cylindrical_ring_ek1 *= TWO;
	cylindrical_ring_ek2 *= TWO;
}


void Correlation::Compute_Tk_shell_spectrum(FluidSF& T)
{
	universal->Compute_shell_spectrum(T.csf.F, T.nlin, 0, shell_ek);
	
	shell_ek *= TWO;
}

void Correlation::Compute_Tk_ring_spectrum(FluidSF& T)
{
	
	universal->Compute_ring_spectrum(T.csf.F, T.nlin, 0, ring_ek);
	
	ring_ek *= TWO;
}

void Correlation::Compute_Tk_cylindrical_ring_spectrum(FluidSF& T)
{
	
	universal->Compute_cylindrical_ring_spectrum(T.csf.F, T.nlin, 0, cylindrical_ring_ek);
	
	cylindrical_ring_ek *= TWO;
	
}


// For pressure

//***********************************  End of nusselt.cc **************************************



