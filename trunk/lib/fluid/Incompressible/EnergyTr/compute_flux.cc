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

/*! \file  compute_flux.cc
 * 
 * @brief  Computes energy flux from inside/outside a sphere of field V/W to inside/outsude
 *			of field V/W.
 *
 *	The Giver field is filled inside/outside the sphere (region A1).  
 *	We do the same for the receiver field (region A2). <BR>
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "EnergyTr.h"

//*********************************************************************************************

// Fluid
void EnergyTr::Compute_flux(FluidVF &U)
{	
	
	flux_self = 0.0;

	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
		Fill_in_sphere(sphere_index, U);	
		
		Nlin_incompress::Compute_nlin(U, Giver);
			//	U.nlin = U.grad Giver<	
		
		flux_self(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, U);	
		// flux_self = -(U.grad U<). U> = U.nlin< . U>
	} 
}


//*********************************************************************************************

// Scalar
//

void EnergyTr::Compute_flux(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Compute_flux_scalar(U, T);
	
	else if (global.program.kind == "RBC")
		Compute_flux_RBC(U, T);
}


void EnergyTr::Compute_flux_scalar(FluidVF& U, FluidSF& T)
{

	flux_SF = 0.0;
	
	Compute_flux(U);
	// flux_self = (U.grad U<). U>	
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {		
		Fill_in_sphere(sphere_index, T);	
					
		Nlin_incompress::Compute_nlin_first_component(U, Giver);												
		// U.nlin1 = U.grad T<	
		
		flux_SF(sphere_index) = -Prod_out_sphere_nlinT(sphere_index, U, T);		
		// T.flux = (U.grad T<). T>		
							
	}
}

//
// RB convection
void EnergyTr::Compute_flux_RBC(FluidVF& U, FluidSF& T)
{
	
	flux_self = 0.0;
	flux_SF = 0.0;
	
	if (global.PHYSICS.Pr_option == "PRZERO") 
		Compute_flux(U);
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {		// fill only Temperature flux
		
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
			
			Fill_in_sphere(sphere_index, T);	
			
			Nlin_incompress::Compute_nlin_first_component(U, Giver);						
			// U.nlin1 = U.grad T<	
			
			flux_SF(sphere_index) = -Prod_out_sphere_nlinT(sphere_index, U, T);		
			// T.flux = -(U.grad T<). T>		
		}
	}
	
	else
		Compute_flux_scalar(U, T);	
}


//*********************************************************************************************

//
// Vector
void EnergyTr::Compute_flux(FluidVF& U, FluidVF& W)
{
	Compute_flux(U);									
	// flux_self = (U.grad U<). U>
	
	flux_VF_Uin_Wout = 0.0;
	flux_VF_Uin_Win = 0.0;
	
	// U< to W>;  U< to W<
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {		
		Fill_in_sphere(sphere_index, U);				
		
		Nlin_incompress::Compute_nlin(W, Giver);												
		// W.nlin = W.grad U<
		
		flux_VF_Uin_Wout(sphere_index) = Prod_out_sphere_nlinV(sphere_index, W, W);	
		// (W.graad U<). W>
		
		flux_VF_Uin_Win(sphere_index) = Prod_in_sphere_nlinV(sphere_index, W, W);		
		// (W.graad U<). W<
	}
	
	// U> to W>
	flux_VF_Uout_Wout = 0.0;
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {		
		Fill_out_sphere(sphere_index, U);										
			// G = U>
		
		Nlin_incompress::Compute_nlin(W, Giver);											
			// W.nlin = W.grad U>
		
		flux_VF_Uout_Wout(sphere_index) = Prod_out_sphere_nlinV(sphere_index, W, W);	
			// (W.graad U>). W>
	}
	
	
	// W< to W>;  W< to U>
	flux_VF_Win_Wout = 0.0;
	flux_VF_Win_Uout = 0.0;
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {	
		Fill_in_sphere(sphere_index, W);										
		// G = W<
		
		Nlin_incompress::Compute_nlin(U, Giver);												
		// U.nlin = U.grad W<
		
		flux_VF_Win_Wout(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, W);		
		// -(U.graad W<). W>
		
		flux_VF_Win_Uout(sphere_index) = Prod_out_sphere_nlinV(sphere_index, U, U);	
		//  (U.graad W<). U>
	}
	
	//
	// Flux for Elsasser vars
	//
	
	if (global.energy_transfer.Elsasser) {
		flux_Elsasser_plus = 0.0;
		flux_Elsasser_minus = 0.0;
		
		
		MHD::UB_to_Elsasser_field(U, W);												
		// U=Zp=(U+B); B=Zm=(U-B);
		
		// Flux: Zp to Zp
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++)  {		
			Fill_in_sphere(sphere_index, U);											
			// G = Zp<
			
			Nlin_incompress::Compute_nlin(W, Giver);												
			// W.nlin = Zm.grad Zp<
			
			flux_Elsasser_plus(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, W, U);	
			// (Zm.graad Zp<). Zp>
		
		}	
		
		
		// Flux: Zm to Zm
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++)  {
			Fill_in_sphere(sphere_index, W);										
			// G = Zm<
			
			Nlin_incompress::Compute_nlin(U, Giver);												
			// U.nlin = Zp.grad Zm<
			
			flux_Elsasser_minus(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, W);	
			// (Zp.graad Zm<). Zm>
		
		}																																			// Flux: Zp to Zp
		MHD::Elsasser_to_UB_field(U, W);													
		// Back to U, B vars
	} // end of if(Elsasser)
}


//*********************************************************************************************
//
// Magnetoconvection

void EnergyTr::Compute_flux(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	flux_SF = 0.0;
	Compute_flux(U, W);
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++)  {		
		Fill_in_sphere(sphere_index, T);	
					
		Nlin_incompress::Compute_nlin_first_component(U, Giver);											
		// U.nlin = U.grad T<	
		
		flux_SF(sphere_index) = -Prod_out_sphere_nlinT(sphere_index, U, T);		
		// T.flux = (U.grad T<). T>		
	}
}
/*********************************************************************************************

 HELICITY FLUX etc..
***********************************************************************************************/

// Fluid
void EnergyTr::Compute_kinetic_helicity_flux(FluidVF& U)
{	
	/*
	flux_hk = 0.0;

	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {	
		Fill_in_sphere(sphere_index, U);	
		
		Nlin_incompress::Compute_nlin(U, Giver);										
		// U.nlin = U.grad U<	
		
		flux_hk(sphere_index) =  - Prod_out_sphere_nlin_vorticity(sphere_index, U, U)/2;	
		// flux_hm = -(U.grad U<). omega>
		
		
		Fill_in_sphere_vorticity(sphere_index, U);
		// G = i k x u(k)
		
		Nlin_incompress::Compute_nlin(U, Giver);										
		// U.nlin = U.grad omega	
		
		flux_hk(sphere_index) +=   -Prod_out_sphere_nlinV(sphere_index, U, U)/2;	
		// Add -(U.grad omega). U>
		
		Fill_in_sphere(sphere_index, U);	
		
		Compute_nlin_vorticity_helper(U);										
		// U.nlin = omega.grad U<	
		
		flux_hk(sphere_index) +=   Prod_out_sphere_nlinV(sphere_index, U, U)/2;	
		// Add -(omega.grad U<). U>
		
							
	}
	 */
}


//*********************************************************************************************

//
// Vector
void EnergyTr::Compute_magnetic_helicity_flux(FluidVF& U, FluidVF& W)
{
	/*
	
	flux_Whk = 0.0;			// magnetic helicity
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) 
	{
		
		Fill_in_sphere(sphere_index, U);								
		// G = U<
		
		Compute_nlin(W);	
		// W.nlin = W.grad U<
		
		
		flux_Whk(sphere_index) = Prod_out_sphere_nlin_vector_potential(sphere_index, W, W)/2;	
		// flux_hm = W.grad U<. a>
		
		
		Fill_in_sphere(sphere_index, W);								
		// G = W<
		
		Compute_nlin(U);	
		// U.nlin = U.grad W<
		
		flux_Whk(sphere_index) += -Prod_out_sphere_nlin_vector_potential(sphere_index, U, W)/2;	
		// add -(U.graad W<). a>
		
		
		
		Fill_in_sphere(sphere_index, W);								
		// G = W<
		
		Compute_nlin_UcrossB(U);
		// U.nlin = u x W<
		
		flux_Whk(sphere_index) += Prod_out_sphere_nlinV(sphere_index, U, W)/2;		
		// add (U x W<). B>
		
	}
	*/
}



//*********************************************************************************************

// Fluid 2D
void EnergyTr::Compute_enstrophy_flux(FluidVF& U)
{	
	/*
	if (N[2] == 1)
	{
		flux_hk = 0.0;
		
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) 
		{
			
			Fill_in_sphere_vorticity(sphere_index, U);
			// G = i k x u(k): 2D
			
			Compute_nlin(U);										
			// U.nlin = U.grad omega	
					
			flux_hk(sphere_index) =   - Prod_out_sphere_nlin_vorticity(sphere_index, U, U);	
			// -(U.grad omega). omega> : 2D
			// Can make this efficient by mult omega_y nlin_y
		}
	}
	 */
}


//*********************************************************************************************

// MHD 2D
void EnergyTr::Compute_magnetic_enstrophy_flux(FluidVF& U, FluidVF& W)
{	
	/*
	if (N[2] == 1)
	{
		flux_Whk = 0.0;
		
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) 
		{
			
			Fill_in_sphere(sphere_index, W);
			// G = W<
			
			Compute_nlin_UcrossB(U);
			// U.nlin = u x W< : 2D
			
			flux_Whk(sphere_index) = Prod_out_sphere_nlin_vector_potential(sphere_index, U, W);	
			// add (U x W<). B>
		}
	}
	 */
}

		

//****************************  End of compute_flux.cc ****************************************






