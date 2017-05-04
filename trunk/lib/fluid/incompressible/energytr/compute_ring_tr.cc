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


/*! \file  compute_ring_tr.cc
 * 
 * @brief  Computes spherical ring-to-ring transfer between V/W to V/W.
 *
 *	The Giver field is filled in a ring  (region A1).  
 *	We do the same for the receiver a ring (region A2). <BR>
 *
 *	We exclude the last shell (Rmax,infty) for both Giver and Receiver rings. This is to
 *	save computer time.
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug	 No known bugs
 */

/**********************************************************************************************
 *
 *	No transfer for the last ring (max-possible inner radius, INF_RADIUS) 				
 * 
 *********************************************************************************************/  

#include "EnergyTr.h"

//*********************************************************************************************

void EnergyTr::Compute_ring_tr(FluidVF& U)
{
	
	ring_to_ring_self = 0.0;
	
	// skip the last shell -- outer rad = infty			
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
		{	
				
					Fill_ring(ring_shell_from_i, sector_from_i, U);
		
					Nlin_incompress::Compute_nlin(U, Giver);												
					// U.nlin = U.grad Um	

					universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
															
					ring_to_ring_self(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = -temp_ring_tr;
										
		}
}

void EnergyTr::Compute_kinetic_helicity_ring_tr(FluidVF& U,  FluidVF& helicalU)
{
  universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());

  ring_to_ring_VF_UtoW = 0.0;
  ring_to_ring_VF_WtoU = 0.0;
  ring_to_ring_VF_UtoU = 0.0;
  
  // skip the last shell -- outer rad = infty
  for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++)
    for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
    {
      
      Fill_ring(ring_shell_from_i, sector_from_i, U);
      
      Nlin_incompress::Compute_nlin(U, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_UtoW(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = -0.5*temp_ring_tr;
      
      /*****************************************************************************************/
      Fill_ring(ring_shell_from_i, sector_from_i, helicalU);
      
      Nlin_incompress::Compute_nlin(U, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_WtoU(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = -0.5*temp_ring_tr;
      
      /*****************************************************************************************/
      Fill_ring(ring_shell_from_i, sector_from_i, U);
      
      Nlin_incompress::Compute_nlin(helicalU, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(helicalU.nlin1, helicalU.nlin2, helicalU.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_UtoU(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = 0.5*temp_ring_tr;
      
    }
}


void EnergyTr::Compute_enstrophy_ring_tr(FluidVF& U, FluidVF& helicalU)
{
	
	universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());

	ring_to_ring_U_to_helicalU = 0.0;
	ring_to_ring_helicalU_to_helicalU = 0.0;	
	// skip the last shell -- outer rad = infty			
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
		{	
				
					Fill_ring(ring_shell_from_i, sector_from_i, U);
		
					Nlin_incompress::Compute_nlin(helicalU, Giver);												
					// U.nlin = U.grad Um	

					universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, temp_ring_tr);
															
					ring_to_ring_U_to_helicalU(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) += temp_ring_tr;
		

					Fill_ring(ring_shell_from_i, sector_from_i, helicalU);
		
					Nlin_incompress::Compute_nlin(U, Giver);												
					// U.nlin = U.grad Um	

					universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, temp_ring_tr);
															
					ring_to_ring_helicalU_to_helicalU (ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = -temp_ring_tr;


		}

}





//*********************************************************************************************
//	SCALAR   
//

void EnergyTr::Compute_ring_tr(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Compute_ring_tr_scalar(U, T);
	
	else if (global.program.kind == "RBC")
		Compute_ring_tr_RBC(U, T);
}


void EnergyTr::Compute_ring_tr_scalar(FluidVF& U, FluidSF& T)
{
	// U to U
	Compute_ring_tr(U);
	
	// T to T
	ring_to_ring_SF = 0.0;
	
	// skip the last shell -- outer rad = infty				
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
		{
						
			Fill_ring(ring_shell_from_i, sector_from_i, T);
		
			Nlin_incompress::Compute_nlin_first_component(U, Giver);										
			// U.nlin1 = U.grad Tm	
	
			universal->Ring_mult_all(U.nlin1, T.csf.F, temp_ring_tr);	
											
			ring_to_ring_SF(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
				= -temp_ring_tr;	
		}	
}


// RBC
void EnergyTr::Compute_ring_tr_RBC(FluidVF& U, FluidSF& T)
{
	
	ring_to_ring_self = 0.0;	
	ring_to_ring_SF = 0.0;
	
	if (global.PHYSICS.Pr_option == "PRZERO") 
		Compute_flux(U);
	
	else if (global.PHYSICS.Pr_option == "PRINFTY")	{	// fill only Temperature flux
		
		// skip the last shell -- outer rad = infty				
		for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
			for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++) {
				
				Fill_ring(ring_shell_from_i, sector_from_i, T);	
				
				Nlin_incompress::Compute_nlin_first_component(U, Giver);									
				// U.nlin1 = U.grad Tm	
				
				universal->Ring_mult_all(U.nlin1, T.csf.F, temp_ring_tr);	
				
				ring_to_ring_SF(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
					= -temp_ring_tr;	
				
			}
	}
	
	else
		Compute_ring_tr_scalar(U, T);
}


//*********************************************************************************************
//	 MHD   
//
void EnergyTr::Compute_ring_tr(FluidVF& U, FluidVF& W)
{
	Compute_ring_tr(U);
	
	
	// W to W
	ring_to_ring_VF_WtoW = 0.0;
	
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++) {	
		
			Fill_ring(ring_shell_from_i, sector_from_i, W);	
			
			Nlin_incompress::Compute_nlin(U, Giver);												
			// U.nlin = U.grad Wm	
			
			universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_ring_tr);
			
			ring_to_ring_VF_WtoW(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
				= -temp_ring_tr;	
			
		}						

	// U to W
	ring_to_ring_VF_UtoW = 0.0;
	
	
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
		{	
				
			Fill_ring(ring_shell_from_i, sector_from_i, U);
	
			Nlin_incompress::Compute_nlin(W, Giver);											
			// W.nlin = W.grad Um	
			
			universal->Ring_mult_all(W.nlin1, W.nlin2, W.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_ring_tr);
											
			ring_to_ring_VF_UtoW(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
				= temp_ring_tr;	
			
		}
	
	
	// Shell_to_shell transfers for Elsasser vars
	//
	if (global.energy_transfer.Elsasser) {
		ring_to_ring_Elsasser_plus = 0.0;
		ring_to_ring_Elsasser_minus = 0.0;
		
		MHD::UB_to_Elsasser_field(U, W);											
		// U=Zp=(U+B); B=Zm=(U-B);
		
		// ring_to_ring: Zp to Zp
		for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
			for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++) {	
		
				Fill_ring(ring_shell_from_i, sector_from_i, U);				
				// G = Zp<									
			
				Nlin_incompress::Compute_nlin(W, Giver);											
				// W.nlin = Zm.grad Zp<	
				
				universal->Ring_mult_all(W.nlin1, W.nlin2, W.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
				
				ring_to_ring_Elsasser_plus(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
					= -temp_ring_tr;	

			}
		
		// ring_to_ring: Zm to Zm	
		for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
			for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++) {	
		
				Fill_ring(ring_shell_from_i, sector_from_i, W);				
				// G = Zm<
			
				Nlin_incompress::Compute_nlin(U, Giver);										
				// U.nlin = Zp.grad Zm<	
		
				universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_ring_tr);
													
				ring_to_ring_Elsasser_minus(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
					= -temp_ring_tr;	
				
			}


		MHD::Elsasser_to_UB_field(U, W);													
		// Back to U, B vars
	}	
	
}

void EnergyTr::Compute_kinetic_helicity_ring_tr(FluidVF& U, FluidVF& W, FluidVF& helicalU, FluidVF& helicalW)
{
  universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
  universal->Compute_vorticity(W.cvf.V1, W.cvf.V2, W.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());
  
  ring_to_ring_VF_UtoW = 0.0;
  ring_to_ring_VF_WtoU = 0.0;
  ring_to_ring_VF_UtoU = 0.0;
  
  ring_to_ring_VF_BtoW = 0.0;
  ring_to_ring_VF_JtoU = 0.0;
  ring_to_ring_VF_BtoU = 0.0;
  
  // skip the last shell -- outer rad = infty
  for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++)
    for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
    {
      
      Fill_ring(ring_shell_from_i, sector_from_i, U);
      
      Nlin_incompress::Compute_nlin(U, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_UtoW(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = -0.5*temp_ring_tr;
      
      /*****************************************************************************************/
      Fill_ring(ring_shell_from_i, sector_from_i, helicalU);
      
      Nlin_incompress::Compute_nlin(U, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_WtoU(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = -0.5*temp_ring_tr;
      
      /*****************************************************************************************/
      Fill_ring(ring_shell_from_i, sector_from_i, U);
      
      Nlin_incompress::Compute_nlin(helicalU, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(helicalU.nlin1, helicalU.nlin2, helicalU.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_UtoU(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = 0.5*temp_ring_tr;
      
      /*****************************************************************************************/
      Fill_ring(ring_shell_from_i, sector_from_i, W);
      
      Nlin_incompress::Compute_nlin(W, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(W.nlin1, W.nlin2, W.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_BtoW(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = 0.5*temp_ring_tr;
      
      /*****************************************************************************************/
      Fill_ring(ring_shell_from_i, sector_from_i, helicalW);
      
      Nlin_incompress::Compute_nlin(W, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(W.nlin1, W.nlin2, W.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_JtoU(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = 0.5*temp_ring_tr;
      
      /*****************************************************************************************/
      Fill_ring(ring_shell_from_i, sector_from_i, W);
      
      Nlin_incompress::Compute_nlin(helicalW, Giver);
      // U.nlin = U.grad Um
      
      universal->Ring_mult_all(helicalW.nlin1, helicalW.nlin2, helicalW.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_ring_tr);
      
      ring_to_ring_VF_BtoU(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) = -0.5*temp_ring_tr;
      
    }
}


//*********************************************************************************************
// Vector + Scalar
//
void EnergyTr::Compute_ring_tr(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	// U/W to U/W
	Compute_ring_tr(U, W);
	
	// T to T
	ring_to_ring_SF = 0.0;
	
	// skip the last shell -- outer rad = infty				
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++) {
						
			Fill_ring(ring_shell_from_i, sector_from_i, T);	
		
			Nlin_incompress::Compute_nlin_first_component(U, Giver);												
			// U.nlin1 = U.grad Tm	
			
			universal->Ring_mult_all(U.nlin1, T.csf.F, temp_ring_tr);	
			
			ring_to_ring_SF(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
				= -temp_ring_tr;

		}
}		


//****************************** End of Compute_ring_tr.cc ************************************


