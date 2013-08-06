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
 * @brief  Computes cylindrical ring-to-ring transfer between V/W to V/W.
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


#include "EnergyTr.h"

//********************************************************************************************

void EnergyTr::Compute_cylindrical_ring_tr(FluidVF& U)
{
	
	cylindrical_ring_to_ring_self = 0.0;
				
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {	
				
			Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, U);	
		
			Nlin_incompress::Compute_nlin(U, Giver);								
			// nlin = U.grad Um	

			universal->Cyl_ring_mult_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_cylindrical_ring_tr); 
			
			cylindrical_ring_to_ring_self(cylindrical_shell_from_i, slab_from_i, Range::all(), Range::all()) = -temp_cylindrical_ring_tr;	
			
		}
			
}

//*********************************************************************************************
//	SCALAR   
//

void EnergyTr::Compute_cylindrical_ring_tr(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Compute_cylindrical_ring_tr_scalar(U, T);
	
	else if (global.program.kind == "RBC")
		Compute_cylindrical_ring_tr_RBC(U, T);
}


void EnergyTr::Compute_cylindrical_ring_tr_scalar(FluidVF& U, FluidSF& T)
{
	// U to U
	Compute_cylindrical_ring_tr(U);
	
	// T to T
	cylindrical_ring_to_ring_SF = 0.0;
		
	// skip the last shell -- outer rad = infty				
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {
			Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, T);	
		
			Nlin_incompress::Compute_nlin_first_component(U, Giver);										
			// U.nlin1 = U.grad Tm	
			
			universal->Cyl_ring_mult_all(U.nlin1, T.csf.F, temp_cylindrical_ring_tr); 
													
			cylindrical_ring_to_ring_SF(cylindrical_shell_from_i, slab_from_i, Range::all(), Range::all()) = -temp_cylindrical_ring_tr;	
		}
					
}


// RBC
void EnergyTr::Compute_cylindrical_ring_tr_RBC(FluidVF& U, FluidSF& T)
{
	
	cylindrical_ring_to_ring_self = 0.0;
	cylindrical_ring_to_ring_SF = 0.0;
	
	if (global.PHYSICS.Pr_option == "PRZERO")
		Compute_cylindrical_ring_tr(U);
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {		// fill only Temperature transfers
		// T to T
		// skip the last shell -- outer rad = infty				
		for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
			for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {
				
				Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, T);	
				
				Nlin_incompress::Compute_nlin_first_component(U, Giver);									
				// U.nlin1 = U.grad Tm	
				
				universal->Cyl_ring_mult_all(U.nlin1, T.csf.F, temp_cylindrical_ring_tr); 
				
				cylindrical_ring_to_ring_SF(cylindrical_shell_from_i, slab_from_i, Range::all(), Range::all()) = -temp_cylindrical_ring_tr;	
			}
	}	
	
	else
		Compute_cylindrical_ring_tr_scalar(U, T);
				
}


//*********************************************************************************************
//	MHD   
//
void EnergyTr::Compute_cylindrical_ring_tr(FluidVF& U, FluidVF& W)
{
	// U to U
	Compute_cylindrical_ring_tr(U);
	
	// W to W
	cylindrical_ring_to_ring_VF_WtoW = 0.0;
		
	
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {	
		
			Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, W);
			
			Nlin_incompress::Compute_nlin(U, Giver);									
			// U.nlin = U.grad Wm	
			
			universal->Cyl_ring_mult_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_cylindrical_ring_tr);
											
			cylindrical_ring_to_ring_VF_WtoW(cylindrical_shell_from_i, slab_from_i, Range::all(), Range::all()) = -temp_cylindrical_ring_tr;	
		}	// of for loop					



	// U to W
	cylindrical_ring_to_ring_VF_UtoW = 0.0;
		
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {	
			Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, U);
			
			Nlin_incompress::Compute_nlin(W, Giver);									
			// W.nlin = W.grad Um	

			universal->Cyl_ring_mult_all(W.nlin1, W.nlin2, W.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_cylindrical_ring_tr);
											
			cylindrical_ring_to_ring_VF_UtoW(cylindrical_shell_from_i, slab_from_i,Range::all(), Range::all()) = temp_cylindrical_ring_tr;

		}  // end of for loop
	
	
	// Shell_to_shell transfers for Elsasser vars
	//
	
	if (global.energy_transfer.Elsasser) {
		cylindrical_ring_to_ring_Elsasser_plus = 0.0;
		cylindrical_ring_to_ring_Elsasser_minus = 0.0;
		
		MHD::UB_to_Elsasser_field(U, W);								
		// U=Zp=(U+B); B=Zm=(U-B);
		
		// ring_to_ring: Zp to Zp
		for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
			for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {	
				Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, U);		
				// G = Zp<
			
				Nlin_incompress::Compute_nlin(W, Giver);										
				// W.nlin = Zm.grad Zp<	

				universal->Cyl_ring_mult_all(W.nlin1, W.nlin2, W.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_cylindrical_ring_tr);	
											
				cylindrical_ring_to_ring_Elsasser_plus(cylindrical_shell_from_i, slab_from_i,Range::all(), Range::all()) = -temp_cylindrical_ring_tr;	
				
			}		// end of for loop
		
		// ring_to_ring: Zm to Zm	
		for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
			for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {	
		
				Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, W);	
				// G = Zm<
			
				Nlin_incompress::Compute_nlin(U, Giver);									
				// nlin = Zp.grad Zm<	
			
				universal->Cyl_ring_mult_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_cylindrical_ring_tr);	
													
				cylindrical_ring_to_ring_Elsasser_minus(cylindrical_shell_from_i, slab_from_i,Range::all(), Range::all()) = -temp_cylindrical_ring_tr;
				
			}


		MHD::Elsasser_to_UB_field(U, W);									
		// Back to U, B vars
	}
	
}

//*********************************************************************************************
//	VECTOR + SCALAR   
//
void EnergyTr::Compute_cylindrical_ring_tr(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	Compute_cylindrical_ring_tr(U, W);
	
	cylindrical_ring_to_ring_SF = 0.0;
		
	
	// skip the last shell -- outer rad = infty				
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++) 
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {
						
			Fill_cylindrical_ring(cylindrical_shell_from_i, slab_from_i, T);
		
			Nlin_incompress::Compute_nlin_first_component(U, Giver);										
			// U.nlin1 = U.grad Tm	

			universal->Cyl_ring_mult_all(U.nlin1, T.csf.F, temp_cylindrical_ring_tr);
											
			cylindrical_ring_to_ring_SF(cylindrical_shell_from_i, slab_from_i,Range::all(), Range::all()) = -temp_cylindrical_ring_tr;	
			
		}		
}		


//**************************** End of Compute_ring_tr.cc **************************************


