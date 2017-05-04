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

/*! \file  compute_shell_tr.cc
 * 
 * @brief  Computes shell-to-shell transfer between V/W to V/W.
 *
 *	The Giver field is filled in shell  (region A1).  
 *	We do the same for the receiver field (region A2). <BR>
 *  
 *  The shell-from-index = 1:no-shells-1.
 *   shell-mult-all() provides us shell-to-shell energy transfer from these shells to
 *	 1:no_shells.  We output shell-to-shell(1:no-shells-1, 1:no-shells).
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 

#include "EnergyTr.h"


//*********************************************************************************************

// Fluid
//
void EnergyTr::Compute_shell_tr(FluidVF& U)
{

	shelltoshell_self = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++)  {	
		Fill_shell(shell_from_index, U);
		
		Nlin_incompress::Compute_nlin(U, Giver);
		// U.nlin = U.grad Um	
		
		universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
		// results(shell_index) in (*temp_shell_tr)(index)
				
		shelltoshell_self(shell_from_index, Range::all()) = -temp_shell_tr;
			
	}
 
}

//
//	SCALAR   
//

void EnergyTr::Compute_shell_tr(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Compute_shell_tr_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Compute_shell_tr_RBC(U, T);
}


void EnergyTr::Compute_shell_tr_scalar(FluidVF& U, FluidSF& T)
{
	// U to U
	Compute_shell_tr(U);
	
	shelltoshell_SF = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) {	
		Fill_shell(shell_from_index, T);	
		
		Nlin_incompress::Compute_nlin_first_component(U, Giver);												
		// U.nlin1 = U.grad Tm	
		
		universal->Shell_mult_all(U.nlin1, T.csf.F, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);	
		
		shelltoshell_SF(shell_from_index, Range::all()) = -temp_shell_tr;	
	}
}

// RB
void EnergyTr::Compute_shell_tr_RBC(FluidVF& U, FluidSF& T)
{
	
	shelltoshell_self = 0.0;
	shelltoshell_SF = 0.0;
	
	if (global.PHYSICS.Pr_option == "PRZERO") 
		Compute_shell_tr(U);
	
	else if (global.PHYSICS.Pr_option == "PRINFTY")	{	// fill only Temperature flux
		shelltoshell_SF = 0.0;
		
		for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) {
			Fill_shell(shell_from_index, T);	
			
			Nlin_incompress::Compute_nlin_first_component(U, Giver);											
			// U.nlin1 = U.grad Tm	
			
			universal->Shell_mult_all(U.nlin1, T.csf.F, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);		
			
			shelltoshell_SF(shell_from_index, Range::all()) = -temp_shell_tr;	
			
		}
	}
	
	else
		Compute_shell_tr_scalar(U, T);
}

//*********************************************************************************************
//	 MHD   
//
void EnergyTr::Compute_shell_tr(FluidVF& U, FluidVF& W)
{
	// U to U
	Compute_shell_tr(U);
	
	// W to W
	shelltoshell_VF_WtoW = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) {
		Fill_shell(shell_from_index, W);	
		
		Nlin_incompress::Compute_nlin(U, Giver);												
		// U.nlin = U.grad Wm	
		
		universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
		
		shelltoshell_VF_WtoW(shell_from_index, Range::all()) = -temp_shell_tr;	
	}


	// U to W
	shelltoshell_VF_UtoW = 0.0;
	

	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) {	
		Fill_shell(shell_from_index, U);	
		
		Nlin_incompress::Compute_nlin(W, Giver);												
		// W.nlin = W.grad Um	
		
		universal->Shell_mult_all(W.nlin1, W.nlin2, W.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
						
		shelltoshell_VF_UtoW(shell_from_index, Range::all()) = temp_shell_tr;
			

	}
	
	
	// Shell_to_shell transfers for Elsasser vars
	//
	if (global.energy_transfer.Elsasser) {
		shelltoshell_Elsasser_plus = 0.0;
		shelltoshell_Elsasser_minus = 0.0;

		MHD::UB_to_Elsasser_field(U, W);											
		// U=Zp=(U+B); B=Zm=(U-B);
		
		// Shell_to_shell: Zp to Zp
		for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) 
		{	
			Fill_shell(shell_from_index, U);										
			// G = Zp<
			
			Nlin_incompress::Compute_nlin(W, Giver);												
			// W.nlin = Zm.grad Zp<	
			
			universal->Shell_mult_all(W.nlin1, W.nlin2, W.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
								
			shelltoshell_Elsasser_plus(shell_from_index, Range::all()) = -temp_shell_tr;	
		}

		
		// Shell_to_shell: Zm to Zm
		for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) 
		{
			Fill_shell(shell_from_index, W);									
			// G = Zm<
			
			Nlin_incompress::Compute_nlin(U, Giver);												
			// U.nlin = Zp.grad Zm<	
			
			universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
			
			shelltoshell_Elsasser_minus(shell_from_index, Range::all()) = -temp_shell_tr;	
		}

		MHD::Elsasser_to_UB_field(U, W);													
		// Back to U, B vars
	}
	
}



//*********************************************************************************************
// Vector + Scalar
//
void EnergyTr::Compute_shell_tr(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	// U/W to U/W
	Compute_shell_tr(U, W);
	
	// T to T
	shelltoshell_SF = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) 
	{	
		Fill_shell(shell_from_index, T);	
		
		Nlin_incompress::Compute_nlin_first_component(U, Giver);											
		// U.nlin1 = U.grad Tm	
		
		universal->Shell_mult_all(U.nlin1, T.csf.F, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);	
		shelltoshell_SF(shell_from_index, Range::all()) = -temp_shell_tr;	
		
	}
		
}	

/*********************************************************************************************

 HELICITY
**********************************************************************************************/

// Fluid
//
void EnergyTr::Compute_kinetic_helicity_shell_tr(FluidVF& U, FluidVF& helicalU)
{
  universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
  
  shelltoshell_VF_UtoW = 0.0;
  shelltoshell_VF_WtoU = 0.0;
  shelltoshell_VF_UtoU = 0.0;


  
  for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++)  {
    Fill_shell(shell_from_index, U);
    
    Nlin_incompress::Compute_nlin(U, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_UtoW(shell_from_index, Range::all()) = -0.5*temp_shell_tr;
    
    /************************************************************************************/
    Fill_shell(shell_from_index, helicalU);
    
    Nlin_incompress::Compute_nlin(U, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_WtoU(shell_from_index, Range::all()) = -0.5*temp_shell_tr;
    
    /************************************************************************************/
    Fill_shell(shell_from_index, U);
    
    Nlin_incompress::Compute_nlin(helicalU, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(helicalU.nlin1, helicalU.nlin2, helicalU.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_UtoU(shell_from_index, Range::all()) = 0.5*temp_shell_tr;
    
    }

}
	


void EnergyTr::Compute_kinetic_helicity_shell_tr(FluidVF& U, FluidVF& W, FluidVF& helicalU, FluidVF& helicalW)
{
  universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
  universal->Compute_vorticity(W.cvf.V1, W.cvf.V2, W.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());
  
  shelltoshell_VF_UtoW = 0.0;
  shelltoshell_VF_WtoU = 0.0;
  shelltoshell_VF_UtoU = 0.0;
  
  shelltoshell_VF_BtoW = 0.0;
  shelltoshell_VF_JtoU = 0.0;
  shelltoshell_VF_BtoU = 0.0;
  
  
  
  for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++)  {
    Fill_shell(shell_from_index, U);
    
    Nlin_incompress::Compute_nlin(U, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_UtoW(shell_from_index, Range::all()) = -0.5*temp_shell_tr;
    
    /************************************************************************************/
    Fill_shell(shell_from_index, helicalU);
    
    Nlin_incompress::Compute_nlin(U, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_WtoU(shell_from_index, Range::all()) = -0.5*temp_shell_tr;
    
    /************************************************************************************/
    Fill_shell(shell_from_index, U);
    
    Nlin_incompress::Compute_nlin(helicalU, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(helicalU.nlin1, helicalU.nlin2, helicalU.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_UtoU(shell_from_index, Range::all()) = 0.5*temp_shell_tr;
    
    /************************************************************************************/
    Fill_shell(shell_from_index, W);
    
    Nlin_incompress::Compute_nlin(W, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(W.nlin1, W.nlin2, W.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_BtoW(shell_from_index, Range::all()) = 0.5*temp_shell_tr;
    
    /************************************************************************************/
    Fill_shell(shell_from_index, helicalW);
    
    Nlin_incompress::Compute_nlin(W, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(W.nlin1, W.nlin2, W.nlin3,  U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_JtoU(shell_from_index, Range::all()) = 0.5*temp_shell_tr;
    
    /************************************************************************************/
    Fill_shell(shell_from_index, W);
    
    Nlin_incompress::Compute_nlin(helicalW, Giver);
    // U.nlin = U.grad Um
    
    universal->Shell_mult_all(helicalW.nlin1, helicalW.nlin2, helicalW.nlin3,  U.cvf.V1, U.cvf.V2, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
    // results(shell_index) in (*temp_shell_tr)(index)
				
    shelltoshell_VF_BtoU(shell_from_index, Range::all()) = -0.5*temp_shell_tr;
    
  }
  
}


//*********************************************************************************************
//	 MHD   
//
void EnergyTr::Compute_magnetic_helicity_shell_tr(FluidVF& U, FluidVF& W)
{
	/*	
	shelltoshell_Whk = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) 
	{
		
		Fill_shell(shell_from_index, U);	
		// G = U<
		
		Compute_nlin(W);								
		// W.nlin = W.grad U<	
	
		Shell_mult_vector_potential_all(W.nlin1, W.nlin2, W.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_shell_tr);
		
		(*shelltoshell_Whk)(shell_from_index, Range::all()) =  temp_shell_tr/2;
		// flux_hm = W.grad U<. a>
		
		
		Fill_shell(shell_from_index, W);	
		// G = W<
		
		Compute_nlin(U);	
		// nlin = U.grad W<
		
		Shell_mult_vector_potential_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_shell_tr);
		
		shelltoshell_Whk(shell_from_index, Range::all()) 
			= shelltoshell_Whk(shell_from_index, Range::all()) - (*temp_shell_tr)/2;
		// flux_hm = -U.grad W<. a>
			
		
		
		Fill_shell(shell_from_index, W);	
		// G = W<
		
		Compute_nlin_UcrossB(U);
		// nlin = u x W<
		
		
		universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
								
		shelltoshell_hk(shell_from_index, Range::all()) 
			= shelltoshell_Whk(shell_from_index, Range::all()) + temp_shell_tr/2;	
		
	}
	 */

}



//*********************************************************************************************
// Fluid: 2D
//
void EnergyTr::Compute_enstrophy_shell_tr(FluidVF& U, FluidVF& helicalU)
{
	if (Ny > 1)
	{
		universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());

		shelltoshell_hk_helicalU_to_helicalU = 0.0;
		shelltoshell_hk_U_to_helicalU = 0.0;
		
		for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) 
		{
			
			Fill_shell(shell_from_index, U);	
			// G = i k x u(k): 2D
			
			Nlin_incompress::Compute_nlin(helicalU, Giver);										
			// U.nlin = U.grad omega	
			
			universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3,global.energy_transfer.shell_to_shell.radii, temp_shell_tr);

//			Shell_mult_vorticity_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, heU.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3,//licalU.cvf.V3, temp_shell_tr);
			
			shelltoshell_hk_U_to_helicalU(shell_from_index, Range::all()) += temp_shell_tr;


			

			Fill_shell(shell_from_index, helicalU);	
			// G = i k x u(k): 2D
			
			Nlin_incompress::Compute_nlin(U, Giver);										
			// U.nlin = U.grad omega	

			universal->Shell_mult_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
			
//			Shell_mult_vorticity_all(U.nlin1, U.nlin2, U.nlin3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, temp_shell_tr);
			
			shelltoshell_hk_helicalU_to_helicalU(shell_from_index, Range::all()) -= temp_shell_tr;
		}
	}
	/*
	if (N[2] == 1)
	{
		shelltoshell_hk = 0.0;
		
		for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) 
		{
			
			Fill_shell_vorticity(shell_from_index, U);	
			// G = i k x u(k): 2D
			
			Compute_nlin(U);										
			// U.nlin = U.grad omega	
			
			Shell_mult_vorticity_all(U.nlin1, U.nlin2, U.nlin3, U.cvf.V1, U.cvf.V2, U.cvf.V3, temp_shell_tr);
			
			shelltoshell_hk(shell_from_index, Range::all()) = - temp_shell_tr;
		}
	}
	
	*/
}

//*********************************************************************************************
//	 MHD   
//
void EnergyTr::Compute_magnetic_enstrophy_shell_tr(FluidVF& U, FluidVF& W)
{	
	/*
	if (N[2] == 1)
	{
		shelltoshell_Whk = 0.0;
		
		for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) 
		{
			
			Fill_shell(shell_from_index, W);	
			// G = W
			
			Compute_nlin_UcrossB(U);
			// nlin = u x W< : 2D
			
			
			Shell_mult_vector_potential_all(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, temp_shell_tr);
			
			shelltoshell_Whk(shell_from_index, Range::all()) =  temp_shell_tr;
		}	
	}
	*/
}




//******************************  End of Compute_shell_tr.cc  *********************************


