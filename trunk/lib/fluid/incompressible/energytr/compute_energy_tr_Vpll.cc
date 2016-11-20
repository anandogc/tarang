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

//****************************************************************************************

// global.field.anisotropy_dirn = dirn of Vpll
void EnergyTr::Compute_flux_Vpll(FluidVF& U, Pressure& P)
{
        
	flux_self = 0.0;
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
		
		Fill_in_sphere_Vpll(sphere_index, U);  // Giver <- U.cvf.V1
					
		Nlin_incompress::Compute_nlin_first_component(U, Giver);
		// U.nlin1 = U.grad Vpll	
		
		flux_self(sphere_index) = -Prod_out_sphere_nlin_Vpll(sphere_index, U, U);
		// T.flux = (U.grad Vpll). Vpll		
	}
    
    if (global.field.anisotropy_dirn == 1)
        universal->Xderiv(U.cvf.V1, global.temp_array.X);
    else if (global.field.anisotropy_dirn == 2)
        universal->Yderiv(U.cvf.V2, global.temp_array.X);
    if (global.field.anisotropy_dirn == 3)
        universal->Zderiv(U.cvf.V3, global.temp_array.X);
    
    P.Compute_pressure(U);

    universal->Shell_mult_all(global.temp_array.X, P.F, global.energy_transfer.flux.radii, sphere_force_x_field);
    //   shell_prod_pressure_Vpll = Re(i kpll upll(k)* conj(p(k)))
}

//*********************************************************************************************

//
// Vector
void EnergyTr::Compute_flux_Vpll(FluidVF& U, FluidVF& W, Pressure& P)
{
	
    Compute_flux_Vpll(U, P);
    
  	flux_VF_Uin_Wout = 0.0;
	flux_VF_Uin_Win = 0.0;
	
	
	// U< to W>;  U< to W<
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
		Fill_in_sphere_Vpll(sphere_index, U);
		
		Nlin_incompress::Compute_nlin_first_component(W, Giver);
		// W.nlin = W.grad U<
		
		flux_VF_Uin_Wout(sphere_index) = Prod_out_sphere_nlin_Vpll(sphere_index, W, W);
		// (W.graad U<). W>
		
		flux_VF_Uin_Win(sphere_index) = Prod_in_sphere_nlin_Vpll(sphere_index, W, W);
		// (W.graad U<). W<
	}
	
	// U> to W>
	flux_VF_Uout_Wout = 0.0;
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
		Fill_out_sphere_Vpll(sphere_index, U);
        // G = U>
		
		Nlin_incompress::Compute_nlin_first_component(W, Giver);
        // W.nlin = W.grad U>
		
		flux_VF_Uout_Wout(sphere_index) = Prod_out_sphere_nlin_Vpll(sphere_index, W, W);
        // (W.graad U>). W>
	}
	
	
	// W< to W>;  W< to U>
	flux_VF_Win_Wout = 0.0;
	flux_VF_Win_Uout = 0.0;
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
		Fill_in_sphere_Vpll(sphere_index, W);
		// G = W<
		
		Nlin_incompress::Compute_nlin_first_component(U, Giver);
		// U.nlin = U.grad W<
		
		flux_VF_Win_Wout(sphere_index) = -Prod_out_sphere_nlin_Vpll(sphere_index, U, W);
		// -(U.graad W<). W>
		
		flux_VF_Win_Uout(sphere_index) = Prod_out_sphere_nlin_Vpll(sphere_index, U, U);
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
			Fill_in_sphere_Vpll(sphere_index, U);
			// G = Zp<
			
			Nlin_incompress::Compute_nlin_first_component(W, Giver);
			// W.nlin = Zm.grad Zp<
			
			flux_Elsasser_plus(sphere_index) = -Prod_out_sphere_nlin_Vpll(sphere_index, W, U);
			// (Zm.graad Zp<). Zp>
            
		}
		
		// Flux: Zm to Zm
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++)  {
			Fill_in_sphere_Vpll(sphere_index, W);
			// G = Zm<
			
			Nlin_incompress::Compute_nlin_first_component(U, Giver);
			// U.nlin = Zp.grad Zm<
			
			flux_Elsasser_minus(sphere_index) = -Prod_out_sphere_nlin_Vpll(sphere_index, U, W);
			// (Zp.graad Zm<). Zm>
            
		}																																			// Flux: Zp to Zp
		MHD::Elsasser_to_UB_field(U, W);
		// Back to U, B vars
	} // end of if(Elsasser)
}




//****************************************************************************************

void EnergyTr::Compute_shell_tr_Vpll(FluidVF& U, Pressure& P)
{
	
	shelltoshell_self = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) {
		Fill_shell_Vpll(shell_from_index, U);
		
		Nlin_incompress::Compute_nlin_first_component(U, Giver);
		// U.nlin1 = U.grad vpll
		
        if (global.field.anisotropy_dirn == 1)
            universal->Shell_mult_all(U.nlin1, U.cvf.V1, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
        
        else if (global.field.anisotropy_dirn == 2)
            universal->Shell_mult_all(U.nlin1, U.cvf.V2, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
        
        else if (global.field.anisotropy_dirn == 3)
            universal->Shell_mult_all(U.nlin1, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
		
		shelltoshell_self(shell_from_index, Range::all()) = -temp_shell_tr;
	}
    
    
}

//*********************************************************************************************
//	 MHD
//
void EnergyTr::Compute_shell_tr_Vpll(FluidVF& U, FluidVF& W, Pressure& P)
{
	// U to U
	Compute_shell_tr(U);
	
	// W to W
	shelltoshell_VF_WtoW = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) {
		Fill_shell_Vpll(shell_from_index, W);
		
		Nlin_incompress::Compute_nlin_first_component(U, Giver);
		// U.nlin = U.grad Wm
        
        if (global.field.anisotropy_dirn == 1)
            universal->Shell_mult_all(U.nlin1, W.cvf.V1, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
        
        else if (global.field.anisotropy_dirn == 2)
            universal->Shell_mult_all(U.nlin1, W.cvf.V2, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
        
        else if (global.field.anisotropy_dirn == 3)
            universal->Shell_mult_all(U.nlin1, W.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
		
		shelltoshell_VF_WtoW(shell_from_index, Range::all()) = -temp_shell_tr;
	}
    
    
	// U to W
	shelltoshell_VF_UtoW = 0.0;
	
    
	for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++) {
		Fill_shell_Vpll(shell_from_index, U);
		
		Nlin_incompress::Compute_nlin_first_component(W, Giver);
		// W.nlin = W.grad Um
        
        if (global.field.anisotropy_dirn == 1)
            universal->Shell_mult_all(W.nlin1, W.cvf.V1, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
        
        else if (global.field.anisotropy_dirn == 2)
            universal->Shell_mult_all(W.nlin1, W.cvf.V2, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
        
        else if (global.field.anisotropy_dirn == 3)
            universal->Shell_mult_all(W.nlin1, W.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
        
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
			Fill_shell_Vpll(shell_from_index, U);
			// G = Zp<
			
			Nlin_incompress::Compute_nlin_first_component(W, Giver);
			// W.nlin = Zm.grad Zp<
            
            if (global.field.anisotropy_dirn == 1)
                universal->Shell_mult_all(W.nlin1, U.cvf.V1, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
            
            else if (global.field.anisotropy_dirn == 2)
                universal->Shell_mult_all(W.nlin1, U.cvf.V2, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Shell_mult_all(W.nlin1, U.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
            
			shelltoshell_Elsasser_plus(shell_from_index, Range::all()) = -temp_shell_tr;
		}
        
		
		// Shell_to_shell: Zm to Zm
		for (int shell_from_index = 1; shell_from_index < global.energy_transfer.shell_to_shell.no_shells; shell_from_index++)
		{
			Fill_shell_Vpll(shell_from_index, W);
			// G = Zm<
			
			Nlin_incompress::Compute_nlin_first_component(U, Giver);
			// U.nlin = Zp.grad Zm<
            
            if (global.field.anisotropy_dirn == 1)
                universal->Shell_mult_all(U.nlin1, W.cvf.V1, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
            
            else if (global.field.anisotropy_dirn == 2)
                universal->Shell_mult_all(U.nlin1, W.cvf.V2, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Shell_mult_all(U.nlin1, W.cvf.V3, global.energy_transfer.shell_to_shell.radii, temp_shell_tr);
			
			shelltoshell_Elsasser_minus(shell_from_index, Range::all()) = -temp_shell_tr;
		}
        
		MHD::Elsasser_to_UB_field(U, W);
		// Back to U, B vars
	}
	
}




//****************************************************************************************


void EnergyTr::Compute_ring_tr_Vpll(FluidVF& U, Pressure& P)
{
    
	// vpll to vpp
	ring_to_ring_self = 0.0;
	
	// skip the last shell -- outer rad = infty
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++)
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
		{
            
			Fill_ring_Vpll(ring_shell_from_i, sector_from_i, U);
            
			Nlin_incompress::Compute_nlin_first_component(U, Giver);
			// U.nlin1 = U.grad Upll
            
            if (global.field.anisotropy_dirn == 1)
                universal->Ring_mult_all(U.nlin1, U.cvf.V1, temp_ring_tr);
            
            else if (global.field.anisotropy_dirn == 2)
                universal->Ring_mult_all(U.nlin1, U.cvf.V2, temp_ring_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Ring_mult_all(U.nlin1, U.cvf.V3, temp_ring_tr);
            
			ring_to_ring_self(ring_shell_from_i, sector_from_i, Range::all(), Range::all())
            = -temp_ring_tr;
		}
    
    if (global.field.anisotropy_dirn == 1)
        universal->Xderiv(U.cvf.V1, global.temp_array.X);
    else if (global.field.anisotropy_dirn == 2)
        universal->Yderiv(U.cvf.V2, global.temp_array.X);
    if (global.field.anisotropy_dirn == 3)
        universal->Zderiv(U.cvf.V3, global.temp_array.X);
	
    P.Compute_pressure(U);
    
    universal->Ring_mult_all(global.temp_array.X, P.F, ring_force_x_field);
}



//****************************************************************************************
//	 MHD
//
void EnergyTr::Compute_ring_tr_Vpll(FluidVF& U, FluidVF& W, Pressure& P)
{

	Compute_ring_tr_Vpll(U, P);

	// W to W
	ring_to_ring_VF_WtoW = 0.0;
	
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++)
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++) {
            
			Fill_ring_Vpll(ring_shell_from_i, sector_from_i, W);
			
			Nlin_incompress::Compute_nlin_first_component(U, Giver);
			// U.nlin = U.grad Wm
			
            if (global.field.anisotropy_dirn == 1)
                universal->Ring_mult_all(U.nlin1, W.cvf.V1, temp_ring_tr);
            
            else if (global.field.anisotropy_dirn == 2)
                universal->Ring_mult_all(U.nlin1, W.cvf.V2, temp_ring_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Ring_mult_all(U.nlin1, W.cvf.V3, temp_ring_tr);
			
			ring_to_ring_VF_WtoW(ring_shell_from_i, sector_from_i, Range::all(), Range::all())
            = -temp_ring_tr;
		}
    
	// U to W
	ring_to_ring_VF_UtoW = 0.0;
	
	
	for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++)
		for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++)
		{
            
			Fill_ring_Vpll(ring_shell_from_i, sector_from_i, U);
            
			Nlin_incompress::Compute_nlin_first_component(W, Giver);
			// W.nlin = W.grad Um
			
            if (global.field.anisotropy_dirn == 1)
                universal->Ring_mult_all(W.nlin1, W.cvf.V1, temp_ring_tr);
            
            else if (global.field.anisotropy_dirn == 2)
                universal->Ring_mult_all(W.nlin1, W.cvf.V2, temp_ring_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Ring_mult_all(W.nlin1, W.cvf.V3, temp_ring_tr);
            
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
                
				Fill_ring_Vpll(ring_shell_from_i, sector_from_i, U);
				// G = Zp<
                
				Nlin_incompress::Compute_nlin_first_component(W, Giver);
				// W.nlin = Zm.grad Zp<
				
                if (global.field.anisotropy_dirn == 1)
                    universal->Ring_mult_all(W.nlin1, U.cvf.V1, temp_ring_tr);
                
                else if (global.field.anisotropy_dirn == 1)
                    universal->Ring_mult_all(W.nlin1, U.cvf.V2, temp_ring_tr);
                
                else if (global.field.anisotropy_dirn == 1)
                    universal->Ring_mult_all(W.nlin1, U.cvf.V3, temp_ring_tr);
				
				ring_to_ring_Elsasser_plus(ring_shell_from_i, sector_from_i, Range::all(), Range::all())
                = -temp_ring_tr;
                
			}
		
		// ring_to_ring: Zm to Zm
		for (int ring_shell_from_i = 1; ring_shell_from_i < global.energy_transfer.ring_to_ring.no_shells; ring_shell_from_i++)
			for (int sector_from_i = 1; sector_from_i <= global.energy_transfer.ring_to_ring.no_sectors; sector_from_i++) {
                
				Fill_ring_Vpll(ring_shell_from_i, sector_from_i, W);
				// G = Zm<
                
				Nlin_incompress::Compute_nlin_first_component(U, Giver);
				// U.nlin = Zp.grad Zm<
                
                if (global.field.anisotropy_dirn == 1)
                    universal->Ring_mult_all(U.nlin1, W.cvf.V1, temp_ring_tr);
                
                if (global.field.anisotropy_dirn == 2)
                    universal->Ring_mult_all(U.nlin1, W.cvf.V2, temp_ring_tr);
                
                if (global.field.anisotropy_dirn == 3)
                    universal->Ring_mult_all(U.nlin1, W.cvf.V3, temp_ring_tr);
                
				ring_to_ring_Elsasser_minus(ring_shell_from_i, sector_from_i, Range::all(), Range::all())
                = -temp_ring_tr;	
				
			}
        
        
		MHD::Elsasser_to_UB_field(U, W);													
		// Back to U, B vars
	}	
	
}


//****************************************************************************************

void EnergyTr::Compute_cylindrical_ring_tr_Vpll(FluidVF& U, Pressure& P)
{
	
	cylindrical_ring_to_ring_self = 0.0;
    
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++)
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {
            
			Fill_cylindrical_ring_Vpll(cylindrical_shell_from_i, slab_from_i, U);
            
			Nlin_incompress::Compute_nlin_first_component(U, Giver);
			// nlin = U.grad Um
            
            if (global.field.anisotropy_dirn == 1)
                universal->Cyl_ring_mult_all(U.nlin1, U.cvf.V1, temp_cylindrical_ring_tr);
            
            else  if (global.field.anisotropy_dirn == 2)
                universal->Cyl_ring_mult_all(U.nlin1, U.cvf.V2, temp_cylindrical_ring_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Cyl_ring_mult_all(U.nlin1, U.cvf.V3, temp_cylindrical_ring_tr);
			
			cylindrical_ring_to_ring_self(cylindrical_shell_from_i, slab_from_i, Range::all(), Range::all()) = -temp_cylindrical_ring_tr;
			
		}
    
    if (global.field.anisotropy_dirn == 1)
        universal->Xderiv(U.cvf.V1, global.temp_array.X);
    else if (global.field.anisotropy_dirn == 2)
        universal->Yderiv(U.cvf.V2, global.temp_array.X);
    if (global.field.anisotropy_dirn == 3)
        universal->Zderiv(U.cvf.V3, global.temp_array.X);
    
    P.Compute_pressure(U);
    
    universal->Cyl_ring_mult_all(global.temp_array.X, P.F, cylindrical_ring_force_x_field);
    
}


//****************************************************************************************
//	MHD

void EnergyTr::Compute_cylindrical_ring_tr_Vpll(FluidVF& U, FluidVF& W, Pressure& P)
{
	// U to U
	Compute_cylindrical_ring_tr_Vpll(U, P);
	
	// W to W
	cylindrical_ring_to_ring_VF_WtoW = 0.0;
    
	
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++)
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {
            
			Fill_cylindrical_ring_Vpll(cylindrical_shell_from_i, slab_from_i, W);
			
			Nlin_incompress::Compute_nlin_first_component(U, Giver);
			// U.nlin1 = U.grad Wm_pll
			
            if (global.field.anisotropy_dirn == 1)
                universal->Cyl_ring_mult_all(U.nlin1, W.cvf.V1, temp_cylindrical_ring_tr);
            
            else if (global.field.anisotropy_dirn == 2)
                universal->Cyl_ring_mult_all(U.nlin1, W.cvf.V2, temp_cylindrical_ring_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Cyl_ring_mult_all(U.nlin1, W.cvf.V3, temp_cylindrical_ring_tr);
            
			cylindrical_ring_to_ring_VF_WtoW(cylindrical_shell_from_i, slab_from_i, Range::all(), Range::all()) = -temp_cylindrical_ring_tr;
		}	// of for loop
    
    
    
	// U to W
	cylindrical_ring_to_ring_VF_UtoW = 0.0;
    
	for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++)
		for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {
			Fill_cylindrical_ring_Vpll(cylindrical_shell_from_i, slab_from_i, U);
			
			Nlin_incompress::Compute_nlin_first_component(W, Giver);
			// W.nlin = W.grad Um
            
            if (global.field.anisotropy_dirn == 1)
                universal->Cyl_ring_mult_all(W.nlin1, W.cvf.V1, temp_cylindrical_ring_tr);
            
            else if (global.field.anisotropy_dirn == 2)
                universal->Cyl_ring_mult_all(W.nlin1, W.cvf.V2, temp_cylindrical_ring_tr);
            
            else if (global.field.anisotropy_dirn == 3)
                universal->Cyl_ring_mult_all(W.nlin1, W.cvf.V3, temp_cylindrical_ring_tr);
            
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
				Fill_cylindrical_ring_Vpll(cylindrical_shell_from_i, slab_from_i, U);
				// G = Zp<
                
				Nlin_incompress::Compute_nlin_first_component(W, Giver);
				// W.nlin = Zm.grad Zp<
                
                if (global.field.anisotropy_dirn == 1)
                    universal->Cyl_ring_mult_all(W.nlin1, U.cvf.V1, temp_cylindrical_ring_tr);
                
                else if (global.field.anisotropy_dirn == 2)
                    universal->Cyl_ring_mult_all(W.nlin1, U.cvf.V2, temp_cylindrical_ring_tr);
                
                else if (global.field.anisotropy_dirn == 3)
                    universal->Cyl_ring_mult_all(W.nlin1, U.cvf.V3, temp_cylindrical_ring_tr);
                
				cylindrical_ring_to_ring_Elsasser_plus(cylindrical_shell_from_i, slab_from_i,Range::all(), Range::all()) = -temp_cylindrical_ring_tr;
				
			}		// end of for loop
		
		// ring_to_ring: Zm to Zm
		for (int cylindrical_shell_from_i = 1; cylindrical_shell_from_i < global.energy_transfer.cylindrical_ring_to_ring.no_shells; cylindrical_shell_from_i++)
			for (int slab_from_i = 1; slab_from_i <= global.energy_transfer.cylindrical_ring_to_ring.no_slabs; slab_from_i++) {
                
				Fill_cylindrical_ring_Vpll(cylindrical_shell_from_i, slab_from_i, W);
				// G = Zm<
                
				Nlin_incompress::Compute_nlin_first_component(U, Giver);
				// nlin = Zp.grad Zm<
                
                if (global.field.anisotropy_dirn == 1)
                    universal->Cyl_ring_mult_all(U.nlin1, W.cvf.V1, temp_cylindrical_ring_tr);
                
                else if (global.field.anisotropy_dirn == 2)
                    universal->Cyl_ring_mult_all(U.nlin1, W.cvf.V2, temp_cylindrical_ring_tr);
                
                else if (global.field.anisotropy_dirn == 3)
                    universal->Cyl_ring_mult_all(U.nlin1, W.cvf.V3, temp_cylindrical_ring_tr);
                
				cylindrical_ring_to_ring_Elsasser_minus(cylindrical_shell_from_i, slab_from_i,Range::all(), Range::all()) = -temp_cylindrical_ring_tr;
			}
        
        
		MHD::Elsasser_to_UB_field(U, W);									
		// Back to U, B vars
	}
	
}


//****************************  End of compute_flux.cc ****************************************






