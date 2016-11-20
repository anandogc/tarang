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
/*! \file EnergyTr.cc
 * 
 * @brief  Class declaration of EnergyTr 
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec 2008
 *
 @ @bug  No known bugs
 */
 
#include "EnergyTr.h"
#include "extern_var_incompress.h"

//*********************************************************************************************

EnergyTr::EnergyTr()						 
{
	
	if (global.energy_transfer.turnon) {

		// in global file.. set sphere radii, shells, etc.
		
		flux_self.resize(global.energy_transfer.flux.no_spheres+1);
		sphere_force_x_field.resize(global.energy_transfer.flux.no_spheres+1);

		flux_hk.resize(global.energy_transfer.flux.no_spheres+1);
		
		shelltoshell_self.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
      
        shelltoshell_VF_UtoW.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
        shelltoshell_VF_WtoU.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
        shelltoshell_VF_UtoU.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
		
		shelltoshell_hk.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1); 
		
		temp_shell_tr.resize(global.energy_transfer.shell_to_shell.no_shells+1);
		

		// For Enstrophy and Helicity transfers in Fluid and MHD
        flux_VF_Uin_Uout.resize(global.energy_transfer.flux.no_spheres+1);
		flux_VF_Uin_Wout.resize(global.energy_transfer.flux.no_spheres+1);
        flux_VF_Win_Uout.resize(global.energy_transfer.flux.no_spheres+1);

		flux_VF_Win_Wout.resize(global.energy_transfer.flux.no_spheres+1);
        flux_VF_Bin_Wout.resize(global.energy_transfer.flux.no_spheres+1);
		flux_VF_Bin_Uout.resize(global.energy_transfer.flux.no_spheres+1);

        flux_VF_Jin_Uout.resize(global.energy_transfer.flux.no_spheres+1);


		
		if (global.energy_transfer.ring_to_ring.turnon) {
			ring_to_ring_self.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
			
			temp_ring_tr.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);

			ring_force_x_field.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
		}
			
		
		if (global.energy_transfer.cylindrical_ring_to_ring.turnon) {
			cylindrical_ring_to_ring_self.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
			
			temp_cylindrical_ring_tr.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
			
			cylindrical_ring_force_x_field.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
		}
		
		if (global.energy_transfer.helicity_flux_switch) {
			
		}
		
		
		// FOR SCALAR
		if (global.program.T_exists) {
			flux_SF.resize(global.energy_transfer.flux.no_spheres+1);
					
			shelltoshell_SF.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1); 
			
			if (global.energy_transfer.ring_to_ring.turnon) 
				ring_to_ring_SF.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
			
			if (global.energy_transfer.cylindrical_ring_to_ring.turnon) 
				cylindrical_ring_to_ring_SF.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
		}
		
		// FOR VECTOR FIELD W
		if (global.program.W_exists) {
			flux_VF_Uin_Win.resize(global.energy_transfer.flux.no_spheres+1);
			flux_VF_Win_Uout.resize(global.energy_transfer.flux.no_spheres+1);
			flux_VF_Uout_Wout.resize(global.energy_transfer.flux.no_spheres+1);
			flux_Elsasser_plus.resize(global.energy_transfer.flux.no_spheres+1);
			flux_Elsasser_minus.resize(global.energy_transfer.flux.no_spheres+1);

			flux_Whk.resize(global.energy_transfer.flux.no_spheres+1);
			
			
			shelltoshell_VF_WtoW.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
			shelltoshell_VF_UtoW.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
            shelltoshell_VF_BtoW.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
            shelltoshell_VF_JtoU.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
            shelltoshell_VF_BtoU.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
          
			shelltoshell_Elsasser_plus.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
			shelltoshell_Elsasser_minus.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
			
			shelltoshell_Whk.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
				
			energy_tr_shell_B0.resize(global.energy_transfer.shell_to_shell.no_shells+1,global.energy_transfer.shell_to_shell.no_shells+1);
			
			
			if (global.energy_transfer.ring_to_ring.turnon) {
				ring_to_ring_VF_WtoW.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
				
				ring_to_ring_VF_UtoW.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
				
				ring_to_ring_Elsasser_plus.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
				
				ring_to_ring_Elsasser_minus.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
				
				energy_tr_ring_B0.resize(global.energy_transfer.ring_to_ring.no_shells+1,global.energy_transfer.ring_to_ring.no_sectors+1);
			}
			
			if (global.energy_transfer.cylindrical_ring_to_ring.turnon)  {	
				cylindrical_ring_to_ring_VF_WtoW.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
				
				cylindrical_ring_to_ring_VF_UtoW .resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
				
				cylindrical_ring_to_ring_Elsasser_plus.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
				
				cylindrical_ring_to_ring_Elsasser_minus.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
				
				energy_tr_cylindrical_ring_B0.resize(global.energy_transfer.cylindrical_ring_to_ring.no_shells+1,global.energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
			}
			
		}
	}
}


//***************************   End  of EnergyTr.cc   *****************************************



