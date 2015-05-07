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

/*! \file  fill_sphere_shell.cc
 * 
 * @brief  Fill sphere or shell in a region.
 * 
 * Sphere_index / shell_index is passed to the function as a parameter.  The function
 *		picks the inner and outer radii and fill the given region. <BR>
 *		Notation shell(n) = (R(n-1),R(n)] (excluding the modes in the inner sphere and 
 *		including the modes in the outer sphere. <BR>
 *		The last shell contains all the modes beyond the maximum inner sphere.
 *
 * @sa void Fill_array_shell(string basis_type, int N[], Array<Complex,2> A, Array<Complex,2> B,
 *		Real inner_radius, Real outer_radius, 	Real kfactor[])
 *
 * @sa universal_ET.cc
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "EnergyTr.h"


//*********************************************************************************************

///  FOR FLUID --- INSIDE SPHERE FILLED
///	 Includes the modes of the outer surface.
void EnergyTr::Fill_in_sphere(int sphere_index, FluidVF& U)						
{
	Real inner_radius = 0;
	Real outer_radius = global.energy_transfer.flux.radii(sphere_index);				

	universal->Fill_array_shell(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius);
	universal->Fill_array_shell(U.cvf.V2, Giver.cvf.V2, inner_radius, outer_radius);
	universal->Fill_array_shell(U.cvf.V3, Giver.cvf.V3, inner_radius, outer_radius);
}


// For Passive Scalar -- INSIDE SPHERE FILLED

void EnergyTr::Fill_in_sphere(int sphere_index, FluidSF& T)						
{
		
	Real inner_radius = 0;
	Real outer_radius = global.energy_transfer.flux.radii(sphere_index);				
			
	universal->Fill_array_shell(T.csf.F, Giver.cvf.V1, inner_radius, outer_radius);			
}



void EnergyTr::Fill_in_sphere_Vpll(int sphere_index, FluidVF& U)
{
    
	Real inner_radius = 0;
	Real outer_radius = global.energy_transfer.flux.radii(sphere_index);
    
    if (global.field.anisotropy_dirn == 1)
        universal->Fill_array_shell(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius);
    
    else if (global.field.anisotropy_dirn == 2)
        universal->Fill_array_shell(U.cvf.V2, Giver.cvf.V1, inner_radius, outer_radius);
    
    else if (global.field.anisotropy_dirn == 3)
        universal->Fill_array_shell(U.cvf.V3, Giver.cvf.V1, inner_radius, outer_radius);
}


/*====================================================================================

void EnergyTr::Fill_out_sphere(): -- OUTSIDE SPHERE FILLED
	 
// Modes on the surface are included outside the sphere

======================================================================================*/ 

// FOR FLUID --- OUTSIDE SPHERE FILLED

void EnergyTr::Fill_out_sphere(int sphere_index, FluidVF& U)
{
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = INF_RADIUS;

	universal->Fill_array_shell(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius);
	universal->Fill_array_shell(U.cvf.V2, Giver.cvf.V2, inner_radius, outer_radius);
	universal->Fill_array_shell(U.cvf.V3, Giver.cvf.V3, inner_radius, outer_radius);
}


// For Passive Scalar 
void EnergyTr::Fill_out_sphere(int sphere_index, FluidSF& T)						
{
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = INF_RADIUS;
			
	universal->Fill_array_shell(T.csf.F, Giver.cvf.V1, inner_radius, outer_radius);	
}



void EnergyTr::Fill_out_sphere_Vpll(int sphere_index, FluidVF& U)
{
    
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = INF_RADIUS;
    
    if (global.field.anisotropy_dirn == 1)
        universal->Fill_array_shell(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius);
    
    else if (global.field.anisotropy_dirn == 2)
        universal->Fill_array_shell(U.cvf.V2, Giver.cvf.V1, inner_radius, outer_radius);
    
    else if (global.field.anisotropy_dirn == 3)
        universal->Fill_array_shell(U.cvf.V3, Giver.cvf.V1, inner_radius, outer_radius);
}


//*********************************************************************************************

/// Fluid -- fill shell
/// Include the modes of the outer surface, and exclude that of inner surface.
void EnergyTr::Fill_shell(int shell_from_index, FluidVF& U)						
{
	Real inner_radius = global.energy_transfer.shell_to_shell.radii(shell_from_index - 1);
	Real outer_radius = global.energy_transfer.shell_to_shell.radii(shell_from_index);
	
	universal->Fill_array_shell(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius);
	universal->Fill_array_shell(U.cvf.V2, Giver.cvf.V2, inner_radius, outer_radius);
	universal->Fill_array_shell(U.cvf.V3, Giver.cvf.V3, inner_radius, outer_radius);
}


// Passive scalar
void EnergyTr::Fill_shell(int shell_from_index, FluidSF& T)						
{
	Real inner_radius = global.energy_transfer.shell_to_shell.radii(shell_from_index - 1);
	Real outer_radius = global.energy_transfer.shell_to_shell.radii(shell_from_index);
			
	universal->Fill_array_shell(T.csf.F, Giver.cvf.V1, inner_radius, outer_radius);
}

void EnergyTr::Fill_shell_Vpll(int shell_from_index, FluidVF& U)
{
	Real inner_radius = global.energy_transfer.shell_to_shell.radii(shell_from_index - 1);
	Real outer_radius = global.energy_transfer.shell_to_shell.radii(shell_from_index);
    
    if (global.field.anisotropy_dirn == 1)
        universal->Fill_array_shell(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius);
    
    else if (global.field.anisotropy_dirn == 2)
        universal->Fill_array_shell(U.cvf.V2, Giver.cvf.V1, inner_radius, outer_radius);
    
    else if (global.field.anisotropy_dirn == 3)
        universal->Fill_array_shell(U.cvf.V3, Giver.cvf.V1, inner_radius, outer_radius);
}


//*********************************************************************************************

// Spherical rings

void EnergyTr::Fill_ring(int ring_shell_from_index, int sector_from_index, FluidVF& U)						
{
	Real inner_radius = global.energy_transfer.ring_to_ring.radii(ring_shell_from_index - 1);
	Real outer_radius = global.energy_transfer.ring_to_ring.radii(ring_shell_from_index);
	
	Real left_angle = global.energy_transfer.ring_to_ring.sector_angles(sector_from_index-1);
	Real right_angle = global.energy_transfer.ring_to_ring.sector_angles(sector_from_index);	
	
	
	universal->Fill_array_ring(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius, left_angle, right_angle);
	universal->Fill_array_ring(U.cvf.V2, Giver.cvf.V2, inner_radius, outer_radius, left_angle, right_angle);
	universal->Fill_array_ring(U.cvf.V3, Giver.cvf.V3, inner_radius, outer_radius, left_angle, right_angle);
}

	// Passive scalar

void EnergyTr::Fill_ring(int ring_shell_from_index, int sector_from_index, FluidSF& T)						
{
	Real inner_radius = global.energy_transfer.ring_to_ring.radii(ring_shell_from_index - 1);
	Real outer_radius = global.energy_transfer.ring_to_ring.radii(ring_shell_from_index);
	
	Real left_angle = global.energy_transfer.ring_to_ring.sector_angles(sector_from_index-1);
	Real right_angle = global.energy_transfer.ring_to_ring.sector_angles(sector_from_index);
	
	universal->Fill_array_ring(T.csf.F, Giver.cvf.V1, inner_radius, outer_radius, left_angle, right_angle);	
}


void EnergyTr::Fill_ring_Vpll(int ring_shell_from_index, int sector_from_index, FluidVF& U)
{
	Real inner_radius = global.energy_transfer.ring_to_ring.radii(ring_shell_from_index - 1);
	Real outer_radius = global.energy_transfer.ring_to_ring.radii(ring_shell_from_index);
	
	Real left_angle = global.energy_transfer.ring_to_ring.sector_angles(sector_from_index-1);
	Real right_angle = global.energy_transfer.ring_to_ring.sector_angles(sector_from_index);
	
    if (global.field.anisotropy_dirn == 1)
        universal->Fill_array_ring(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius, left_angle, right_angle);
    
    else if (global.field.anisotropy_dirn == 2)
        universal->Fill_array_ring(U.cvf.V2, Giver.cvf.V1, inner_radius, outer_radius, left_angle, right_angle);
    
    else if (global.field.anisotropy_dirn == 3)
        universal->Fill_array_ring(U.cvf.V3, Giver.cvf.V1, inner_radius, outer_radius, left_angle, right_angle);
}



//****************************************************************************************

// Cylindrical rings

void EnergyTr::Fill_cylindrical_ring(int cylindrical_shell_from_index, int slab_from_index, FluidVF& U)			
{	
	
	Real inner_radius = global.energy_transfer.cylindrical_ring_to_ring.radii(cylindrical_shell_from_index - 1);
	Real outer_radius = global.energy_transfer.cylindrical_ring_to_ring.radii(cylindrical_shell_from_index);
	
	Real h_lower = global.energy_transfer.cylindrical_ring_to_ring.kpll_array(slab_from_index-1);
	Real h_upper = global.energy_transfer.cylindrical_ring_to_ring.kpll_array(slab_from_index);
	
	universal->Fill_array_cylindrical_ring(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius, h_lower, h_upper);
	universal->Fill_array_cylindrical_ring(U.cvf.V2, Giver.cvf.V2, inner_radius, outer_radius, h_lower, h_upper);
	universal->Fill_array_cylindrical_ring(U.cvf.V3, Giver.cvf.V3, inner_radius, outer_radius, h_lower, h_upper);
}


// Scalar

void EnergyTr::Fill_cylindrical_ring(int cylindrical_shell_from_index, int slab_from_index, FluidSF& T)			
{	
	Real inner_radius = global.energy_transfer.cylindrical_ring_to_ring.radii(cylindrical_shell_from_index - 1);
	Real outer_radius = global.energy_transfer.cylindrical_ring_to_ring.radii(cylindrical_shell_from_index);
	
	Real h_lower = global.energy_transfer.cylindrical_ring_to_ring.kpll_array(slab_from_index-1);
	Real h_upper = global.energy_transfer.cylindrical_ring_to_ring.kpll_array(slab_from_index);
	
	universal->Fill_array_cylindrical_ring(T.csf.F, Giver.cvf.V1, inner_radius, outer_radius, h_lower, h_upper);
}

void EnergyTr::Fill_cylindrical_ring_Vpll(int cylindrical_shell_from_index, int slab_from_index, FluidVF& U)
{
	Real inner_radius = global.energy_transfer.cylindrical_ring_to_ring.radii(cylindrical_shell_from_index - 1);
	Real outer_radius = global.energy_transfer.cylindrical_ring_to_ring.radii(cylindrical_shell_from_index);
	
	Real h_lower = global.energy_transfer.cylindrical_ring_to_ring.kpll_array(slab_from_index-1);
	Real h_upper = global.energy_transfer.cylindrical_ring_to_ring.kpll_array(slab_from_index);
	
    if (global.field.anisotropy_dirn == 1)
        universal->Fill_array_cylindrical_ring(U.cvf.V1, Giver.cvf.V1, inner_radius, outer_radius, h_lower, h_upper);
    
    else if (global.field.anisotropy_dirn == 2)
        universal->Fill_array_cylindrical_ring(U.cvf.V2, Giver.cvf.V1, inner_radius, outer_radius, h_lower, h_upper);
    
    else if (global.field.anisotropy_dirn == 3)
        universal->Fill_array_cylindrical_ring(U.cvf.V3, Giver.cvf.V1, inner_radius, outer_radius, h_lower, h_upper);
}


//*****************************  End of fill_sphere_shell.cc **********************************




