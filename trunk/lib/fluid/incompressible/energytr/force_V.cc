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

#include "EnergyTr.h"


//*********************************************************************************************
// Force*V  or Force*T.F (for the flux diagnostics). Therefore we use global.energy_transfer.flux.radii for these computations.
//*********************************************************************************************

void EnergyTr::Power_supply_within_sphere(FluidVF& U)
{
	if (U.force_switch)
		universal->Shell_mult_all(U.Force1, U.Force2, U.Force3, U.cvf.V1,  U.cvf.V2,  U.cvf.V3, global.energy_transfer.flux.radii, sphere_force_x_field);
												
}


void EnergyTr::Power_supply_within_sphere(FluidSF& T)
{
	if (T.force_switch)
		universal->Shell_mult_all(T.Force, T.csf.F, global.energy_transfer.flux.radii, sphere_force_x_field);
	
}



/**********************************************************************************************

						FOR THE RINGS
	
****************************************************************************************/

void EnergyTr::Power_supply_ring(FluidVF& U)
{
	if (U.force_switch)
		universal->Ring_mult_all(U.Force1, U.Force2, U.Force3, U.cvf.V1,  U.cvf.V2,  U.cvf.V3, ring_force_x_field);
}

void EnergyTr::Power_supply_ring(FluidSF& T)
{
	if (T.force_switch)
		universal->Ring_mult_all(T.Force, T.csf.F, ring_force_x_field);
}

/**********************************************************************************************

						FOR THE CYLINDRICAL RINGS
	
***********************************************************************************************/


void EnergyTr::Power_supply_cylindrical_ring(FluidVF& U)
{
	if (U.force_switch)
		universal->Cyl_ring_mult_all(U.Force1, U.Force2, U.Force3, U.cvf.V1,  U.cvf.V2,  U.cvf.V3, cylindrical_ring_force_x_field);
}

void EnergyTr::Power_supply_cylindrical_ring(FluidSF& T)
{
	if (T.force_switch)
		universal->Cyl_ring_mult_all(T.Force, T.csf.F, cylindrical_ring_force_x_field);
}



//***********************************  End of nusselt.cc **************************************



