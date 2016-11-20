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


/*! \file  compute_force_TG.cc
 * 
 * @brief  Set up force using Taylor Green theory.
 *
 * @note Fx = F*sin(k0 x) cos(k0 y) cos(k0 z) <BR>
 *		Fy = -F*cos(k0 x)sin(k0 y) cos(k0 z) <BR>
 *		Fz = 0
 *
 *  @note is_force_field_para_read is used to read only once.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "FORCE.h"

//*********************************************************************************************


void FORCE::Compute_force_Taylor_Green(FluidVF& U)
{

	if (!global.force.configuration_done)	{
		int k0 = ((int) global.force.double_para(0));
		Real force_amp = global.force.double_para(1);
		global.force.configuration_done = true; // To read only once.
		
		Setup_Taylor_Green_force_field(U, k0, force_amp);
	}	
	
}


//*********************************************************************************************
//
// Scalar
//

void FORCE::Compute_force_Taylor_Green(FluidVF& U, FluidSF& T)
{
	if (T.force_switch == 1)
		T.Force = 0.0;
	
	Compute_force_Taylor_Green(U);
}


//*********************************************************************************************
//
// force both V and W
//

void FORCE::Compute_force_Taylor_Green(FluidVF& U, FluidVF& W)
{
	
	if (!global.force.configuration_done) {
		int k0 = int(global.force.double_para(0));
		Real force_amp = global.force.double_para(1);
		Real forceW_amp = global.force.double_para(2);
		
		Setup_Taylor_Green_force_field(U, k0, force_amp);
		
		if (W.force_switch == 1)
			Setup_Taylor_Green_force_field(W, k0, forceW_amp);
		
		global.force.configuration_done = true;
	}
		
}


//*********************************************************************************************	
//
//  Force both V and W
//	
	
void FORCE::Compute_force_Taylor_Green(FluidVF& U, FluidVF& W, FluidSF& T)
{
	if (T.force_switch == 1)
		T.Force = 0.0;
	
	Compute_force_Taylor_Green(U, W);
}	


//****************************** End of compute_force_TG.cc ***********************************


