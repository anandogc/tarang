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
/*! \file  compute_force_ABC.cc
 * 
 * @brief  Set up force using ABC model.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FOUR basis; not for SCFT basis.
 *
 * @note \f$ F_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ F_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ F_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 *  @note is_force_field_para_read is used to read only once.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "FORCE.h"


//*********************************************************************************************

void FORCE::Compute_force_ABC(FluidVF& U)
{
	
	if (!global.force.configuration_done) {
		int k0 = ((int) global.force.double_para(0));
		Real A = global.force.double_para(1);
		Real B = global.force.double_para(2);
		Real C = global.force.double_para(3);
		Real force_amp = global.force.double_para(4);
		
		Setup_ABC_force_field(U, k0, force_amp, A, B, C);
	
		global.force.configuration_done = true;
	}	
	
}


//*********************************************************************************************
//
// Scalar
//

void FORCE::Compute_force_ABC(FluidVF& U, FluidSF& T)
{

	if (!global.force.configuration_done) {
		T.Force = 0.0;
		
		Compute_force_ABC(U);
		global.force.configuration_done = true;	
	}
}


//*********************************************************************************************
//
// Force both V and W
//

void FORCE::Compute_force_ABC(FluidVF&U, FluidVF& W)
{	

	if (!global.force.configuration_done) {

		int k0 = ((int) global.force.double_para(1));
		Real A = global.force.double_para(2);
		Real B = global.force.double_para(3);
		Real C = global.force.double_para(4);
		Real force_amp = global.force.double_para(5);
		Real forceW_amp = global.force.double_para(6);
		
		Setup_ABC_force_field(U, k0, force_amp, A, B, C);
		Setup_ABC_force_field(W, k0, forceW_amp, A, B, C);
		
		global.force.configuration_done = true;
	}	
	
}


//*********************************************************************************************
//
// Vector + scalar
//

void FORCE::Compute_force_ABC(FluidVF& U, FluidVF& W, FluidSF& T)
{	

	if (!global.force.configuration_done) {
		T.Force = 0.0;
		Compute_force_ABC(U, W);
		
		global.force.configuration_done = true;
	}

}	

//*****************************   End of compute_force_ABC.cc *********************************



