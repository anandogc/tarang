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

/*! \file  compute_force_dynamo_6mode.cc
 * 
 * @brief  Set up force for 6 mode dynamo model of Verma et al.
 *
 *  @note is_force_field_para_read is used to read only once.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug No known bugs
 */

#include "FORCE.h"

//*********************************************************************************************


void FORCE::Compute_force_DYNAMO_SIX_MODE(FluidVF& U, FluidVF& W)
{
	if (U.force_switch == 1) {
		
		if (!global.force.configuration_done) {
			int k0 = ((int) global.force.double_para(0));
			Real amp101 = global.force.double_para(1);
			Real amp011 = global.force.double_para(2);
			Real amp112 = global.force.double_para(3);
			Real h = global.force.double_para(4);
			
			Setup_SIX_MODE_force_field(U, k0, amp101, amp011, amp112, h); 
			// force only u field		
			
			global.force.configuration_done = true;
		}
	}		
}



//****************************** End of compute_force_dynamo_6mode.cc *************************


