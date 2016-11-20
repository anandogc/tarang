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

/*! \file  dynamo_6mode.cc
 * 
 * @brief Initial conditions for dynamo 6 model of Verma et al.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FOUR basis; not for SCFT basis.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug No known bugs
 */

#include "FluidIO.h"

//*********************************************************************************************

void FluidIO::Init_cond_DYNAMO_SIX_MODE(FluidVF& U, FluidVF& W)
{
	int k0 = int(global.io.double_para(0));
	Real amp101 = global.io.double_para(1);
	Real amp011 = global.io.double_para(2);
	Real amp112 = global.io.double_para(3);
	Real h = global.io.double_para(4);
	
	Real ampW101 = global.io.double_para(5);
	Real ampW011 = global.io.double_para(6);
	Real ampW112 = global.io.double_para(7);
	Real hW = global.io.double_para(8);
	
	Setup_SIX_MODE_field(U, k0, amp101, amp011, amp112, h);
	Setup_SIX_MODE_field(W, k0, ampW101, ampW011, ampW112, hW);
}

//******************************* End of Init_cond_dynamo_6mode.cc ****************************


