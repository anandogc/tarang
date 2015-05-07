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


/*! \file  init_cond_ABC.cc
 * 
 * @brief Initial conditions as ABC flow.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FFF basis; not for SCFT basis.
 *
 * @note \f$ V_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ V_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ V_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "FluidIO.h"


//*********************************************************************************************

void FluidIO::Init_cond_ABC(FluidVF& U)
{
	int k0 = (int) global.io.double_para(0);
	Real A = global.io.double_para(1);
	Real B = global.io.double_para(2);
	Real C = global.io.double_para(3);
	Real amp = global.io.double_para(4);
	
	Setup_ABC_field(U, k0, amp, A, B, C);
}

//*********************************************************************************************
// Scalar 
//

void FluidIO::Init_cond_ABC(FluidVF& U, FluidSF& T)
{
	Init_cond_ABC(U);
		
	T.csf.F = 0.0;
}


//*********************************************************************************************
// Vector
//

void FluidIO::Init_cond_ABC(FluidVF& U, FluidVF& W)
{
	int k0 = int(global.io.double_para(0));
	Real A = global.io.double_para(1);
	Real B = global.io.double_para(2);
	Real C = global.io.double_para(3);
	Real amp = global.io.double_para(4);
	Real ampW = global.io.double_para(5);
	
	Setup_ABC_field(U, k0, amp, A, B, C);
	Setup_ABC_field(W, k0, ampW, A, B, C);

}

//*********************************************************************************************
//	Vector+scalar
//

void FluidIO::Init_cond_ABC(FluidVF& U, FluidVF& W, FluidSF& T)
{
	Init_cond_ABC(U, W);
		
	T.csf.F = 0.0;
}



//******************************** End of Init_cond_ABC.cc ************************************

