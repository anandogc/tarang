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


/*! \file  init_cond_TG.cc
 * 
 * @brief Initial conditions as Taylor Green flow (TG).
 *
 * @note The parameters are read from parameter file. 
 *
 *		Vx = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Vy = -amp*cos(k0 x)sin(k0 y) cos(k0 z)
 *		Vz = 0
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "FluidIO.h"


//*********************************************************************************************


void FluidIO::Init_cond_Taylor_Green(FluidVF& U)
{
	
	int k0 = int(global.io.double_para(0));
	Real amp = global.io.double_para(1);
	
	Setup_Taylor_Green_field(U, k0, amp);
}

//
// Scalar 
//


void FluidIO::Init_cond_Taylor_Green(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR") 
		Init_cond_Taylor_Green_scalar(U, T);
	
	else if (global.program.kind == "RBC") 
		Init_cond_Taylor_Green_RBC(U, T);
}


void FluidIO::Init_cond_Taylor_Green_scalar(FluidVF& U, FluidSF& T)
{	
	Init_cond_Taylor_Green(U);
	
	T.csf.F = 0.0;
}


void  FluidIO::Init_cond_Taylor_Green_RBC(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO") {
		Init_cond_Taylor_Green(U);	
		
		U.Zero_Prandtl_number_compute_temperature(T);
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
		cout << "TG forcing not allowed for PRINFTY option" << endl;
		exit(0);
	}
	
	else
		Init_cond_Taylor_Green(U, T);
}

//
// Vector
//

void FluidIO::Init_cond_Taylor_Green(FluidVF& U, FluidVF& W)
{

	int k0 = int(global.io.double_para(0));
	Real amp = global.io.double_para(1);
	Real ampW = global.io.double_para(2);
	
	Setup_Taylor_Green_field(U, k0, amp);
	Setup_Taylor_Green_field(W, k0, ampW);		
}

//
//	Vector+scalar
//

void FluidIO::Init_cond_Taylor_Green(FluidVF& U, FluidVF& W, FluidSF& T)
{	
	Init_cond_Taylor_Green(U, W);
	
	T.csf.F = 0.0;

}



//******************************** End of Init_cond_TG.cc *************************************







