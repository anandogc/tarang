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

//
// Vector
//

void FluidIO::Init_cond_dynamo_full_velocity_field(FluidVF& U, FluidVF& W)
{
	Init_cond_complex_field(U);  // Read the velocity field..
	
// Init_cond_double_para(1,2,3) = totalEb,	kkmin, kkmax
// Distribute totalEb among the modes in 
	
	Real totalEb = global.io.double_para(0);
	Real Kmin = global.io.double_para(1);
	Real Kmax = global.io.double_para(2);
	
	int total_no_modes = (4*M_PI/3)*(my_pow(Kmax,3) - my_pow(Kmin,3));
	Real ekW = totalEb/total_no_modes;   // Energy per mode
	
	Real Kmag, ampW, phase1W, phase2W, phase3W;
	int index, maxN3;
	
    for (int lx=0; lx<global.field.maxlx; lx++)
	    for (int ly=0; ly<global.field.maxly; ly++)
	        for (int lz=0; lz<global.field.maxlz; lz++) {
                
				Kmag = universal->Kmagnitude(lx, ly, lz);
				
				if ((Kmag > Kmin) && (Kmag <= Kmax)) {
					index = (int) ceil(Kmag);
					
					ampW = sqrt(2*ekW);
					
					phase1W = 2*M_PI * SPECrand.random();
					phase2W = 2*M_PI * SPECrand.random();
					phase3W = 2*M_PI * SPECrand.random();
					
					Put_vector_amp_phase_comp_conj(W, lx, ly, lz, ampW, phase1W, phase2W, phase3W);
				}	
			}
	
	if (my_id == master_id) 
		W.cvf.V1(0,0,0) = W.cvf.V2(0,0,0) = W.cvf.V3(0,0,0) = 0.0;
	
}



//******************************** End of Init_cond_TG.cc *************************************







