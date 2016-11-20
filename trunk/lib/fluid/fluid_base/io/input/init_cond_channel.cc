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

/*! \file  init_cond_RB_misc.cc
 * 
 * @brief Initial conditions as RB convection: Lorenz equations.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FFF basis; not for SCFT basis.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 */

#include "FluidIO.h"


//**********************************************************************************************


/** @brief Set the modes for Rayleigh Taylor instability
 * 
 * @return  rho = 0 for x=0:0.5, and 1 for x=0.5:1.
 *		in T:  1 for x=0:0.5, and 0 for x=0.5:1.
 *		we model by T(x) = ( 1-tanh((x-0.5)*slope) )/2 with slope = large no.
 *		theta = T(x) -1+x.
 */
void  FluidIO::Init_cond_channel_flow(FluidVF& U)
{
 /*   int rx;
    Real x;
	
	U.cvf.V1 = 0.0;
	U.cvf.V2 = 0.0;
	U.cvf.V3 = 0.0;
    
	for (int lx=0; lx<global.field.maxlx; lx++) {
        rx = universal->Get_rx_real_space(lx);
        x = cos((Nx-rx)*M_PI/Nx);  // theta = (Nx-rx)*pi/Nx
		real(U.rvf.V2r(Range::all(),Range(0,Nz/2-1),lx)) = 1-x*x;
        imag(U.rvf.V2r(Range::all(),Range(0,Nz/2-1),lx)) = 1-x*x;
	}

    universal->Forward_transform(U.rvf.V2r, U.cvf.V2);
    //global.io.double_para.resize(1);
    //global.io.double_para(0) = 0.01;
    // Init_cond_energy_helicity_spectrum(U); */
 
}


//******************************** End of Init_cond_RB_misc.cc ********************************
