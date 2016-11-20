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

/*! \file  compute_force_decay.cc
 * 
 * @brief  F = 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "FORCE.h"

//*********************************************************************************************


void FORCE::Compute_force_decay(FluidVF& U)
{
	if (U.force_switch) {
		U.Force1 = 0.0; 
		U.Force2 = 0.0; 
		U.Force3 = 0.0;
	}	
}

//*********************************************************************************************

void FORCE::Compute_force_decay(FluidVF& U, FluidSF& T)
{
	Compute_force_decay(U);
	
	if (T.force_switch) {
		T.Force = 0.0;
	}	
}

//*********************************************************************************************

void FORCE::Compute_force_decay(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	Compute_force_decay(U);
	
	T1.Force = 0.0;	
	T2.Force = 0.0;	
}


//*********************************************************************************************


void FORCE::Compute_force_decay(FluidVF& U, FluidVF& W)
{
	Compute_force_decay(U);

	if (W.force_switch) {
		W.Force1 = 0.0; 
		W.Force2 = 0.0; 
		W.Force3 = 0.0;
	}	
}


//*********************************************************************************************

void FORCE::Compute_force_decay(FluidVF &U, FluidVF& W, FluidSF& T)
{
	Compute_force_decay(U, W);
	
	if (T.force_switch) {
		T.Force = 0.0;
	}
}


//*********************************************************************************************


void FORCE::Compute_force_pressure_grad(FluidVF& U)
{
	if (U.force_switch) {
		U.Force1 = 0.0; 
		U.Force2 = 0.0;
		U.Force3 = 0.0;
	}	
}

//************************ End of compute_force_decay.cc **************************************


