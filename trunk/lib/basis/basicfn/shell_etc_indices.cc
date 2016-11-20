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


/*! \file	basis_basicfn.cc
 * 
 * @brief basic functions definitions that are common to all basis functions.
 *
 * @sa basis_basicfn.h
 *
 * @version 4.0 Parallel version
 * @author  M. K. Verma
 * @date	Sept 2008
 * @bug		The Transpose for 2D is not working at present.. it must transport 
 *				real array.
 */ 


#include "shell_etc_indices.h"

//
//
//*********************************************************************************************


int Get_shell_index(Real kkmag, Array<Real, 1> shell_radius_array)
{
	
	return first(shell_radius_array >= kkmag);
}

//
//*********************************************************************************************

int Get_sector_index(Real theta, Array<Real, 1> sector_angle_array)
{
	int index;
	
	if (theta > M_PI/2)
		theta = M_PI - theta;		// 0 <= theta <= M_PI/2
	
	return (abs(theta) < MYEPS)	? 1 :  first(sector_angle_array >= theta);
}

//
//*********************************************************************************************

void Compute_ring_index
(
	Real kkmag, 
	Real theta, 
	Array<Real, 1> shell_radius_array, 
	Array<Real, 1> sector_angle_array, 
	int& shell_index, 
	int& sector_index
)
{	
	shell_index = first(shell_radius_array >= kkmag);

	if (theta > M_PI/2)
		theta = M_PI - theta;		// 0 <= theta <= M_PI/2
	
	sector_index = ((abs(theta) < MYEPS) ? 1 : first(sector_angle_array >= theta));
}


//
//*********************************************************************************************

int Get_slab_index(Real kkpll, Real kkperp, Array<Real,1> cylinder_kpll_array)
{
	
	return (abs(kkpll) < MYEPS) ? 1 : first(cylinder_kpll_array >= abs(kkpll));
	
}

//
//*********************************************************************************************

void Compute_cylindrical_ring_index
(
	Real kkpll, 
	Real kkperp, 
	Array<Real,1> cylinder_shell_radius_array,
	Array<Real,1> cylinder_kpll_array, 
	int& shell_index, 
	int& slab_index
)	
{

	shell_index = first(cylinder_shell_radius_array >= kkperp);
	
	slab_index = (abs(kkpll) < MYEPS) ? 1 : first(cylinder_kpll_array >= abs(kkpll));
	
}


//***********************************  End of basis_basicfn.cc ********************************



