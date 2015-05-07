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

/*! \files	struct_fn.h
 * 
 * @brief Structure function calculations
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */ 


#ifndef _STRUCT_FN_H
#define _STRUCT_FN_H

#include "Rvf.h"
#include "Rsf.h"

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif


using namespace blitz;

//*********************************************************************************************

/*! @brief 3D: Compute \f$ dr(i) = dj(i)*xfactor(i) \f$ for a pair of points
 * 
 *	@param	TinyVector<int,2> j1, j2		Grid index of the two points
 *
 *	@return TinyVector<Real,2> dr			dr vector between the two point.
 *  @note	For periodic direction: dr(i)=(j2(i)-j1(i))*xfactor(i) 
 *				if diff <= N[i]/2;
 *				Otherwise dr(i) =(j2(i)-j1(i)-N[i]/2)*xfactor(i) 
 */
void Compute_dr
(
	string basis_type,
	int N[], 
	TinyVector<int,3> j1, TinyVector<int,3> j2, 
	TinyVector<Real,3> &dr,
	Real xfactor[]
);


//******************************************************************************

/** @brief 3D: Compute the structure function of vector field 
 *					\f$ \vec{v}(\vec{x}) \f$.
 *
 * @param	Vx, Vy, Vz  Velocity field; Real values stored in complex arrays.
 * @param	N[]  The size of the arrays A.
 * @param	structurefn_q_min
 * @param	structurefn_q_max
 * @param	xfactor[i] = min_dx[i]; 
 *			we ignore \f$ 2 \pi \f$ which is the overall factor.
 *
 * @return	\f$ St(r, q, 0) = |\Delta v_{||}|^q \f$  with r=ceil(|r2-r1|).
 * @return	\f$ St(r, q, 1) = |\Delta v_{\perp}|^q \f$  
 *
 * @note	Loop over all points P & Q (Q>=P); Find r between the two points, 
 *			and compute St(r,q,0-1) by summing over all r's.  
 * @note	When P=Q, the real and imag parts are two different points in R.
 * @note	The final St(r,q,0-1) is divided by total no of points with sep r.
 */

//
/// Same as above but for scalar in 3D
/* void Compute_local_structure_function
(
	string basis_type,
	int N[], 
	Array<Complex,3> F	
); */


#endif

//*********************************   End of struct_fn.h  *************************************



