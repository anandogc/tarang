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

/*! \file array_basic_inline.h
 * 
 * @brief Some basic Array operations: Array_real_mult, Output_asreal, 
 *			Model_initial_shell_spectrum
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug  No know bug
 */ 


#ifndef _FIELD_ARRAY_BASIC_INLINE_H
#define _FIELD_ARRAY_BASIC_INLINE_H

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif


#include "../../basis_basicfn/basis_basicfn_inline.h"
#include "../../basis_basicfn/basis_basicfn.h"


using namespace blitz;

//******************************************************************************

/// return \f$ |\vec{r}_2 - \vec{r}_1| \f$ for 2D.
inline Real Get_abs_dr
(
	TinyVector<Real,2> r1, 
	TinyVector<Real,2> r2
) 
{
	return sqrt(pow2(r2(0)-r1(0)) + pow2(r2(1)-r1(1)));
}


/// return \f$ |\vec{r}_2 - \vec{r}_1| \f$ for 3D.
inline Real Get_abs_dr
(	
	TinyVector<Real,3> r1, 
	TinyVector<Real,3> r2, 
	TinyVector<Real,3> r3
) 
{
	return sqrt(pow2(r2(0)-r1(0)) + pow2(r2(1)-r1(1)) + pow2(r2(2)-r1(2)));
}


//******************************************************************************

/** @brief Maximum radius inside the real sphere.
 *	Along the periodic dirn, max(r)=N[i]/2. \\
 *	Along the nonperiodic dirn, max(r)=N[i]-1. 
 *
 *	@return  min of max(r) along the three axis.
 */
inline Real Max_radius_inside_real_space(string basis_type, int N[], Real Delta_x[])
{
	if (N[2] > 1) // 3D case
		if (basis_type == "FOUR")	{
			Real xmag = min( (N[1]/2)*Delta_x[1], (N[2]/2)*Delta_x[2] );  
			xmag = min(xmag, (N[3]/2)*Delta_x[3]);
			
			return xmag;	
		}
		
		else if (basis_type == "SCFT")	{
			Real xmag = min( (N[1]-1)*Delta_x[1], (N[2]/2)*Delta_x[2] );  
			xmag = min(xmag, (N[3]/2)*Delta_x[3]);
			
			return xmag;	
		}
		
		else
			return 0; // for -Wall
	
	else if (N[2] == 1) // 2D case
		if (basis_type == "FOUR")
			return min( (N[1]/2)*Delta_x[1], (N[3]/2)*Delta_x[3] ); 
		
		else if (basis_type == "SCFT")
			return  min( (N[1]-1)*Delta_x[1], (N[3]/2)*Delta_x[3] ); 
		
		else 
			return 0;  // for -Wall
}


//******************************************************************************

inline Real Min_radius_outside_real_space(string basis_type, int N[], Real Delta_x[])
{
	if (basis_type == "FOUR")
		return  sqrt( pow2((N[1]/2)*Delta_x[1]) + pow2((N[2]/2)*Delta_x[2]) 
						+ pow2((N[3]/2)*Delta_x[3]) );
						
	else if (basis_type == "SCFT")
		return  sqrt( pow2((N[1]-1)*Delta_x[1]) + pow2((N[2]/2)*Delta_x[2]) 
						+ pow2((N[3]/2)*Delta_x[3]) );	
						
	else
		return 0; // for -Wall																	
}



#endif

//*********************************  End of array_basic_inline.h ******************************



