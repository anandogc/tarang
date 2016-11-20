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

/*! \file  PlainCvf.h
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 */


#ifndef _PLAIN_CVF
#define _PLAIN_CVF

#include "basis.h"

//*********************************************************************************************

//! Complex vector field
/*!
 * D dimensional complex vector field has D components with D array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, ...N_D/2+1] \f$.   <BR>
 * These arrays can also store real values, and their dimension are \f$[N_1, N_2, ...N_D] \f$. 
 *
 */

//*********************************************************************************************

class PlainCVF 
{ 
	public:
	 
	//!  \f$ V_x(N_1, N_2, N_3/2+1) \f$.
	Array<Complex,3> V1;							
	
	//!  \f$ V_y(N_1, N_2, N_3/2+1) \f$.
	Array<Complex,3> V2;							
	
	//!  \f$ V_z(N_1, N_2, N_3/2+1) \f$.
	Array<Complex,3> V3;																
								
	
//*********************************************************************************************
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi etc. 
	 * Initialize the arrays.
	 */
	PlainCVF();
	
	void Initialize();
};
#endif


//******************************** End of CVF.h ***********************************************	


