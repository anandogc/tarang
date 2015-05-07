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

/*! \file  PlainCsf.h
 * 
 * @brief  Class declaration of Csf, a Complex Scalar Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 */

#ifndef _PLAINCSF_H
#define _PLAINCSF_H

#include "basis.h"

//*********************************************************************************************

//! Complex Scalar field
/*!
 * D dimensional complex scalar field with D array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, ...N_D/2+1] \f$.   <BR>
 * This array can also store real values, and its dimension is \f$[N_1, N_2, ...N_D] \f$. 
 *
 */

//*********************************************************************************************

class PlainCSF
{ 

public:
  
	//!  \f$F(N_1, N_2, N_3/2+1) \f$.
	Array<Complex,3> F;				
																			
	
//********************************************************************************************* 

public:

	/*! A constructor; 
	 *
	 *  Allocation for F, Ek etc. 
	 * Initialize the arrays.
	 */
	PlainCSF();


void Initialize();

};
#endif

//********************************* End of CSF class decl *************************************
