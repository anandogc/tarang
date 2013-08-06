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

/*! \file  Rsf.h
 * 
 * @brief  Class declaration of Rsf, a Real Scalar Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * 
 * @bug  No known bugs
 */
 

#ifndef _RSF
#define _RSF

#include "basis.h"
class CSF;

//*********************************************************************************************

//! Real Scalar field
/*!
 * 3 dimensional real scalar field with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, N_3/2+1] \f$.   <BR>
 * This array contains two real values for every complex value.
 *  The dimension of the real array is \f$[N_1, N_2, ...N_3] \f$. 
 *
 * Inverse transform: Complex -> Real.  <BR>
 * Forward transform: Real -> Complex. <BR>
 * The implemented basis functions: <BR>
 *  FOUR: Fourier along all directions. <BR>
 *  SCFT: Sin/Cos along x, and Fourier along all perpendicular directions. <BR>
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */

//*********************************************************************************************

class RSF 
{
private:
	string field_name;

public:

	//!  \f$F(local_N1, N_2, N_3/2+1) \f$.
	Array<DP,3> Fr;
	

//*********************************************************************************************
public:

	/*! A constructor; 
	 *
	 *  Allocation for Fr, Ek etc. 
	 * Initialize the arrays.
	 */
	RSF(string field_name="");

	//*****************************************************************************************
	
	/*! @brief Inplace Forward transform of CSF; \f$ \mathcal{F}(Fr) -> Fr \f$.
	*
	*  @param	 temp_r			Complex 3D array, a temporary array.
	*
	*  @return FOUR: FOURIER transform
	*  @return SCFT: SFT(Fr) -> Fr.
	*/
//	void Forward_transform();
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}(Fr) -> F	\f$.
	 *
	 *  @param	 Fr				Fr that is to be transformed. 
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(Fr) -> F  (Ftr unaffected).
	 */	
	void Forward_transform(CSF csf);
	
	
	//*****************************************************************************************
	
	/*! @brief Inplace Inverse transform of CSF; \f$ \mathcal{F}^{-1}(Fr) -> Fr \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(Fr) -> Fr.
	 */	
//	void Inverse_transform();
	
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}^{-1}(F) -> Fr	\f$.
	 *
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(F) -> Fr  (F unaffected).
	 */	
	void Inverse_transform(CSF csf);

	//*****************************************************************************************

	// Read/Write functions
	void Write_real_field();
	void Read_real_field();
}; 

#endif

//**************************  End of RSF class declaration ************************************
