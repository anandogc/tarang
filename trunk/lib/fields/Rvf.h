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

/*! \file  Rvf.h
 * 
 * @brief  Class declaration of Rvf, a Real Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug  No known bugs
 */
 
 
#ifndef _RVF
#define _RVF

#include "basis.h"

class CVF;

//*********************************************************************************************

//! Real vector field
/*!
 * 3 dimensional real vector field has 3 components with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, N_3/2+1] \f$.   <BR>
 * These arrays hold two real values for each complex entry.
 *  The real dimension are \f$[N_1, N_2, N_3] \f$. 
 *
 * Inverse transform: Complex -> Real.  <BR>
 * Forward transform: Real -> Complex. <BR>
 * The implemented basis functions: <BR>
 *  FOUR: Fourier along all directions. <BR>
 *  SCFT: Sin/Cos along x, and Fourier along all perpendicular directions. <BR>
 * 
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */

//*********************************************************************************************

class RVF 
{
private: 
	string field_name;

public:
	//!  \f$ V_x(local_N1, N_2, N_3/2+1) \f$.
	Array<Real,3> V1r;
	
	//!  \f$ V_y(local_N1, N_2, N_3/2+1) \f$.
	Array<Real,3> V2r;
	
	//!  \f$ V_z(local_N1, N_2, N_3/2+1) \f$.
	Array<Real,3> V3r;
	
	
//*********************************************************************************************			
 
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */ 
	RVF(string field_name=""); 

	//*****************************************************************************************
	
	/*! @brief Inplace Forward transform of CVF; \f$ \mathcal{F}(\vec{Vr}) -> \vec{Vr} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: FOURIER transform
	 *  @return SCFT: SFT(V1r) -> V1r; CFT(V2r) -> V2r; CFT(V3r) -> V3r;
	 */
//	void Forward_transform();
	
	//
	//
	
	/*! @brief  Forward transform in transpose order of CVF; 
	 *				\f$ \mathcal{F}(\vec{Vr}) -> \vec{V}	\f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(V1r) -> V1; CFT(V2r) -> V2; CFT(V3r) -> V3;  (Vr unaffected).
	 */	
	void Forward_transform(CVF cvf);
		
	//*****************************************************************************************
																									
	/*! @brief Inplace Inverse transform of CVF; \f$ \mathcal{F}(\vec{Vr}) -> \vec{Vr} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(V1r) -> V1r; ICFT(V2r) -> V2r; ICFT(V3r) -> V3r;
	 */													
//	void Inverse_transform();

	//
	//
													
	/*! @brief  Inverse transform in transpose order of CVF; 
	 *				\f$ \mathcal{F}(\vec{V}) -> \vec{Vr}	\f$.
	 *
	 *  @param	 V1				x component of V that is to be transformed. 
	 *  @param	 V2				y component of V that is to be transformed. 
	 *  @param	 V3r			z component of V that is to be transformed. 
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(V1) -> V1r.
	 *  @return SCFT: ICFT(V2) -> V2r.
	 *  @return SCFT: ICFT(V3) -> V3r.
	 */																
	void Inverse_transform (CVF cvf);
												
	
	//*****************************************************************************************

	void Read_real_field();
	void Write_real_field();
	void Write_real_field_slice(unsigned int slice_file_counter);
	
};

#endif

//**************************  End of RVF.h  ***************************************************


