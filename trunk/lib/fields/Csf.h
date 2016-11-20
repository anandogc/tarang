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

/*! \file  Csf.h
 * 
 * @brief  Class declaration of Csf, a Complex Scalar Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bugs Out of date functions output_split_arrays, CS_iutput_split_arrays.
 */

#ifndef _CSF
#define _CSF

#include "plainCsf.h" 

class RSF;

//*********************************************************************************************

//! Complex Scalar field
/*!
 * 3 dimensional complex scalar field with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, N_3/2+1] \f$.   <BR>
 * This array can also store real values, and its dimension is \f$[N_1, N_2, N_3] \f$. 
 *
 * Inverse transform: Complex -> Real.  <BR>
 * Forward transform: Real -> Complex. <BR>
 * The implemented basis functions: <BR>
 *  FOUR: Fourier along all directions. <BR>
 *  SCFT: Sin/Cos along x, and Fourier along all perpendicular directions. <BR>
 * 
 *  The modal energy \f$ C(\vec{k}) \f$ is defined in each basis function.  
 *  The dissipation rate for the mode is \f$ 2 K^2 C(\vec{k}) \f$.
 * 
 *  Entropy = = \f$ \sum p(k) log(1/p(k)) \f$ where probability  \f$ p(k) =  E(k)/E_{total)\f$ 
 *					with \f$ E(k) = |A(k)|^2/2 \f$.
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */



class CSF : public PlainCSF
{ 
private:
	string field_name;

public:

	//!  Total energy of Csf.													
	Real		total_energy;						
	
	//!  Total dissipation rate of Csf (without \f$ \kappa \f$).
	Real		total_k2energy;					
	
	//!  Entropy of Csf
	Real		entropy;																				
	
												
//*********************************************************************************************					
	
public:

	/*! A constructor; 
	 *
	 *  Allocation for F, Ek etc. 
	 * Initialize the arrays.
	 */
	CSF(string field_name="");


	void Copy_to(CSF& to);
	void Copy_to(PlainCSF& to);
	void Copy_from(CSF& from);
	void Copy_from(PlainCSF& from);
	
	
	/*! @brief Inplace Forward transform of CSF; \f$ \mathcal{F}(F) -> F \f$.
	*
	*  @param	 temp_r			Complex 3D array, a temporary array.
	*
	*  @return FOUR: FOURIER transform
	*  @return SCFT: SFT(F) -> F.
	*/
//	void Forward_transform();
	
	
	/*! @brief Inplace Inverse transform of CSF; \f$ \mathcal{F}^{-1}(F) -> F \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(F) -> F.
	 */	
//	void Inverse_transform();
	
	
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}(F) -> F	\f$.
	 *
	 *  @param	 Ftr			Ftr that is to be transformed. 
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(Ftr) -> F  (Ftr unaffected).
	 */	
	void Forward_transform(RSF rsf);
	
	
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}^{-1}(F) -> Ftr	\f$.
	 *
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(F) -> Ftr  (F unaffected).
	 */		
	void Inverse_transform(RSF rsf);


	//***************************************************************************************** 
	
	/// F(k) = F(k)/K^2.  F(0)=0.
	void Divide_ksqr();										
  	
	//***************************************************************************************** 
	
	/// Compute total energy and dissipation of F.
	void Compute_total_energy();
	void Compute_total_k2energy();
	void Compute_total_kn_energy(int n, Real &result);
	
	/// Compute entropy of F.
	void Compute_entropy();
	
	
	//*****************************************************************************************	
	
	/// 3D: Return modal energy for grid index (i1,i2,i3).
	Real Modal_energy(int i1, int i2, int i3);
	
	//*****************************************************************************************	
	
	/// Dealiase F.
	void Dealias();
	
	//*****************************************************************************************	

	void Read_complex_field();
	void Read_reduced_complex_field();
	void Write_complex_field();
	void Write_reduced_complex_field();

	
}; 

#endif

//******************************** End of CSF.h ***********************************************



