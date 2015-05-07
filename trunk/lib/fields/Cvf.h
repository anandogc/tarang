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

/*! \file  Cvf.h
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008 
 * @bugs  No known bugs
 */

#ifndef _CVF
#define _CVF

#include "plainCvf.h"

class RVF;


//*********************************************************************************************

//! Complex vector field
/*!
 * 3 dimensional complex vector field has 3 components with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2,N_3/2+1] \f$.   <BR>
 * These arrays can also store real values, and their dimension are \f$[N_1, N_2, N_3] \f$. 
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
 *  Helicity1 = \f$ H1(K) = \vec{K} . [\vec{Vr} \times \vec{Vi}] \f$. <BR>
 *  Helicity2 = \f$ H2(K) = \vec{K} . [(\vec{Vr} \times \vec{Vi})] /K^2 \f$. <BR>
 * 
 *  Entropy = \f$ \sum p(k) log(1/p(k)) \f$ where probability  
 *				\f$ p(k) =  E(k)/E_{total) \f$ with \f$ E(k) = |A(k)|^2/2 \f$.
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */

//*********************************************************************************************

class CVF: public PlainCVF 
{

private:
	string field_name;
 
public:		

	
	//*****************************************************************************************
	
	// Nondynamic allocation now
	
	//!  Total energy of Cvf.
	Real		total_E1, total_E2, total_E3;
	Real		total_energy;						
	
	//!  Total dissipation rate of Cvf (without \f$ \nu \f$).
	Real		total_k2energy;	
	Real		total_k2Hc;	
	
	//!  \f$ \sum \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K})) \f$				
	Real		total_helicity1;						
	
	//!  \f$ \sum \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K}))/K^2 \f$
	Real		total_helicity2;						
	
	//!  Dissipation of H1. 
	Real		total_k2H1;				
	
	//!  Dissipation of H2. 
	Real		total_k2H2;				
	
	//!  Entropy of Cvf
	Real		entropy;	
	
																							
	//*****************************************************************************************
	
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */
	CVF(string field_name="");
 
	/// Copy curr Vi's to VF to.
	void Copy_to(CVF& to); 
	void Copy_to(PlainCVF& to); 
	void Copy_from(CVF& from);
	void Copy_from(PlainCVF& from);

	//*****************************************************************************************	
	
	/*! @brief Inplace Forward transform of CVF; \f$ \mathcal{F}(\vec{V}) -> \vec{V} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: FOURIER transform
	 *  @return SCFT: SFT(V1) -> V1; CFT(V2) -> V2; CFT(V3) -> V3;
	 */
//	void Forward_transform();
	

	/*! @brief  Forward transform in transpose order of CVF; 
	 *			\f$ \mathcal{F}(\vec{Vtr}) -> \vec{V}	\f$.
	 *
	 *  @param	 V1tr			x component of Vtr that is to be transformed. 
	 *  @param	 V2tr			y component of Vtr that is to be transformed. 
	 *  @param	 V3tr			z component of Vtr that is to be transformed. 
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(V1tr) -> V1; CFT(V2tr) -> V2; CFT(V3tr) -> V3;  (Vtr unaffected).
	 */	
	void Forward_transform(RVF rvf);
													
	//*****************************************************************************************													
													
	/*! @brief Inplace Inverse transform of CVF; \f$ \mathcal{F}(\vec{V}) -> \vec{V} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(V1) -> V1; ICFT(V2) -> V2; ICFT(V3) -> V3;
	 */													
//	void Inverse_transform();
	
	
	
	/*! @brief  Inverse transform in transpose order of CVF; 
	 *			\f$ \mathcal{F}(\vec{V}) -> \vec{Vtr}	\f$.
	 *
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(V1) -> V1tr.
	 *  @return SCFT: ICFT(V2) -> V2tr.
	 *  @return SCFT: ICFT(V3) -> V3tr.
	 */																								
	void Inverse_transform(RVF rvf);
													
	//*****************************************************************************************   
	
	void Divide_ksqr();
	
	/// Compute total energy and dissipation of \f$ \vec{V} \f$. 
	void Compute_total_energy();
	void Compute_total_k2energy();
	void Compute_total_kn_energy(int n, Real &result);
	
	void Compute_total_k2Hc(Array<Complex,3> W1, Array<Complex,3> W2, Array<Complex,3> W3);
	
	/// Compute entropy of \f$ \vec{V} \f$.
	void Compute_entropy();
	
	//*****************************************************************************************
	
	
	/// Return modal energy for grid index (i1,i2,i3).
	Real Modal_energy(int i1, int i2, int i3);
	
	/// vorticity = modal vorticity for grid index (i1,i2,i3).
	void Compute_Modal_vorticity(int i1, int i2, int i3, TinyVector<Complex,3> &vorticity);
	
	void Compute_Modal_vorticity_y_component(int l1, int l2, int l3, Complex &vort_y);
	
	//*****************************************************************************************
	
	/// Return modal helicity for the grid index (i1,i2,i3).
	Real Modal_helicity(int i1, int i2, int i3);
	
	/// Compute total helicity1, helicity2, and their dissipation for Cvf.
	void Compute_total_helicity();
	
	
	//*****************************************************************************************	
	
	 /// Dealiase \f$ \vec{V} \f$ using 2/3 rule.
	 void Dealias();

	//*****************************************************************************************																																			
	void Read_complex_field();
	void Read_reduced_complex_field();
 
	/// Output \f$ \vec{V} \f$  to file_out.	
	void Write_complex_field();
	
	void Write_reduced_complex_field();
	
	
	//*****************************************************************************************


	
};

#endif

//******************************** End of CVF.h ***********************************************	






