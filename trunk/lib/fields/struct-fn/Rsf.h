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

#include "universal/universal_inline.h" 
#include "universal/universal_basic.h" 
#include "universal/universal_tr.h" 
#include "universal/universal_energy.h" 
#include "universal/universal_ET.h" 
#include "universal/universal_misc.h" 

#include "array_fn/array_basic_inline.h"
#include "array_fn/array_basic.h"
// #include "struct_fn.h"
// #include "planar_struct_fn.h"


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
public:

	//!  \f$F(local_N1, N_2, N_3/2+1) \f$.
	Array<Complex,3> *Fr;
	

	//! max grid coordinate along the anisotropy direction.
//	Real RS_structure_fn_rpll_max;
	
	int Qi1_start, Qi2_start, Qi3_start;
	int Qi1_end, Qi2_end, Qi3_end;
	
//	int dist_farthest_proc;  // for i1 idex.
	
	//! structure function
	/*
	 *	St(r, q, m) with  r=[0, rmax], q=[q_min, q_max], m=0:1
	 *  rmax = maximum radius that fits inside real field.
	 *  m=0:  for \f$ \Delta u_{||} \f$
	 *  m=1:  for \f$ \Delta u_{\perp} \f$
	 */
	Array<Real,2> *RS_St;
	Array<Real,1> *RS_St_count;
	
	Array<Real,2> *RS_St_final;
	Array<Real,1> *RS_St_count_final;
	
	Array<Real,3> *RS_st_planar;
	

//*********************************************************************************************
public:

	/*! A constructor; 
	 *
	 *  Allocation for Fr, Ek etc. 
	 * Initialize the arrays.
	 */
	RSF();

	//*****************************************************************************************
	
	/*! @brief Inplace Forward transform of CSF; \f$ \mathcal{F}(Fr) -> Fr \f$.
	*
	*  @param	 temp_r			Complex 3D array, a temporary array.
	*
	*  @return FOUR: FOURIER transform
	*  @return SCFT: SFT(Fr) -> Fr.
	*/
	void RS_Forward_transform();
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}(Fr) -> F	\f$.
	 *
	 *  @param	 Fr				Fr that is to be transformed. 
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(Fr) -> F  (Ftr unaffected).
	 */	
	void RS_Forward_transform_transpose_order(Array<Complex,3> F);
	
	
	//*****************************************************************************************
	
	/*! @brief Inplace Inverse transform of CSF; \f$ \mathcal{F}^{-1}(Fr) -> Fr \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(Fr) -> Fr.
	 */	
	void RS_Inverse_transform();
	
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}^{-1}(F) -> Fr	\f$.
	 *
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(F) -> Fr  (F unaffected).
	 */	
	void RS_Inverse_transform_transpose_order(Array<Complex,3> F);

	//*****************************************************************************************
	

	//*****************************************************************************************	
	
	void Compute_dSt
	(
		TinyVector<int,3> j1, TinyVector<int,3> j2,
		Real F1, Real F2
	);
	
	void Compute_structure_function_unit
	(
		int otherproc,
		Array<Complex,3> otherprocF
	);
	
	void Compute_local_structure_function();
	
	void RS_Compute_structure_function();
	
	void RS_structure_fn_Q_ends(int Pi1, int Pi2, int Pi3);
	
	//*****************************************************************************************
	
	void Compute_planar_structure_function_unit_anis1(ofstream& file_out, int pll_index);
	
	void Compute_planar_structure_function_unit_anis2
	(
		int pll_index,
		int otherproc, 
		Array<Complex,2> otherprocF
	);
	
	
	void Compute_planar_structure_function_unit_anis3
	(
		int pll_index,
		int otherproc, 
		Array<Real,2> otherprocF,
		int imag_switch
	);
	
	void Compute_local_planar_structure_function_anis2(int pll_index);
	
	void Compute_local_planar_structure_function_anis3
	(
		int pll_index,
		int imag_switch
	);
	
	void RS_Compute_planar_structure_function_anis1(ofstream& file_out, int pll_index);
	
	void RS_Compute_planar_structure_function_anis2(int pll_index);
	
	void RS_Compute_planar_structure_function_anis3(int pll_index, int imag_switch);
	
	//*****************************************************************************************
	
	/// Output F to fileout as real.
	void RS_Output(ofstream& fileout);
	
	void RS_Output_transpose_order(ofstream& fileout);
	
	void RS_input(ifstream& file_in);
	
}; 

#endif

//**************************  End of RSF class declaration ************************************



