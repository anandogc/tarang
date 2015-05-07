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

/*

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

 */ 
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

public:
	//!  \f$ V_x(local_N1, N_2, N_3/2+1) \f$.
	Array<Complex,3> *V1r;	
	
	//!  \f$ V_y(local_N1, N_2, N_3/2+1) \f$.
	Array<Complex,3> *V2r;			
	
	//!  \f$ V_z(local_N1, N_2, N_3/2+1) \f$.
	Array<Complex,3> *V3r;		
	

	
	//! max grid coordinate along the anisotropy direction.
//	int RV_structurefn_pll_ind_max;
	
	
	int Qi1_start, Qi2_start, Qi3_start;
	int Qi1_end, Qi2_end, Qi3_end;
	
	int dist_farthest_proc;  // for i1 idex.
	
	
	//! structure function
	/*
	 *	CV_St(r, q, m) with  r=[0, rmax], q=[q_min, q_max], m=0:1
	 *  rmax = maximum radius that fits inside real field.
	 *  m=0:  for \f$ \Delta u_{||} \f$
	 *  m=1:  for \f$ \Delta u_{\perp} \f$
	 */
	Array<Real,3> *RV_St;
	Array<Real,1> *RV_St_count;
	
	Array<Real,3> *RV_St_final;
	Array<Real,1> *RV_St_count_final;
	
	
	
//*********************************************************************************************			
 
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */
	RVF(); 

	//*****************************************************************************************
	
	/*! @brief Inplace Forward transform of CVF; \f$ \mathcal{F}(\vec{Vr}) -> \vec{Vr} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: FOURIER transform
	 *  @return SCFT: SFT(V1r) -> V1r; CFT(V2r) -> V2r; CFT(V3r) -> V3r;
	 */
	void RV_Forward_transform(Array<Complex,3> temp_r);
	
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
	void RV_Forward_transform_transpose_order
	(
		Array<Complex,3> V1, Array<Complex,3> V2, Array<Complex,3> V3
	);
		
	//*****************************************************************************************
																									
	/*! @brief Inplace Inverse transform of CVF; \f$ \mathcal{F}(\vec{Vr}) -> \vec{Vr} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(V1r) -> V1r; ICFT(V2r) -> V2r; ICFT(V3r) -> V3r;
	 */													
	void RV_Inverse_transform();

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
	void RV_Inverse_transform_transpose_order
	(
		Array<Complex,3> V1, Array<Complex,3> V2, Array<Complex,3> V3
	);
												
	//*****************************************************************************************
												
	
	//*****************************************************************************************	

	// Compute structure function
	
	void Compute_dSt
	(
		 TinyVector<int,3> j1, TinyVector<int,3> j2,
		 TinyVector<Real,3> V1, TinyVector<Real,3> V2
	);
	
	void Compute_structure_function_unit
	(
		 int otherproc, 
		 Array<Complex,3> otherV1,
		 Array<Complex,3> otherV2,
		 Array<Complex,3> otherV3
	);
	
	void Compute_local_structure_function();
	
	void RV_Compute_structure_function();
	
	void RV_structure_fn_Q_ends(int Pi1, int Pi2, int Pi3);
	
	//*****************************************************************************************
	
	void Compute_planar_structure_function_unit_anis1(ofstream& file_out, int pll_index);
	
	void Compute_planar_structure_function_unit_anis2
	(
		 int pll_index,
		 int otherproc, 
		 Array<Complex,2> otherprocV1,
		 Array<Complex,2> otherprocV2,
		 Array<Complex,2> otherprocV3
	);
	 
	
	void Compute_planar_structure_function_unit_anis3
	(
		 int pll_index,
		 int otherproc, 
		 Array<Real,2> otherprocV1,
		 Array<Real,2> otherprocV2,
		 Array<Real,2> otherprocV3,
		 int imag_switch
	);
	 
	void Compute_local_planar_structure_function_anis2(int pll_index);
	 
	void Compute_local_planar_structure_function_anis3
	(
		int pll_index,
		int imag_switch
	);
	 
	void RV_Compute_planar_structure_function_anis1(ofstream& file_out, int pll_index);
	
	void RV_Compute_planar_structure_function_anis2(int pll_index);
	
	void RV_Compute_planar_structure_function_anis3(int pll_index, int imag_switch);
	
	
	//*****************************************************************************************
	/// Output \f$ \vec{Vr} \f$ as real to fileout.
	void RV_Output(ofstream& fileout);
	
	void RV_Output_transpose_order(ofstream& fileout);
	
	void RV_input(ifstream& file_in);
	
	//void RV_Output_hdf5(DataSet* dataset1, DataSet* dataset2, DataSet* dataset3,
	//	       	 Array<Complex,3> temp_array);
	
};

#endif

//**************************  End of RVF.h  ***************************************************


