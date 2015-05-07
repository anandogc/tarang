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

#ifndef _H_FluidSF
#define _H_FluidSF

#include "fields.h" 

class FluidSF
{
 
public:
	CSF csf;
	RSF rsf;

    Real diffusion_coefficient;
    Real hyper_diffusion_coefficient;
    int hyper_diffusion_exponent;

	bool hyper_diffusion_switch;	
	//!   \f$ nlin(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> nlin;						
	
	//!  Force \f$ F_x(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> Force;	
	bool force_switch;
	
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */
	FluidSF
	(
		Real diffusion_coefficient, 
		Real hyper_diffusion_coefficient, 
		int hyper_diffusion_exponent,
		bool force_switch,
		string field_name
	);
 
	void Inverse_transform();
	void Forward_transform();
	
	void Copy_field_to(CSF& T);
	void Copy_field_to(PlainCSF& T);
	void Copy_field_from(CSF& T);
	void Copy_field_from(PlainCSF& T);
	
	void Satisfy_strong_reality_condition_field();
	void Satisfy_strong_reality_condition_force_field();
    
    void Satisfy_weak_reality_condition_field();
	void Satisfy_weak_reality_condition_force_field();
    
    void Test_reality_condition_field();
	void Test_reality_condition_force_field();
    
    void Dealias_force_field();
	
	void Mult_field_exp_ksqr_dt(Real a);
	void Mult_nlin_exp_ksqr_dt(Real a);
	
	void Add_nlin_factor_dt(Real factor);
	void Add_field_nlin_factor_dt(PlainCSF& W, Real factor);
	void Add_complex_conj(int kx, int ky, int kz, Complex localG);
	
	Real Get_Tk(int kx, int ky, int kz);
    
    Real Get_dt();
    

//Added during compilation
//================================================================================
	void Assign_field_and_comp_conj(int kx, int ky, int kz, Complex localG);
	void Assign_random_complex_scalar(int kx, int ky, int kz, Real rand_range);
	void Assign_random_real_scalar(int kx, int ky, int kz, Real rand_range);
	void Assign_field(int kx, int ky, int kz, Real localG);
//================================================================================
};
	
#endif

//******************************** End of CVF.h ***********************************************	






