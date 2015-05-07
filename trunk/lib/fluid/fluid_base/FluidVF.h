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

#ifndef _H_FluidVF
#define _H_FluidVF

#include "PlainFluidVF.h" 
#include "FluidSF.h"

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


class FluidVF:  public PlainFluidVF
{
 
public:
	//!  Force along x \f$ F_x(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> nlin1;						
	
	//!  Force along y \f$ F_y(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> nlin2;						
	
	//!  Force along z \f$ F_x(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> nlin3;	
	
	
	//!  Force along x \f$ F_x(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> Force1;						
	
	//!  Force along y \f$ F_y(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> Force2;						
	
	//!  Force along z \f$ F_x(local_{N1}, N_2, N_3/2+1) \f$.
	Array<Complex,3> Force3;
	
	
	//! Dissipation coefficient appearing before laplacian \f$ \nu \f$.
	Real		dissipation_coefficient;
	bool hyper_diffusion_switch;	
	//! Hyper_dissipation coefficient appearing before \f$\nabla^exponent \f$: 
	bool	hyper_dissipation_switch;
	Real		hyper_dissipation_coefficient;
	int		hyper_dissipation_exponent;
	
	bool force_switch;
                                                                                        
	
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */
	FluidVF
	(
		Real dissipation_coefficient, 
		Real hyper_dissipation_coefficient, 
		int hyper_dissipation_exponent,
		bool force_switch,
		string field_name
	);
	

	void Compute_divergence_field(Array<Complex,3> div, Real &total_abs_div, bool print_switch);
    void Compute_divergence_nlin(Array<Complex,3> div);
	
	Real Get_mag_V0();
	
	void Satisfy_strong_reality_condition_field();
	void Satisfy_strong_reality_condition_force_field();
    
    void Satisfy_weak_reality_condition_field();
	void Satisfy_weak_reality_condition_force_field();
    
    void Test_reality_condition_field();
	void Test_reality_condition_force_field();
    
    void Dealias_force_field();
	
	void Zero_Prandtl_number_compute_temperature(FluidSF& T);
	void Infinite_Prandtl_number_compute_velocity(FluidSF& T);
	
	void Copy_field_to(CVF& W);
	void Copy_field_to(PlainCVF& W);
	void Copy_field_from(CVF& W);
	void Copy_field_from(PlainCVF& W);
	
	void Add_nlin_factor_dt(Real factor);
	void Add_field_nlin_factor_dt(PlainCVF& Y, Real factor);
	
	void Mult_field_exp_ksqr_dt(Real a);
	void Mult_nlin_exp_ksqr_dt(Real a);
	
	void Add_complex_conj(int kx, int ky, int kz, Complex Vx, Complex Vy, Complex Vz);
	void Assign_field_and_comp_conj(int kx, int ky, int kz, Complex Vx, Complex Vy, Complex Vz);
	void Assign_random_complex_vector(int kx, int ky, int kz, Real rand_range);
	void Assign_random_real_vector(int kx, int ky, int kz, Real rand_range);
	
	void Assign_field(int kx, int ky, int kz, Real Vx, Real Vy, Real Vz);
	
	Real Get_Tk(int kx, int ky, int kz);

	void Zero_modes_RB_slip(FluidSF& T);
	
	void Get_local_max_real_space(Real local_ux[]);
	
	Real Get_dt();
	Real Get_dt(FluidSF& T);
	Real Get_dt(FluidVF& W);
	Real Get_dt(FluidVF& W, FluidSF& T);
};
	
#endif

//******************************** End of CVF.h ***********************************************	
