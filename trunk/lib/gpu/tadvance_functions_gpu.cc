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


/*! \file  compute_rhs.cc
 * 
 * @brief Compute rhs.
 * @sa void Time_advance_incompress::Compute_rhs() 
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "gpu.h"

//*********************************************************************************************
// Add force to nlin

void Add_force_fluid_gpu()
{
	
	if (U_force_switch_gpu) {
		U_nlin1_gpu -= U_Force1_gpu;
		U_nlin2_gpu -= U_Force2_gpu;
		U_nlin3_gpu -= U_Forc31_gpu;
	}
}


//******************************************************************************

// nlin = nlin + grad(p)
void Add_pressure_gradient_gpu()
{
	Add_Xderiv(P_gpu, U_nlin1_gpu);
	
	Add_Yderiv(P_gpu, U_nlin2_gpu);
	
	Add_Zderiv(P_gpu, U_nlin3_gpu);
}


//*********************************************************************************************

/** @brief Compute rhs for fluid simulation
 * 
 * @note Status before entering the fn
 *		nlin[i] contains Dj T[Vr[j]*Vr[i]] - f[i];
 *		F.p contains pressure;
 * 
 * @return (rhs) = nlin[i] =  -Dj T[Vr[j]*Vr[i]] + f[i] - Di [p(k)] (-grad(pressure)) 	
 */
void Compute_rhs_fluid_gpu()
{	

	Add_pressure_gradient_gpu();				// nlin = nlin + grad(p)

	U_nlin1_gpu = -U_nlin1_gpu;
	U_nlin2_gpu = -U_nlin2_gpu;
	U_nlin3_gpu = -U_nlin3_gpu;
}

//*********************************************************************************************
//*********************************************************************************************
/** @brief Single time step
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) *\exp(-K^2*\nu*dt*a) 
 *								+ c dt \vec{N}(t+dt) *\exp(-K^2*\nu*dt*b) ] 
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$
 */
void Single_time_step(Real a, Real b, Real c)
{	
	
	if ((abs(a) < MYEPS) && (abs(b) < MYEPS))
		Add_U_nlin_factor_dt(c);
	
	else if ((abs(a) < MYEPS) && (abs(b) > MYEPS)) {
		Mult_U_nlin_exp_ksqr_dt(b);
		Add_U_nlin_factor_dt(c);
	}
	
	else if ((abs(b) < MYEPS) && (abs(a) > MYEPS)) {
		Mult_U_exp_ksqr_dt(a);
		Add_U_nlin_factor_dt(c);	
	}
	
	else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS)) {
		Add_U_nlin_factor_dt(c);									// V = V + c*dt*nlin
		Mult_U_exp_ksqr_dt(a);
	}
	
	else if ((abs(a) > MYEPS) && (abs(b) > MYEPS)) {
		Mult_U_exp_ksqr_dt((a-b));
		Add_U_nlin_factor_dt(c);	
		Mult_U_exp_ksqr_dt(b);
	}
}
//*********************************************************************************************
/** @brief Compute_force_TO_rhs()
 *
 * @return Computes force and nlin; add them; compute pressure; combine them to compute rhs.
 *
 */
void Compute_force_TO_rhs_fluid_gpu()
{
	
	Compute_force_fluid_gpu();
	
	Compute_nlin_fluid_gpu();
	
	Add_force_fluid_gpu();							// nlin = nlin - f
	
	Compute_pressure_gpu();							// Compute pressure using V(t+dt/2)
	
	Compute_rhs_fluid_gpu();
}


/**********************************************************************************************
 
 Function 1:  V = V + nlin*dt;
 Function 2: for CVF: Y = Y+ factor*dt*U.nlin;
 
 ***********************************************************************************************/

void Add_U_nlin_factor_dt(Real factor)
{
	if (abs(factor) > MYEPS2) {
		U_V1_gpu += (Tdt_gpu)*(U_nlin1_gpu);
		U_V2_gpu += (Tdt_gpu)*(U_nlin2_gpu);
		U_V3_gpu += (Tdt_gpu)*(U_nlin3_gpu);
	}
	
}

//*********************************************************************************************

/** @brief Multiply vector field by \f$ \exp(-\nu K^2 dt) \f$ or
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 *
 *  @param dt
 *
 * @return \f$ V(k) = V(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1,
 *			\f$ V(k) = V(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */
void Mult_U_exp_ksqr_dt(Real a)
{
	if (abs(a) > MYEPS) {
		if (!hyper_dissipation_switch_gpu) {
			Array_exp_ksqr(U_V1_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu);
			Array_exp_ksqr(U_V2_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu);
			Array_exp_ksqr(U_V3_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu);
		}
		
		else {
			Array_exp_ksqr(U_V1_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu -hyper_dissipation_coefficient_gpu*a*Tdt_gpu, hyper_dissipation_exponent_gpu);
			
			Array_exp_ksqr(U_V2_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu -hyper_dissipation_coefficient_gpu*a*Tdt_gpu, hyper_dissipation_exponent_gpu);
			
			Array_exp_ksqr(U_V3_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu -hyper_dissipation_coefficient_gpu*a*Tdt_gpu, hyper_dissipation_exponent_gpu);
		}
	}
	
}



//*********************************************************************************************
/** @brief Multiply nonlinear field by \f$ \exp(-\nu K^2 dt) \f$ or
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 *
 *  @param dt
 *
 * @return \f$ N(k) = N(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1,
 *			\f$ N(k) = N(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */

void Mult_U_nlin_exp_ksqr_dt(Real a)
{
	if (abs(a) > MYEPS) {
		if (!hyper_dissipation_switch_gpu) {
			Array_exp_ksqr(U_nlin1_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu);
			Array_exp_ksqr(U_nlin2_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu);
			Array_exp_ksqr(U_nlin3_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu);
		}
		
		else {
			Array_exp_ksqr(U_nlin1_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu, -hyper_dissipation_coefficient_gpu*a*Tdt_gpu, hyper_dissipation_exponent_gpu);
			 
			Array_exp_ksqr(U_nlin2_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu, -hyper_dissipation_coefficient_gpu*a*Tdt_gpu, hyper_dissipation_exponent_gpu);
			
			Array_exp_ksqr(U_nlin3_gpu, -dissipation_coefficient_gpu*a*Tdt_gpu, -hyper_dissipation_coefficient_gpu*a*Tdt_gpu, hyper_dissipation_exponent_gpu);
		}
	}
}


void Compute_pressure_gpu()
{

	Compute_divergence_U_nlin(P_gpu);  // div(nlin) -> F  (the proc uses temp arrays)
	
	Array_divide_ksqr(P_gpu);			// F = F/k^2
    
}



void Compute_divergence_U_nlin(Array<Complex,3> div)
{
    Real total_abs_div;  // not reqd for this function.
    Compute_divergence(U_nlin1_gpu, U_nlin2_gpu, U_nlin3_gpu, div, "nlin", total_abs_div, false);
    // PS: global.temp_array.X is used in Compute_divergence_field.. so
    // be careful.
    // Used in pressure computation.
}

//*********************************************************************************************


/** @brief Compute divergence of VF \f$ \vec{V} \f$.
 *
 *	@note	Divergence is stored in CVF of class NLIN.  We use the same class to store
 *			pressure also.  Hence divergence function should be called in the beginning
 *			or end of loop when *F is free.
 *
 *  @return  \f$ *F = \mathcal{F}(D_i V_i) \f$.
 */
void Compute_divergence_U(Array<Complex,3> div, Real &total_abs_div, bool print_switch)
{
    Compute_divergence(U_V1_gpu, U_V2_gpu, U_V3_gpu, div, "field", total_abs_div, print_switch);
    // PS: global.temp_array.X is used in Compute_divergence_field.. so
    // be careful.  We call div using global.temp_array.X2
}




//**********************************   End of compute_rhs.cc  *********************************

