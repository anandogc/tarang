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


/*! \file  IncFluid.h
 * 
 *	@brief  Class declaration of IncFluid, Incompressible fluid (V field).
 * 
 *	For fluid, scalar, MHD, this class is for velocity. <BR>
 *	For scalar simulation, we need IncFluid V & IncSF T. <BR>
 *	For MHD simulation, we need IncFluid V & IncVF W.
 *
 *	@author  M. K. Verma
 *	@version 4.0 MPI
 *	@date Sept 2008
 *
 * @bug  No known bugs
 */
 
//*********************************************************************************************

#ifndef _H_Time_advance
#define _H_Time_advance

#include "fluid_base.h"   
#include "Nlin.h"
#include "Pressure.h"
#include "FORCE.h"

//! @brief Incompressible fluid
/*!
 *  Inherits IncVF, which itself inherits CVF, RVF, NLIN, and EnergyTr. <BR>
 *  Inherits Time that contains time vars. <BR>
 * 
 *	This is fluid 
 *
 *	@sa IncSF.h
 *	@sa Nlin.h
 *	@sa EnergyTr.h
 */
 
//*********************************************************************************************	

	// class nlin;

class Time_advance_incompress
{

public:	
	
	void Add_force(FluidVF& U);
	void Add_force(FluidVF& U, FluidSF& T);
	void Add_force_scalar(FluidVF& U, FluidSF& T);
	void Add_force_RBC(FluidVF& U, FluidSF& T);
	void Add_force(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Add_force(FluidVF& U, FluidVF& W);
	void Add_force(FluidVF& U, FluidVF& W, FluidSF& T);
	void Add_force(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	
	void Add_pressure_gradient(FluidVF& U, Pressure& P);
			
	void Compute_rhs(FluidVF& U, Pressure& P);       // nlin[i]=-ik[i]*pressure(k)+nlin[i] 
	void Compute_rhs(FluidVF& U, FluidSF& T, Pressure& P); 
	void Compute_rhs_scalar(FluidVF& U, FluidSF& T, Pressure& P);
	void Compute_rhs_RBC(FluidVF& U, FluidSF& T, Pressure& P); 
	void Compute_rhs(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P);
	void Compute_rhs(FluidVF& U, FluidVF& W, Pressure& P);  
	void Compute_rhs(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P); 
    void Compute_rhs(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P);
	
    void Single_time_step(FluidVF& U, Real a, Real b, Real c);    
	void Single_time_step(FluidSF& T, Real a, Real b, Real c);
	void Single_time_step(FluidVF& U, FluidSF& T, Real a, Real b, Real c);
	void Single_time_step_scalar(FluidVF& U, FluidSF& T, Real a, Real b, Real c);
	void Single_time_step_RBC(FluidVF& U, FluidSF& T, Real a, Real b, Real c); 
	void Single_time_step(FluidVF& U, FluidSF& T1, FluidSF& T2, Real a, Real b, Real c);
	void Single_time_step(FluidVF& U, FluidVF& W, Real a, Real b, Real c);
	void Single_time_step(FluidVF& U, FluidVF& W, FluidSF& T, Real a, Real b, Real c);
	void Single_time_step(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Real a, Real b, Real c);
	
	void Compute_force_TO_rhs(FluidVF& U, Pressure& P, FORCE &Force);
	void Compute_force_TO_rhs(FluidVF& U,  FluidSF& T, Pressure& P, FORCE &Force);
	void Compute_force_TO_rhs_scalar(FluidVF& U, FluidSF& T, Pressure& P, FORCE &Force);
	void Compute_force_TO_rhs_RBC(FluidVF& U, FluidSF& T, Pressure& P, FORCE &Force);
	void Compute_force_TO_rhs(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force);
	void Compute_force_TO_rhs(FluidVF& U, FluidVF& W, Pressure& P, FORCE &Force);
	void Compute_force_TO_rhs(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P,FORCE &Force);
    void Compute_force_TO_rhs(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P, FORCE& Force);
	
    void Compute_force_TO_rhs(FluidSF& T, FORCE& Force);
	
	
	void Compute_RK4_parts(FluidVF& U, PlainCVF& tot_Vrhs, Real b, Real factor);
	void Compute_RK4_parts(FluidSF& T, PlainCSF& tot_Srhs, Real b, Real factor);
	void Compute_RK4_parts(FluidVF& U, FluidSF& T, PlainCVF& tot_Vrhs, PlainCSF& tot_Srhs, Real b, Real factor);
	void Compute_RK4_parts_scalar(FluidVF& U, FluidSF& T, PlainCVF& tot_Vrhs, PlainCSF& tot_Srhs, Real b, Real factor);
	void Compute_RK4_parts_RBC(FluidVF& U, FluidSF& T, PlainCVF& tot_Vrhs, PlainCSF& tot_Srhs, Real b, Real factor);
	void Compute_RK4_parts(FluidVF& U, FluidSF& T1, FluidSF& T2, PlainCVF& tot_Vrhs, PlainCSF& tot_S1rhs, PlainCSF& tot_S2rhs, Real b, Real factor);
	void Compute_RK4_parts(FluidVF& U, FluidVF& W, PlainCVF& tot_Vrhs, PlainCVF& tot_Wrhs, Real b, Real factor);
	void Compute_RK4_parts(FluidVF& U, FluidVF& W, FluidSF& T, PlainCVF& tot_Vrhs, PlainCVF& tot_Wrhs, PlainCSF& tot_Srhs, Real b, Real factor);
	void Compute_RK4_parts(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, PlainCVF& tot_Vrhs,PlainCVF& tot_Wrhs,PlainCSF& tot_Srhs, PlainCSF& tot_Crhs, Real b, Real factor);
	
	void Make_field_incompressible(FluidVF& U);
	
	// time advance
	void Time_advance_step(FluidVF& U, Pressure& P, FORCE& Force);		 
	void Time_advance_step(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force);
	void Time_advance_step(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force);
	void Time_advance_step(FluidVF& U, FluidVF& W, Pressure& P, FORCE& Force);	
	void Time_advance_step(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P, FORCE& Force);
	void Time_advance_step(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C,Pressure& P, FORCE& Force);
	void Time_advance_step(FluidSF& T, FORCE& Force);
    
	void Euler(FluidVF& U, Pressure& P, FORCE& Force);
	void RK2(FluidVF& U, Pressure& P, FORCE& Force);
	void RK4(FluidVF& U, Pressure& P, FORCE& Force);
	
	void Euler(FluidVF& U,FluidSF& T, Pressure& P, FORCE& Force);
	void RK2(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force);
	void RK4(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force);
	
	void Euler(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force);
	void RK2(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force);
	void RK4(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force);
	
	void Euler(FluidVF& U, FluidVF& W, Pressure& P, FORCE& Force);
	void RK2(FluidVF& U, FluidVF& W, Pressure& P, FORCE& Force);
	void RK4(FluidVF& U, FluidVF& W, Pressure& P, FORCE& Force);
	
	void Euler(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P, FORCE& Force);
	void RK2(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P, FORCE& Force);
	void RK4(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P, FORCE& Force);
	
	void Euler(FluidVF& U,  FluidVF& W, FluidSF&  T,  FluidSF& C, Pressure& P, FORCE& Force);
	void RK2(FluidVF& U,  FluidVF& W, FluidSF&  T,  FluidSF& C, Pressure& P, FORCE& Force);
	void RK4(FluidVF& U,  FluidVF& W, FluidSF& T,  FluidSF& C, Pressure& P, FORCE& Force);
    
	// For chebyshev
    void BDF1(FluidVF& U, Pressure& P, FORCE& Force);
	void BDF1(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force);
	
	void Adam_Bashforth(FluidVF& U, Pressure& P, FORCE& Force);
	void Adam_Bashforth(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force);
    
	// For GP
    void Euler(FluidSF& T, FORCE& Force);
	void RK2(FluidSF& T, FORCE& Force);
	void RK4(FluidSF& T, FORCE& Force);
    
    void Compute_derivative_at_boundary(Array<Complex,1> v, Complex &v_prime_plus, Complex &v_prime_minus);
    void Compute_derivative_at_boundary(Array<Real,1> v, Real &v_prime_plus, Real &v_prime_minus);
    
//    void Add_ubydt(FluidVF& U);
    void Compute_homgeneous_soln_influence_matrix(FluidVF& U, Pressure& P);
	void Compute_pressure_Ux_chebyshev(FluidVF& U, Pressure& P);
    
    void Tridiagonal_solver(Array<Real,1> f, Array<Real,1> y, Real lambda, Real d0, int odd_switch);
    void Helmholtz_real(Array<Real,1> f, Array<Real,1> u, Real lambda, Real value_plus1, Real value_minus1);
    void Helmholtz_complex(Array<Complex, 1> f, Array<Complex, 1> u, Real lambda, Complex value_plus1, Complex value_minus1);
    
    void Helmholtz_real_full_array(Array<Real,3> R, Array<Real,3> u, bool pressure_switch, Real lambda_supplement,  Real value_plus1, Real value_minus1);
    
    void Helmholtz_complex_full_array(Array<Complex,3> R, Array<Complex,3> u, Real lambda_supplement, Complex value_plus1, Complex value_minus1);
};

#endif

//========================= Class declaration of IncFluid ends ============================== 
 
