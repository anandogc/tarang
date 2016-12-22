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

#include "Time_advance.h"

//*********************************************************************************************
// Add force to nlin

void Time_advance_incompress::Add_force(FluidVF& U)
{
	
	if (U.force_switch) {
		U.nlin1 = U.nlin1 - U.Force1;	
		U.nlin2 = U.nlin2 - U.Force2;
		U.nlin3 = U.nlin3 - U.Force3;
	}
}

	// Passive scalar or RB convection

void Time_advance_incompress::Add_force(FluidVF& U, FluidSF& T)
{
	
	if (global.program.kind == "INC_SCALAR")
		Add_force_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Add_force_RBC(U, T);						
}

void Time_advance_incompress::Add_force_scalar(FluidVF& U, FluidSF& T)
{
	
	Add_force(U);
	
	if (T.force_switch)
		T.nlin = T.nlin - T.Force;						
}

void Time_advance_incompress::Add_force_RBC(FluidVF& U, FluidSF& T)
{
	if (global.program.basis_type == "ChFF")
		Add_force_scalar(U,T);
	
	else {
		if (global.PHYSICS.Pr_option == "PRZERO")
			Add_force(U);
		
		else if (global.PHYSICS.Pr_option == "PRINFTY") 
			T.nlin = T.nlin - T.Force;
		
		else
			Add_force_scalar(U,T);
	}
}

	// MHD

void Time_advance_incompress::Add_force(FluidVF& U, FluidVF& W)
{
	Add_force(U);
	
	if (W.force_switch) {
		W.nlin1 = W.nlin1 - W.Force1;	
		W.nlin2 = W.nlin2 - W.Force2;
		W.nlin3 = W.nlin3 - W.Force3;	
	}	
}

	// MHD + Scalar

void Time_advance_incompress::Add_force(FluidVF& U, FluidVF& W, FluidSF& T)
{
	Add_force(U, W);
	
	if (T.force_switch)
		T.nlin = T.nlin - T.Force;						
}	


void Time_advance_incompress::Add_force(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	Add_force(U, W, T);
	
	if (C.force_switch)
		C.nlin = C.nlin - C.Force;
}

void Time_advance_incompress::Add_force(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	
	Add_force(U);
	
	T1.nlin = T1.nlin - T1.Force;	
	T2.nlin = T2.nlin - T2.Force;	
}

//******************************************************************************

/*void Time_advance_incompress::Add_ubydt(FluidVF& U)
{
	if (U.force_switch) {
		U.nlin1 = U.nlin1 - U.cvf.V1/global.time.dt;
		U.nlin2 = U.nlin2 - U.cvf.V2/global.time.dt;
		U.nlin3 = U.nlin3 - U.cvf.V3/global.time.dt;
	}
} */

//******************************************************************************

// nlin = nlin + grad(p)
void Time_advance_incompress::Add_pressure_gradient(FluidVF& U, Pressure& P)       
{
	
	global.program.sincostr_switch = global.program.sincostr_switch_divergence;
	universal->Add_Xderiv(P.F, U.nlin1);
	
	global.program.sincostr_switch = global.program.sincostr_switch_divergence;
	universal->Add_Yderiv(P.F, U.nlin2);
	
	global.program.sincostr_switch = global.program.sincostr_switch_divergence;	   	
	universal->Add_Zderiv(P.F, U.nlin3 );
    

    if (global.program.kind == "KEPLERIAN") {
        Real omega_keplerian = global.force.double_para(0);
        Real q_keplerian = global.force.double_para(1);
        
        Real q_omega_t = q_keplerian*omega_keplerian*global.time.previous;
        
        universal->Yderiv(P.F, global.temp_array.X);
        U.nlin1 = U.nlin1 + Complex(q_omega_t,0)*global.temp_array.X;
    }
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
void Time_advance_incompress::Compute_rhs(FluidVF& U, Pressure& P)       
{	

	Add_pressure_gradient(U, P);				// nlin = nlin + grad(p)

	U.nlin1 = -U.nlin1;  
	U.nlin2 = -U.nlin2;
	U.nlin3 = -U.nlin3;
}


//*********************************************************************************************

/** @brief Compute rhs for fluid simulation with scalar
 * 
 * @note Status before entering the fn
 *		nlin[i] contains Dj T[Vr[j]*Vr[i]] - f[i];
 *		F.p contains pressure;
 *		T.nlin = T[U.grad T]	
 * 
 * @return (rhs) = nlin[i] =  -Dj T[Vr[j]*Vr[i]] + f[i] - Di [p(k)] (-grad(pressure)) 
 *			(T.nlin) = T.nlin = -T.nlin	
 */

void Time_advance_incompress::Compute_rhs(FluidVF& U, FluidSF& T, Pressure& P)
{
	if (global.program.kind == "INC_SCALAR") 
		Compute_rhs_scalar(U, T, P);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Compute_rhs_RBC(U, T, P);
	
}


void Time_advance_incompress::Compute_rhs_scalar(FluidVF& U, FluidSF& T, Pressure& P)       
{
	Compute_rhs(U, P);				// Fluid:  U.nlin[i] = -U.nlin[i] -FT[grad(p)]
	
	T.nlin = -T.nlin;			// For scalar: rhs = -T.nlin  
}


void Time_advance_incompress::Compute_rhs_RBC(FluidVF& U, FluidSF& T, Pressure& P)       
{
	
	if (global.PHYSICS.Pr_option == "PRZERO")
		Compute_rhs(U, P);
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") 
		T.nlin = -T.nlin;
	
	else 
		Compute_rhs_scalar(U, T, P);
}


void Time_advance_incompress::Compute_rhs(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P)       
{
	Compute_rhs(U, P);				// Fluid:  U.nlin[i] = -U.nlin[i] -FT[grad(p)]
	T1.nlin = -T1.nlin;				// For scalars: rhs = -T.nlin 
	T2.nlin = -T2.nlin;	
}

//*********************************************************************************************

/** @brief Compute rhs for fluid simulation with a vector
 * 
 * @note Status before entering the fn
 *		nlin[i] contains Dj T[Vr[j]*Vr[i]] - f[i];
 *		F.p contains pressure;
 *		W.nlin = T[U.grad B - B.grad U]	
 * 
 * @return (rhs) = nlin[i] =  -Dj T[Vr[j]*Vr[i]] + f[i] - Di [p(k)] (-grad(pressure)) 
 *			(W.nlin) = W.nlin[i] = -W.nlin[i]
 */
void Time_advance_incompress::Compute_rhs(FluidVF& U, FluidVF& W, Pressure& P)       
{

	Compute_rhs(U, P);					// Fluid:  U.nlin[i] = -U.nlin[i] -FT[grad(p)]
	   
	W.nlin1 = -W.nlin1;			// VF-W rhs[i] = -nlin[i]
	W.nlin2 = -W.nlin2;
	W.nlin3 = -W.nlin3;
}


//*********************************************************************************************

void Time_advance_incompress::Compute_rhs(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P)       
{
	Compute_rhs(U, W, P);
	
	T.nlin = -T.nlin;			// For scalar: rhs = -T.nlin     
}



void Time_advance_incompress::Compute_rhs(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P)
{
	Compute_rhs(U, W, T, P);
	
	C.nlin = -C.nlin;			// For scalar: rhs = -T.nlin
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
void Time_advance_incompress::Single_time_step(FluidVF& U, Real a, Real b, Real c)
{	
	
	if ((abs(a) < MYEPS) && (abs(b) < MYEPS)) 
		U.Add_nlin_factor_dt(c);
	
	else if ((abs(a) < MYEPS) && (abs(b) > MYEPS)) {
		U.Mult_nlin_exp_ksqr_dt(b);
		U.Add_nlin_factor_dt(c);
	}
	
	else if ((abs(b) < MYEPS) && (abs(a) > MYEPS)) {
		U.Mult_field_exp_ksqr_dt(a);
		U.Add_nlin_factor_dt(c);	
	}
	
	else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS)) {
		U.Add_nlin_factor_dt(c);									// V = V + c*dt*nlin
		U.Mult_field_exp_ksqr_dt(a);
	}
	
	else if ((abs(a) > MYEPS) && (abs(b) > MYEPS)) {
		U.Mult_field_exp_ksqr_dt((a-b));
		U.Add_nlin_factor_dt(c);	
		U.Mult_field_exp_ksqr_dt(b);
	}
}


//*********************************************************************************************
void Time_advance_incompress::Single_time_step(FluidSF& T, Real a, Real b, Real c)
{	
	
	if ((abs(a) < MYEPS) && (abs(b) < MYEPS)) 
		T.Add_nlin_factor_dt(c);
	
	else if ((abs(a) < MYEPS) && (abs(b) > MYEPS)) {
		T.Mult_nlin_exp_ksqr_dt(b);
		T.Add_nlin_factor_dt(c);
	}
	
	else if ((abs(b) < MYEPS) && (abs(a) > MYEPS))	{
		T.Mult_field_exp_ksqr_dt(a);
		T.Add_nlin_factor_dt(c);
	}
	
	else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS)) {
		T.Add_nlin_factor_dt(c);								
		T.Mult_field_exp_ksqr_dt(a);
	}
	
	else if ((abs(a) > MYEPS) && (abs(b) > MYEPS)) {
		T.Mult_field_exp_ksqr_dt((a-b));
		T.Add_nlin_factor_dt(c);	
		T.Mult_field_exp_ksqr_dt(b);
	}
	
}


//*********************************************************************************************	
//	Passive scalar  
void Time_advance_incompress::Single_time_step(FluidVF& U, FluidSF& T, Real a, Real b, Real c)
{	
	if (global.program.kind == "INC_SCALAR") 
		Single_time_step_scalar(U, T, a, b, c);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Single_time_step_RBC(U, T, a, b, c);
	
}	


//*********************************************************************************************
void Time_advance_incompress::Single_time_step_scalar(FluidVF& U, FluidSF& T, Real a, Real b, Real c)
{	
	Single_time_step(U, a, b, c);
	Single_time_step(T, a, b, c);	
}


//*********************************************************************************************
void Time_advance_incompress::Single_time_step_RBC(FluidVF& U, FluidSF& T, Real a, Real b, Real c)
{
	
	if (global.PHYSICS.Pr_option == "PRZERO") {	
		Single_time_step(U, a, b, c);						// V advance
		U.Zero_Prandtl_number_compute_temperature(T);
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
		Single_time_step(T, a, b, c);	// Time advance T
		U.Infinite_Prandtl_number_compute_velocity(T); 
	}
	
	else 
		Single_time_step_scalar(U, T, a, b, c);
} 

//*********************************************************************************************
void Time_advance_incompress::Single_time_step(FluidVF& U, FluidSF& T1, FluidSF& T2, Real a, Real b, Real c)
{	
	Single_time_step(U, a, b, c);
	Single_time_step(T1, a, b, c);	
	Single_time_step(T2, a, b, c);
}


//*********************************************************************************************
//	MHD		

void Time_advance_incompress::Single_time_step(FluidVF& U, FluidVF& W, Real a, Real b, Real c)
{  
	Single_time_step(U, a, b, c);
	Single_time_step(W, a, b, c);
}  


//*********************************************************************************************
//	MHD + passive scalar  

void Time_advance_incompress::Single_time_step(FluidVF& U, FluidVF& W, FluidSF& T, Real a, Real b, Real c)
{
	Single_time_step(U, a, b, c);
	Single_time_step(W, a, b, c);
	Single_time_step(T, a, b, c);
}

// MHD Astro
void Time_advance_incompress::Single_time_step(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Real a, Real b, Real c)
{
	Single_time_step(U, a, b, c);
	Single_time_step(W, a, b, c);
	Single_time_step(T, a, b, c);
	Single_time_step(C, a, b, c);
}

//*********************************************************************************************
/** @brief Compute_force_TO_rhs()
 *
 * @return Computes force and nlin; add them; compute pressure; combine them to compute rhs.
 *
 */
void Time_advance_incompress::Compute_force_TO_rhs(FluidVF& U, Pressure& P, FORCE& Force)
{

	Force.Compute_force(U);
	
	Nlin_incompress::Compute_nlin(U);					// Compute nlin using V(t+dt/2)

	Add_force(U);											// nlin = nlin - f

	P.Compute_pressure(U);									// Compute pressure using V(t+dt/2)
	Compute_rhs(U, P);
	
}

void Time_advance_incompress::Compute_force_TO_rhs(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force)
{
	if (global.program.kind == "INC_SCALAR") 
		Compute_force_TO_rhs_scalar(U, T, P, Force);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Compute_force_TO_rhs_RBC(U, T, P, Force);
}

void Time_advance_incompress::Compute_force_TO_rhs_scalar(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force)
{
	Force.Compute_force(U, T);
	
	Nlin_incompress::Compute_nlin(U, T);					// Compute nlin using V(t+dt/2)
	
	Add_force(U, T);										// nlin = nlin - f
	
	P.Compute_pressure(U);									// Compute pressure using V(t+dt/2)
	
	Compute_rhs(U, T, P);		
}

void Time_advance_incompress::Compute_force_TO_rhs_RBC(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force)
{
	Force.Compute_force(U, T);
	
	Nlin_incompress::Compute_nlin(U, T);						// Compute nlin using V(t+dt/2)
    
	Add_force(U, T);											// nlin = nlin - f

	if (global.PHYSICS.Pr_option != "PRINFTY")
		P.Compute_pressure(U);									// Compute pressure using V(t+dt/2)

	Compute_rhs(U, T, P);
}


void Time_advance_incompress::Compute_force_TO_rhs(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force)
{
	Force.Compute_force(U, T1, T2);
	
	Nlin_incompress::Compute_nlin(U, T1, T2);					// Compute nlin using V(t+dt/2)
	
	Add_force(U, T1, T2);										// nlin = nlin - f
	
	P.Compute_pressure(U);									// Compute pressure using V(t+dt/2)
	
	Compute_rhs(U, T1, T2, P);		
}


void Time_advance_incompress::Compute_force_TO_rhs(FluidVF& U, FluidVF& W, Pressure& P, FORCE& Force)
{
	Force.Compute_force(U, W);
	
	Nlin_incompress::Compute_nlin(U, W);							// Compute nlin using V(t+dt/2)

	Add_force(U, W);											// nlin = nlin - f
	
	P.Compute_pressure(U);										// Compute pressure using V(t+dt/2)
    
	Compute_rhs(U, W, P);
}


void Time_advance_incompress::Compute_force_TO_rhs(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P, FORCE& Force)
{
	
	Force.Compute_force(U, W, T);
	
	Nlin_incompress::Compute_nlin(U, W, T);					// Compute nlin using V(t+dt/2)
	
	Add_force(U, W, T);										// nlin = nlin - f
	
	P.Compute_pressure(U);										// Compute pressure using V(t+dt/2)
	
	Compute_rhs(U, W, T, P);
}


void Time_advance_incompress::Compute_force_TO_rhs(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P, FORCE& Force)
{
	Force.Compute_force(U, W, T, C);

	Nlin_incompress::Compute_nlin(U, W, T, C);					// Compute nlin using V(t+dt/2)

	Add_force(U, W, T, C);										// nlin = nlin - f

	P.Compute_pressure(U);										// Compute pressure using V(t+dt/2)

	Compute_rhs(U, W, T, C, P);
}

// GP
void Time_advance_incompress::Compute_force_TO_rhs(FluidSF& T, FORCE& Force)
{
//	Force.Compute_force(T);
	
	Nlin_incompress::Compute_nlin(T);					// nlin in rhs
	
    // T.nlin = T.nlin + T.Force;  // THis is the final (rhs)
}

//*********************************************************************************************
//*********************************************************************************************
/** @brief Compute_RK4_parts(PlainCVF& tot_Vrhs, Real dt, Real b, Real factor) conmputes Ci's for 
 *		for computing fields at t+dt.
 *
 * @return tot_Vrhs = (U+factor*dt*RHS(t))*exp(-nu k^2 b*dt)
 *			RHS is contained in *nlin
 */
void Time_advance_incompress::Compute_RK4_parts(FluidVF& U, PlainCVF& tot_Vrhs, Real b, Real factor)
{	
	U.Mult_nlin_exp_ksqr_dt(b);				// U.nlin = U.nlin x exp(-nu k^2 b*dt)
	U.Add_field_nlin_factor_dt(tot_Vrhs, factor);		// tot_Vrhs += factor*dt*U.nlin
}

void Time_advance_incompress::Compute_RK4_parts(FluidSF& T, PlainCSF& tot_Srhs, Real b, Real factor)
{	
	T.Mult_nlin_exp_ksqr_dt(b);					// T.nlin = T.nlin x exp(kappa k^2 b*dt)
	T.Add_field_nlin_factor_dt(tot_Srhs, factor);		// tot_Srhs += factor*dt*S.nlin
}


void Time_advance_incompress::Compute_RK4_parts(FluidVF& U, FluidSF& T, PlainCVF& tot_Vrhs, PlainCSF& tot_Srhs, 
                                                Real b, Real factor)
{	
	if (global.program.kind == "INC_SCALAR") 
		Compute_RK4_parts_scalar(U, T,tot_Vrhs, tot_Srhs, b, factor);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Compute_RK4_parts_RBC(U, T,tot_Vrhs, tot_Srhs, b, factor);
}


void Time_advance_incompress::Compute_RK4_parts_scalar(FluidVF& U, FluidSF& T, PlainCVF& tot_Vrhs, 
                                                       PlainCSF& tot_Srhs, Real b, Real factor)
{	
	Compute_RK4_parts(U, tot_Vrhs, b, factor);
	Compute_RK4_parts(T, tot_Srhs, b, factor);
}


void Time_advance_incompress::Compute_RK4_parts_RBC(FluidVF& U, FluidSF& T, PlainCVF& tot_Vrhs, 
                                                    PlainCSF& tot_Srhs, Real b, Real factor)
{	
	if (global.PHYSICS.Pr_option == "PRZERO") {
		Compute_RK4_parts(U, tot_Vrhs, b, factor);
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY")	{
		Compute_RK4_parts(T, tot_Srhs, b, factor);	
	}
	
	else 
		Compute_RK4_parts_scalar(U, T, tot_Vrhs, tot_Srhs, b, factor);
}

void Time_advance_incompress::Compute_RK4_parts(FluidVF& U, FluidSF& T1, FluidSF& T2, PlainCVF& tot_Vrhs, PlainCSF& tot_S1rhs, PlainCSF& tot_S2rhs, Real b, Real factor)
{	
	Compute_RK4_parts(U, tot_Vrhs, b, factor);
	Compute_RK4_parts(T1, tot_S1rhs, b, factor);
	Compute_RK4_parts(T2, tot_S2rhs, b, factor);
}


void Time_advance_incompress::Compute_RK4_parts(FluidVF& U, FluidVF& W, PlainCVF& tot_Vrhs, PlainCVF& tot_Wrhs, Real b, Real factor)
{	
	Compute_RK4_parts(U, tot_Vrhs, b, factor);
	Compute_RK4_parts(W, tot_Wrhs, b, factor);
}


void Time_advance_incompress::Compute_RK4_parts(FluidVF& U, FluidVF& W, FluidSF& T, PlainCVF& tot_Vrhs,PlainCVF& tot_Wrhs,PlainCSF& tot_Srhs, Real b, Real factor)
{	
	Compute_RK4_parts(U, tot_Vrhs, b, factor);
	Compute_RK4_parts(W, tot_Wrhs, b, factor);
	Compute_RK4_parts(T, tot_Srhs, b, factor);
}

void Time_advance_incompress::Compute_RK4_parts(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, PlainCVF& tot_Vrhs,PlainCVF& tot_Wrhs,PlainCSF& tot_Srhs, PlainCSF& tot_Crhs, Real b, Real factor)
{
	Compute_RK4_parts(U, tot_Vrhs, b, factor);
	Compute_RK4_parts(W, tot_Wrhs, b, factor);
	Compute_RK4_parts(T, tot_Srhs, b, factor);
	Compute_RK4_parts(C, tot_Crhs, b, factor);
}


void Time_advance_incompress::Make_field_incompressible(FluidVF& U)
{
    Real total_abs_div;
    
    universal->Compute_divergence(U.cvf.V1, U.cvf.V2, U.cvf.V3, global.temp_array.X2, "field", total_abs_div, false);
    
    Real omega_keplerian = global.force.double_para(0);
    Real q_keplerian = global.force.double_para(1);
    
    Real q_omega_t = q_keplerian*omega_keplerian*global.time.now;
    
    Real Kx,Ky,Kz, mu_x, Ksqr, mu_sqr;
    
    Complex Ux000 = U.cvf.V1(0,0,0);
    Complex Uy000 = U.cvf.V2(0,0,0);
    Complex Uz000 = U.cvf.V3(0,0,0);
    
    for (int lx=0; lx<local_Nx; lx++)
        for (int ly=0; ly<Ny; ly++)
            for (int lz=0; lz<=Nz/2; lz++) {
                Kx = (universal->Get_kx(lx))*kfactor[1];
                Ky = (universal->Get_ky(ly))*kfactor[2];
                Kz = (universal->Get_kz(lz))*kfactor[3];
                Ksqr = pow2(Kx)+pow2(Ky)+pow2(Kz);
                mu_sqr = Ksqr+pow2(q_omega_t*Ky)+2*q_omega_t*Kx*Ky;
                mu_x = Kx+q_omega_t*Ky;
                
                U.cvf.V1(lx,ly,lz) = U.cvf.V1(lx,ly,lz) + global.temp_array.X2(lx,ly,lz)*Complex(0,mu_x/mu_sqr);
                U.cvf.V2(lx,ly,lz) = U.cvf.V2(lx,ly,lz) + global.temp_array.X2(lx,ly,lz)*Complex(0,Ky/mu_sqr);
                U.cvf.V3(lx,ly,lz) = U.cvf.V3(lx,ly,lz) + global.temp_array.X2(lx,ly,lz)*Complex(0,Kz/mu_sqr);
            }
    
    if (master) {
        U.cvf.V1(0,0,0) = Ux000;
        U.cvf.V2(0,0,0) = Uy000;
        U.cvf.V3(0,0,0) = Uz000;
    }
}


//**********************************   End of compute_rhs.cc  *********************************

