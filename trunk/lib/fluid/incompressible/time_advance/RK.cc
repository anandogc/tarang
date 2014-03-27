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


/*! \file  time_advance.cc
 * 
 * @brief Time advances velocity field in Incompressible NS by one unit
 *
 * Status before entering the fn   <BR>
 *		nlin = FT[U.grad U]; U.F = p(k);  <BR>
 *		(computed by Compute_nlin and compute_pressure)  <BR>
 *	
 *	Procedure (EULER)	 <BR>											
 *		1. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()     <BR>          
 *		2. V(t+dt) = V(t) + dt*nlin in fn single_time_step() <BR>
 *		
 *	Procedure (RK2)	 <BR>
 *		1. Save V[i] in Vcopy[i] <BR>
 *		2. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()  <BR>
 *		3. V(t+dt) = V(t) + (dt/2)*nlin using fn single_time_step(): MID POINT <BR>
 *		4. compute_nlin(V(t+dt/2)) <BR>
 *		5. compute p(t+dt/2) using compute_pressure() <BR>
 *		6. V[i] = Vcopy[i] (copy back) <BR>
 *		7. V(t+dt) = V(t) + dt*nlin(t+dt/2) using fn single_time_step()  <BR>
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug No known bugs
 */                
   
#include "Time_advance.h"

//*********************************************************************************************
 
void Time_advance_incompress::Euler(FluidVF& U, Pressure& P, FORCE& Force)
{	
    Compute_force_TO_rhs(U,P,Force);

	Single_time_step(U, 1, 1, 1);
}
//*********************************************************************************************		

void Time_advance_incompress::RK2(FluidVF& U, Pressure& P, FORCE& Force)	
{	
	// Allocated once and for all- Vcopy location (static) 
	static PlainCVF Vcopy;
	
	U.Copy_field_to(Vcopy);								// Vcopy = U
	
	Compute_force_TO_rhs(U,P,Force);
	
	Single_time_step(U, 0.5, 0.5, 0.5);					// Goto the mid point
	
	Compute_force_TO_rhs(U,P,Force);
	
	U.Copy_field_from(Vcopy);								// U = Vcopy
	
	Single_time_step(U, 1, 0.5, 1);						
	// Time-step by Tdt now using the mid-point slopes 
}

	
//*********************************************************************************************

void Time_advance_incompress::RK4(FluidVF& U, Pressure& P, FORCE& Force)
{
	
	// Allocated once and for all- Vcopy location (static)
	static PlainCVF Vcopy;									
	static PlainCVF tot_Vrhs;
										
	U.Copy_field_to(Vcopy);
	tot_Vrhs.Initialize();
	
	// Step 1
	
	Compute_force_TO_rhs(U,P,Force);   
	Single_time_step(U, 0.5, 0.5, 0.5);					// u1: Goto the mid point

	// tot_Vrhs = (dt/6)*RHS(t)*exp(-nu k^2 b*dt) with b=1.0
	//	RHS is contained in *nlin
	Compute_RK4_parts(U, tot_Vrhs, 1.0, 1.0/6);
	
	// Step 2
	
	Compute_force_TO_rhs(U,P,Force);
	
	U.cvf.Copy_from(Vcopy);
	
	Single_time_step(U, 0.5, 0, 0.5);				// u2: Goto the mid point
	
	Compute_RK4_parts(U, tot_Vrhs, 0.5, 1.0/3);
	
	// Step 3
	
	Compute_force_TO_rhs(U,P,Force);
	
	U.cvf.Copy_from(Vcopy);
	
	Single_time_step(U, 1, 0.5, 1);		// u3: Go to the end point
	
	Compute_RK4_parts(U, tot_Vrhs, 0.5, 1.0/3);
	
	// Step 4
	
	Compute_force_TO_rhs(U,P,Force);
									
	Compute_RK4_parts(U, tot_Vrhs, 0, 1.0/6);	// tot_rhs += rhs(t+dt, u3)
	
	// Final result
	U.cvf.Copy_from(Vcopy);
	U.Mult_field_exp_ksqr_dt(1.0);
		
	U.cvf.V1 = U.cvf.V1 + (tot_Vrhs.V1);
	U.cvf.V2 = U.cvf.V2 + (tot_Vrhs.V2);
	U.cvf.V3 = U.cvf.V3 + (tot_Vrhs.V3);
}

//**********************************   End of time_advance.cc  ********************************
   

