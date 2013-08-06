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
/*! \file  time_advance_SF.cc
 * 
 * @brief Time advances velocity and passive scalar fields in Incompressible flow by one unit
 *
 * Status before entering the fn   <BR>
 *		nlin = FT[U.grad U]; U.F = p(k);  <BR>
 *		T.nlin = FT[U.grad T]
 *		(computed by Compute_nlin and compute_pressure)  <BR>
 *	
 *	Procedure (EULER)	 <BR>											
 *		1. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()     <BR>          
 *		2. V,T(t+dt) = V,T(t) + dt*nlin(V,T) in fn single_time_step() <BR>
 *		
 *	Procedure (RK2)	 <BR>
 *		1. Save V[i] in Vcopy[i] <BR>
 *		2. Compute_rhs(T): nlin[i] = -nlin[i] - grad(pressure) and T.nlin = -T.nlin <BR>
 *		3. V,T(t+dt) = V,T(t) + (dt/2)*nlin(V,T) using fn single_time_step(): MID POINT <BR>
 *		4. Compute_nlin(V(t+dt/2),T(t+dt/2)) <BR>
 *		5. compute p(t+dt/2) using compute_pressure() <BR>
 *		6. V,T[i] = V,Tcopy[i] (copy back) <BR>
 *		7. V,T(t+dt) = V,T(t) + dt*nlinU,T(t+dt/2) using fn single_time_step()  <BR>
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
                
   
#include "Time_advance.h"

//********************************************************************************************* 
void Time_advance_incompress::Euler(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force)
{
	
	Compute_force_TO_rhs(U,T1, T2,P,Force);
	
	Single_time_step(U, T1, T2, 1, 1, 1);  
}	


		
//*********************************************************************************************						
	
void Time_advance_incompress::RK2(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force)		
{
	// Allocated once and for all- Vcopy,Scopy  (static)
	static PlainCVF Vcopy;							
	static PlainCSF S1copy; 
	static PlainCSF S2copy;

	
	U.cvf.Copy_to(Vcopy);	
	T1.csf.Copy_to(S1copy);		// Vcopy[i] <- V[i]; Scopy <- T.csf.F     
	T2.csf.Copy_to(S2copy);	
	
	Compute_force_TO_rhs(U,T1,T2,P,Force);  
	  
	Single_time_step(U, T1, T2, 0.5, 0.5, 0.5);			// Goto the mid point
	
	Compute_force_TO_rhs(U,T1,T2,P,Force);

	U.cvf.Copy_from(Vcopy);		
	T1.csf.Copy_from(S1copy);	// Copy back into V[i],T
	T2.csf.Copy_from(S2copy);
	
	Single_time_step(U, T1, T2, 1, 0.5, 1);							
	// Time-step by global.time.Tdt now using the mid-point slopes 
}
	
//*********************************************************************************************			
	
void Time_advance_incompress::RK4(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force)	
{
	
	static PlainCVF Vcopy;													
	static PlainCSF S1copy;
	static PlainCSF S2copy;
	
	static PlainCVF tot_Vrhs;								
	static PlainCSF tot_S1rhs;
	static PlainCSF tot_S2rhs;
							
	
	U.cvf.Copy_to(Vcopy);	
	T1.csf.Copy_to(S1copy);		// Vcopy[i] <- V[i]; Scopy <- T.csf.F    
	T2.csf.Copy_to(S2copy);
	
	tot_Vrhs.Initialize();  
	tot_S1rhs.Initialize();
	tot_S2rhs.Initialize();
	
	// Step 1
	
	Compute_force_TO_rhs(U,T1,T2,P,Force);   
	Single_time_step(U, T1, T2, 0.5, 0.5, 0.5);				// u1: Goto the mid point

	// tot_Vrhs += (global.time.Tdt/6)*RHS(t)*exp(-nu k^2 b*dt) with b=1.0
	//	RHS is contained in *nlin
	Compute_RK4_parts(U, T1, T2, tot_Vrhs, tot_S1rhs, tot_S2rhs, 1.0, 1.0/6); 
	
	// Step 2
	
	Compute_force_TO_rhs(U,T1,T2,P,Force);
	
	U.cvf.Copy_from(Vcopy);		
	T1.csf.Copy_from(S1copy);	// Copy back into V[i],T
	T2.csf.Copy_from(S2copy);
	
	Single_time_step(U, T1, T2, 0.5, 0, 0.5);					// u2: Goto the mid point
	
	Compute_RK4_parts(U, T1, T2, tot_Vrhs, tot_S1rhs, tot_S2rhs, 0.5, 1.0/3);
	
	// Step 3
	
	Compute_force_TO_rhs(U,T1,T2,P,Force);
	
	U.cvf.Copy_from(Vcopy);		
	T1.csf.Copy_from(S1copy);	// Copy back into V[i],T
	T2.csf.Copy_from(S2copy);
	
	Single_time_step(U, T1, T2, 1, 0.5, 1);					// u3: Go to the end point
	
	Compute_RK4_parts(U, T1, T2, tot_Vrhs,  tot_S1rhs, tot_S2rhs, 0.5, 1.0/3);
	
	// Step 4
	
	Compute_force_TO_rhs(U,T1,T2,P,Force);
									
	Compute_RK4_parts(U, T1, T2, tot_Vrhs, tot_S1rhs, tot_S2rhs, 0, 1.0/6);
	
	// Final result: Assemble the parts
	
	U.cvf.Copy_from(Vcopy);		
	T1.csf.Copy_from(S1copy);	// Copy back into V[i],T
	T2.csf.Copy_from(S2copy);
	
	U.Mult_field_exp_ksqr_dt(1.0);  
	T1.Mult_field_exp_ksqr_dt(1.0);
	T2.Mult_field_exp_ksqr_dt(1.0);
	
		// Finally add tot_Vrsh and tot_S1rhs, tot_S2rhs
	U.cvf.V1 = U.cvf.V1 + (tot_Vrhs.V1);
	U.cvf.V2 = U.cvf.V2 + (tot_Vrhs.V2);
	U.cvf.V3 = U.cvf.V3 + (tot_Vrhs.V3);
	
	T1.csf.F = T1.csf.F + (tot_S1rhs.F);
	T2.csf.F = T2.csf.F + (tot_S2rhs.F);
}	

//**********************************   End of time_advance_SF.cc  *****************************


