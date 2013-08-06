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

/*! \file  time_advance_VF_MHD.cc
 * 
 * @brief Time advances velocity and magnetic fields in MHD flow by one unit
 *
 * Status before entering the fn   <BR>
 *		nlin = FT[U.grad U-B.grad B]; U.F = p(k);  <BR>
 *		W.nlin = FT[U.grad B - B.grad U]
 *		(computed by Compute_nlin and compute_pressure)  <BR>
 *	
 *	Procedure (EULER)	 <BR>											
 *		1. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()     <BR>          
 *		2. V,W(t+dt) = V,W(t) + dt*nlin(V,W) in fn single_time_step() <BR>
 *		
 *	Procedure (RK2)	 <BR>
 *		1. Save V[i] in Vcopy[i] <BR>
 *		2. Compute_rhs(T): nlin[i] = -nlin[i] - grad(pressure) and T.nlin = -T.nlin <BR>
 *		3. V,T(t+dt) = V,T(t) + (dt/2)*nlin(V,T) using fn single_time_step(): MID POINT <BR>
 *		4. Compute_nlin(V(t+dt/2),T(t+dt/2)) <BR>
 *		5. compute p(t+dt/2) using compute_pressure() <BR>
 *		6. V,T[i] = V,Tcopy[i] (copy back) <BR>
 *		7. V,T(t+dt) = V,T(t) + dt*nlinU,T(t+dt/2) using fn single_time_step()  <BR> *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 */
          
   
#include "Time_advance.h"

//***************************************************************************


void Time_advance_incompress::Euler(FluidVF &U, FluidVF& W, Pressure &P, FORCE& Force)
{
    global.time.previous = global.time.now;
    global.time.keplerian = global.time.now;
    
    Compute_force_TO_rhs(U,W,P,Force);
    
	Single_time_step(U, W, 1, 1, 1);
	// Single-time-step both U and W
}

	
void Time_advance_incompress::RK2(FluidVF &U, FluidVF& W, Pressure &P, FORCE& Force)
{
	if (global.program.kind != "KEPLERIAN") {
        // Allocated once and for all- Vcopy, Wcopy (static)
        static PlainCVF Vcopy;
        static PlainCVF Wcopy;
        
        U.cvf.Copy_to(Vcopy);
        W.cvf.Copy_to(Wcopy);
        
        Compute_force_TO_rhs(U,W,P,Force);
        Single_time_step(U, W, 0.5, 0.5, 0.5);
        // Goto the mid point
        
        global.time.keplerian += global.time.dt/2;
        Compute_force_TO_rhs(U,W,P,Force);
        
        U.cvf.Copy_from(Vcopy);
        W.cvf.Copy_from(Wcopy);
        
        Single_time_step(U, W, 1, 0.5, 1);
        // Time-step by global.time.Tdt now using the mid-point slopes
    }
    
    else { // Keplerian
        // Allocated once and for all- Vcopy, Wcopy (static)
        static PlainCVF Vcopy;		
        static PlainCVF Wcopy;     
        
        U.cvf.Copy_to(Vcopy);
        W.cvf.Copy_to(Wcopy);
     
        global.time.previous = global.time.now;
        global.time.keplerian = global.time.now;
        
        Compute_force_TO_rhs(U,W,P,Force);
        Single_time_step(U, W, 0.5, 0.5, 0.5);			
        // Goto the mid point

        // Compute nlin at t=t0+dt/2
        global.time.keplerian += global.time.dt/2;
        Compute_force_TO_rhs(U,W,P,Force);
        
        U.cvf.Copy_from(Vcopy);	
        W.cvf.Copy_from(Wcopy);	

        global.time.keplerian = global.time.previous; // t=initial time for exp() calc
        Single_time_step(U, W, 1, 0.5, 1);							
        // Time-step by global.time.Tdt now using the mid-point slopes
    } 
}
	
	


void Time_advance_incompress::RK4(FluidVF &U, FluidVF& W, Pressure &P, FORCE& Force)
{
	if (global.program.kind != "KEPLERIAN") {
        static PlainCVF Vcopy;
        static PlainCVF Wcopy;
                            
        static PlainCVF tot_Vrhs;								
        static PlainCVF tot_Wrhs;
                            
        U.cvf.Copy_to(Vcopy);
        W.cvf.Copy_to(Wcopy);
        
        tot_Vrhs.Initialize(); 
        tot_Wrhs.Initialize();
        
        Compute_force_TO_rhs(U,W,P,Force);       
        Single_time_step(U, W, 0.5, 0.5, 0.5);				// u1: Goto the mid point

        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 1.0, 1.0/6); 
        
        // Step 2
        
        Compute_force_TO_rhs(U,W,P,Force);
        
        U.cvf.Copy_from(Vcopy);	
        W.cvf.Copy_from(Wcopy);
        
        Single_time_step(U, W, 0.5, 0, 0.5);					// u2: Goto the mid pt 
        
        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 0.5, 1.0/3);
        
        // Step 3
        
        Compute_force_TO_rhs(U,W,P,Force);
        
        U.cvf.Copy_from(Vcopy);	
        W.cvf.Copy_from(Wcopy);
        
        Single_time_step(U, W, 1, 0.5, 1);					// u3 : Goto the mid pt
        
        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 0.5, 1.0/3);
        
        // Step 4
        
        Compute_force_TO_rhs(U,W,P,Force);
                        
        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 0, 1.0/6);
        
        // Final result
        
        U.cvf.Copy_from(Vcopy);	
        W.cvf.Copy_from(Wcopy); 
        
        U.Mult_field_exp_ksqr_dt(1.0);
        W.Mult_field_exp_ksqr_dt(1.0);
                    
        U.cvf.V1 = U.cvf.V1 + (tot_Vrhs.V1);
        U.cvf.V2 = U.cvf.V2 + (tot_Vrhs.V2);
        U.cvf.V3 = U.cvf.V3 + (tot_Vrhs.V3);
        
        W.cvf.V1 = W.cvf.V1 + (tot_Wrhs.V1);
        W.cvf.V2 = W.cvf.V2 + (tot_Wrhs.V2);		 
        W.cvf.V3 = W.cvf.V3 + (tot_Wrhs.V3);
    }
    
    else {
        static PlainCVF Vcopy;
        static PlainCVF Wcopy;
        
        static PlainCVF tot_Vrhs;
        static PlainCVF tot_Wrhs;
        
        U.cvf.Copy_to(Vcopy);
        W.cvf.Copy_to(Wcopy);
        
        tot_Vrhs.Initialize();
        tot_Wrhs.Initialize();
        
        global.time.previous = global.time.now;
        global.time.keplerian = global.time.now;
        
        Compute_force_TO_rhs(U,W,P,Force);
        Single_time_step(U, W, 0.5, 0.5, 0.5);				// u1: Goto the mid point
        
        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 1.0, 1.0/6);
        
        // Step 2
        global.time.keplerian += global.time.dt/2; // for both steps 2 and 3
        
        Compute_force_TO_rhs(U,W,P,Force);
        
        U.cvf.Copy_from(Vcopy);
        W.cvf.Copy_from(Wcopy);
        
        Single_time_step(U, W, 0.5, 0, 0.5);					// u2: Goto the mid pt
        
        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 0.5, 1.0/3);
        
        // Step 3
        
        Compute_force_TO_rhs(U,W,P,Force);
        
        U.cvf.Copy_from(Vcopy);
        W.cvf.Copy_from(Wcopy);
        
        Single_time_step(U, W, 1, 0.5, 1);					// u3 : Goto the mid pt
        
        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 0.5, 1.0/3);
        
        // Step 4
        
        global.time.keplerian  = global.time.now; // t0+dt
        Compute_force_TO_rhs(U,W,P,Force);
        
        Compute_RK4_parts(U, W, tot_Vrhs, tot_Wrhs, 0, 1.0/6);
        
        // Final result
        global.time.keplerian = global.time.previous; // t=initial time for exp() calc
        
        U.cvf.Copy_from(Vcopy);
        W.cvf.Copy_from(Wcopy);
        
        U.Mult_field_exp_ksqr_dt(1.0);
        W.Mult_field_exp_ksqr_dt(1.0);
        
        U.cvf.V1 = U.cvf.V1 + (tot_Vrhs.V1);
        U.cvf.V2 = U.cvf.V2 + (tot_Vrhs.V2);
        U.cvf.V3 = U.cvf.V3 + (tot_Vrhs.V3);
        
        W.cvf.V1 = W.cvf.V1 + (tot_Wrhs.V1);
        W.cvf.V2 = W.cvf.V2 + (tot_Wrhs.V2);
        W.cvf.V3 = W.cvf.V3 + (tot_Wrhs.V3);
    }
}


//**********************************   End of time_advance_VF_MHD.cc  *************************



