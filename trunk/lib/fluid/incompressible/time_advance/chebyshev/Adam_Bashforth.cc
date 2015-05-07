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

// R = nlin-F
// S = R-u^n/dt
//****************************************************************************************

void Time_advance_incompress::Adam_Bashforth(FluidVF& U, Pressure& P, FORCE& Force)	
{
	Real lambda_supplement;
	
	static Array<Complex,3> R(shape_complex_array);
	static Array<Complex,3> nlinV1_prev(shape_complex_array);
	static Array<Complex,3> nlinV2_prev(shape_complex_array);
	static Array<Complex,3> nlinV3_prev(shape_complex_array);
    
	Force.Compute_force(U);
	
	Nlin_incompress::Compute_nlin(U);
	
	Add_force(U);
	// S = nlin = nlin-f
	
	if (abs(global.time.now - global.time.init) < MYEPS) {  // init: nlin_prev = nlin_now
		nlinV1_prev = U.nlin1;  U.nlin1 *= TWO;
		nlinV2_prev = U.nlin2;  U.nlin2 *= TWO;
		nlinV1_prev = U.nlin3;  U.nlin3 *= TWO;
	}
	else {
		R = U.nlin1;  U.nlin1 = THREE*(U.nlin1) - nlinV1_prev;  nlinV1_prev = R;
		R = U.nlin2;  U.nlin2 = THREE*(U.nlin2) - nlinV2_prev;  nlinV2_prev = R;
		R = U.nlin3;  U.nlin3 = THREE*(U.nlin3) - nlinV3_prev;  nlinV3_prev = R;
	}
	
	
	Add_pressure_gradient(U, P);
    
	// Pressure computation
    // rhs = div(3N^n -N^n-1+\grad p^n)
	U.Compute_divergence_nlin(global.temp_array.X2);
    
    Helmholtz_complex_full_array(global.temp_array.X2, P.F, 0, Complex(0,0), Complex(0,0));
    // P.F = p_particular; lambda_supplement=0 (third arg).
    
	// V1 computation
	// For Pr=infty, we U/dt is set to zero; nlin=-force for this case
	universal->Subtract_Laplacian(U.dissipation_coefficient, U.cvf.V1, U.nlin1);  // U.nlin1 -= nu*Laplcian(U.cvf.V1)
    
    U.nlin1 -=  U.cvf.V1/global.time.dt;
    lambda_supplement = 2/(global.time.dt*U.dissipation_coefficient);
	
	universal->Add_Xderiv(P.F, U.nlin1);
	U.nlin1 /= (-U.dissipation_coefficient);
	
    Helmholtz_complex_full_array(U.nlin1, U.cvf.V1, lambda_supplement, Complex(0,0), Complex(0,0));
    // X2 = u_particular
    
	Compute_pressure_Ux_chebyshev(U, P);
    // P.F has pressure now; U.cvf.V1 = Vx
	
    // V2 soln
	if (!global.program.two_dimension)  {
		universal->Subtract_Laplacian(U.dissipation_coefficient, U.cvf.V2, U.nlin2);
		
        U.nlin2 -= U.cvf.V2/global.time.dt; // S computed
        lambda_supplement = 1/(global.time.dt*U.dissipation_coefficient);
		
		universal->Add_Yderiv(P.F, U.nlin2);
		U.nlin2 /= (-U.dissipation_coefficient);
		
		Helmholtz_complex_full_array(U.nlin2, U.cvf.V2, lambda_supplement, Complex(0,0), Complex(0,0));
	}
    
    // V3 soln
	universal->Subtract_Laplacian(U.dissipation_coefficient, U.cvf.V3, U.nlin3);
	
    U.nlin3 -= U.cvf.V3/global.time.dt;  // S computed
    lambda_supplement = 1/(global.time.dt*U.dissipation_coefficient);

	universal->Add_Zderiv(P.F, U.nlin3);
	U.nlin3 /= (-U.dissipation_coefficient);
    
    Helmholtz_complex_full_array(U.nlin3, U.cvf.V3, lambda_supplement, Complex(0,0), Complex(0,0));
}


//**********************************   End of time_advance.cc  ********************************
   

