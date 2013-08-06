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

	// r = h + u/dt
	
//*********************************************************************************************		

void Time_advance_incompress::BDF1(FluidVF& U, Pressure& P, FORCE& Force)	
{	
/*	int kx,ky,kz;
	DP Kperpsqr, lambda, sqrKx;
    complx vx_prime_plus, vx_prime_minus;
    complx delta_plus, delta_minus;
	
	Force.Compute_force(U);
	
	Nlin_incompress::Compute_nlin(U);	
	
	Add_force(U);
	// h = nlin = nlin-f

	// Compute particular solution for pressure. div(r) = div(h)
	U.Compute_divergence_nlin(global.temp_array.X2);
    
    // Make sure Xr is not used anywhere else
    Transform::Transpose_array_SLAB(global.temp_array.X2, global.temp_array.Xr, 'X', 'Y', 'Z', 'Z');
    
    Transform::Transpose_array_SLAB(U.cvf.V1, U.rvf.V1r, 'X', 'Y', 'Z', 'Z');
    
    Transform::Transpose_array_SLAB(U.cvf.V2, U.rvf.V2r, 'X', 'Y', 'Z', 'Z');
    
    Transform::Transpose_array_SLAB(U.cvf.V3, U.rvf.V3r, 'X', 'Y', 'Z', 'Z');
    
    Transform::Transpose_array_SLAB(U.nlin1, global.temp_array.nlin1TR, 'X', 'Y', 'Z', 'Z');
    
    Transform::Transpose_array_SLAB(U.nlin2, global.temp_array.nlin2TR, 'X', 'Y', 'Z', 'Z');
    
    Transform::Transpose_array_SLAB(U.nlin3, global.temp_array.nlin3TR, 'X', 'Y', 'Z', 'Z');
	
    
	for(int ly=0; ly<local_Ny; ly++) 
		for (int lz=0; lz<global.field.maxlz; lz++) {
			global.temp_array.in_helm_complex = global.temp_array.Xr(ly, Range::all(), lz);
			
			ky = universal->Get_ky(ly);
			kz = universal->Get_kz(lz);
			
			Kperpsqr = pow2(ky*kfactor[2]) + pow2(kz*kfactor[3]);
            lambda = Kperpsqr;
        
			// Compute p_particular p_1d
            Helmholtz_complex(global.temp_array.in_helm_complex, global.temp_array.out_helm_complex, lambda, 0, 0);	// out = p_p
            
			global.temp_array.pressure_1d = global.temp_array.out_helm_complex; 
            
            // Compute particular soln for V1_1d
            Xderiv1d_Chebyshev(global.temp_array.pressure_1d, global.temp_array.pressure_derivative_1d);
		
            global.temp_array.nlin1_1d = global.temp_array.nlin1TR(ly,Range::all(),lz) -U.rvf.V1r(ly,Range::all(),lz)/global.time.dt + global.temp_array.pressure_derivative_1d;
            
            global.temp_array.in_helm_complex = global.temp_array.nlin1_1d/U.dissipation_coefficient; 
            
            if (Kperpsqr > MYEPS) {
                lambda = Kperpsqr + 1/(global.time.dt*U.dissipation_coefficient);
                Helmholtz_complex(global.temp_array.in_helm_complex, global.temp_array.out_helm_complex, lambda, 0, 0);	
                global.temp_array.V1_1d = global.temp_array.out_helm_complex;
            }
            else 
                global.temp_array.V1_1d = 0;
            
            // Vx_prime_plus,minus contain derivative at the boundary.
            Compute_derivative_at_boundary(global.temp_array.out_helm_complex, vx_prime_plus, vx_prime_minus);
            
            // complete p, v1 soln
            Compute_delta_pm(ly, lz, vx_prime_plus, vx_prime_minus, delta_plus, delta_minus);
            
            global.temp_array.pressure_1d = global.temp_array.pressure_1d + delta_plus*global.temp_array.pressure_plus(ly,Range::all(),lz) + delta_minus*global.temp_array.pressure_minus(ly,Range::all(),lz);
            
            global.temp_array.V1_1d = global.temp_array.V1_1d + delta_plus*global.temp_array.vx_plus(ly,Range::all(),lz) + delta_minus*global.temp_array.vx_minus(ly,Range::all(),lz);
            
            global.temp_array.pressureTR(ly,Range::all(),lz) = global.temp_array.pressure_1d;
            U.rvf.V1r(ly,Range::all(),lz) = global.temp_array.V1_1d;
            
            // V2 soln
            Yderiv1d_Chebyshev(ly, global.temp_array.pressure_1d, global.temp_array.pressure_derivative_1d); 
            
            global.temp_array.nlin2_1d = global.temp_array.nlin2TR(ly,Range::all(),lz) -U.rvf.V2r(ly,Range::all(),lz)/global.time.dt + global.temp_array.pressure_derivative_1d; 
            
            global.temp_array.in_helm_complex = global.temp_array.nlin2_1d/U.dissipation_coefficient;
            lambda = Kperpsqr + 1/(global.time.dt*U.dissipation_coefficient);
            Helmholtz_complex(global.temp_array.in_helm_complex, global.temp_array.out_helm_complex, lambda, 0, 0);	
            global.temp_array.V2_1d = global.temp_array.out_helm_complex;
            
             U.rvf.V2r(ly,Range::all(),lz) = global.temp_array.V2_1d;
            
            // V3 soln
            Zderiv1d_Chebyshev(lz, global.temp_array.pressure_1d, global.temp_array.pressure_derivative_1d); 
            
            global.temp_array.nlin3_1d = global.temp_array.nlin3TR(ly,Range::all(),lz) -U.rvf.V3r(ly,Range::all(),lz)/global.time.dt + global.temp_array.pressure_derivative_1d; 
            
            global.temp_array.in_helm_complex = global.temp_array.nlin3_1d/U.dissipation_coefficient;
            lambda = Kperpsqr + 1/(global.time.dt*U.dissipation_coefficient);
            Helmholtz_complex(global.temp_array.in_helm_complex, global.temp_array.out_helm_complex, lambda, 0, 0);	
            global.temp_array.V3_1d = global.temp_array.out_helm_complex;
            
             U.rvf.V3r(ly,Range::all(),lz) = global.temp_array.V3_1d; 
        }
    
	Transform::Transpose_array_SLAB(global.temp_array.pressureTR, P.F, 'Y', 'X', 'Z', 'Z');
    Transform::Transpose_array_SLAB(U.rvf.V1r, U.cvf.V1, 'Y', 'X', 'Z', 'Z');
    Transform::Transpose_array_SLAB(U.rvf.V2r, U.cvf.V2, 'Y', 'X', 'Z', 'Z');
    Transform::Transpose_array_SLAB(U.rvf.V3r, U.cvf.V3, 'Y', 'X', 'Z', 'Z'); 
 
 */
	
}


//**********************************   End of time_advance.cc  ********************************
   

