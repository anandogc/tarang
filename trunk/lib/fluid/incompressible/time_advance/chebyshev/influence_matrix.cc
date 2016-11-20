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



//****************************************************************************************

void Print_large_Fourier_elements_real(Array<Real,3> A)
{
    for (int ly=0; ly<local_Ny; ly++)
        for (int lz=0; lz<local_Nz; lz++)
            for (int lx=0; lx<Nx; lx++){
				if (abs(A(ly, lz, lx)) > MYEPS2) {
					cout << "my_id = " << my_id <<  " vect(k) = ( " << universal->Get_ky(ly) << "," << universal->Get_kz(lz) << "," << universal->Get_kx(lx)  <<");  A(k) = " << A(ly,lz,lx) << '\n';
				}
			}
    cout << endl;
}

void Time_advance_incompress::Compute_homgeneous_soln_influence_matrix(FluidVF& U, Pressure& P)
{
	Real lambda_supplement;
    
	global.temp_array.influence_matrix = 0.0;
	
    // Compute p+
    global.temp_array.Xreal = 0.0;
    Helmholtz_real_full_array(global.temp_array.Xreal, global.temp_array.pressure_plus, true, 0, 1.0, 0);
    
//	cout << "P+ = " << global.temp_array.pressure_plus.shape() << endl;  Print_large_Fourier_elements_real(global.temp_array.pressure_plus);
	
	
    // compute u+
    universal->Xderiv(global.temp_array.pressure_plus, global.temp_array.Xreal);
    
    global.temp_array.Xreal = global.temp_array.Xreal/(-U.dissipation_coefficient);
    
    lambda_supplement = 1/(global.time.dt*U.dissipation_coefficient);
    Helmholtz_real_full_array(global.temp_array.Xreal,  global.temp_array.vx_plus, false, lambda_supplement, 0, 0);
       
    // Compute p-
    global.temp_array.Xreal = 0.0;
    Helmholtz_real_full_array(global.temp_array.Xreal, global.temp_array.pressure_minus, true, 0, 0, 1.0);
	// p- = global.temp_array.pressure_minus
    
    // compute u-
    universal->Xderiv(global.temp_array.pressure_minus, global.temp_array.Xreal);
    
    global.temp_array.Xreal = global.temp_array.Xreal/(-U.dissipation_coefficient);
    
    lambda_supplement = 1/(global.time.dt*U.dissipation_coefficient);
    

    Helmholtz_real_full_array(global.temp_array.Xreal,  global.temp_array.vx_minus, false, lambda_supplement,0, 0);
    
    //Compute influence matrix
    global.temp_array.influence_matrix = 0.0;

    
    Real factor1, factor2;
    
    // vx_plus = v
    // T'_k(1) = n^2;  T'_k(-1) = (-1)^k n^2
   for (int kx=0; kx<=(Nx-2); kx=kx+2) {
        factor1 = kx*kx;
        factor2 = (kx+1)*(kx+1);
        
        global.temp_array.influence_matrix(Range::all(),Range::all(),0,0) += factor1*global.temp_array.vx_plus(Range::all(),Range::all(),kx) + factor2*global.temp_array.vx_plus(Range::all(),Range::all(),kx+1);
        
        global.temp_array.influence_matrix(Range::all(),Range::all(),0,1) += factor1*global.temp_array.vx_minus(Range::all(),Range::all(),kx) + factor2*global.temp_array.vx_minus(Range::all(),Range::all(),kx+1);
        
        global.temp_array.influence_matrix(Range::all(),Range::all(),1,0) += (-factor1)*global.temp_array.vx_plus(Range::all(),Range::all(),kx) + factor2*global.temp_array.vx_plus(Range::all(),Range::all(),kx+1);
        
        
        global.temp_array.influence_matrix(Range::all(),Range::all(),1,1) += (-factor1)*global.temp_array.vx_minus(Range::all(),Range::all(),kx) + factor2*global.temp_array.vx_minus(Range::all(),Range::all(),kx+1);
    } 
    
    global.temp_array.influence_matrix *= kfactor[1];
	
/*	for (int lz=0; lz<local_Nz; lz++)
		for(int ly=0; ly<local_Ny; ly++)
			cout << "lz,ly,influence matrix = " << lz << " " << ly << global.temp_array.influence_matrix(ly,lz,Range::all(),Range::all()) << endl; */
	
/*	cout << "vx+(0,0), vx-(0,0) " << global.temp_array.vx_plus(0,0,Range::all()) << global.temp_array.vx_minus(0,0,Range::all()) << endl;
	
	Real v_prime_plus, v_prime_minus;
	
	Compute_derivative_at_boundary(global.temp_array.vx_plus(0,0,Range::all()) , v_prime_plus, v_prime_minus);
	cout << "vx+(+1), vx+(-1) = " << v_prime_plus << " " << v_prime_minus << endl;
	
	Compute_derivative_at_boundary(global.temp_array.vx_minus(0,0,Range::all()) , v_prime_plus, v_prime_minus);
	cout << "vx-(+1), vx-(-1) = " << v_prime_plus << " " << v_prime_minus << endl; */
}

//******************************************************************************
void Time_advance_incompress::Compute_pressure_Ux_chebyshev(FluidVF& U, Pressure& P)
{
    Real Vplus_top, Vminus_top, Vplus_bottom, Vminus_bottom;
    Complex delta_plus,delta_minus;
    Complex vx_prime_plus,vx_prime_minus;
	Real determinant;
    
    for(int ly=0; ly<local_Ny; ly++)
		for (int lz=0; lz<local_Nz; lz++) {
            Vplus_top = global.temp_array.influence_matrix(ly,lz,0,0);
            Vminus_top = global.temp_array.influence_matrix(ly,lz,0,1);
            Vplus_bottom = global.temp_array.influence_matrix(ly,lz,1,0);
            Vminus_bottom = global.temp_array.influence_matrix(ly,lz,1,1);
            
            determinant = Vplus_top*Vminus_bottom - Vminus_top*Vplus_bottom;
            
            global.temp_array.V1_x = U.cvf.V1(ly,lz,Range::all());
            
            Compute_derivative_at_boundary(global.temp_array.V1_x, vx_prime_plus, vx_prime_minus);
            
            if (abs(determinant) > MYEPS) {
                delta_plus  = (Vminus_top*vx_prime_minus - Vminus_bottom*vx_prime_plus)/determinant;
                
                delta_minus = (Vplus_bottom*vx_prime_plus - Vplus_top*vx_prime_minus)/determinant;
            }
            else {
             //   cout << "Division by 0 in Compute_delta_pm() for lz, ly = " << lz << " " << ly << endl;
              //  cout << "lz, ly, V's, det = " << lz << " " << ly << " " << Vplus_top << " " << Vminus_top << " "  << Vplus_bottom << " " << Vminus_bottom << " " << determinant << endl;
            }
            
            // Pl check...
            P.F(ly,lz,Range::all()) +=  (delta_plus*global.temp_array.pressure_plus(ly,lz,Range::all()) + delta_minus*global.temp_array.pressure_minus(ly,lz,Range::all()));
            
            U.cvf.V1(ly,lz,Range::all()) += (delta_plus*global.temp_array.vx_plus(ly,lz,Range::all()) + delta_minus*global.temp_array.vx_minus(ly,lz,Range::all()));
        }
	
	for (int lx=0; lx<(U.cvf.V1).extent(2); lx++) {
		universal->Assign_spectral_field(lx,0,0, P.F,Complex(0,0));
		universal->Assign_spectral_field(lx,0,0, U.cvf.V1,Complex(0,0));
	}
}


//****************************************************************************************

void Time_advance_incompress::Compute_derivative_at_boundary(Array<Complex,1> v,Complex &v_prime_plus, Complex &v_prime_minus)
{
    v_prime_plus = 0.0;
    v_prime_minus = 0.0;
    
    Real factor1, factor2;
    
    // Nx even
    for (int kx=0; kx<=(Nx-2); kx=kx+2) {
        factor1 = kx*kx;
        factor2 = (kx+1)*(kx+1);
        
        v_prime_plus += (factor1*v(kx) + factor2*v(kx+1));
        v_prime_minus += (-factor1*v(kx) + factor2*v(kx+1));
    }
    
    v_prime_plus = kfactor[1]*v_prime_plus;
    v_prime_minus = kfactor[1]*v_prime_minus;
    
}

 //****************************************************************************************


 void Time_advance_incompress::Compute_derivative_at_boundary(Array<Real,1> v, Real &v_prime_plus, Real &v_prime_minus)
 {
     v_prime_plus = 0.0;
     v_prime_minus = 0.0;
     
     Real factor1, factor2;
 
     // Nx even
     for (int kx=0; kx<=(Nx-2); kx=kx+2) {
         factor1 = kx*kx;
         factor2 = (kx+1)*(kx+1);
         
         v_prime_plus += (factor1*v(kx) + factor2*v(kx+1));
         v_prime_minus += (-factor1*v(kx) + factor2*v(kx+1));
     }
     
     v_prime_plus = kfactor[1]*v_prime_plus;
     v_prime_minus = kfactor[1]*v_prime_minus;
 }

     



//**********************************   End of time_advance.cc  ********************************
   

