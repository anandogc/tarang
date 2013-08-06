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

/*! \file four_basic.cc
 * 
 * @sa four_basic.h
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date	22/08/2008
 * @bug		No known bugs
 */ 

#include "CFFF_slab.h"
#include "CFFF_slab_inline.h"


/**********************************************************************************************

	-- Number of modes included inside shell (inner_radius, outer_radius]
	-- also count the complex conj
	
***********************************************************************************************/


int  CFFF_SLAB::Get_number_modes_in_shell(DP inner_radius, DP outer_radius)
{
    
/*	int kx_max, ky_max, kz_max;
	
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (N[2] > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0.0;
    
    kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	int count = 0;
	DP  Kmag;
	
	for (int kx = -kx_max; kx <= kx_max; kx++)
		for (int ky = -ky_max; ky <= ky_max; ky++)  
			for (int kz = 0; kz <= kz_max; kz++) {
                Kmag = sqrt(pow2(kx*kfactor[1]) + pow2(ky*kfactor[2]) + pow2(kz*kfactor[3]));
				
				if (( Kmag > inner_radius) && ( Kmag <= outer_radius))	{
					if (kz == 0)
						count++;
					else 
						count = count + 2;
				}		
			}
	
	return count;	*/

}

//*********************************************************************************************

void  CFFF_SLAB::Print_large_Fourier_elements(Array<complx,3> A)
{
	for (int lx=0; lx<local_Nx; lx++) 
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<Nz; lz++) {
				if (abs(A(lx, ly, lz)) > MYEPS2) {
					cout << "my_id = " << my_id <<  " vect(k) = ( " << Get_kx(lx) << "," << Get_ky(ly) << "," << Get_kz(lz) <<");  Array(k) = " << A(lx, ly, lz) << '\n';
				}
			}
    cout << endl;
}

/**********************************************************************************************

       		Replaces A(k) by A(k)*K^2.

***********************************************************************************************/


void  CFFF_SLAB::Array_mult_ksqr(Array<complx,3> A)
{
/*	firstIndex ly;
	secondIndex lz;
  
	DP Kxsqr;    // Kx^2
	
	for (int lx = 0; lx < local_Nx; lx++) 	{
		Kxsqr = pow2(Get_kx(lx)*kfactor[1]);	
			
		A(lx, Range(0,N[2]/2), Range::all()) 
        = A(lx, Range(0,N[2]/2), Range::all())* ( Kxsqr + pow2(ly*kfactor[2]) + pow2(lz*kfactor[3]));
  
		if (N[2] > 1)
			A(lx, Range(N[2]/2+1,N[2]-1), Range::all())  
            = A(lx, Range(N[2]/2+1,N[2]-1), Range::all())* ( Kxsqr + pow2((ly+1-N[2]/2)*kfactor[2])+ pow2(lz*kfactor[3]));
	}  */
}


/**********************************************************************************************

	Replaces A(k) by A(k)/K^2 with K^2 = sum_i [ki*kfactor(i)]^2 ; A(0)=0

***********************************************************************************************/



void  CFFF_SLAB::Array_divide_ksqr(Array<complx,3> A)
{
/*	firstIndex ly;
	secondIndex lz;
  
	DP Kxsqr;
	
	for (int lx = 0; lx < local_Nx; lx++) {
		Kxsqr = pow2(Get_kx(lx)*kfactor[1]);	
		
		A(lx, Range(0,N[2]/2), Range::all()) 
        = A(lx, Range(0,N[2]/2), Range::all()) /( Kxsqr + pow2(ly*kfactor[2])+pow2(lz*kfactor[3]));
  
		if (N[2] > 1)
			A(lx, Range(N[2]/2+1,N[2]-1), Range::all())  
            = A(lx, Range(N[2]/2+1,N[2]-1), Range::all()) /( Kxsqr + pow2((ly+1-N[2]/2)*kfactor[2]) + pow2(lz*kfactor[3]));
	}
	   
	// To avoid division by zero
	if (my_id == master_id)   
		A(0,0,0) = 0.0;   */
}


/**********************************************************************************************

	Replaces A(k) by A(k)*exp(factor*K^2).
 
    In the solver, factor = a

***********************************************************************************************/


void  CFFF_SLAB::Array_exp_ksqr(Array<complx,3> A, DP factor)
{	
    DP Kxsqr, Kysqr;
    DP Ksqr;
    
    DP dt = global.time.dt;
    DP mu = 0;
    DP kappa = global.field.diss_coefficients[0];
    DP sinx,cosx,re,im;
    DP diss_factor;
    
 //   cout << "asqr, before exp = " << sum(sqr(abs(A))) << endl;
    
    for (int lx = 0; lx < local_Nx; lx++) {
        Kxsqr = pow2(Get_kx(lx)*kfactor[1]);	
        for (int ly=0; ly<Ny; ly++) {
            Kysqr = pow2(Get_ky(ly)*kfactor[2]);
            for (int lz=0; lz<Nz; lz++) {
                Ksqr = Kxsqr + Kysqr + pow2(Get_kz(lz)*kfactor[3]);
                
                diss_factor = exp(-factor*kappa*dt*Ksqr*Ksqr);
                cosx = cos(factor*(Ksqr/2-mu)*dt);
                sinx = sin(factor*(Ksqr/2-mu)*dt);
                re = real(A(lx, ly, lz));
                im = imag(A(lx, ly, lz));
                real(A(lx, ly, lz)) = (re*cosx - im*sinx)*diss_factor;
                imag(A(lx, ly, lz)) = (im*cosx + re*sinx)*diss_factor;
            }
        }
    }
    
  //   cout << "asqr, after exp = " << sum(sqr(abs(A))) << endl;
}



/**********************************************************************************************

	Replaces A(k) by A(k)*[ exp( factor*K^2+hyper_factor*K^4 )]

***********************************************************************************************/

void  CFFF_SLAB::Array_exp_ksqr(Array<complx,3> A, DP factor, DP hyper_factor, int hyper_exponent)
{	
/*	if (global.program.kind != "KEPLERIAN") {
        DP Kxsqr, Kysqr;
        DP Ksqr;
        DP Kpownm2;	// K^{q-2} where q = hyper_exponent
        
        
        for (int lx = 0; lx < local_Nx; lx++) {
            Kxsqr = pow2(Get_kx(lx)*kfactor[1]);
            
            for (int ly=0; ly<=N[2]/2; ly++) {
                Kysqr = pow2(ly*kfactor[2]);
                for (int lz=0; lz<=N[3]/2; lz++) {
                    Ksqr = Kxsqr + Kysqr + pow2(lz*kfactor[3]);
                    
                    if (hyper_exponent == 4) 
                        Kpownm2 = Ksqr;
                    else
                        Kpownm2 = my_pow(Ksqr,(hyper_exponent-2)/2);
                        
                    A(lx, ly, lz) *= exp((factor+hyper_factor*Kpownm2)* Ksqr); 
                }
            }
            
            if (N[2] >1)
                for (int ly=N[2]/2+1; ly<N[2]; ly++) {
                    Kysqr = pow2((ly-N[2])*kfactor[2]);
                    for (int lz=0; lz<=N[3]/2; lz++) {
                        Ksqr = Kxsqr + Kysqr + pow2(lz*kfactor[3]);
                        
                        if (hyper_exponent == 4) 
                            Kpownm2 = Ksqr;
                        else
                            Kpownm2 = my_pow(Ksqr,(hyper_exponent-2)/2);
                        
                        A(lx, ly, lz) *= exp((factor+hyper_factor*Kpownm2)* Ksqr);
                    }
                }
        }
    }
    
    else {
        DP Kx, Ky, Kxsqr, Kysqr;
        DP Ksqr, mu_sqr;
        DP Kpownm2;	// K^{q-2} where q = hyper_exponent
        
        DP omega_keplerian = global.force.double_para(0);
        DP q_keplerian = global.force.double_para(1);
        
        DP q_omega_t = q_keplerian*omega_keplerian*global.time.keplerian;
        
        for (int lx = 0; lx < local_Nx; lx++) {
            Kx = Get_kx(lx)*kfactor[1];
            Kxsqr = pow2(Kx);
            
            for (int ly=0; ly<Ny; ly++) {
                Ky = Get_ky(ly)*kfactor[2];
                Kysqr = pow2(Ky);
                for (int lz=0; lz<=N[3]/2; lz++) {
                    Ksqr = Kxsqr + Kysqr + pow2(lz*kfactor[3]);
                    mu_sqr = Ksqr+pow2(q_omega_t*Ky)+2*q_omega_t*Kx*Ky;
                    
                    if (hyper_exponent == 4)
                        Kpownm2 = mu_sqr;
                    else
                        Kpownm2 = my_pow(mu_sqr,(hyper_exponent-2)/2);
                    
                    A(lx, ly, lz) *= exp((factor+hyper_factor*Kpownm2)* mu_sqr);
                }
            }
        }
    } */
	
}
 


/**********************************************************************************************

	Replaces A(k) by A(k)*(V0.K)^2 / K^2.

***********************************************************************************************/


void  CFFF_SLAB::Array_mult_V0_khat_sqr(Array<complx,3> A, TinyVector<DP,3> V0)
{
/*	firstIndex ly;
	secondIndex lz;
  
	DP V0x = V0(0);
	DP V0y = V0(1);
	DP V0z = V0(2);
	
	DP Kx;
	
	for (int lx = 0; lx < local_Nx; lx++) {
		Kx = Get_kx(lx)*kfactor[1];	
		A(lx, Range(0,N[2]/2), Range::all()) 
        = A(lx, Range(0,N[2]/2), Range::all())* sqr(V0x*Kx + V0y*ly*kfactor[2] + V0z*lz*kfactor[3]);
			
		if (N[2] > 1)
			A(lx, Range(N[2]/2+1,N[2]-1), Range::all())  
            = A(lx, Range(N[2]/2+1,N[2]-1), Range::all())* sqr(V0x*Kx + V0y*(ly+1-N[2]/2)*kfactor[2] + V0z*lz*kfactor[3]);
	}  
	
	Array_divide_ksqr(A);*/
}

/**********************************************************************************************
 
 COmment
 
 ***********************************************************************************************/

void  CFFF_SLAB::Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz)
{}

void  CFFF_SLAB::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
{
/*	DP Kx = kx*kfactor[1];
	DP Ky = ky*kfactor[2];
	DP Kz = kz*kfactor[3];
	
	if (kz != 0) 
		Vz = (complx(0,Kx)*Vx + complx(0,Ky)*Vy)/complx(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) { // 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = (complx(0,Kx)*Vx)/complx(0,-Ky); 
		}
		else {	// k = (kx,0,0); input fields are (Vy, Vz)
			Vz = Vy;
			Vy = Vx;
			Vx = complx(0,0);
		}
	} */
	
}

/*****************************************************************************************
 
Dealias
 step 1: A(Range(Nx/3+1,2*Nx/3-1),:,:) = 0;
 step 2: A(Range(0,Nx/3),Range(Ny/3+1,2*Ny/3-1),:) = 0;
         A(Range(2*Nx/3,toEnd),Range(Ny/3+1,2*Ny/3-1),:) = 0;
 
 step 3: A(Range(0,Nx/3),Range(0,Ny/3),Range(Nz/3+1,Nz/2)) = 0;
         A(Range(0,Nx/3),Range(2Ny/3,end),Range(Nz/3+1,Nz/2)) = 0;
 
         A(Range(2*Nx/3,toEnd),Range(0,Ny/3),Range(Nz/3+1,Nz/2)) = 0;
         A(Range(2*Nx/3,toEnd),Range(2Ny/3,end),Range(Nz/3+1,Nz/2)) = 0;
 
 *****************************************************************************************/

void  CFFF_SLAB::Dealias(Array<complx,3> A)
{	
/*	int index_Nxby3 = Get_lx(Nx/3);     // Only when kx=Nx/3 in the proc
    int index_2Nxby3 = Get_lx(2*Nx/3);  // Only when kx=2*Nx/3 in the proc
    
    int lx_min_for_zeroing,lx_max_for_zeroing;
	
	// Array index of the full A
	int local_ix_min = local_Nx_start;
	int local_ix_max = local_Nx_start+local_Nx-1;
    
    Range lx_range;
	
	if (Nx >= 3) {
		if (local_ix_min <= Nx/3) {
            if (local_ix_max <= Nx/3)
                ; // Not to be zeroed; Do nothing
            else if (local_ix_max < 2*Nx/3) // kx=Nx/3 but not 2*Nx/3 inside proc
                A(Range(index_Nxby3+1,toEnd),Range::all(),Range::all()) = 0.0;
            else // Both kx=Nx/3 but not 2*Nx/3 inside proc
                A(Range(index_Nxby3+1,index_2Nxby3-1),Range::all(),Range::all()) = 0.0;
        }
        
        else if (local_ix_min < 2*Nx/3) {
            if (local_ix_max < 2*Nx/3)
                A(Range::all(),Range::all(),Range::all()) = 0.0;
            else // kx=2*Nx/3 inside proc
                A(Range(0,index_2Nxby3-1),Range::all(),Range::all()) = 0.0;
        }
    }
          
        
    if (Ny >= 3) {
        if (local_ix_min <= Nx/3) {
            if (local_ix_max <= Nx/3)  
                A(Range(0,toEnd), Range(Ny/3+1,2*Ny/3-1), Range::all())= 0.0;
            else  // kx=Nx/3  inside proc
                A(Range(0,index_Nxby3), Range(Ny/3+1,2*Ny/3-1), Range::all())= 0.0;
        }
        
        if (local_ix_max >= 2*Nx/3) {
            if (local_ix_min < 2*Nx/3)  // kx=2*Nx/3 inside proc
                A(Range(index_2Nxby3,toEnd), Range(Ny/3+1,2*Ny/3-1), Range::all())= 0.0;
            else    
                A(Range(0,toEnd), Range(Ny/3+1,2*Ny/3-1), Range::all())= 0.0;
        }
    }
    
    
    if (Nz >= 3) {
        if (local_ix_min <= Nx/3) {
            if (local_ix_max <= Nx/3)
                lx_range = Range(0,local_Nx-1);
            else    // kx=Nx/3 inside the proc
                lx_range = Range(0,index_Nxby3);
            
            if (Ny==1)
                A(lx_range,0,Range(Nz/3+1,Nz/2))= 0.0;
            else {
                A(lx_range, Range(0,Ny/3), Range(Nz/3+1,Nz/2))= 0.0;
                A(lx_range, Range(2*Ny/3,toEnd), Range(Nz/3+1,Nz/2))= 0.0;
            }
        }
        
        if (local_ix_max >= 2*Nx/3) {
            if (local_ix_min < 2*Nx/3) // kx=2*Nx/3 inside the proc (assumes Nx>=2)
                lx_range = Range(index_2Nxby3,local_Nx-1);
            else
                lx_range = Range(0,local_Nx-1);
            
            if (Ny==1)
                A(lx_range,0,Range(Nz/3+1,Nz/2))= 0.0;
            else {
                A(lx_range, Range(0,Ny/3), Range(Nz/3+1,Nz/2))= 0.0;
                A(lx_range, Range(2*Ny/3,toEnd), Range(Nz/3+1,Nz/2))= 0.0;
            }
        }
    }*/
}


// Data resides till outer_radius in k-space
bool  CFFF_SLAB::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
{
/*	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > Nx/3) || (ky_max > Ny/3) || (kz_max > Nz/3))
		return true;
	else
		return false;
 */
}


/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
 
 // For kz =0 and kz = N[3]/2 planes
 // f(-kx, -ky, 0) = conj(f(kx, ky, 0))- do for kz=N[3]/2
 // Implement:  f(kx, ky, 0 or Nz/2) = conj(f(-kx, -ky, 0 or Nz/2))
  // ky=Ny/2 : set A(kx,ky,Nz/2)=0; quite efficient without losing many data points
 
 ***********************************************************************************************/

void  CFFF_SLAB::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
}

void  CFFF_SLAB::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{
}


//*********

void  CFFF_SLAB::Test_reality_condition_in_Array(Array<complx,3> A)
{
}


//********************************************************************************************* 


/** @brief Compute divergence of VF \f$ \vec{V} \f$.
 *
 *	@note	Divergence is stored in div
 *
 *  @return  \f$ *F = \mathcal{F}(D_i A_i) \f$. 
 */
void  CFFF_SLAB::Compute_divergence(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, Array<complx,3> div, string field_or_nlin, DP &total_abs_div, bool print_switch)
{
    Xderiv(Ax, div);			

	Yderiv(Ay, global.temp_array.X);
	div = div + global.temp_array.X;
	
	Zderiv(Az, global.temp_array.X);
	div = div + global.temp_array.X;
    
//    Print_large_Fourier_elements
	
  /*  if (field_or_nlin == "field") {
        total_abs_div = sqrt(Array_sqr(div));
        
        if ((print_switch) && (total_abs_div > MYEPS2)) {
            cout << "NON-ZERO DIVERGENCE for the following modes:"  << endl;
            Print_large_Fourier_elements(div);
        }
    } */
}


//********************************************************************************************* 


// kz=0 is already filled.

void  CFFF_SLAB::Fill_Vz(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
/*    if (global.io.input_vx_vy_switch && global.field.incompressible) {
		int kx, ky, kz;
		complx vz;
        
		for (int lx=0; lx<local_Nx; lx++)
			for (int ly=0; ly<Ny; ly++)
				for (int lz=1; lz<=(Nz/2); lz++) {
					kx = Get_kx(lx);
					ky = Get_ky(ly);
					kz = lz;
                    
					Last_component(kx, ky, kz, Ax(lx,ly,lz), Ay(lx,ly,lz), vz);
                    
					Az(lx,ly,lz) = vz;
				}
	} */
}

void  CFFF_SLAB::Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
    // Do nothing
}

void  CFFF_SLAB::Zero_modes(Array<complx,3> F)
{
    // Do nothing
}

//*********************************************************************************************

void  CFFF_SLAB::Copy_array_complex_to_real(Array<complx,3> A, Array<DP,3> Acopy)
{
     // Dummy function: Do nothing... Note that Acopy is not allocated.
}

//*********************************************************************************************

void  CFFF_SLAB::Copy_array_real_to_complex(Array<DP,3> Acopy, Array<complx,3> A)
{
     // Dummy function: Do nothing... Note that Acopy is not allocated.
}


//******************************************************************************

void  CFFF_SLAB::Copy_realarray_complex_to_real(Array<complx,3> Ar, Array<DP,3> Ar_copy)
{
     // Dummy function: Do nothing... Note that Ar_copy is not allocated.
}

//******************************************************************************

void  CFFF_SLAB::Copy_realarray_real_to_complex(Array<DP,3> Ar_copy, Array<complx,3> Ar)
{
     // Dummy function: Do nothing... Note that Ar_copy is not allocated.
}

//******************************************************************************

// A(Nx-global.io.N_in_reduced[1]/2,:,:) contains kx=Nx_red/2 data.  Must be set to zero.
// Same for Ny
void  CFFF_SLAB::Adjust_array_after_reading_reduced(Array<complx,3> A)
{
  /*  int lx=Get_lx(-global.io.N_in_reduced[1]/2);
    
    if ((lx >= 0) && (lx < local_Nx))
        A(lx,Range::all(),Range::all()) = 0.0;
    
    if (Ny > 1)
		A(Range::all(),Ny-global.io.N_in_reduced[2]/2,Range::all()) = 0.0;
   */
}

//*****************************  End of four_basic.cc **************************








