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
 * @date    22/08/2008
 * @bug     No known bugs
 */ 

#include "FFFW_slab.h"
#include "FFFW_slab_inline.h"




//*********************************************************************************************


/**********************************************************************************************
 
 Replaces A(k) by A(k)*exp(factor*K^2).
  
 ***********************************************************************************************/


void FFFW_SLAB::Array_exp_ksqr(Array<Complex,3> A, Real factor)
{
    if (global.program.kind != "KEPLERIAN") {
        Universal::Array_exp_ksqr(A, factor);
    }
    
    else { // Keplerian
        //TODO - Hybrid
        Real Kx, Ky, Kxsqr, Kysqr;
        Real Ksqr, mu_sqr;
        
        Real omega_keplerian = global.force.double_para(0);
        Real q_keplerian = global.force.double_para(1);
        
        Real q_omega_t = q_keplerian*omega_keplerian*global.time.keplerian;
        
        for (int lx = 0; lx < local_Nx; lx++) {
            Kx = Get_kx(lx)*kfactor[1];
            Kxsqr = pow2(Kx);
            
            for (int ly=0; ly<Ny; ly++) {
                Ky = Get_ky(ly)*kfactor[2];
                Kysqr = pow2(Ky);
                for (int lz=0; lz<=N[3]/2; lz++) {
                    Ksqr = Kxsqr + Kysqr + pow2(lz*kfactor[3]);
                    mu_sqr = Ksqr+pow2(q_omega_t*Ky)+2*q_omega_t*Kx*Ky;
                    A(lx, ly, lz) *= exp(factor*Ksqr);
                }
            }
        }
    }
}



/**********************************************************************************************
 
 Replaces A(k) by A(k)*[ exp( factor*K^2+hyper_factor*K^4 )]
 
 ***********************************************************************************************/

void FFFW_SLAB::Array_exp_ksqr(Array<Complex,3> A, Real factor, Real hyper_factor, int hyper_exponent)
{
    
    if (global.program.kind != "KEPLERIAN") {
        Universal::Array_exp_ksqr(A, factor, hyper_factor, hyper_exponent);
    }
    
    /*  else {
     Real Kx, Ky, Kxsqr, Kysqr;
     Real Ksqr, mu_sqr;
     Real Kpownm2;    // K^{q-2} where q = hyper_exponent
     
     Real omega_keplerian = global.force.double_para(0);
     Real q_keplerian = global.force.double_para(1);
     
     Real q_omega_t = q_keplerian*omega_keplerian*global.time.keplerian;
     
     for (int ly=0; ly<=global.field.maxly; ly++) {
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
 
 COmment
 
 ***********************************************************************************************/

void FFFW_SLAB::Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz)
{}

void FFFW_SLAB::Last_component(int kx, int ky, int kz, Complex &Vx, Complex &Vy, Complex &Vz)
{
    Real Kx = kx*kfactor[1];
    Real Ky = ky*kfactor[2];
    Real Kz = kz*kfactor[3];
    
    if (kz != 0) 
        Vz = (Complex(0,Kx)*Vx + Complex(0,Ky)*Vy)/Complex(0,-Kz);
        // works for both 2D (Ky=0) and 3D.
    
    else {
        if (ky != 0) { // 3D: input fields are (Vx, Vz); compute Vy
            Vz = Vy;
            Vy = (Complex(0,Kx)*Vx)/Complex(0,-Ky); 
        }
        else {  // k = (kx,0,0); input fields are (Vy, Vz)
            Vz = Vy;
            Vy = Vx;
            Vx = Complex(0,0);
        }
    }
    
}

/*****************************************************************************************
 
Dealias
 A(Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,2*Nx/3-1) = 0;
 
 *****************************************************************************************/

void FFFW_SLAB::Dealias(Array<Complex,3> A)
{   

    Array<int,1> Ax_filter(Nx);
    
    Ax_filter = 0;
    
    Ax_filter(Range(Nx/3+1,2*Nx/3-1)) = 1;
    
    int first_x = first(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
    int last_x = last(Ax_filter(Range(my_id*local_Nx,(my_id+1)*local_Nx-1)) == 1);
    
    if (Ny>1)
        A(Range(first_x,last_x), Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,toEnd)) = 0.0;
    else
        A(Range(first_x,last_x), 0, Range(Nz/3+1,toEnd)) = 0.0;
}


// Data resides till outer_radius in k-space
bool FFFW_SLAB::Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius)
{

    int kx_max = (int) ceil(outer_radius/kfactor[1]);
    int ky_max = (int) ceil(outer_radius/kfactor[2]);
    int kz_max = (int) ceil(outer_radius/kfactor[3]);
    
    if ((kx_max > Nx/3) || (ky_max > Ny/3) || (kz_max > Nz/3))
        return true;
    else
        return false;
}


/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
 
 // For kz =0 and kz = N[3]/2 planes
 // f(-ky, 0, -kx) = conj(f(ky,0,kx))- do for kz=N[3]/2
 // Implement:  f(ky, 0 or Nz/2, kx) = conj(f(-ky, 0 or Nz/2, -kx))
  // ky=Ny/2 : set A(kx,ky,Nz/2)=0; quite efficient without losing many data points
 
 ***********************************************************************************************/

void FFFW_SLAB::Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A)
{
    
    int array_index_minus_kx, array_index_minus_ky;
    
    // for kz=0
    Get_XY_plane(A, global.temp_array.plane_xy, 0);
    
    // For a given (minusky, minuskx), locate (ky,ky) and then subst.
    // A(minusky, 0, minuskx) = conj(A(ky,0,kx))
    for (int lx=0; lx<local_Nx; lx++) {
        if (Get_kx(lx) != Ny/2) {  // Do not apply for ky=Ny/2
            array_index_minus_kx =  Get_ix(-Get_kx(lx));  // minusky = -Get_ky(ly);
            
            for (int ly=Ny/2+1; ly<Ny; ly++) {
                array_index_minus_ky = -Get_kx(lx);        // kx<0; minuskx = -Get_kx(lx) > 0
                A(lx, ly, 0) = conj(global.temp_array.plane_xy(array_index_minus_kx,array_index_minus_ky));
            }
            // for (ky=0,kz=0) line, change only the A(-kx,0,0)
            if (Get_kx(lx) < 0)
                A(lx,0,0) = conj(global.temp_array.plane_xy(array_index_minus_kx,0));
        }
    }
        
    // for kz=Nz/2; 
    A(Range::all(),Range::all(),Nz/2) = 0.0;
}

void FFFW_SLAB::Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A)
{

    // for kz=Nz/2
    A(Range::all(),Range::all(),Nz/2) = 0.0;
}


//*********

void FFFW_SLAB::Test_reality_condition_in_Array(Array<Complex,3> A)
{
    //TODO - Hybrid
    int array_index_minus_kx, array_index_minus_ky;
    
    // for kz=0
    Get_XY_plane(A, global.temp_array.plane_xy, 0);
    
    // For a given (minusky, minuskx), locate (ky,ky) and then subst.
    // A(minusky, 0, minuskx) = conj(A(ky,0,kx))
    for (int lx=0; lx<local_Nx; lx++) {
        if (Get_kx(lx) != Nx/2) {  // Do not apply for ky=Ny/2
            array_index_minus_kx =  Get_ix(-Get_kx(lx));  // minusky = -Get_ky(ly);
            
            for (int ly=Ny/2+1; ly<Ny; ly++) {
                array_index_minus_ky = -Get_ky(ly);        // kx<0; minuskx = -Get_kx(lx) > 0
                
                if (abs(A(lx,ly,0)-conj(global.temp_array.plane_xy(array_index_minus_kx,array_index_minus_ky))) > MYEPS2)
                    cout << "Reality condition voilated for (kx,ky,kz)=(" <<array_index_minus_kx <<  "," << Get_kx(lx) << "," << 0 << ")" << endl;
            }
            // for (ky=0,kz=0) line, change only the A(-kx,0,0)
            if (Get_kx(lx) < 0) {
                if (abs(A(lx,0,0)-conj(global.temp_array.plane_xy(array_index_minus_kx,array_index_minus_ky))) > MYEPS2)
                    cout << "Reality condition voilated for (kx,ky,kz)=(" << 0 <<  "," << Get_kx(lx) << "," << 0 << ")" << endl;
            }
        }
    }
    
    // for kz=Nz/2 plane
    for (int lx=0; lx<local_Nx; lx++)
        for (int ly=0; ly<Ny; ly++)
            if (abs(A(lx,ly,Nz/2)) > MYEPS)
                cout << "Reality condition voilated for (kx,ky,kz)=(" << Get_kx(lx) <<  "," << Get_ky(ly) << "," << Nz/2 << ")" << endl;
}


//********************************************************************************************* 


/** @brief Compute divergence of VF \f$ \vec{V} \f$.
 *
 *  @note   Divergence is stored in div
 *
 *  @return  \f$ *F = \mathcal{F}(D_i A_i) \f$. 
 */
void FFFW_SLAB::Compute_divergence(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, Array<Complex,3> div, string field_or_nlin, Real &total_abs_div, bool print_switch)
{
    Xderiv(Ax, div);            

    Add_Yderiv(Ay, div);
    
    Add_Zderiv(Az, div);
    
    
    if (global.program.kind == "KEPLERIAN") {
        Real omega_keplerian = global.force.double_para(0);
        Real q_keplerian = global.force.double_para(1);
        Real q_omega_t;
        
        if (field_or_nlin == "nlin") // for pressure computation
            q_omega_t = q_keplerian*omega_keplerian*global.time.keplerian;
        else if (field_or_nlin == "field") // for div(field computations
            q_omega_t = q_keplerian*omega_keplerian*global.time.now;
        
        Yderiv(Ax, global.temp_array.X);
        div = div + Complex(q_omega_t,0)*global.temp_array.X;
    }
    
    if (field_or_nlin == "field") {
        total_abs_div = sqrt(Array_sqr(div));
        
        if ((print_switch) && (total_abs_div > MYEPS2)) {
            cout << "NON-ZERO DIVERGENCE for the following modes:"  << endl;
            Print_large_Fourier_elements(div);
        }
    }
}

void FFFW_SLAB::Get_XY_plane(Array<Complex,3> A, Array<Complex,2> plane_xy, int kz)
{
    if (numprocs == 1) {
        plane_xy= A(Range::all(), Range::all(), kz);
        return;
    }
    
    else {
        int data_size, full_data_size;

        global.temp_array.plane_xy_inproc = A(Range::all(),Range::all(),kz);
    
        data_size = 2*local_Nx*Ny;
        
        full_data_size = 2*Nx*Ny;
        
        MPI_Gather(reinterpret_cast<Real*>(global.temp_array.plane_xy_inproc.data()), data_size, MPI_Real, reinterpret_cast<Real*>(plane_xy.data()), data_size, MPI_Real, 0, MPI_COMM_WORLD);
        
        MPI_Bcast(reinterpret_cast<Real*>(global.temp_array.plane_xy.data()), full_data_size, MPI_Real, 0, MPI_COMM_WORLD);
    } 
}


int FFFW_SLAB::Read(Array<Complex,3> A, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
    return BasicIO::Read(A.data(), plan, file_name, dataset_name);
}

int FFFW_SLAB::Read(Array<Real,3> Ar, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
    return BasicIO::Read(Ar.data(), plan, file_name, dataset_name);
}


int FFFW_SLAB::Write(Array<Complex,3> A, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
    return BasicIO::Write(A.data(), plan, folder_name, file_name, dataset_name);
}

int FFFW_SLAB::Write(Array<Real,3> Ar, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
    return BasicIO::Write(Ar.data(), plan, folder_name, file_name, dataset_name);  
}
//*****************************  End of four_basic.cc **************************








