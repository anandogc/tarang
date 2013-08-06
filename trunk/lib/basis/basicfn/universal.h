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


/*! \file  universal.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  A. G. Chatterjee
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug	No known bugs
 */


#ifndef _UNIVERAL_H
#define _UNIVERAL_H

#include "def_vars.h"

//*********************************************************************************************	


class Universal
{				
public:
    
    // basic
	virtual int Get_number_modes_in_shell(DP inner_radius, DP outer_radius);
	virtual void Print_large_Fourier_elements(Array<complx,3> A);
	virtual void Array_mult_ksqr(Array<complx,3> A);
	virtual void Array_divide_ksqr(Array<complx,3> A);
	virtual void Array_exp_ksqr(Array<complx,3> A, DP factor);
	virtual void Array_exp_ksqr(Array<complx,3> A, DP factor, DP hyper_factor, int hyper_exponent);
	virtual void Array_mult_V0_khat_sqr(Array<complx,3> A, TinyVector<DP,3> V0);
    virtual void Fill_Vz(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
    virtual void Compute_divergence(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, Array<complx,3> div, string field_or_nlin, DP &total_abs_div, bool print_switch);
	virtual void Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
    virtual void Zero_modes(Array<complx,3> F);
    

	// Local functions...
    virtual void Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz) = 0;
    virtual void Last_component(int kx, int ky, int kz, complx& Vx, complx& Vy, complx& Vz) = 0;
    virtual void Dealias(Array<complx,3> A) = 0;
    virtual bool Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius) = 0;
    virtual void Satisfy_strong_reality_condition_in_Array(Array<complx,3> A) = 0;
    virtual void Satisfy_weak_reality_condition_in_Array(Array<complx,3> A) = 0;
    virtual void Test_reality_condition_in_Array(Array<complx,3> A) = 0;

    
    // transform
    
	virtual void Forward_transform(Array<DP,3> Ar, Array<complx,3> A) = 0;
	virtual void Inverse_transform(Array<complx,3> A, Array<DP,3> Ar) = 0;
    
	virtual void Xderiv(Array<complx,3> A, Array<complx,3> B) = 0;
    virtual void Yderiv(Array<complx,3> A, Array<complx,3> B) = 0;
    virtual void Zderiv(Array<complx,3> A, Array<complx,3> B) = 0;

    virtual void Add_Xderiv(Array<complx,3> A, Array<complx,3> B) = 0;
	virtual void Add_Yderiv(Array<complx,3> A, Array<complx,3> B) = 0;
	virtual void Add_Zderiv(Array<complx,3> A, Array<complx,3> B) = 0;
	
	virtual void Xderiv(Array<DP,3> A, Array<DP,3> B) = 0;
	
	virtual void Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B) = 0;
	virtual void Subtract_Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B) = 0;
	// Energy
    
    //virtual DP Get_local_energy_XZ_plane(Array<complx,3> A, int ny) = 0;
    
 //    virtual DP Get_local_energy_XZ_plane
	// (
	// 	Array<complx,3> A, Array<complx,3> B, 
	// 	int ny
	// ) = 0;
	
	virtual DP Get_local_energy_real_space(Array<DP,3> Ar) = 0;
	virtual DP Get_local_energy(Array<complx,3> A) = 0;
	virtual DP Get_local_energy_real_space(Array<DP,3> Ar, Array<DP,3> Br) = 0;
	virtual DP Get_local_energy(Array<complx,3> A, Array<complx,3> B) = 0;
	
	virtual void Compute_local_helicity
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 DP &local_helicity1, DP &local_helicity2,
	 DP &local_dissipation_H1, DP &local_dissipation_H2
	 ) = 0;
    
	virtual void Compute_total_helicity
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 DP &total_helicity1, DP &total_helicity2,
	 DP &total_dissipation_H1, DP &total_dissipation_H2
     ) = 0;
    
	virtual void Compute_local_shell_spectrum_helicity
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3,
	 Array<DP,1> local_H1k_count
	 ) = 0;
    
	virtual void Compute_shell_spectrum_helicity
	(
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3
	 ) = 0;
	
	virtual void Compute_local_ring_spectrum_helicity
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 Array<DP,2> local_H1k1, Array<DP,2> local_H1k2, Array<DP,2> local_H1k3
	 ) = 0;
	
	
	virtual void Compute_ring_spectrum_helicity
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 Array<DP,2> local_H1k1, Array<DP,2> local_H1k2, Array<DP,2> local_H1k3
	 ) = 0;
	
	virtual void Compute_local_cylindrical_ring_spectrum_helicity
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 Array<DP,2> local_H1k1, Array<DP,2> local_H1k2
	 ) = 0;
    
	virtual void Compute_cylindrical_ring_spectrum_helicity
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
     Array<DP,2> local_H1k1, Array<DP,2> local_H1k2
	 ) = 0;



	
	virtual DP Get_total_energy_real_space(Array<DP,3> Ar);
    virtual DP Get_total_energy(Array<complx,3> A);
	virtual DP Get_total_energy_real_space(Array<DP,3> Ar, Array<DP,3> Br);
    virtual DP Get_total_energy(Array<complx,3> A, Array<complx,3> B);
    
    virtual DP Get_local_Sn(Array<complx,3> A, DP n);
    virtual DP Get_total_Sn(Array<complx,3> A, DP n);
	
    virtual DP Get_local_Sn(Array<complx,3> A, Array<complx,3> B, DP n);
    virtual DP Get_total_Sn(Array<complx,3> A, Array<complx,3> B, DP n);
    
	virtual void Compute_local_shell_spectrum
    (  
        Array<complx,3> A, 
        int n, 
        Array<DP,1> local_Sk, 
        Array<DP,1> local_Sk_count
    );

	virtual void Compute_shell_spectrum
    (  
        Array<complx,3> A, 
        int n, 
        Array<DP,1> local_Sk
    );

	virtual void Compute_local_shell_spectrum
    (
        Array<complx,3> A,  Array<complx,3> B,
        int n, 
        Array<DP,1> local_Sk, 
        Array<DP,1> local_Sk_count
    );

	virtual void Compute_shell_spectrum
    (  
     Array<complx,3> A, Array<complx,3> B, 
        int n, 
        Array<DP,1> Sk
    );
    
    virtual DP Get_local_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
    virtual DP Get_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
    
    virtual DP Get_local_entropy_scalar(Array<complx,3> A);
    virtual DP Get_entropy_scalar(Array<complx,3> A);
    
	virtual void Compute_local_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        int n, Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
    );
    
	virtual void Compute_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        int n, Array<DP,2> E1k, Array<DP,2> E2k, Array<DP,2> E3k
    );
    
	virtual void Compute_local_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        int n, Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
    );
    
	virtual void Compute_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        int n, Array<DP,2> E1k, Array<DP,2> E2k, Array<DP,2> E3k
     );
    
	virtual void Compute_local_ring_spectrum
    (
        Array<complx,3> F, 
        int n, Array<DP,2> local_Sk
    );
    
	virtual void Compute_ring_spectrum
    ( 
        Array<complx,3> F, 
        int n, Array<DP,2> Sk
    );
    
   virtual void Compute_local_ring_spectrum
    (
        Array<complx,3> F, Array<complx,3> G, 
        int n, Array<DP,2> local_Sk
    );	
    
   virtual void Compute_ring_spectrum
    (
        Array<complx,3> F, Array<complx,3> G,  
        int n, Array<DP,2> Sk
    );
    
	virtual void Compute_local_cylindrical_ring_spectrum
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 int n, Array<DP,2> local_S1k, Array<DP,2> local_S2k
	 );
    
	virtual void Compute_cylindrical_ring_spectrum
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 int n, Array<DP,2> S1k, Array<DP,2> S2k
	 );
	
	
    
	virtual void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        int n, Array<DP,2> local_S1k, Array<DP,2> local_S2k
    );
    
   virtual void Compute_cylindrical_ring_spectrum
    ( 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        int n, Array<DP,2> S1k, Array<DP,2> S2k
     );
    
    virtual void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> F,
        int n, Array<DP,2> local_S1k, Array<DP,2> local_S2k
    );
    
   virtual void Compute_cylindrical_ring_spectrum
    ( 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> F,  
        int n, Array<DP,2> S1k, Array<DP,2> S2k
     );
    
	virtual void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> F, 
        int n, Array<DP,2> local_Sk
    );
    
	virtual void Compute_cylindrical_ring_spectrum
    (
        Array<complx,3> F, 
        int n, Array<DP,2> Sk
    );
    
	virtual void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> F, Array<complx,3> G,
        int n, Array<DP,2> local_Sk
    );
    
	virtual void Compute_cylindrical_ring_spectrum
    (
        Array<complx,3> F, Array<complx,3> G,
        int n, Array<DP,2> Sk
    );
    

    
	virtual void Compute_local_imag_shell_spectrum_B0
    ( 
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,1> local_Sk
    );
    
	virtual void Compute_imag_shell_spectrum_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,1> Sk
    );
    
	virtual void Compute_local_imag_ring_spectrum_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_Sk
    );
    
	virtual void Compute_imag_ring_spectrum_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> Sk
    );
    
	virtual void Compute_local_imag_cylindrical_ring_spectrum_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_Sk
    );
    
	virtual void Compute_imag_cylindrical_ring_spectrum_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        Array<DP,2> Sk
    );
    
    // EnergyTr
    
    virtual DP Local_shell_mult_single
    ( 
        Array<complx,3> A, Array<complx,3> B, 
        DP inner_radius, DP outer_radius
    );
    
    virtual DP Shell_mult_single
    ( 
        Array<complx,3> A, Array<complx,3> B, 
        DP inner_radius, DP outer_radius
    );
    
	virtual void Local_shell_mult_all
    ( 
        Array<complx,3> A, Array<complx,3> B, 
		Array<DP, 1> shell_radius_array,
        Array<DP,1> local_result
    );
    
	virtual void Shell_mult_all
    ( 
        Array<complx,3> A, Array<complx,3> B, 
		Array<DP, 1> shell_radius_array,
        Array<DP,1> result
    );
    
	virtual void Local_shell_mult_all
    (	 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
		Array<DP, 1> shell_radius_array,
        Array<DP,1> local_result
    );
    
	virtual void Shell_mult_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
		Array<DP, 1> shell_radius_array,
        Array<DP,1> result
    );
    
	virtual void Local_ring_mult_all
	( 
        Array<complx,3> A, Array<complx,3> B, 
        Array<DP,2> local_result
    );
    
    
	virtual void Ring_mult_all
    ( 
        Array<complx,3> A, Array<complx,3> B, 
        Array<DP,2> result
    );
    
	virtual void Local_ring_mult_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_result
    );
    
	virtual void Ring_mult_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        Array<DP,2> result
    );
    
	virtual void Local_cyl_ring_mult_all
    (
        Array<complx,3> A, Array<complx,3> B, 
        Array<DP,2> local_result
    );
    
	virtual void Cyl_ring_mult_all
    (
        Array<complx,3> A, Array<complx,3> B, 
        Array<DP,2> result
    );
    
	virtual void Local_cyl_ring_mult_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        Array<DP,2> local_result
    );
    
	virtual void Cyl_ring_mult_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> result
    );
    
	virtual void Local_shell_mult_all_imagVW_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,1> local_result
    );
    
	virtual void Shell_mult_all_imagVW_B0
    ( 
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        Array<DP,1> result
    );
    
	virtual void Local_ring_mult_all_imagVW_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_result
    );
    
	virtual void Ring_mult_all_imagVW_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        Array<DP,2> result
    );
    
	virtual void Local_cyl_ring_mult_all_imagVW_B0
    ( 
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_result
    );
    
	virtual void Cyl_ring_mult_all_imagVW_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> result
    );
    
    virtual DP Local_shell_mult_vorticity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        DP inner_radius, DP outer_radius
    );
    
	virtual  DP Shell_mult_vorticity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        DP inner_radius, DP outer_radius
    );
    
    virtual DP Local_shell_mult_vector_potential
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        DP inner_radius, DP outer_radius
    );
    
    virtual DP Shell_mult_vector_potential
    ( 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        DP inner_radius, DP outer_radius
    );
    
	virtual void Local_shell_mult_vorticity_all
    ( 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        Array<DP,1> result
    );
    
	virtual void Shell_mult_vorticity_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        Array<DP,1> result
    );
    
	virtual void Local_shell_mult_vector_potential_all
    (
         Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
         Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
         Array<DP,1> result
    );
    
	virtual void Shell_mult_vector_potential_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        Array<DP,1> result
	);
	
	virtual void Fill_array_shell(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius);
	
	virtual void Fill_array_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP left_angle, DP right_angle);
	
	virtual void Fill_array_cylindrical_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP h_lower, DP h_upper);
	
		// Inline functions							 
	virtual int Get_kx(int lx) = 0;
	virtual int Get_lx(int kx) = 0;
	virtual int Get_ix(int kx) = 0;
	
	virtual int Get_ky(int ly) = 0;
	virtual int Get_ly(int ky) = 0;
	virtual int Get_kz(int lz) = 0;
	virtual int Get_lz(int kz) = 0;
	virtual int Get_iz(int kz) = 0;
	
	virtual bool Probe_in_me(int kx, int ky, int kz) = 0;
	
	virtual complx Get_spectral_field(int kx, int ky, int kz, Array<complx,3> A) = 0;
	virtual TinyVector<complx,3> Get_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az) = 0;
//	virtual DP Get_spectral_field(int kx, int ky, int kz, Array<complx,3> A) = 0;
//	virtual TinyVector<DP,3> Get_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az) = 0;
	
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> A,complx field) = 0;
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V) = 0;
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> A, DP field) = 0;
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V) = 0;
	
    virtual void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> A,complx field) = 0;
    virtual void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V) = 0;
    virtual void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> A, DP field) = 0;
    virtual void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V) = 0;

	virtual complx Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A) = 0;
	virtual TinyVector<complx,3> Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az) = 0;
//	virtual DP Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A) = 0;
//	virtual TinyVector<DP,3> Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az) = 0;
	
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A,complx field) = 0;
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V) = 0;
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, DP field) = 0;
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V) = 0;

    virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A,complx field) = 0;
    virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V) = 0;
    virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, DP field) = 0;
    virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V) = 0;
	
	
	virtual	int Get_lx_real_space(int rx) = 0;
	virtual int Get_ly_real_space(int ry) = 0;
	virtual int Get_lz_real_space(int rz) = 0;
	
	virtual	int Get_rx_real_space(int lx) = 0;
	virtual int Get_ry_real_space(int ly) = 0;
	virtual int Get_rz_real_space(int lz) = 0;
	
	virtual bool Probe_in_me_real_space(int rx, int ry, int rz) = 0;
	virtual DP Get_real_field(int rx, int ry, int rz, Array<DP,3> A) = 0;
	virtual TinyVector<DP,3> Get_real_field(int rx, int ry, int rz, Array<DP,3> Ax, Array<DP,3> Ay, Array<DP,3> Az) = 0;
	virtual void Assign_real_field(int rx, int ry, int rz, Array<DP,3> A, DP field) = 0;
	virtual void Assign_real_field(int rx, int ry, int rz, Array<DP,3> Ax, Array<DP,3> Ay, Array<DP,3> Az, TinyVector<DP,3> V) = 0;
	
	
	virtual void Wavenumber(int lx, int ly, int lz, TinyVector<DP,3> & K) = 0;
	virtual void Wavenumber(int lx, int ly, int lz, TinyVector<complx,3> & K) = 0;
	
	virtual DP Kmagnitude(int lx, int ly, int lz) = 0;
	
    virtual int Min_radius_outside() = 0;
	virtual int Max_radius_inside() = 0;
	
	virtual DP Approx_number_modes_in_shell(int radius) = 0;
	
	virtual DP Multiplicity_factor(int lx, int ly, int lz) = 0;
	
	virtual DP Modal_energy(int lx, int ly, int lz, Array<complx,3> A) = 0;
	virtual DP Get_Modal_helicity
	(
		int lx, int ly, int lz,
		Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
	) = 0;
	
	virtual void Compute_Modal_vorticity
	(
		 int lx, int ly, int lz, 
		 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
		 TinyVector<complx,3> &vorticity
	) = 0;
	
	virtual void Compute_Modal_vorticity_y_component
	(
		 int lx, int ly, int lz, 
		 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
		 complx &vort_y
	 ) = 0;
	
	virtual DP AnisKpll(int lx, int ly, int lz) = 0;
	virtual DP AnisKperp(int lx, int ly, int lz) = 0;
	virtual DP AnisKh1(int lx, int ly, int lz) = 0;
	virtual DP AnisKh2(int lx, int ly, int lz) = 0;
	virtual DP Anis_min_Kpll() = 0;
	virtual DP Anis_max_Kpll()  = 0;
	
	virtual int Anis_max_Krho_radius_inside() = 0;
	virtual DP Get_max_polar_angle()  = 0;
	
	virtual DP AnisKvect_polar_angle(int lx, int ly, int lz) = 0;
	virtual DP AnisKvect_azimuthal_angle(int lx, int ly, int lz) = 0;
    
};

#endif


//******************************** End of four_pencil.h  **************************************


