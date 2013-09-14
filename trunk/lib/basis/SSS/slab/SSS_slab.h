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


/*! \file  universal_basic.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug	No known bugs
 */


#ifndef _SSS_SLAB_H
#define _SSS_SLAB_H

#include "def_vars.h"
#include "basicfn_inline.h"
#include "spectral_transform.h"
#include "universal.h"
#include "ArrayOps.h"

using namespace blitz;

//*********************************************************************************************	


class SSS_SLAB: public Universal
{				
public:
	
	SSS_SLAB();
	
	#include "universal_fn_names.h"
	
	void Print_large_Fourier_elements(Array<complx,3> A, string array_name);
	
	void Zero_modes(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
	
	void Zero_modes(Array<complx,3> F);
	
	void Array_mult_ksqr(Array<complx,3> A);
	void Array_divide_ksqr(Array<complx,3> A);
	void Array_exp_ksqr(Array<complx,3> A, DP factor);
	void Array_exp_ksqr(Array<complx,3> A, DP factor, DP hyper_factor, int hyper_exponent);
	void Array_mult_V0_khat_sqr(Array<complx,3> A, TinyVector<DP,3> V0);
    void Fill_Vz(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
	
	// energy
	DP Get_local_Sn(Array<complx,3> A, DP n);
	
    DP Get_local_Sn(Array<complx,3> A, Array<complx,3> B, DP n);
    
	void Compute_local_shell_spectrum
    (  
        Array<complx,3> A, 
        int n, 
        Array<DP,1> local_Sk, 
        Array<DP,1> local_Sk_count
    );


	void Compute_local_shell_spectrum
    (
        Array<complx,3> A,  Array<complx,3> B,
        int n, 
        Array<DP,1> local_Sk, 
        Array<DP,1> local_Sk_count
    );

    DP Get_local_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
    
    DP Get_local_entropy_scalar(Array<complx,3> A);
    
	void Compute_local_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        int n, Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
    );
    
    
	void Compute_local_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        int n, Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
    );
    
    
	void Compute_local_ring_spectrum
    (
        Array<complx,3> F, 
        int n, Array<DP,2> local_Sk
    );

   void Compute_local_ring_spectrum
    (
        Array<complx,3> F, Array<complx,3> G, 
        int n, Array<DP,2> local_Sk
    );	
    
	void Compute_local_cylindrical_ring_spectrum
    (
	 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	 int n, Array<DP,2> local_S1k, Array<DP,2> local_S2k
	 );
	
    
	void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        int n, Array<DP,2> local_S1k, Array<DP,2> local_S2k
    );
    
    
    void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<complx,3> F,
        int n, Array<DP,2> local_S1k, Array<DP,2> local_S2k
    );
    
    
	void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> F, 
        int n, Array<DP,2> local_Sk
    );

	void Compute_local_cylindrical_ring_spectrum
    (
        Array<complx,3> F, Array<complx,3> G,
        int n, Array<DP,2> local_Sk
    );

    
	void Compute_local_imag_shell_spectrum_B0
    ( 
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,1> local_Sk
    );
    
	void Compute_local_imag_ring_spectrum_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_Sk
    );
    
    
	void Compute_local_imag_cylindrical_ring_spectrum_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_Sk
    );
	
	// Energy Tr
	
	DP Local_shell_mult_single
    ( 
        Array<complx,3> A, Array<complx,3> B, 
        DP inner_radius, DP outer_radius
    );
    
	void Local_shell_mult_all
    ( 
        Array<complx,3> A, Array<complx,3> B, 
		Array<DP, 1> shell_radius_array,
        Array<DP,1> local_result
    );
    
    
	void Local_shell_mult_all
    (	 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
		Array<DP, 1> shell_radius_array,
        Array<DP,1> local_result
    );
    
	void Local_ring_mult_all
	( 
        Array<complx,3> A, Array<complx,3> B, 
        Array<DP,2> local_result
    );
    
    
	void Local_ring_mult_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_result
    );
    
	void Local_cyl_ring_mult_all
    (
        Array<complx,3> A, Array<complx,3> B, 
        Array<DP,2> local_result
    );
    
	void Local_cyl_ring_mult_all
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        Array<DP,2> local_result
    );
    

	void Local_shell_mult_all_imagVW_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,1> local_result
    );

    
	void Local_ring_mult_all_imagVW_B0
    (
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_result
    );
    
    
	void Local_cyl_ring_mult_all_imagVW_B0
    ( 
        TinyVector<DP,3> B0,
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
        Array<DP,2> local_result
    );
    
    DP Local_shell_mult_vorticity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        DP inner_radius, DP outer_radius
    );
    
    DP Local_shell_mult_vector_potential
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
        DP inner_radius, DP outer_radius
    );
    
	void Local_shell_mult_vorticity_all
    ( 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
        Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
        Array<DP,1> result
    );
    
    
	void Local_shell_mult_vector_potential_all
    (
         Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
         Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
         Array<DP,1> result
    );
	
	void Fill_array_shell(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius);
	
	void Fill_array_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP left_angle, DP right_angle);
	
	void Fill_array_cylindrical_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP h_lower, DP h_upper);
};

#endif


//******************************** End of scft-slab.h  **************************************


