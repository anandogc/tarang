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


#ifndef _SSS_PENCIL_H
#define _SSS_PENCIL_H

#include "def_vars.h"
#include "basicfn_inline.h"
#include "fftk.h"
#include "universal.h"
#include "ArrayOps.h"


//*********************************************************************************************	


class SSS_PENCIL:public Universal
{				
public:
	
	SSS_PENCIL();
	
	#include "universal_fn_names.h"
	
	void Print_large_Fourier_elements(Array<Complex,3> A, string array_name);
	
	void Zero_modes(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	void Zero_modes(Array<Complex,3> F);
	
	void Assign_sub_array(Range y_range, Range z_range, Range x_range, Array<Real,3> A, Real value);
	
	void Array_mult_ksqr(Array<Complex,3> A);
	void Array_divide_ksqr(Array<Complex,3> A);
	void Array_exp_ksqr(Array<Complex,3> A, Real factor);
	void Array_exp_ksqr(Array<Complex,3> A, Real factor, Real hyper_factor, int hyper_exponent);
	void Array_mult_V0_khat_sqr(Array<Complex,3> A, TinyVector<Real,3> V0);
    void Fill_Vz(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	
	// energy
	Real Get_local_Sn(Array<Complex,3> A, Real n);
	
    Real Get_local_Sn(Array<Complex,3> A, Array<Complex,3> B, Real n);
    
	void Compute_local_shell_spectrum
    (  
        Array<Complex,3> A, 
        int n, 
        Array<Real,1> local_Sk, 
        Array<Real,1> local_Sk_count
    );


	void Compute_local_shell_spectrum
    (
        Array<Complex,3> A,  Array<Complex,3> B,
        int n, 
        Array<Real,1> local_Sk, 
        Array<Real,1> local_Sk_count
    );

    Real Get_local_entropy(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
    
    Real Get_local_entropy_scalar(Array<Complex,3> A);
    
	void Compute_local_ring_spectrum
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        int n, Array<Real,2> local_E1k, Array<Real,2> local_E2k, Array<Real,2> local_E3k
    );
    
    
	void Compute_local_ring_spectrum
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
        int n, Array<Real,2> local_E1k, Array<Real,2> local_E2k, Array<Real,2> local_E3k
    );
    
    
	void Compute_local_ring_spectrum
    (
        Array<Complex,3> F, 
        int n, Array<Real,2> local_Sk
    );

   void Compute_local_ring_spectrum
    (
        Array<Complex,3> F, Array<Complex,3> G, 
        int n, Array<Real,2> local_Sk
    );	
    
	void Compute_local_cylindrical_ring_spectrum
    (
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 int n, Array<Real,2> local_S1k, Array<Real,2> local_S2k
	 );
	
    
	void Compute_local_cylindrical_ring_spectrum
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
        int n, Array<Real,2> local_S1k, Array<Real,2> local_S2k
    );
    
    
    void Compute_local_cylindrical_ring_spectrum
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Complex,3> F,
        int n, Array<Real,2> local_S1k, Array<Real,2> local_S2k
    );
    
    
	void Compute_local_cylindrical_ring_spectrum
    (
        Array<Complex,3> F, 
        int n, Array<Real,2> local_Sk
    );

	void Compute_local_cylindrical_ring_spectrum
    (
        Array<Complex,3> F, Array<Complex,3> G,
        int n, Array<Real,2> local_Sk
    );

    
	void Compute_local_imag_shell_spectrum_B0
    ( 
        TinyVector<Real,3> B0,
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
        Array<Real,1> local_Sk
    );
    
	void Compute_local_imag_ring_spectrum_B0
    (
        TinyVector<Real,3> B0,
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
        Array<Real,2> local_Sk
    );
    
    
	void Compute_local_imag_cylindrical_ring_spectrum_B0
    (
        TinyVector<Real,3> B0,
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
        Array<Real,2> local_Sk
    );
	
	// Energy Tr
	
	Real Local_shell_mult_single
    ( 
        Array<Complex,3> A, Array<Complex,3> B, 
        Real inner_radius, Real outer_radius
    );
    
	void Local_shell_mult_all
    ( 
        Array<Complex,3> A, Array<Complex,3> B, 
		Array<Real, 1> shell_radius_array,
        Array<Real,1> local_result
    );
    
    
	void Local_shell_mult_all
    (	 
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		Array<Real, 1> shell_radius_array,
        Array<Real,1> local_result
    );
    
	void Local_ring_mult_all
	( 
        Array<Complex,3> A, Array<Complex,3> B, 
        Array<Real,2> local_result
    );
    
    
	void Local_ring_mult_all
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
        Array<Real,2> local_result
    );
    
	void Local_cyl_ring_mult_all
    (
        Array<Complex,3> A, Array<Complex,3> B, 
        Array<Real,2> local_result
    );
    
	void Local_cyl_ring_mult_all
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
        Array<Real,2> local_result
    );
    

	void Local_shell_mult_all_imagVW_B0
    (
        TinyVector<Real,3> B0,
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
        Array<Real,1> local_result
    );

    
	void Local_ring_mult_all_imagVW_B0
    (
        TinyVector<Real,3> B0,
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
        Array<Real,2> local_result
    );
    
    
	void Local_cyl_ring_mult_all_imagVW_B0
    ( 
        TinyVector<Real,3> B0,
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
        Array<Real,2> local_result
    );
    
    Real Local_shell_mult_vorticity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
        Real inner_radius, Real outer_radius
    );
    
    Real Local_shell_mult_vector_potential
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
        Real inner_radius, Real outer_radius
    );
    
	void Local_shell_mult_vorticity_all
    ( 
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
        Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
        Array<Real,1> result
    );
    
    
	void Local_shell_mult_vector_potential_all
    (
         Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
         Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
         Array<Real,1> result
    );
	
	void Fill_array_shell(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius);
	
	void Fill_array_ring(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius, Real left_angle, Real right_angle);
	
	void Fill_array_cylindrical_ring(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius, Real h_lower, Real h_upper);
    
};

#endif


//******************************** End of scft-slab.h  **************************************


