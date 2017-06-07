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
#include "BasicIO.h"

//*********************************************************************************************	


class Universal
{				
public:
	//HDF5 IO plans
	h5::Plan H5_real;
	h5::Plan H5_full;
	h5::Plan H5_kz0_full;
	h5::Plan H5_in_reduced;
	h5::Plan H5_out_reduced;
	h5::Plan H5_in_kz0_reduced;
	h5::Plan H5_out_kz0_reduced;
	vector<h5::Plan> H5_slices;

	// basic

	virtual int Get_number_modes_in_shell(Real inner_radius, Real outer_radius);
	virtual void Print_large_Fourier_elements(Array<Complex,3> A, string array_name="Array");
	virtual void Array_mult_ksqr(Array<Complex,3> A);
	virtual void Array_divide_ksqr(Array<Complex,3> A);
	virtual void Array_exp_ksqr(Array<Complex,3> A, Real factor);
	virtual void Array_exp_ksqr(Array<Complex,3> A, Real factor, Real hyper_factor, int hyper_exponent);
	virtual void Array_mult_V0_khat_sqr(Array<Complex,3> A, TinyVector<Real,3> V0);
	virtual void Fill_Vz(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	virtual void Compute_divergence(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, Array<Complex,3> div, string field_or_nlin, Real &total_abs_div, bool print_switch);
	virtual void Zero_modes(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	virtual void Zero_modes(Array<Complex,3> F);
	

	// Local functions...
	virtual void Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz) = 0;
	virtual void Last_component(int kx, int ky, int kz, Complex& Vx, Complex& Vy, Complex& Vz) = 0;
	virtual void Dealias(Array<Complex,3> A) = 0;
	virtual bool Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius) = 0;
	virtual void Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A) = 0;
	virtual void Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A) = 0;
	virtual void Test_reality_condition_in_Array(Array<Complex,3> A) = 0;

	
	// transform
	
	virtual void Forward_transform(Array<Real,3> Ar, Array<Complex,3> A) = 0;
	virtual void Inverse_transform(Array<Complex,3> A, Array<Real,3> Ar) = 0;
	
	virtual void Xderiv(Array<Complex,3> A, Array<Complex,3> B) = 0;
	virtual void Yderiv(Array<Complex,3> A, Array<Complex,3> B) = 0;
	virtual void Zderiv(Array<Complex,3> A, Array<Complex,3> B) = 0;

	virtual void Add_Xderiv(Array<Complex,3> A, Array<Complex,3> B) = 0;
	virtual void Add_Yderiv(Array<Complex,3> A, Array<Complex,3> B) = 0;
	virtual void Add_Zderiv(Array<Complex,3> A, Array<Complex,3> B) = 0;
	
	virtual void Xderiv(Array<Real,3> A, Array<Real,3> B) = 0;
	
	virtual void Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B) = 0;
	virtual void Subtract_Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B) = 0;
	// Energy
	
	//virtual Real Get_local_energy_XZ_plane(Array<Complex,3> A, int ny) = 0;
	
 //    virtual Real Get_local_energy_XZ_plane
	// (
	// 	Array<Complex,3> A, Array<Complex,3> B, 
	// 	int ny
	// ) = 0;
	
	virtual Real Get_local_energy_real_space(Array<Real,3> Ar) = 0;
	virtual Real Get_local_energy(Array<Complex,3> A) = 0;
	virtual Real Get_local_energy_real_space(Array<Real,3> Ar, Array<Real,3> Br) = 0;
	virtual Real Get_local_energy(Array<Complex,3> A, Array<Complex,3> B) = 0;
	
	virtual void Compute_vorticity
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Real inner_radius, Real outer_radius
	);

	virtual void Compute_local_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Real &local_helicity1, Real &local_helicity2,
	 Real &local_dissipation_H1, Real &local_dissipation_H2
	 ) = 0;
	
	virtual void Compute_total_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Real &total_helicity1, Real &total_helicity2,
	 Real &total_dissipation_H1, Real &total_dissipation_H2
	 ) = 0;
	
	virtual void Compute_local_shell_spectrum_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Real,1> local_H1k1, Array<Real,1> local_H1k2, Array<Real,1> local_H1k3,
	 Array<Real,1> local_H1k_count
	 ) = 0;
	
	virtual void Compute_shell_spectrum_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Real,1> H1k1, Array<Real,1> H1k2, Array<Real,1> H1k3
	 ) = 0;
	
  
  virtual void Compute_local_shell_spectrum_helicity2
  (
   Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
   Array<Real,1> local_H1k1, Array<Real,1> local_H1k2, Array<Real,1> local_H1k3,
   Array<Real,1> local_H1k_count
   ) = 0;
  
  virtual void Compute_shell_spectrum_helicity2
  (
   Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
   Array<Real,1> H1k1, Array<Real,1> H1k2, Array<Real,1> H1k3
   ) = 0;
  virtual
  
  void Compute_local_ring_spectrum_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Real,2> local_H1k1, Array<Real,2> local_H1k2, Array<Real,2> local_H1k3
	 ) = 0;
	
	
	virtual void Compute_ring_spectrum_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Real,2> local_H1k1, Array<Real,2> local_H1k2, Array<Real,2> local_H1k3
	 ) = 0;
	
	virtual void Compute_local_cylindrical_ring_spectrum_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Real,2> local_H1k1, Array<Real,2> local_H1k2
	 ) = 0;
	
	virtual void Compute_cylindrical_ring_spectrum_helicity
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Real,2> local_H1k1, Array<Real,2> local_H1k2
	 ) = 0;



	
	virtual Real Get_total_energy_real_space(Array<Real,3> Ar);
	virtual Real Get_total_energy(Array<Complex,3> A);
	virtual Real Get_total_energy_real_space(Array<Real,3> Ar, Array<Real,3> Br);
	virtual Real Get_total_energy(Array<Complex,3> A, Array<Complex,3> B);
	
	virtual Real Get_local_Sn(Array<Complex,3> A, Real n);
	virtual Real Get_total_Sn(Array<Complex,3> A, Real n);
	
	virtual Real Get_local_Sn(Array<Complex,3> A, Array<Complex,3> B, Real n);
	virtual Real Get_total_Sn(Array<Complex,3> A, Array<Complex,3> B, Real n);
	
	virtual void Compute_local_shell_spectrum
	(  
		Array<Complex,3> A, 
		int n, 
		Array<Real,1> local_Sk, 
		Array<Real,1> local_Sk_count
	);

	virtual void Compute_shell_spectrum
	(  
		Array<Complex,3> A, 
		int n, 
		Array<Real,1> local_Sk
	);

	virtual void Compute_local_shell_spectrum
	(
		Array<Complex,3> A,  Array<Complex,3> B,
		int n, 
		Array<Real,1> local_Sk, 
		Array<Real,1> local_Sk_count
	);

	virtual void Compute_shell_spectrum
	(  
	 Array<Complex,3> A, Array<Complex,3> B, 
		int n, 
		Array<Real,1> Sk
	);
	
	virtual Real Get_local_entropy(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	virtual Real Get_entropy(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	
	virtual Real Get_local_entropy_scalar(Array<Complex,3> A);
	virtual Real Get_entropy_scalar(Array<Complex,3> A);
	
	virtual void Compute_local_ring_spectrum
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		int n, Array<Real,2> local_E1k, Array<Real,2> local_E2k, Array<Real,2> local_E3k
	);
	
	virtual void Compute_ring_spectrum
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		int n, Array<Real,2> E1k, Array<Real,2> E2k, Array<Real,2> E3k
	);
	
	virtual void Compute_local_ring_spectrum
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		int n, Array<Real,2> local_E1k, Array<Real,2> local_E2k, Array<Real,2> local_E3k
	);
	
	virtual void Compute_ring_spectrum
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		int n, Array<Real,2> E1k, Array<Real,2> E2k, Array<Real,2> E3k
	 );
	
	virtual void Compute_local_ring_spectrum
	(
		Array<Complex,3> F, 
		int n, Array<Real,2> local_Sk
	);
	
	virtual void Compute_ring_spectrum
	( 
		Array<Complex,3> F, 
		int n, Array<Real,2> Sk
	);
	
   virtual void Compute_local_ring_spectrum
	(
		Array<Complex,3> F, Array<Complex,3> G, 
		int n, Array<Real,2> local_Sk
	);	
	
	virtual void Compute_ring_spectrum
	(
		Array<Complex,3> F, Array<Complex,3> G,  
		int n, Array<Real,2> Sk
	);

	virtual void Compute_local_helical_ring_spectrum
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Complex,3> helicalAx, Array<Complex,3> helicalAy, Array<Complex,3> helicalAz,
	 int n,
	 Array<Real,2> local_hk1, Array<Real,2> local_hk2, Array<Real,2> local_hk3
	);

	virtual void Compute_helical_ring_spectrum
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 Array<Complex,3> helicalAx, Array<Complex,3> helicalAy, Array<Complex,3> helicalAz,
	 int n,
	 Array<Real,2> hk1, Array<Real,2> hk2, Array<Real,2> hk3
	);
	
	virtual void Compute_local_cylindrical_ring_spectrum
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 int n, Array<Real,2> local_S1k, Array<Real,2> local_S2k
	);
	
	virtual void Compute_cylindrical_ring_spectrum
	(
	 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	 int n, Array<Real,2> S1k, Array<Real,2> S2k
	);
	
	
	
	virtual void Compute_local_cylindrical_ring_spectrum
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		int n, Array<Real,2> local_S1k, Array<Real,2> local_S2k
	);
	
   virtual void Compute_cylindrical_ring_spectrum
	( 
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		int n, Array<Real,2> S1k, Array<Real,2> S2k
	 );
	
	virtual void Compute_local_cylindrical_ring_spectrum
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		Array<Complex,3> F,
		int n, Array<Real,2> local_S1k, Array<Real,2> local_S2k
	);
	
   virtual void Compute_cylindrical_ring_spectrum
	( 
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		Array<Complex,3> F,  
		int n, Array<Real,2> S1k, Array<Real,2> S2k
	 );
	
	virtual void Compute_local_cylindrical_ring_spectrum
	(
		Array<Complex,3> F, 
		int n, Array<Real,2> local_Sk
	);
	
	virtual void Compute_cylindrical_ring_spectrum
	(
		Array<Complex,3> F, 
		int n, Array<Real,2> Sk
	);
	
	virtual void Compute_local_cylindrical_ring_spectrum
	(
		Array<Complex,3> F, Array<Complex,3> G,
		int n, Array<Real,2> local_Sk
	);
	
	virtual void Compute_cylindrical_ring_spectrum
	(
		Array<Complex,3> F, Array<Complex,3> G,
		int n, Array<Real,2> Sk
	);
	

	
	virtual void Compute_local_imag_shell_spectrum_B0
	( 
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,1> local_Sk
	);
	
	virtual void Compute_imag_shell_spectrum_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,1> Sk
	);
	
	virtual void Compute_local_imag_ring_spectrum_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> local_Sk
	);
	
	virtual void Compute_imag_ring_spectrum_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> Sk
	);
	
	virtual void Compute_local_imag_cylindrical_ring_spectrum_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> local_Sk
	);
	
	virtual void Compute_imag_cylindrical_ring_spectrum_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Array<Real,2> Sk
	);
	
	// EnergyTr
	
	virtual Real Local_shell_mult_single
	( 
		Array<Complex,3> A, Array<Complex,3> B, 
		Real inner_radius, Real outer_radius
	);
	
	virtual Real Shell_mult_single
	( 
		Array<Complex,3> A, Array<Complex,3> B, 
		Real inner_radius, Real outer_radius
	);
	
	virtual void Local_shell_mult_all
	( 
		Array<Complex,3> A, Array<Complex,3> B, 
		Array<Real, 1> shell_radius_array,
		Array<Real,1> local_result
	);
	
	virtual void Shell_mult_all
	( 
		Array<Complex,3> A, Array<Complex,3> B, 
		Array<Real, 1> shell_radius_array,
		Array<Real,1> result
	);
	
	virtual void Local_shell_mult_all
	(	 
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		Array<Real, 1> shell_radius_array,
		Array<Real,1> local_result
	);
	
	virtual void Shell_mult_all
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		Array<Real, 1> shell_radius_array,
		Array<Real,1> result
	);
	
	virtual void Local_ring_mult_all
	( 
		Array<Complex,3> A, Array<Complex,3> B, 
		Array<Real,2> local_result
	);
	
	
	virtual void Ring_mult_all
	( 
		Array<Complex,3> A, Array<Complex,3> B, 
		Array<Real,2> result
	);
	
	virtual void Local_ring_mult_all
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> local_result
	);
	
	virtual void Ring_mult_all
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Array<Real,2> result
	);
	
	virtual void Local_cyl_ring_mult_all
	(
		Array<Complex,3> A, Array<Complex,3> B, 
		Array<Real,2> local_result
	);
	
	virtual void Cyl_ring_mult_all
	(
		Array<Complex,3> A, Array<Complex,3> B, 
		Array<Real,2> result
	);
	
	virtual void Local_cyl_ring_mult_all
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Array<Real,2> local_result
	);
	
	virtual void Cyl_ring_mult_all
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> result
	);
	
	virtual void Local_shell_mult_all_imagVW_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,1> local_result
	);
	
	virtual void Shell_mult_all_imagVW_B0
	( 
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Array<Real,1> result
	);
	
	virtual void Local_ring_mult_all_imagVW_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> local_result
	);
	
	virtual void Ring_mult_all_imagVW_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Array<Real,2> result
	);
	
	virtual void Local_cyl_ring_mult_all_imagVW_B0
	( 
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> local_result
	);
	
	virtual void Cyl_ring_mult_all_imagVW_B0
	(
		TinyVector<Real,3> B0,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
		Array<Real,2> result
	);
	
	virtual Real Local_shell_mult_vorticity
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Real inner_radius, Real outer_radius
	);
	
	virtual  Real Shell_mult_vorticity
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Real inner_radius, Real outer_radius
	);
	
	virtual Real Local_shell_mult_vector_potential
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Real inner_radius, Real outer_radius
	);
	
	virtual Real Shell_mult_vector_potential
	( 
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		Real inner_radius, Real outer_radius
	);
	
	virtual void Local_shell_mult_vorticity_all
	( 
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		Array<Real,1> result
	);
	
	virtual void Shell_mult_vorticity_all
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		Array<Real,1> result
	);
	
	virtual void Local_shell_mult_vector_potential_all
	(
		 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
		 Array<Real,1> result
	);
	
	virtual void Shell_mult_vector_potential_all
	(
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
		Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
		Array<Real,1> result
	);
	
	virtual void Fill_array_shell(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius);
	
	virtual void Fill_array_ring(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius, Real left_angle, Real right_angle);
	
	virtual void Fill_array_cylindrical_ring(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius, Real h_lower, Real h_upper);
	
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
	
	virtual Complex Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> A) = 0;
	virtual TinyVector<Complex,3> Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az) = 0;
//	virtual Real Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> A) = 0;
//	virtual TinyVector<Real,3> Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az) = 0;
	
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A,Complex field) = 0;
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V) = 0;
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Real field) = 0;
	virtual void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V) = 0;
	
	virtual void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A,Complex field) = 0;
	virtual void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V) = 0;
	virtual void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Real field) = 0;
	virtual void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V) = 0;

	virtual Complex Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A) = 0;
	virtual TinyVector<Complex,3> Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az) = 0;
//	virtual Real Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A) = 0;
//	virtual TinyVector<Real,3> Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az) = 0;
	
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Complex field) = 0;
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V) = 0;
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Real field) = 0;
	virtual void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V) = 0;

	virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Complex field) = 0;
	virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V) = 0;
	virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Real field) = 0;
	virtual void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V) = 0;
	
	virtual void Craya_to_helical(Complex U1, Complex U2, Complex &U_plus, Complex &U_minus) = 0;
	virtual void Helical_to_Craya(Complex U_plus, Complex U_minus, Complex &U1, Complex &U2) = 0;
	
	virtual void Craya_to_cartesian(int lx, int ly, int lz, Complex U1, Complex U2, Complex &vpll, Complex &vh1, Complex &vh2) = 0;
	virtual void Cartesian_to_Craya(int lx, int ly, int lz, Complex vpll, Complex vh1, Complex vh2, Complex &U1, Complex &U2) = 0;
	// shubhadeep //
	virtual void Helical_to_cartesian(int lx, int ly, int lz, Complex U_plus, Complex U_minus, Complex &vpll, Complex &vh1, Complex &vh2) = 0;
	virtual void Cartesian_to_helical(int lx, int ly, int lz, Complex vpll, Complex vh1, Complex vh2, Complex &U_plus, Complex &U_minus) = 0;
	// Shubhadeep //
	virtual	int Get_lx_real_space(int rx) = 0;
	virtual int Get_ly_real_space(int ry) = 0;
	virtual int Get_lz_real_space(int rz) = 0;
	
	virtual	int Get_rx_real_space(int lx) = 0;
	virtual int Get_ry_real_space(int ly) = 0;
	virtual int Get_rz_real_space(int lz) = 0;
	
	virtual bool Probe_in_me_real_space(int rx, int ry, int rz) = 0;
	virtual Real Get_real_field(int rx, int ry, int rz, Array<Real,3> A) = 0;
	virtual TinyVector<Real,3> Get_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az) = 0;
	virtual void Assign_real_field(int rx, int ry, int rz, Array<Real,3> A, Real field) = 0;
	virtual void Assign_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az, TinyVector<Real,3> V) = 0;
	
	
	virtual void Wavenumber(int lx, int ly, int lz, TinyVector<Real,3> & K) = 0;
	virtual void Wavenumber(int lx, int ly, int lz, TinyVector<Complex,3> & K) = 0;
	
	virtual Real Kmagnitude(int lx, int ly, int lz) = 0;
	
	virtual int Min_radius_outside() = 0;
	virtual int Max_radius_inside() = 0;
	
	virtual Real Approx_number_modes_in_shell(int radius) = 0;
	
	virtual Real Multiplicity_factor(int lx, int ly, int lz) = 0;
	
	virtual Real Modal_energy(int lx, int ly, int lz, Array<Complex,3> A) = 0;
	virtual Real Get_Modal_helicity
	(
		int lx, int ly, int lz,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az
	) = 0;
	
	virtual void Compute_Modal_vorticity
	(
		 int lx, int ly, int lz, 
		 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		 TinyVector<Complex,3> &vorticity
	) = 0;
	
	virtual void Compute_Modal_vorticity_y_component
	(
		 int lx, int ly, int lz, 
		 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		 Complex &vort_y
	 ) = 0;
	
	virtual Real AnisKpll(int lx, int ly, int lz) = 0;
	virtual Real AnisKperp(int lx, int ly, int lz) = 0;
	virtual Real AnisKh1(int lx, int ly, int lz) = 0;
	virtual Real AnisKh2(int lx, int ly, int lz) = 0;
	virtual Real Anis_min_Kpll() = 0;
	virtual Real Anis_max_Kpll()  = 0;
	
	virtual int Anis_max_Krho_radius_inside() = 0;
	virtual Real Get_max_polar_angle()  = 0;
	
	virtual Real AnisKvect_polar_angle(int lx, int ly, int lz) = 0;
	virtual Real AnisKvect_azimuthal_angle(int lx, int ly, int lz) = 0;


	virtual int Read(Array<Complex,3> A, h5::Plan plan, string file_name, string dataset_name="") = 0;
	virtual int Read(Array<Real,3> Ar, h5::Plan plan, string file_name, string dataset_name="") = 0;

	virtual int Write(Array<Complex,3> A, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name="") = 0;
	virtual int Write(Array<Real,3> Ar, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name="") = 0;
};

#endif


//******************************** End of four_pencil.h  **************************************


