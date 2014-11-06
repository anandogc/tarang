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

/*! \file scft_energy.cc
 * 
 * @sa scft_energy.h
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 22/08/2008
 * @bug Line 1299: Helicity -- do we need factor 2?
 */

#include "SSS_pencil.h"
#include "SSS_pencil_inline.h"
#include "shell_etc_indices.h"

/**********************************************************************************************

	Computes total energy |f(m,k)|^2  without k=0.
		Factor 2 because of negative k modes are not stored.
		
	for ky plane Eky = 2*sum(A^2) - sum(A(m,ky,0)^2) - sum(A(0,ky,kz)^2) + A(0,ky,0)^2
	
	for ky = 0 plane: subtract (1/2)A(0)^2 (origin) from the above
		
***********************************************************************************************/


DP SSS_PENCIL::Get_local_energy_real_space(Array<DP,3> Ar)
{
	
	return Array_sqr(Ar);
}


DP SSS_PENCIL::Get_local_energy(Array<complx,3> A)  
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP total =  4*Array_sqr(Ar);
	
	// subtractions |A(:,ky=0,:)|^2;
	if (my_y_pcoord == 0) {
		total -= 2*Array_sqr(Ar(Range::all(), 0, Range::all()));
		total += Array_sqr(Ar(Range::all(), 0, 0));
	}
	
	// subtractions |A(:, :, kz=0)|^2
	total -= 2*Array_sqr(Ar(Range::all(), Range::all(), 0));
	
	
	if (my_x_pcoord == 0) {
		total -= 2*Array_sqr(Ar(0, Range::all(), Range::all()));
		total += Array_sqr(Ar(0, Range::all(), 0));
	}
	
	if ((my_x_pcoord == 0) && (my_y_pcoord == 0)) {
		total +=  Array_sqr(Ar(0, 0, Range::all()));
		total -= my_pow(Ar(0,0,0),2)/2;
	}
	
	return total;
}

/**********************************************************************************************

		Computes total A(k)*conj(B(k)) without k=0.
	
		
***********************************************************************************************/

DP SSS_PENCIL::Get_local_energy_real_space(Array<DP,3> Ar, Array<DP,3> Br)
{
	return mydot(Ar, Br)/ (DP(Nx)*DP(Ny)*DP(Nz));
}

DP SSS_PENCIL::Get_local_energy(Array<complx,3> A, Array<complx,3> B)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP total =  4*mydot(Ar, Br);
	
	// subtractions |A(:,ky=0,:)|^2;
	if (my_y_pcoord == 0) {
		total -= 2*mydot(Ar(Range::all(),0,Range::all()), Br(Range::all(),0,Range::all()));
		total += mydot(Ar(Range::all(),0,0), Br(Range::all(),0,0));
	}
	
	// subtractions |A(:, :, kz=0)|^2
	total -= 2*mydot(Ar(Range::all(),Range::all(),0), Br(Range::all(),Range::all(),0));
	
	
	if (my_x_pcoord == 0) {
		total -= 2*mydot(Ar(0,Range::all(),Range::all()), Br(0,Range::all(),Range::all()));
		total += mydot(Ar(0,Range::all(),0), Br(0,Range::all(),0));
	}
	
	if ((my_x_pcoord == 0) && (my_y_pcoord == 0)) {
		total +=  mydot(Ar(0,0,Range::all()), Br(0,0,Range::all()));
		total -=  Ar(0,0,0)*Br(0,0,0)/2;
	}
	
	return total;
}

/**********************************************************************************************

	Compute Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2
	Multiplication factor = 1 for kz>0 because of energy spectrum details.
 
	not for 2D

***********************************************************************************************/


void SSS_PENCIL::Compute_local_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &local_helicity1, DP &local_helicity2, 
	DP &local_k2H1, DP &local_k2H2
)
{
	local_helicity1 = local_helicity2 = 0;
	local_k2H1 = local_k2H2 = 0;
}

//

void SSS_PENCIL::Compute_total_helicity
( 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &local_k2H1, DP &local_k2H2
)
{
	total_helicity1 = total_helicity2 = 0;
	local_k2H1 = local_k2H2 = 0;
	
}



/**********************************************************************************************

	Compute helicity spectrum
	Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2

	Not for 2D
***********************************************************************************************/

void SSS_PENCIL::Compute_local_shell_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3, 
	Array<DP,1> local_H1k_count
)
{
	local_H1k1 = 0;
	local_H1k2 = 0;
	local_H1k3 = 0;
}

//

void SSS_PENCIL::Compute_shell_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3
)	
{
	H1k1 = 0;
	H1k2 = 0;
	H1k3 = 0;
}



//*********************************************************************************************
//
// helicity spectrum
// Not for 2D
//

void SSS_PENCIL::Compute_local_ring_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,2> local_H1k1, Array<DP,2> local_H1k2, Array<DP,2> local_H1k3
)												
{

	local_H1k1 = 0.0;
	local_H1k2 = 0.0;
	local_H1k3 = 0.0;
}


void SSS_PENCIL::Compute_ring_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,2> H1k1, Array<DP,2> H1k2, Array<DP,2> H1k3
)		
{
	H1k1 = 0.0;
	H1k2 = 0.0;
	H1k3 = 0.0;
}


//*********************************************************************************************

void SSS_PENCIL::Compute_local_cylindrical_ring_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
 Array<DP,2> local_H1k1,  Array<DP,2> local_H1k2
)
{
	local_H1k1 = 0.0;
	local_H1k2 = 0.0;
}

//
//
void SSS_PENCIL::Compute_cylindrical_ring_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	 Array<DP,2> H1k1,  Array<DP,2> H1k2
)
{
	H1k1 = 0.0;
	H1k2 = 0.0;
	
}


/**********************************************************************************************
 
 Computes sum[ (1/2) K^n |A(k)|^2];   k=0 excluded
 
 ***********************************************************************************************/


//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
DP SSS_PENCIL::Get_local_Sn(Array<complx,3> A, DP n)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Sn = 0.0;
	DP Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) { 
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz) *my_pow(Kmag,n)* Vsqr(Ar(lx,ly,lz));
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if ((my_id == master_id) && (n==0))
		Sn -= Multiplicity_factor(0, 0, 0) * pow2(abs(A(0,0,0)));
	
	return Sn;
}

//*******

DP SSS_PENCIL::Get_local_Sn(Array<complx,3> A, Array<complx,3> B, DP n)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP Sn = 0.0;
	DP Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) { 
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz) *my_pow(Kmag,n)* Ar(lx,ly,lz)*Br(lx,ly,lz);
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if ((my_id == master_id) && (n==0))
		Sn -= Multiplicity_factor(0, 0, 0) * real(A(0,0,0)*B(0,0,0));
	
	return Sn;
}
//*************


void SSS_PENCIL::Compute_local_shell_spectrum
(
 Array<complx,3> A,
 int n,
 Array<DP,1> local_Sk,
 Array<DP,1> local_Sk_count
 )
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_Sk_count = 0.0;
	local_Sk = 0.0;
	
	DP Kmag;
	int index;
	DP factor;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) { 
				Kmag = Kmagnitude(lx, ly, lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) {
					factor = Multiplicity_factor(lx,ly,lz);
					
					local_Sk(index) += factor* my_pow(Kmag,n)* Vsqr(Ar(lx,ly,lz));
					local_Sk_count(index) += 2*factor;
					// Mult by 2 because factor is offset by 2 in Multiply_factor
				}
			}
}


//*********************************************************************************************
//
// Re(A. conj(B))
//


void SSS_PENCIL::Compute_local_shell_spectrum
(
 Array<complx,3> A,  Array<complx,3> B,
 int n,
 Array<DP,1> local_Sk,
 Array<DP,1> local_Sk_count
 )
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_Sk_count=0;
	local_Sk = 0.0;
	
	DP Kmag;
	int index;
	DP factor;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) { 
				Kmag = Kmagnitude(lx, ly, lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) {
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(index) += factor* my_pow(Kmag,n)* Ar(lx,ly,lz)*Br(lx,ly,lz);
					
					local_Sk_count(index) += 2*factor;
				}
			}
	
}
//*******************


DP SSS_PENCIL::Get_local_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
	
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP modal_energy, prob, local_entropy;
	
	DP total_energy = Get_total_energy(Ax) + Get_total_energy(Ay) + Get_total_energy(Az);
	
	MPI_Bcast( &total_energy, 1, MPI_DP, master_id, MPI_COMM_WORLD);
	
	local_entropy = 0.0;
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(Axr(lx,ly,lz),Ayr(lx,ly,lz),Azr(lx,ly,lz));
				
				prob = modal_energy/ total_energy;
				
				if (prob > MYEPS)
					local_entropy += (-prob * log(prob)/log(2.0));
			}
	
	return local_entropy;
}


//
// Scalar
//


DP SSS_PENCIL::Get_local_entropy_scalar(Array<complx,3> A)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	DP modal_energy, prob;
	
	DP local_entropy = 0.0;
	
	DP total_energy = Get_total_energy(A);
	
	MPI_Bcast( &total_energy, 1, MPI_DP, master_id, MPI_COMM_WORLD);
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) { 
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(Ar(lx,ly,lz));
				
				prob = modal_energy / total_energy;
				
				if (prob > MYEPS)
					local_entropy += (-prob*log(prob)/log(2.0));
			}
	
	return local_entropy;
	
}


void SSS_PENCIL::Compute_local_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 int n,
 Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
 )
{
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_E1k = 0.0;
	local_E2k = 0.0;
	local_E3k = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx,ly,lz);
					local_E1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(Axr(lx,ly,lz));
					local_E2k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(Ayr(lx,ly,lz));
					local_E3k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(Azr(lx,ly,lz));
				}
			}
	
}





//*********************************************************************************************
//
//  A.B
//


void SSS_PENCIL::Compute_local_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 int n, Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
 )
{
	
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Array<DP,3> Bxr=Array<DP,3>(reinterpret_cast<DP*>(Bx.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Byr=Array<DP,3>(reinterpret_cast<DP*>(By.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Bzr=Array<DP,3>(reinterpret_cast<DP*>(Bz.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_E1k = 0.0;
	local_E2k = 0.0;
	local_E3k = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
					local_E1k(shell_index, sector_index) += factor*my_pow(Kmag,n)*  Axr(lx,ly,lz)* Bxr(lx,ly,lz);
					local_E2k(shell_index, sector_index) += factor*my_pow(Kmag,n)*  Ayr(lx,ly,lz)* Byr(lx,ly,lz);
					local_E3k(shell_index, sector_index) += factor*my_pow(Kmag,n)*  Azr(lx,ly,lz)* Bzr(lx,ly,lz);
				}
			}
}



//*********************************************************************************************
// scalar
//

void SSS_PENCIL::Compute_local_ring_spectrum
(
 Array<complx,3> F,
 int n,
 Array<DP,2> local_Sk
 )
{
	Array<DP,3> Fr=Array<DP,3>(reinterpret_cast<DP*>(F.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_Sk = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)* Vsqr(Fr(lx,ly,lz));
				}
			}
}


//
// F*G
//

void SSS_PENCIL::Compute_local_ring_spectrum
(
 Array<complx,3> F,
 Array<complx,3> G,
 int n,
 Array<DP,2> local_Sk
 )
{
	Array<DP,3> Fr=Array<DP,3>(reinterpret_cast<DP*>(F.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Gr=Array<DP,3>(reinterpret_cast<DP*>(G.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_Sk = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)*Fr(lx,ly,lz)*Gr(lx,ly,lz);
				}
			}
	
	
}



/**********************************************************************************************
 
 CYLINDERICAL RING SPECTRUM
 
 Not for 2D
 
 ***********************************************************************************************/

void SSS_PENCIL::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 int n,
 Array<DP,2> local_S1k, Array<DP,2> local_S2k
 )
{
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	DP anisV1;
	
	DP Kmag, Kpll, Kperp;
	DP V2sqr;
	DP factor;
	int shell_index, slab_index;
	TinyVector<DP,3> V;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				
				V = Axr(lx,ly,lz), Ayr(lx,ly,lz), Azr(lx,ly,lz);
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					anisV1 = V(global.field.anisotropy_dirn-1);
					
					local_S1k(shell_index, slab_index) += factor*my_pow(Kmag,n)* Vsqr(anisV1);
					
					V2sqr = Vsqr(V)- Vsqr(anisV1);
					local_S2k(shell_index, slab_index) += factor*my_pow(Kmag,n)*V2sqr;
				}
			}
	
}


//*********************************************************************************************

void SSS_PENCIL::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 int n,
 Array<DP,2> local_S1k, Array<DP,2> local_S2k
 )
{
	
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Array<DP,3> Bxr=Array<DP,3>(reinterpret_cast<DP*>(Bx.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Byr=Array<DP,3>(reinterpret_cast<DP*>(By.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Bzr=Array<DP,3>(reinterpret_cast<DP*>(Bz.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	DP anisV1, anisW1;
	TinyVector<DP,3> V, Vrho, W, Wrho;
	
	DP Kmag, Kpll, Kperp;
	DP factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					V = Axr(lx,ly,lz), Ayr(lx,ly,lz), Azr(lx,ly,lz);
					W = Bxr(lx,ly,lz), Byr(lx,ly,lz), Bzr(lx,ly,lz);
					
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					if (global.field.anisotropy_dirn == 1) {
						anisV1 = V(0); anisW1 = W(0);
						
						Vrho = ZERO, V(1), V(2);
						Wrho = ZERO, W(1), W(2);
					}
					
					if (global.field.anisotropy_dirn == 2) {
						anisV1 = V(1); anisW1 = W(1);
						
						Vrho = V(0), ZERO, V(2);
						Wrho = W(0), ZERO, W(2);
					}
					
					if (global.field.anisotropy_dirn == 3) {
						anisV1 = V(2); anisW1 = W(2);
						
						Vrho = V(0), V(1), ZERO;
						Wrho = W(0), W(1), ZERO;
					}
					
					local_S1k(shell_index, slab_index) += factor*my_pow(Kmag,n)* anisV1*anisW1;
					local_S2k(shell_index, slab_index) += factor*my_pow(Kmag,n)* dot(Vrho,Wrho);
				}
			}
	
}



//*********************************************************************************************

void SSS_PENCIL::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> F,
 int n,
 Array<DP,2> local_S1k, Array<DP,2> local_S2k
 )
{
	
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Array<DP,3> Fr=Array<DP,3>(reinterpret_cast<DP*>(F.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	DP anisV1, Fval;
	
	DP Kmag, Kpll, Kperp;
	DP VTperp;
	DP factor;
	int shell_index, slab_index;
	TinyVector<DP,3> V, Vless;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				
				V = Axr(lx,ly,lz), Ayr(lx,ly,lz), Azr(lx,ly,lz);
				Fval = Fr(lx,ly,lz);
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					anisV1 = V(global.field.anisotropy_dirn-1);
					Vless = V - anisV1;
					
					local_S1k(shell_index, slab_index) += factor*my_pow(Kmag,n)* anisV1*Fval;
					
					VTperp = sum(Vless) * Fval;
					local_S2k(shell_index, slab_index) += factor*my_pow(Kmag,n)*VTperp;
				}
			}
	
	
}




//*********************************************************************************************

void SSS_PENCIL::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> F,
 int n,
 Array<DP,2> local_Sk
 )
{
	Array<DP,3> Fr=Array<DP,3>(reinterpret_cast<DP*>(F.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_Sk = 0.0;
	
	DP Kmag, Kperp, Kpll;
	DP factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)* Vsqr(Fr(lx,ly,lz));
				}
			}
	
}



//*********************************************************************************************

void SSS_PENCIL::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> F,
 Array<complx,3> G,
 int n,
 Array<DP,2> local_Sk
 )
{
	Array<DP,3> Fr=Array<DP,3>(reinterpret_cast<DP*>(F.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<DP,3> Gr=Array<DP,3>(reinterpret_cast<DP*>(G.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	local_Sk = 0.0;
	
	DP Kmag, Kperp, Kpll;
	DP factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)* Fr(lx,ly,lz)*Gr(lx,ly,lz);
				}
			}
	
}



/**********************************************************************************************
 
 (B0.k) Im[vectA. conj(vectB)]
 
 ***********************************************************************************************/

void SSS_PENCIL::Compute_local_imag_shell_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,1> local_Sk
 )
{
	
	/*
	 Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	 Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	 Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	 
	 Array<DP,3> Bxr=Array<DP,3>(reinterpret_cast<DP*>(Bx.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	 Array<DP,3> Byr=Array<DP,3>(reinterpret_cast<DP*>(By.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	 Array<DP,3> Bzr=Array<DP,3>(reinterpret_cast<DP*>(Bz.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	 
	 DP Kmag;
	 int shell_index;
	 TinyVector<DP,3> K;
	 TinyVector<DP,3> V, W;
	 
	 
	 local_Sk = 0.0;
	 
	 int	Kmax = Min_radius_outside();
	 
	 for (int ly=0; ly<local_Ny; ly++)
	 for (int lz=0; lz<local_Nz; lz++)
	 for (int lx=0; lx<local_Nx; lx++) {
	 Kmag = Kmagnitude(lx, ly, lz);
	 
	 shell_index = (int) ceil(Kmag);
	 
	 if (shell_index <= Kmax) {
	 Wavenumber(lx, ly, lz, K);
	 V = Axr(lx,ly,lz), Ayr(lx,ly,lz), Azr(lx,ly,lz);
	 W = Bxr(lx,ly,lz), Byr(lx,ly,lz), Bzr(lx,ly,lz);
	 
	 local_Sk(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
	 
	 // factor multiplied by 2. recall the defn of energy spectrum and factor
	 }
	 }*/
	
	// To work
	
}



//******************************** for rings **************************************************

void SSS_PENCIL::Compute_local_imag_ring_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,2> local_Sk
 )
{
	/*	DP Kmag, theta;
	 int shell_index, sector_index;
	 TinyVector<DP,3> K;
	 TinyVector<complx,3> V, W;
	 
	 local_Sk = 0.0;
	 
	 int	Kmax = Max_radius_inside();
	 
	 for (int ly=0; ly<Ax.extent(0); ly++)
	 for (int lz=0; lz<Ax.extent(1); lz++)
	 for (int lx=0; lx<Ax.extent(2); lx++) {
	 Kmag = Kmagnitude(lx, ly, lz);
	 shell_index = (int) ceil(Kmag);
	 
	 if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
	 theta = AnisKvect_polar_angle(lx, ly, lz);
	 
	 sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
	 
	 Wavenumber(lx, ly, lz, K);
	 V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
	 W = Bx(lx,ly,lz), By(lx,ly,lz), Bz(lx,ly,lz);
	 
	 local_Sk(shell_index, sector_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
	 }
	 }
	 */
	
	
}

//******************************** for cylinder rings *****************************************

// Not for 2D
void SSS_PENCIL::Compute_local_imag_cylindrical_ring_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,2> local_Sk
 )
{
	/*
	 DP Kpll, Kperp;
	 int shell_index, slab_index;
	 TinyVector<DP,3> K;
	 TinyVector<complx,3> V, W;
	 
	 local_Sk = 0.0;
	 
	 int	Kperp_max = Anis_max_Krho_radius_inside();
	 
	 for (int ly=0; ly<Ax.extent(0); ly++)
	 for (int lz=0; lz<Ax.extent(1); lz++)
	 for (int lx=0; lx<Ax.extent(2); lx++) {
	 Kpll = AnisKpll(lx, ly, lz);
	 
	 Kperp = AnisKperp(lx, ly, lz);
	 
	 shell_index = (int) ceil(Kperp);
	 
	 if (shell_index <= Kperp_max) {
	 slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
	 
	 Wavenumber(lx, ly, lz, K);
	 V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
	 W = Bx(lx,ly,lz), By(lx,ly,lz), Bz(lx,ly,lz);
	 
	 local_Sk(shell_index, slab_index) +=  2* Multiplicity_factor(lx, ly, lz)* dot(B0,K)* mydot_imag(V,W);
	 
	 }
	 }
	 */
}





//*****************************  End of scft_energy.cc ****************************************













