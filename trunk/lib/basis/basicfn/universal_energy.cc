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

#include "universal.h"
#include "Global_extern_vars.h"
#include "basicfn_inline.h"
#include "shell_etc_indices.h"


/**********************************************************************************************

	-- Number of modes included inside shell (inner_radius, outer_radius]
	-- also count the complex conj
	
***********************************************************************************************/


DP Universal::Get_total_energy_real_space(Array<DP,3> Ar)
{
	
	DP local_total = Get_local_energy_real_space(Ar);
	
	DP total=0;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total;
}



DP Universal::Get_total_energy(Array<complx,3> A)
{
	DP local_total = Get_local_energy(A);
	
	DP total=0;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total;
}

DP Universal::Get_total_energy_real_space(Array<DP,3> Ar, Array<DP,3> Br)
{
	
	DP local_total = Get_local_energy_real_space(Ar, Br);
	
	DP total=0;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total;
	
}


DP Universal::Get_total_energy(Array<complx,3> A, Array<complx,3> B)
{
	
	DP local_total = Get_local_energy(A, B);
	
	DP total=0;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total;
	
}

/**********************************************************************************************
 
 Computes sum[ (1/2) K^n |A(k)|^2];   k=0 excluded
 
 ***********************************************************************************************/


//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
DP Universal::Get_local_Sn(Array<complx,3> A, DP n)
{
	DP Sn = 0.0;
	DP Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	
	for (int ly=0; ly<A.extent(0); ly++)
        for (int lz=0; lz<A.extent(1); lz++)
            for (int lx=0; lx<A.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz) * my_pow(Kmag,n)* Vsqr(A(ly,lz,lx));
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if ((my_id == master_id) && (n==0))
		Sn -= Multiplicity_factor(0, 0, 0) * pow2(abs(A(0,0,0)));
	
	return Sn;
}

//

DP Universal::Get_total_Sn(Array<complx,3> A, DP n)
{
	DP local_total;
	local_total = Get_local_Sn(A, n);
	
	DP total=0;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total;
}

//**********************************************************************************************
//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
DP Universal::Get_local_Sn(Array<complx,3> A, Array<complx,3> B, DP n)
{
	DP Sn = 0.0;
	DP Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	
	for (int ly=0; ly<A.extent(0); ly++)
        for (int lz=0; lz<A.extent(1); lz++)
            for (int lx=0; lx<A.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz) *my_pow(Kmag,n)* real_cprod(A(ly,lz,lx), B(ly,lz,lx));
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if ((my_id == master_id) && (n==0))
		Sn -= Multiplicity_factor(0, 0, 0) * real(A(0,0,0)*B(0,0,0));
	
	return Sn;
}

//

DP Universal::Get_total_Sn(Array<complx,3> A, Array<complx,3> B, DP n)
{
	DP local_total;
	local_total = Get_local_Sn(A, B, n);
	
	DP total=0;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total;
}



/**********************************************************************************************
 
 Computes Sk(K) = sum[ K'^n|A(K')|^2/2 ] with K-1 < K'<= K.
 Sk(0) for n=0 is the mean energy.
 
 Sk(A, B, k) = sum[ K'^n Real(A(K')*conj(B(K'))/2 ]
 
 For the partial shell: avg* volume of the shell.
 The volume is 2*pi*K^2.
 
 ***********************************************************************************************/


void Universal::Compute_local_shell_spectrum
(
 Array<complx,3> A,
 int n,
 Array<DP,1> local_Sk,
 Array<DP,1> local_Sk_count
 )
{
	local_Sk_count = 0.0;
	local_Sk = 0.0;
	
	DP Kmag;
	int index;
	DP factor;
	
    int	Kmax = Min_radius_outside();
	
    for (int ly=0; ly<A.extent(0); ly++)
        for (int lz=0; lz<A.extent(1); lz++)
            for (int lx=0; lx<A.extent(2); lx++) {
                Kmag = Kmagnitude(lx, ly, lz);
                index = (int) ceil(Kmag);
                
                if (index <= Kmax) {
                    factor = Multiplicity_factor(lx,ly,lz);
					
                    local_Sk(index) += factor* my_pow(Kmag,n)* Vsqr(A(ly,lz,lx));
                    local_Sk_count(index) += 2*factor;
                    // Mult by 2 because factor is offset by 2 in Multiply_factor
                }
            }
}


void Universal::Compute_shell_spectrum
(
 Array<complx,3> A,
 int n,
 Array<DP,1> Sk
 )
{
	static Array<DP,1> local_Sk(Sk.length());
	static Array<DP,1> local_Sk_count(Sk.length());
	
	local_Sk_count = 0.0;
	local_Sk = 0.0;
	
	Compute_local_shell_spectrum(A, n, local_Sk, local_Sk_count);
	
	static Array<DP,1> Sk_count(Sk.length());
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk_count.data()), reinterpret_cast<DP*>(Sk_count.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	// For semi-filled shells
	
	if (my_id == master_id)  {
		int Kmax_inside = Max_radius_inside();
		int	Kmax = Min_radius_outside();
		
		for (int index = Kmax_inside+1; index <= Kmax; index++)
			if (Sk_count(index) >= 1)
				Sk(index) = Sk(index) * Approx_number_modes_in_shell(index)/ Sk_count(index);
	}
}

//*********************************************************************************************
//
// Re(A. conj(B))
//


void Universal::Compute_local_shell_spectrum
(
 Array<complx,3> A,  Array<complx,3> B,
 int n,
 Array<DP,1> local_Sk,
 Array<DP,1> local_Sk_count
 )
{
	local_Sk_count=0;
	local_Sk = 0.0;
	
	DP Kmag;
	int index;
	DP factor;
	
	int	Kmax = Min_radius_outside();
	
	for (int ly=0; ly<A.extent(0); ly++)
        for (int lz=0; lz<A.extent(1); lz++)
            for (int lx=0; lx<A.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) {
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(index) += factor* my_pow(Kmag,n)*real_cprod(A(ly,lz,lx),B(ly,lz,lx));
					
					local_Sk_count(index) += 2*factor;
				}
			}
	
}

//
//
void Universal::Compute_shell_spectrum
(
 Array<complx,3> A, Array<complx,3> B,
 int n,
 Array<DP,1> Sk
 )
{
	static Array<DP,1> local_Sk(Sk.length());
	static Array<DP,1> local_Sk_count(Sk.length());
	
	local_Sk_count = 0.0;
	local_Sk = 0.0;
	Compute_local_shell_spectrum(A, B, n, local_Sk, local_Sk_count);
	
	static Array<DP,1> Sk_count(Sk.length());
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk_count.data()), reinterpret_cast<DP*>(Sk_count.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	// For semi-filled shells
	
	if (my_id == master_id) {
		int Kmax_inside = Max_radius_inside();
		int	Kmax = Min_radius_outside();
		
		for (int index =  Kmax_inside+1; index <= Kmax; index++)
			if (Sk_count(index) >= 1)
				Sk(index) = Sk(index) *Approx_number_modes_in_shell(index)/ Sk_count(index);
	} 
}

//**************************************************************************************
/**********************************************************************************************
 
 Compute entropy = sum(pi log(1/pi) where pi = Ei/Total E
 
 Skip the mean mode k= (0,0,0)
 
 ***********************************************************************************************/



DP Universal::Get_local_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
	
	DP modal_energy, prob, local_entropy;
	
	DP total_energy = Get_total_energy(Ax) + Get_total_energy(Ay) + Get_total_energy(Az);
	
	MPI_Bcast( &total_energy, 1, MPI_DP, master_id, MPI_COMM_WORLD);
	
	local_entropy = 0.0;
	
    for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
                
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx));
				
				prob = modal_energy/ total_energy;
				
				if (prob > MYEPS)
					local_entropy += (-prob * log(prob)/log(2.0));
			}
	
	return local_entropy;
}


DP Universal::Get_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
	
	DP local_entropy = Get_local_entropy(Ax, Ay, Az);
	
	DP entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	if (my_id == master_id)
		return entropy;
	
	else
		return 0;
}


//
// Scalar
//


DP Universal::Get_local_entropy_scalar(Array<complx,3> A)
{
	
	DP modal_energy, prob;
	
	DP local_entropy = 0.0;
	
	DP total_energy = Get_total_energy(A);
	
	MPI_Bcast( &total_energy, 1, MPI_DP, master_id, MPI_COMM_WORLD);
	
	for (int ly=0; ly<A.extent(0); ly++)
        for (int lz=0; lz<A.extent(1); lz++)
            for (int lx=0; lx<A.extent(2); lx++) {
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(A(ly,lz,lx));
				
				prob = modal_energy / total_energy;
				
				if (prob > MYEPS)
					local_entropy += (-prob*log(prob)/log(2.0));
			}
	
	return local_entropy;
	
}


DP Universal::Get_entropy_scalar(Array<complx,3> A)
{
	
	DP local_entropy = Get_local_entropy_scalar(A);
	
	DP entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	if (my_id == master_id)
		return entropy;
	
	else
		return 0;
}




/**********************************************************************************************
 
 Compute ring spectrum: Craya basis.
 2D: Sk(k, theta)
 3d: Sk1(k, theta) along horzontal direction; S2k(k, theta)
 in the perpendicular dirn (polar).
 
 ***********************************************************************************************/


void Universal::Compute_local_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 int n,
 Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
 )
{
	local_E1k = 0.0;
	local_E2k = 0.0;
	local_E3k = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx,ly,lz);
					local_E1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(Ax(ly,lz,lx));
					local_E2k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(Ay(ly,lz,lx));
					local_E3k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(Az(ly,lz,lx));
				}
			}
	
}


//*********************************************************************************************

void Universal::Compute_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 int n,
 Array<DP,2> E1k, Array<DP,2> E2k, Array<DP,2> E3k
 )
{
	
	static Array<DP,2> local_E1k( (E1k.length())(0), (E1k.length())(1));
	static Array<DP,2> local_E2k( (E1k.length())(0), (E1k.length())(1));
	static Array<DP,2> local_E3k( (E1k.length())(0), (E1k.length())(1));
	
	Compute_local_ring_spectrum(Ax, Ay, Az, n,local_E1k, local_E2k, local_E3k);
	
	int data_size = (E1k.length())(0) * (E1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_E1k.data()), reinterpret_cast<DP*>(E1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_E2k.data()), reinterpret_cast<DP*>(E2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_E3k.data()), reinterpret_cast<DP*>(E3k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}



//*********************************************************************************************
//
//  A.B
//


void Universal::Compute_local_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 int n, Array<DP,2> local_E1k, Array<DP,2> local_E2k, Array<DP,2> local_E3k
 )
{
	
	local_E1k = 0.0;
	local_E2k = 0.0;
	local_E3k = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
					local_E1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(Ax(ly, lz, lx), Bx(ly, lz, lx));
					local_E2k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(Ay(ly, lz, lx), By(ly, lz, lx));
					local_E3k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(Az(ly, lz, lx), Bz(ly, lz, lx));
				}
			}
}


//*********************************************************************************************

void Universal::Compute_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 int n, Array<DP,2> E1k, Array<DP,2> E2k, Array<DP,2> E3k
 )
{
	
	static Array<DP,2> local_E1k( (E1k.length())(0), (E1k.length())(1));
	static Array<DP,2> local_E2k( (E1k.length())(0), (E1k.length())(1));
	static Array<DP,2> local_E3k( (E1k.length())(0), (E1k.length())(1));
	
	Compute_local_ring_spectrum(Ax, Ay, Az, Bx, By, Bz, n, local_E1k, local_E2k, local_E3k);
	
	int data_size = (E1k.length())(0) * (E1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_E1k.data()), reinterpret_cast<DP*>(E1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_E2k.data()), reinterpret_cast<DP*>(E2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_E3k.data()), reinterpret_cast<DP*>(E3k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//*********************************************************************************************
// scalar
//

void Universal::Compute_local_ring_spectrum
(
 Array<complx,3> F,
 int n,
 Array<DP,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int ly=0; ly<F.extent(0); ly++)
        for (int lz=0; lz<F.extent(1); lz++)
            for (int lx=0; lx<F.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)* Vsqr(F(ly,lz,lx));
				}
			}
}



void Universal::Compute_ring_spectrum
(
 Array<complx,3> F,
 int n,
 Array<DP,2> Sk
 )
{
	
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_ring_spectrum(F, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()),reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//
// F*G
//

void Universal::Compute_local_ring_spectrum
(
 Array<complx,3> F,
 Array<complx,3> G,
 int n,
 Array<DP,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int ly=0; ly<F.extent(0); ly++)
        for (int lz=0; lz<F.extent(1); lz++)
            for (int lx=0; lx<F.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)*real_cprod(F(ly,lz,lx),G(ly,lz,lx));
				}
			}
	
	
}


void Universal::Compute_ring_spectrum
(
 Array<complx,3> F,
 Array<complx,3> G,
 int n,
 Array<DP,2> Sk
 )
{
	
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_ring_spectrum(F, G, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}




/**********************************************************************************************
 
 CYLINDERICAL RING SPECTRUM
 
 Not for 2D
 
 ***********************************************************************************************/

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 int n,
 Array<DP,2> local_S1k, Array<DP,2> local_S2k
 )
{
	
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	complx anisV1;
	
	DP Kmag, Kpll, Kperp;
	DP V2sqr;
	DP factor;
	int shell_index, slab_index;
	TinyVector<complx,3> V;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				
				V = Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx);
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

//
//
void Universal::Compute_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 int n,
 Array<DP,2> S1k, Array<DP,2> S2k
 )
{
	
	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, n, local_S1k, local_S2k);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S1k.data()), reinterpret_cast<DP*>(S1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S2k.data()), reinterpret_cast<DP*>(S2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}

//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 int n,
 Array<DP,2> local_S1k, Array<DP,2> local_S2k
 )
{
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	complx anisV1, anisW1;
	TinyVector<complx,3> V, Vrho, W, Wrho;
	
	DP Kmag, Kpll, Kperp;
	DP factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					V = Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx);
					W = Bx(ly, lz, lx), By(ly, lz, lx), Bz(ly, lz, lx);
					
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					if (global.field.anisotropy_dirn == 1) {
						anisV1 = V(0); anisW1 = W(0);
						
						Vrho = complx(0.0,0.0), V(1), V(2);
						Wrho = complx(0.0,0.0), W(1), W(2);
					}
					
					if (global.field.anisotropy_dirn == 2) {
						anisV1 = V(1); anisW1 = W(1);
						
						Vrho = V(0), complx(0.0,0.0), V(2);
						Wrho = W(0), complx(0.0,0.0), W(2);
					}
					
					if (global.field.anisotropy_dirn == 3) {
						anisV1 = V(2); anisW1 = W(2);
						
						Vrho = V(0), V(1), complx(0.0,0.0);
						Wrho = W(0), W(1), complx(0.0,0.0);
					}
					
					local_S1k(shell_index, slab_index) += factor*my_pow(Kmag,n)* real_cprod(anisV1, anisW1);
					local_S2k(shell_index, slab_index) += factor*my_pow(Kmag,n)* mydot(Vrho,Wrho);
				}
			}
	
}



//
//
void Universal::Compute_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 int n,
 Array<DP,2> S1k, Array<DP,2> S2k
 )
{
	
	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, Bx, By, Bz, n, local_S1k, local_S2k);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S1k.data()), reinterpret_cast<DP*>(S1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S2k.data()), reinterpret_cast<DP*>(S2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> F,
 int n,
 Array<DP,2> local_S1k, Array<DP,2> local_S2k
 )
{
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	complx anisV1, Fval;
	
	DP Kmag, Kpll, Kperp;
	DP VTperp;
	DP factor;
	int shell_index, slab_index;
	TinyVector<complx,3> V, Vless;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				
				V = Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx);
				Fval = F(ly,lz,lx);
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					anisV1 = V(global.field.anisotropy_dirn-1);
					Vless = V - anisV1;
					
					local_S1k(shell_index, slab_index) += factor*my_pow(Kmag,n)* real_cprod(anisV1, Fval);
					
					VTperp = real_cprod(sum(Vless), Fval);
					local_S2k(shell_index, slab_index) += factor*my_pow(Kmag,n)*VTperp;
				}
			}

	
}



//
//
void Universal::Compute_cylindrical_ring_spectrum
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> F,
 int n,
 Array<DP,2> S1k, Array<DP,2> S2k
 )
{
	
	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, F, n, local_S1k, local_S2k);
	
	int data_size = S1k.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S1k.data()), reinterpret_cast<DP*>(S1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S2k.data()), reinterpret_cast<DP*>(S2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}




//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> F,
 int n,
 Array<DP,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	DP Kmag, Kperp, Kpll;
	DP factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<F.extent(0); ly++)
        for (int lz=0; lz<F.extent(1); lz++)
            for (int lx=0; lx<F.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)* Vsqr(F(ly,lz,lx));
				}
			}
	
}

//
//
void Universal::Compute_cylindrical_ring_spectrum
(
 Array<complx,3> F,
 int n,
 Array<DP,2> Sk
 )
{
	
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(F, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<complx,3> F,
 Array<complx,3> G,
 int n,
 Array<DP,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	DP Kmag, Kperp, Kpll;
	DP factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<F.extent(0); ly++)
        for (int lz=0; lz<F.extent(1); lz++)
            for (int lx=0; lx<F.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)*real_cprod(F(ly,lz,lx),G(ly,lz,lx));
				}
			}
	
}

//
//
void Universal::Compute_cylindrical_ring_spectrum
(
 Array<complx,3> F,
 Array<complx,3> G,
 int n,
 Array<DP,2> Sk
 )
{
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(F, G, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


/**********************************************************************************************
 
 (B0.k) Im[vectA. conj(vectB)]
 
 ***********************************************************************************************/

void Universal::Compute_local_imag_shell_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,1> local_Sk
 )
{
	DP Kmag;
	int shell_index;
	TinyVector<DP,3> K;
	TinyVector<complx,3> V, W;
	
	
	local_Sk = 0.0;
	
	int	Kmax = Min_radius_outside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				shell_index = (int) ceil(Kmag);
				
				if (shell_index <= Kmax) {
					Wavenumber(lx, ly, lz, K);
					V = Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx);
					W = Bx(ly, lz, lx), By(ly, lz, lx), Bz(ly, lz, lx);
					
					local_Sk(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
					
					// factor multiplied by 2. recall the defn of energy spectrum and factor
				}
			}
	
}


//
//
void Universal::Compute_imag_shell_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,1> Sk
 )
{
	
	static Array<DP,1> local_Sk(Sk.length());
	
	local_Sk = 0.0;
	
	Compute_local_imag_shell_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz,local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}

//******************************** for rings **************************************************

void Universal::Compute_local_imag_ring_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,2> local_Sk
 )
{
	DP Kmag, theta;
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
					V = Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx);
					W = Bx(ly, lz, lx), By(ly, lz, lx), Bz(ly, lz, lx);
					
					local_Sk(shell_index, sector_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
				}
			}
	
}

//
//
void Universal::Compute_imag_ring_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,2> Sk
 )
{
	
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	local_Sk = 0.0;
	
	Compute_local_imag_ring_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz,local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//******************************** for cylinder rings *****************************************

// Not for 2D
void Universal::Compute_local_imag_cylindrical_ring_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,2> local_Sk
 )
{
	
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
					V = Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx);
					W = Bx(ly, lz, lx), By(ly, lz, lx), Bz(ly, lz, lx);
					
					local_Sk(shell_index, slab_index) +=  2* Multiplicity_factor(lx, ly, lz)* dot(B0,K)* mydot_imag(V,W);
					
				}
			}
	
}

//
//
void Universal::Compute_imag_cylindrical_ring_spectrum_B0
(
 TinyVector<DP,3> B0,
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP,2> Sk
 )
{
	
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	local_Sk = 0.0;
	
	Compute_local_imag_cylindrical_ring_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz, local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}



//*****************************  End of four_basic.cc **************************








