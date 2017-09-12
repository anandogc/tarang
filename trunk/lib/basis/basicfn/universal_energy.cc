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


Real Universal::Get_total_energy_real_space(Array<Real,3> Ar)
{
	Real local_energy = Get_local_energy_real_space(Ar);
	Real total_energy=0;

	MPI_Reduce(&local_energy, &total_energy, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total_energy;
}



Real Universal::Get_total_energy(Array<Complex,3> A)
{
	Real local_energy = Get_local_energy(A);
	Real total_energy=0;

	MPI_Reduce(&local_energy, &total_energy, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total_energy;
}

Real Universal::Get_total_energy_real_space(Array<Real,3> Ar, Array<Real,3> Br)
{
	
	Real local_energy = Get_local_energy_real_space(Ar, Br);
	Real total_energy=0;

	MPI_Reduce(&local_energy, &total_energy, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total_energy;
	
}


Real Universal::Get_total_energy(Array<Complex,3> A, Array<Complex,3> B)
{
	
	Real local_energy = Get_local_energy(A, B);
	Real total_energy=0;

	MPI_Reduce(&local_energy, &total_energy, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total_energy;
	
}

/**********************************************************************************************
 
 Computes sum[ (1/2) K^n |A(k)|^2];   k=0 excluded
 
 ***********************************************************************************************/


//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
Real Universal::Get_local_Sn(Array<Complex,3> A, Real n)
{
	Real Sn = 0.0;
	Real Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz)*my_pow(Kmag,n)*Vsqr(A(lx,ly,lz));
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if (master && (n==0))
		Sn -= Multiplicity_factor(0,0,0)*pow2(abs(A(0,0,0)));
	
	return Sn;
}

//

Real Universal::Get_total_Sn(Array<Complex,3> A, Real n)
{
	Real local_Sn;
	local_Sn = Get_local_Sn(A, n);
	
	Real total_Sn=0;
	MPI_Reduce(&local_Sn, &total_Sn, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total_Sn;
}

//**********************************************************************************************
//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
Real Universal::Get_local_Sn(Array<Complex,3> A, Array<Complex,3> B, Real n)
{
	Real Sn = 0.0;
	Real Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz)*my_pow(Kmag,n)*real_cprod(A(lx,ly,lz), B(lx,ly,lz));
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if (master && (n==0))
		Sn -= Multiplicity_factor(0,0,0) * real(A(0,0,0)*B(0,0,0));
	
	return Sn;
}

//

Real Universal::Get_total_Sn(Array<Complex,3> A, Array<Complex,3> B, Real n)
{
	Real local_Sn;
	local_Sn = Get_local_Sn(A, B, n);
	
	Real total_Sn=0;
	MPI_Reduce(&local_Sn, &total_Sn, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	return total_Sn;
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
 Array<Complex,3> A,
 int n,
 Array<Real,1> local_Sk,
 Array<Real,1> local_Sk_count
 )
{
	local_Sk_count = 0.0;
	local_Sk = 0.0;
	
	Real Kmag;
	int index;
	Real factor;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) {
					factor = Multiplicity_factor(lx,ly,lz);
					
                    //local_Sk(index) += factor*my_pow(Kmag,n)*Vsqr(A(lx,ly,lz));//old spectrum formula
                    //Abhishek implemented Stepanov et al. PRE 2014
                    if (Ny==1)
                      local_Sk(index) += factor*my_pow(Kmag,n)*Vsqr(A(lx,ly,lz))*(Kmag);
                    else
                      local_Sk(index) += factor*my_pow(Kmag,n)*Vsqr(A(lx,ly,lz))*pow2(Kmag);

                  local_Sk_count(index) += 2*factor;
					// Mult by 2 because factor is offset by 2 in Multiply_factor
				}
			}
}


void Universal::Compute_shell_spectrum
(
 Array<Complex,3> A,
 int n,
 Array<Real,1> Sk
 )
{
	static Array<Real,1> local_Sk(Sk.shape());
	static Array<Real,1> local_Sk_count(Sk.shape());
	
	local_Sk_count = 0.0;
	local_Sk = 0.0;
	
	Compute_local_shell_spectrum(A, n, local_Sk, local_Sk_count);
	
	static Array<Real,1> Sk_count(Sk.shape());
	
	int data_size = Sk.size();
	
	MPI_Allreduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, MPI_COMM_WORLD);
	
	MPI_Allreduce(reinterpret_cast<Real*>(local_Sk_count.data()), reinterpret_cast<Real*>(Sk_count.data()), data_size, MPI_Real, MPI_SUM, MPI_COMM_WORLD);
	
    //Abhishek implemented Stepanov et al. PRE 2014
    int dimension_factor;
    if (Ny==1)
      dimension_factor = 2;
    else
      dimension_factor = 4;
    for (int i=0; i<=Sk_count.extent(0); i++){
      if (Sk_count(i) <= 0)
        Sk(i) = 0;
      else
        Sk(i) = Sk(i)*(dimension_factor*M_PI)/(Sk_count(i)*(kfactor[1]*kfactor[2]*kfactor[3]));
    }
  
  // For semi-filled shells

}

//*********************************************************************************************
//
// Re(A. conj(B))
//


void Universal::Compute_local_shell_spectrum
(
 Array<Complex,3> A,  Array<Complex,3> B,
 int n,
 Array<Real,1> local_Sk,
 Array<Real,1> local_Sk_count
 )
{
	local_Sk_count=0;
	local_Sk = 0.0;
	
	Real Kmag;
	int index;
	Real factor;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) {
					factor = Multiplicity_factor(lx,ly,lz);
					
					local_Sk(index) += factor* my_pow(Kmag,n)*real_cprod(A(lx,ly,lz),B(lx,ly,lz));
					
					local_Sk_count(index) += 2*factor;
				}
			}
	
}

//
//
void Universal::Compute_shell_spectrum
(
 Array<Complex,3> A, Array<Complex,3> B,
 int n,
 Array<Real,1> Sk
 )
{
	static Array<Real,1> local_Sk(Sk.shape());
	static Array<Real,1> local_Sk_count(Sk.shape());
	
	local_Sk_count = 0.0;
	local_Sk = 0.0;
	Compute_local_shell_spectrum(A, B, n, local_Sk, local_Sk_count);
	
	static Array<Real,1> Sk_count(Sk.shape());
	
	int data_size = Sk.size();
	
	MPI_Allreduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, MPI_COMM_WORLD);
	
	MPI_Allreduce(reinterpret_cast<Real*>(local_Sk_count.data()), reinterpret_cast<Real*>(Sk_count.data()), data_size, MPI_Real, MPI_SUM, MPI_COMM_WORLD);
	
	// For semi-filled shells
	
	if (master) {
		int Kmax_inside = Max_radius_inside();
		int	Kmax = Min_radius_outside();
		
		for (int index =  Kmax_inside+1; index <= Kmax; index++)
			if (Sk_count(index) >= 1)
				Sk(index) = Sk(index)*Approx_number_modes_in_shell(index)/Sk_count(index);
	} 
}

//**************************************************************************************
/**********************************************************************************************
 
 Compute entropy = sum(pi log(1/pi) where pi = Ei/Total E
 
 Skip the mean mode k= (0,0,0)
 
 ***********************************************************************************************/



Real Universal::Get_local_entropy(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	
	Real modal_energy, prob, local_entropy;
	
	Real total_energy = Get_total_energy(Ax) + Get_total_energy(Ay) + Get_total_energy(Az);
	
	MPI_Bcast( &total_energy, 1, MPI_Real, master_id, MPI_COMM_WORLD);
	
	local_entropy = 0.0;
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz));
				
				prob = modal_energy/ total_energy;
				
				if (prob > MYEPS)
					local_entropy += (-prob * log(prob)/log(TWO));
			}
	
	return local_entropy;
}


Real Universal::Get_entropy(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	
	Real local_entropy = Get_local_entropy(Ax, Ay, Az);
	
	Real entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	if (master)
		return entropy;
	
	else
		return 0;
}


//
// Scalar
//


Real Universal::Get_local_entropy_scalar(Array<Complex,3> A)
{
	
	Real modal_energy, prob;
	
	Real local_entropy = 0.0;
	
	Real total_energy = Get_total_energy(A);
	
	MPI_Bcast( &total_energy, 1, MPI_Real, master_id, MPI_COMM_WORLD);
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(A(lx,ly,lz));
				
				prob = modal_energy / total_energy;
				
				if (prob > MYEPS)
					local_entropy += (-prob*log(prob)/log(TWO));
			}
	
	return local_entropy;
	
}


Real Universal::Get_entropy_scalar(Array<Complex,3> A)
{
	
	Real local_entropy = Get_local_entropy_scalar(A);
	
	Real entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	if (master)
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
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 int n,
 Array<Real,2> local_E1k, Array<Real,2> local_E2k, Array<Real,2> local_E3k
 )
{
	local_E1k = 0.0;
	local_E2k = 0.0;
	local_E3k = 0.0;
	
	Real Kmag, theta;
	Real factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				
				Kmag = Kmagnitude(lx,ly,lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					
					theta = AnisKvect_polar_angle(lx,ly,lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx,ly,lz);
					local_E1k(shell_index, sector_index) += factor*my_pow(Kmag,n)*Vsqr(Ax(lx,ly,lz));
					local_E2k(shell_index, sector_index) += factor*my_pow(Kmag,n)*Vsqr(Ay(lx,ly,lz));
					local_E3k(shell_index, sector_index) += factor*my_pow(Kmag,n)*Vsqr(Az(lx,ly,lz));
				}
			}
	
}


//*********************************************************************************************

void Universal::Compute_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 int n,
 Array<Real,2> E1k, Array<Real,2> E2k, Array<Real,2> E3k
 )
{
	
	static Array<Real,2> local_E1k(E1k.shape());
	static Array<Real,2> local_E2k(E1k.shape());
	static Array<Real,2> local_E3k(E1k.shape());
	
	Compute_local_ring_spectrum(Ax, Ay, Az, n,local_E1k, local_E2k, local_E3k);
	
	int data_size = E1k.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_E1k.data()), reinterpret_cast<Real*>(E1k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_E2k.data()), reinterpret_cast<Real*>(E2k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_E3k.data()), reinterpret_cast<Real*>(E3k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}



void Universal::Compute_local_helical_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> helicalAx, Array<Complex,3> helicalAy, Array<Complex,3> helicalAz,
 int n,
 Array<Real,2> local_hk1, Array<Real,2> local_hk2, Array<Real,2> local_hk3
 )
{
	local_hk1 = 0.0;
	local_hk2 = 0.0;
	local_hk3 = 0.0;
	
	Real Kmag, theta;
	Real factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				
				Kmag = Kmagnitude(lx,ly,lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					
					theta = AnisKvect_polar_angle(lx,ly,lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx,ly,lz);
					local_hk1(shell_index, sector_index) += (1/4.0)*factor*my_pow(Kmag,n)*real(conj(Ax(lx,ly,lz))*helicalAx(lx,ly,lz) + conj(Ax(lx,ly,lz))*helicalAx(lx,ly,lz));
					local_hk2(shell_index, sector_index) += (1/4.0)*factor*my_pow(Kmag,n)*real(conj(Ay(lx,ly,lz))*helicalAy(lx,ly,lz) + conj(Ay(lx,ly,lz))*helicalAy(lx,ly,lz));
					local_hk3(shell_index, sector_index) += (1/4.0)*factor*my_pow(Kmag,n)*real(conj(Az(lx,ly,lz))*helicalAz(lx,ly,lz) + conj(Az(lx,ly,lz))*helicalAz(lx,ly,lz));
				}
			}
	
}


void Universal::Compute_helical_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> helicalAx, Array<Complex,3> helicalAy, Array<Complex,3> helicalAz,
 int n,
 Array<Real,2> hk1, Array<Real,2> hk2, Array<Real,2> hk3
 )
{
	
	static Array<Real,2> local_hk1(hk1.shape());
	static Array<Real,2> local_hk2(hk2.shape());
	static Array<Real,2> local_hk3(hk3.shape());
	
	Compute_local_helical_ring_spectrum(Ax, Ay, Az, helicalAx, helicalAy, helicalAz, n, local_hk1, local_hk2, local_hk3);
	
	int data_size = hk1.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_hk1.data()), reinterpret_cast<Real*>(hk1.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_hk2.data()), reinterpret_cast<Real*>(hk2.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_hk3.data()), reinterpret_cast<Real*>(hk3.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}

//*********************************************************************************************
//
//  A.B
//


void Universal::Compute_local_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 int n, Array<Real,2> local_E1k, Array<Real,2> local_E2k, Array<Real,2> local_E3k
 )
{
	
	local_E1k = 0.0;
	local_E2k = 0.0;
	local_E3k = 0.0;
	
	Real Kmag, theta;
	Real factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				
				Kmag = Kmagnitude(lx,ly,lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx,ly,lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx,ly,lz);
					local_E1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(Ax(lx,ly,lz), Bx(lx,ly,lz));
					local_E2k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(Ay(lx,ly,lz), By(lx,ly,lz));
					local_E3k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(Az(lx,ly,lz), Bz(lx,ly,lz));
				}
			}
}


//*********************************************************************************************

void Universal::Compute_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 int n, Array<Real,2> E1k, Array<Real,2> E2k, Array<Real,2> E3k
 )
{
	
	static Array<Real,2> local_E1k(E1k.shape());
	static Array<Real,2> local_E2k(E1k.shape());
	static Array<Real,2> local_E3k(E1k.shape());
	
	Compute_local_ring_spectrum(Ax, Ay, Az, Bx, By, Bz, n, local_E1k, local_E2k, local_E3k);
	
	int data_size = E1k.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_E1k.data()), reinterpret_cast<Real*>(E1k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_E2k.data()), reinterpret_cast<Real*>(E2k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_E3k.data()), reinterpret_cast<Real*>(E3k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//*********************************************************************************************
// scalar
//

void Universal::Compute_local_ring_spectrum
(
 Array<Complex,3> F,
 int n,
 Array<Real,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	Real Kmag, theta;
	Real factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx,ly,lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx,ly,lz);
					
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)* Vsqr(F(lx,ly,lz));
				}
			}
}



void Universal::Compute_ring_spectrum
(
 Array<Complex,3> F,
 int n,
 Array<Real,2> Sk
 )
{
	
	static Array<Real,2> local_Sk(Sk.shape());
	
	Compute_local_ring_spectrum(F, n, local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_Sk.data()),reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//
// F*G
//

void Universal::Compute_local_ring_spectrum
(
 Array<Complex,3> F,
 Array<Complex,3> G,
 int n,
 Array<Real,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	Real Kmag, theta;
	Real factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx,ly,lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx,ly,lz);
					
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)*real_cprod(F(lx,ly,lz),G(lx,ly,lz));
				}
			}
	
	
}


void Universal::Compute_ring_spectrum
(
 Array<Complex,3> F,
 Array<Complex,3> G,
 int n,
 Array<Real,2> Sk
 )
{
	
	static Array<Real,2> local_Sk(Sk.shape());
	
	Compute_local_ring_spectrum(F, G, n, local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}




/**********************************************************************************************
 
 CYLINDERICAL RING SPECTRUM
 
 Not for 2D
 
 ***********************************************************************************************/

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 int n,
 Array<Real,2> local_S1k, Array<Real,2> local_S2k
 )
{
	
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	Complex anisV1;
	
	Real Kmag, Kpll, Kperp;
	Real V2sqr;
	Real factor;
	int shell_index, slab_index;
	TinyVector<Complex,3> V;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				
				V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
				Kmag = Kmagnitude(lx,ly,lz);
				
				Kperp = AnisKperp(lx,ly,lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx,ly,lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx,ly,lz);
					
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
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 int n,
 Array<Real,2> S1k, Array<Real,2> S2k
 )
{
	
	static Array<Real,2> local_S1k(S1k.shape());
	static Array<Real,2> local_S2k(S1k.shape());
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, n, local_S1k, local_S2k);
	
	int data_size = S1k.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_S1k.data()), reinterpret_cast<Real*>(S1k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_S2k.data()), reinterpret_cast<Real*>(S2k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}

//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 int n,
 Array<Real,2> local_S1k, Array<Real,2> local_S2k
 )
{
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	Complex anisV1, anisW1;
	TinyVector<Complex,3> V, Vrho, W, Wrho;
	
	Real Kmag, Kpll, Kperp;
	Real factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				
				Kmag = Kmagnitude(lx,ly,lz);
				Kperp = AnisKperp(lx,ly,lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
					W = Bx(lx,ly,lz), By(lx,ly,lz), Bz(lx,ly,lz);
					
					Kpll = AnisKpll(lx,ly,lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx,ly,lz);
					
					if (global.field.anisotropy_dirn == 1) {
						anisV1 = V(0); anisW1 = W(0);
						
						Vrho = Complex(0.0,0.0), V(1), V(2);
						Wrho = Complex(0.0,0.0), W(1), W(2);
					}
					
					if (global.field.anisotropy_dirn == 2) {
						anisV1 = V(1); anisW1 = W(1);
						
						Vrho = V(0), Complex(0.0,0.0), V(2);
						Wrho = W(0), Complex(0.0,0.0), W(2);
					}
					
					if (global.field.anisotropy_dirn == 3) {
						anisV1 = V(2); anisW1 = W(2);
						
						Vrho = V(0), V(1), Complex(0.0,0.0);
						Wrho = W(0), W(1), Complex(0.0,0.0);
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
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 int n,
 Array<Real,2> S1k, Array<Real,2> S2k
 )
{
	
	static Array<Real,2> local_S1k(S1k.shape());
	static Array<Real,2> local_S2k(S1k.shape());
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, Bx, By, Bz, n, local_S1k, local_S2k);
	
	int data_size = S1k.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_S1k.data()), reinterpret_cast<Real*>(S1k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_S2k.data()), reinterpret_cast<Real*>(S2k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> F,
 int n,
 Array<Real,2> local_S1k, Array<Real,2> local_S2k
 )
{
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;
	
	Complex anisV1, Fval;
	
	Real Kmag, Kpll, Kperp;
	Real VTperp;
	Real factor;
	int shell_index, slab_index;
	TinyVector<Complex,3> V, Vless;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				
				V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
				Fval = F(lx,ly,lz);
				Kmag = Kmagnitude(lx,ly,lz);
				
				Kperp = AnisKperp(lx,ly,lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx,ly,lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx,ly,lz);
					
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
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> F,
 int n,
 Array<Real,2> S1k, Array<Real,2> S2k
 )
{
	
	static Array<Real,2> local_S1k(S1k.shape());
	static Array<Real,2> local_S2k(S1k.shape());
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, F, n, local_S1k, local_S2k);
	
	int data_size = S1k.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_S1k.data()), reinterpret_cast<Real*>(S1k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<Real*>(local_S2k.data()), reinterpret_cast<Real*>(S2k.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}




//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<Complex,3> F,
 int n,
 Array<Real,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	Real Kmag, Kperp, Kpll;
	Real factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				
				Kperp = AnisKperp(lx,ly,lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx,ly,lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx,ly,lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)* Vsqr(F(lx,ly,lz));
				}
			}
	
}

//
//
void Universal::Compute_cylindrical_ring_spectrum
(
 Array<Complex,3> F,
 int n,
 Array<Real,2> Sk
 )
{
	
	static Array<Real,2> local_Sk(Sk.shape());
	
	Compute_local_cylindrical_ring_spectrum(F, n, local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//*********************************************************************************************

void Universal::Compute_local_cylindrical_ring_spectrum
(
 Array<Complex,3> F,
 Array<Complex,3> G,
 int n,
 Array<Real,2> local_Sk
 )
{
	
	local_Sk = 0.0;
	
	Real Kmag, Kperp, Kpll;
	Real factor;
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				
				Kperp = AnisKperp(lx,ly,lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx,ly,lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx,ly,lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)*real_cprod(F(lx,ly,lz),G(lx,ly,lz));
				}
			}
	
}

//
//
void Universal::Compute_cylindrical_ring_spectrum
(
 Array<Complex,3> F,
 Array<Complex,3> G,
 int n,
 Array<Real,2> Sk
 )
{
	static Array<Real,2> local_Sk(Sk.shape());
	
	Compute_local_cylindrical_ring_spectrum(F, G, n, local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


/**********************************************************************************************
 
 (B0.k) Im[vectA. conj(vectB)]
 
 ***********************************************************************************************/

void Universal::Compute_local_imag_shell_spectrum_B0
(
 TinyVector<Real,3> B0,
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 Array<Real,1> local_Sk
 )
{
	Real Kmag;
	int shell_index;
	TinyVector<Real,3> K;
	TinyVector<Complex,3> V, W;
	
	
	local_Sk = 0.0;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				
				shell_index = (int) ceil(Kmag);
				
				if (shell_index <= Kmax) {
					Wavenumber(lx,ly,lz, K);
					V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
					W = Bx(lx,ly,lz), By(lx,ly,lz), Bz(lx,ly,lz);
					
					local_Sk(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
					
					// factor multiplied by 2. recall the defn of energy spectrum and factor
				}
			}
	
}


//
//
void Universal::Compute_imag_shell_spectrum_B0
(
 TinyVector<Real,3> B0,
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 Array<Real,1> Sk
 )
{
	
	static Array<Real,1> local_Sk(Sk.shape());
	
	local_Sk = 0.0;
	
	Compute_local_imag_shell_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz,local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}

//******************************** for rings **************************************************

void Universal::Compute_local_imag_ring_spectrum_B0
(
 TinyVector<Real,3> B0,
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 Array<Real,2> local_Sk
 )
{
	Real Kmag, theta;
	int shell_index, sector_index;
	TinyVector<Real,3> K;
	TinyVector<Complex,3> V, W;
	
	local_Sk = 0.0;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kmag = Kmagnitude(lx,ly,lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx,ly,lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					Wavenumber(lx,ly,lz, K);
					V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
					W = Bx(lx,ly,lz), By(lx,ly,lz), Bz(lx,ly,lz);
					
					local_Sk(shell_index, sector_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
				}
			}
	
}

//
//
void Universal::Compute_imag_ring_spectrum_B0
(
 TinyVector<Real,3> B0,
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 Array<Real,2> Sk
 )
{
	
	static Array<Real,2> local_Sk(Sk.shape());
	
	local_Sk = 0.0;
	
	Compute_local_imag_ring_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz,local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//******************************** for cylinder rings *****************************************

// Not for 2D
void Universal::Compute_local_imag_cylindrical_ring_spectrum_B0
(
 TinyVector<Real,3> B0,
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 Array<Real,2> local_Sk
 )
{
	
	Real Kpll, Kperp;
	int shell_index, slab_index;
	TinyVector<Real,3> K;
	TinyVector<Complex,3> V, W;
	
	local_Sk = 0.0;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
			for (int lz=0; lz<maxlz; lz++) {
				Kpll = AnisKpll(lx,ly,lz);
				
				Kperp = AnisKperp(lx,ly,lz);
				
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					Wavenumber(lx,ly,lz, K);
					V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
					W = Bx(lx,ly,lz), By(lx,ly,lz), Bz(lx,ly,lz);
					
					local_Sk(shell_index, slab_index) +=  2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
					
				}
			}
	
}

//
//
void Universal::Compute_imag_cylindrical_ring_spectrum_B0
(
 TinyVector<Real,3> B0,
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
 Array<Real,2> Sk
 )
{
	
	static Array<Real,2> local_Sk(Sk.shape());
	
	local_Sk = 0.0;
	
	Compute_local_imag_cylindrical_ring_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz, local_Sk);
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<Real*>(local_Sk.data()), reinterpret_cast<Real*>(Sk.data()), data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}



//*****************************  End of four_basic.cc **************************








