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


/*! \file four_energy.cc
 * 
 * @sa four_energy.h
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 22/08/2008
 * @bug Line 1299: Helicity -- do we need factor 2?
 */  

#include "CFFF_slab.h"
#include "CFFF_slab_inline.h"
#include "shell_etc_indices.h"


/**********************************************************************************************

       		Computes A^2/2 & Re(A.B*)/2 except k=0

***********************************************************************************************/


DP CFFF_SLAB::Get_local_energy_XZ_plane(Array<complx,3> A, int ny)
{
	return 0;
}

DP CFFF_SLAB::Get_local_energy_real_space(Array<complx,3> Ar)
{	
	return Array_sqr(Ar);
}


DP CFFF_SLAB::Get_total_energy_real_space(Array<complx,3> Ar)
{
	
	DP local_total = Get_local_energy_real_space(Ar);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	
	
	if (my_id == master_id) 
		return total/(2* DP(Nx) * DP(Ny) * DP(Nz));
	
	else
		return 0;
}



DP CFFF_SLAB::Get_local_energy(Array<complx,3> A)
{
	
	DP total = Array_sqr(A);
	
	return total; 
	
}

//


DP CFFF_SLAB::Get_total_energy(Array<complx,3> A)
{

	DP local_total = Get_local_energy(A);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total - pow2(abs(A(0,0,0)))/2;
		
	else
		return 0;
}

//*********************************************************************************************

DP CFFF_SLAB::Get_local_energy_XZ_plane(Array<complx,3> A, Array<complx,3> B, int ny)
{
	return 0.0;
}

DP CFFF_SLAB::Get_local_energy(Array<complx,3> A, Array<complx,3> B)
{
	DP total = mydot(A, B);
	
	return total;
}


// 

DP CFFF_SLAB::Get_total_energy(Array<complx,3> A, Array<complx,3> B)
{

	DP local_total = Get_local_energy(A, B);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total;
		
	else
		return 0;	
}


/**********************************************************************************************

       		Computes sum[ (1/2) K^n |A(k)|^2];   k=0 excluded

***********************************************************************************************/


//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
DP CFFF_SLAB::Get_local_Sn(Array<complx,3> A, DP n)
{
/*	DP Sn = 0.0;
	DP Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	  
	for (int lx=0; lx<local_Nx; lx++)			
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz) * my_pow(Kmag,n)* Vsqr(A(lx,ly,lz));
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if ((my_id == master_id) && (n==0))
			Sn -= Multiplicity_factor(0, 0, 0) * pow2(abs(A(0,0,0)));		

	return Sn; */
}

//

DP CFFF_SLAB::Get_total_Sn(Array<complx,3> A, DP n)
{
/*	DP local_total;
	local_total = Get_local_Sn(A, n);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total;
		
	else
		return 0;	*/
}

//**********************************************************************************************
//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
DP CFFF_SLAB::Get_local_Sn(Array<complx,3> A, Array<complx,3> B, DP n)
{
/*	DP Sn = 0.0;
	DP Kmag;	// Kmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)			
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax)
					Sn += Multiplicity_factor(lx,ly,lz) *my_pow(Kmag,n)* real_cprod(A(lx,ly,lz), B(lx,ly,lz));
			}
	
	// The above sum adds for k=0 mode for n=0 since my_pow(0,0) = 1.
	// Subtract the contribution from the origin only for n=0.
	if ((my_id == master_id) && (n==0))
		Sn -= Multiplicity_factor(0, 0, 0) * real(A(0,0,0)*B(0,0,0));
	
	return Sn; */
}

//

DP CFFF_SLAB::Get_total_Sn(Array<complx,3> A, Array<complx,3> B, DP n)
{
/*	DP local_total;
	local_total = Get_local_Sn(A, B, n);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	
	
	if (my_id == master_id) 
		return total;
	
	else
		return 0;	*/
}



/**********************************************************************************************
 
		Computes Sk(K) = sum[ K'^n|A(K')|^2/2 ] with K-1 < K'<= K.
             Sk(0) for n=0 is the mean energy.
			  
			  Sk(A, B, k) = sum[ K'^n Real(A(K')*conj(B(K'))/2 ] 
		
			  For the partial shell: avg* volume of the shell.
			  The volume is 2*pi*K^2.

***********************************************************************************************/


void CFFF_SLAB::Compute_local_shell_spectrum
(  
	Array<complx,3> A, 
	int n, 
	Array<DP,1> local_Sk, 
	Array<DP,1> local_Sk_count
)
{
/*	local_Sk_count = 0.0;
	local_Sk = 0.0;

	DP Kmag;													
	int index;
	DP factor;
	
    int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) {
					factor = Multiplicity_factor(lx,ly,lz);
				
					local_Sk(index) += factor* my_pow(Kmag,n)* Vsqr(A(lx,ly,lz));
					local_Sk_count(index) += 2*factor;
					// Mult by 2 because factor is offset by 2 in Multiply_factor 
				}
			} */
}


void CFFF_SLAB::Compute_shell_spectrum
(  
    Array<complx,3> A, 
    int n, 
    Array<DP,1> Sk
)
{


/*	static Array<DP,1> local_Sk(Sk.length());
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
	} */
}
			

//*********************************************************************************************			
//
// Re(A. conj(B))
//


void CFFF_SLAB::Compute_local_shell_spectrum
(
	Array<complx,3> A,  Array<complx,3> B,
	int n, 
	Array<DP,1> local_Sk, 
	Array<DP,1> local_Sk_count
)
{
/*	local_Sk_count=0;
	local_Sk = 0.0;		
	
	DP Kmag;																
	int index;
	DP factor;

	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) {
					factor = Multiplicity_factor(lx, ly, lz);
				
					local_Sk(index) += factor* my_pow(Kmag,n)*real_cprod(A(lx,ly,lz),B(lx,ly,lz));

					local_Sk_count(index) += 2*factor;
				}	
			} */

}

//
//
void CFFF_SLAB::Compute_shell_spectrum
(  
	Array<complx,3> A, Array<complx,3> B, 
	int n, 
	Array<DP,1> Sk
)
{
/*	static Array<DP,1> local_Sk(Sk.length());
	static Array<DP,1> local_Sk_count(Sk.length());
	
	local_Sk_count = 0.0; 
	local_Sk = 0.0;	
	Compute_local_shell_spectrum(A, B, n, local_Sk, local_Sk_count);
	
	static Array<DP,1> Sk_count(Sk.length());
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), 
						reinterpret_cast<DP*>(Sk.data()),
						data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk_count.data()), 
						reinterpret_cast<DP*>(Sk_count.data()), 
						data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 

	// For semi-filled shells
	
	if (my_id == master_id) {
		int Kmax_inside = Max_radius_inside(); 	
		int	Kmax = Min_radius_outside();

		for (int index =  Kmax_inside+1; index <= Kmax; index++) 
			if (Sk_count(index) >= 1) 
				Sk(index) = Sk(index) *Approx_number_modes_in_shell(index)/ Sk_count(index); 
	} */
}



/**********************************************************************************************

	Compute Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2
	Multiplication factor = 1 for kz>0 because of energy spectrum details.
 
	not for 2D

***********************************************************************************************/



void CFFF_SLAB::Compute_local_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &local_helicity1, DP &local_helicity2, 
	DP &local_k2H1, DP &local_k2H2
)
{

/*	DP modal_helicity, Kmag, Ksqr, factor;

	local_helicity1 = local_helicity2 = 0.0;
	local_k2H1 =  0.0;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)
		for (int ly=0; ly<Ny; ly++)
			for (int lz=0; lz<=Nz/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax) {
					factor = 2*Multiplicity_factor(lx, ly, lz);
					// factor multiplied by 2 because of the defn  Hk = K . (Vr x Vi).
					// recall the defn of energy spectrum that contains 1/2.
					
					modal_helicity = factor*Get_Modal_helicity(lx, ly, lz, Ax, Ay, Az);
					local_helicity1 += modal_helicity;
					
                    Ksqr = pow2(Kmag);
					if (Ksqr > MYEPS)
						local_helicity2 += modal_helicity / Ksqr;
						
					local_k2H1 += Ksqr * modal_helicity;
				}	
			}	
			
	local_k2H1 *= 2.0;
	local_k2H2 = 2.0 * local_helicity1;		*/
}

//

void CFFF_SLAB::Compute_total_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &total_k2H1, DP &total_k2H2
)
{
/*	DP local_helicity1, local_helicity2;
	DP local_k2H1, local_k2H2;
	
	Compute_local_helicity(Ax, Ay, Az,local_helicity1, local_helicity2,local_k2H1, local_k2H2);
									
	MPI_Reduce(&local_helicity1, &total_helicity1, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
								
	MPI_Reduce(&local_helicity2, &total_helicity2, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
								
	MPI_Reduce(&local_k2H1, &total_k2H1, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
									
	MPI_Reduce(&local_k2H2, &total_k2H2, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);*/
}


/**********************************************************************************************

	Compute helicity spectrum
	Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2
 
	Not for 2D

***********************************************************************************************/



void CFFF_SLAB::Compute_local_shell_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3, 
	Array<DP,1> local_H1k_count
)												
{

/*	local_H1k_count = 0.0;
	local_H1k1 = 0.0;	
	local_H1k2 = 0.0;
	local_H1k3 = 0.0;
	
	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, K;
	DP Kmag;															
	int index;
	DP factor;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax) 	{
					factor = 2*Multiplicity_factor(lx, ly, lz);
					// factor multiplied by 2 because of the defn  Hk = K . (Vr x Vi).
					// recall the defn of energy spectrum that contains 1/2.
					
					Vreal = real(Ax(lx, ly, lz)), real(Ay(lx, ly, lz)), real(Az(lx, ly, lz));
					Vimag = imag(Ax(lx, ly, lz)), imag(Ay(lx, ly, lz)), imag(Az(lx, ly, lz));
			
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber(lx, ly, lz, K);
					
					// modal_helicity = factor * dot(K, VrcrossVi);	
					local_H1k1(index) += factor* (K(0)*VrcrossVi(0));
					local_H1k2(index) += factor* (K(1)*VrcrossVi(1));
					local_H1k3(index) += factor* (K(2)*VrcrossVi(2));
					
					local_H1k_count(index) = local_H1k_count(index) + 2*factor;
				}	
			}  */

}


void CFFF_SLAB::Compute_shell_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3 
)	
{
/*
	static Array<DP,1> local_H1k1(H1k1.length());
	static Array<DP,1> local_H1k2(H1k1.length());
	static Array<DP,1> local_H1k3(H1k1.length());
	
	static Array<DP,1> local_H1k_count(H1k1.length());
	
	local_H1k1 = 0.0;
	local_H1k2 = 0.0;
	local_H1k3 = 0.0;
	
	local_H1k_count = 0.0; 
		
	
	Compute_local_shell_spectrum_helicity(Ax, Ay, Az, local_H1k1, local_H1k2, local_H1k3, local_H1k_count);
	
	static Array<DP,1> H1k1_count(H1k1.length());	
	int data_size = H1k1.size();
				
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k1.data()), reinterpret_cast<DP*>(H1k1.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
					
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k2.data()), reinterpret_cast<DP*>(H1k2.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k3.data()), reinterpret_cast<DP*>(H1k3.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);								  								  
					
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k_count.data()), reinterpret_cast<DP*>(H1k1_count.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 

	// The shells near the edges do not complete half sphere, so normalize the shells.
	
	if (my_id == master_id) {
		int Kmax_inside = Max_radius_inside(); 	
		int	Kmax = Min_radius_outside();

		for (int index = Kmax_inside+1; index <= Kmax; index++) 
			if (H1k1_count(index) >= 1) {
                
				H1k1(index) = H1k1(index)* Approx_number_modes_in_shell(index)/H1k1_count(index); 
										
				H1k2(index) = H1k2(index)* Approx_number_modes_in_shell(index)/H1k1_count(index); 
										
				H1k3(index) = H1k3(index)* Approx_number_modes_in_shell(index)/H1k1_count(index); 
			}
			
	} */
	
}



/**********************************************************************************************

	Compute entropy = sum(pi log(1/pi) where pi = Ei/Total E
	
	Skip the mean mode k= (0,0,0)

***********************************************************************************************/



DP CFFF_SLAB::Get_local_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{
	
/*	DP modal_energy, prob, local_entropy;
	
	DP total_energy = Get_total_energy(Ax) + Get_total_energy(Ay) + Get_total_energy(Az);
	
	MPI_Bcast( &total_energy, 1, MPI_DP, master_id, MPI_COMM_WORLD); 
	
	local_entropy = 0.0;
			
	for (int lx=0; lx<local_Nx; lx++)									
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
                
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz));
								
				prob = modal_energy/ total_energy;		
				
				if (prob > MYEPS) 
					local_entropy += (-prob * log(prob)/log(2.0));
			}
			
	return local_entropy;	*/	
}											


DP CFFF_SLAB::Get_entropy(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
{

/*	DP local_entropy = Get_local_entropy(Ax, Ay, Az);
	
	DP entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return entropy;
		
	else
		return 0;	*/
}


//
// Scalar
// 


DP CFFF_SLAB::Get_local_entropy_scalar(Array<complx,3> A)
{

/*	DP modal_energy, prob;
	
	DP local_entropy = 0.0;
	
	DP total_energy = Get_total_energy(A);
	
	MPI_Bcast( &total_energy, 1, MPI_DP, master_id, MPI_COMM_WORLD);
			
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=10; lz<=Nz/2; lz++) {
			
				modal_energy = Multiplicity_factor(lx,ly,lz)* Vsqr(A(lx,ly,lz));
				
				prob = modal_energy / total_energy;		
				
				if (prob > MYEPS) 
					local_entropy += (-prob*log(prob)/log(2.0));
			}
			
	return local_entropy;	*/	

}				


DP CFFF_SLAB::Get_entropy_scalar(Array<complx,3> A)
{

/*	DP local_entropy = Get_local_entropy_scalar(A);
	
	DP entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return entropy;
	
	else
		return 0;	*/	
}



/**********************************************************************************************

	Compute ring spectrum: Craya basis.
	2D: Sk(k, theta)
	3d: Sk1(k, theta) along horzontal direction; S2k(k, theta) 
			in the perpendicular dirn (polar).

***********************************************************************************************/


void CFFF_SLAB::Compute_local_ring_spectrum
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	int n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k
)											
{
/*	local_S1k = 0.0;
	local_S2k = 0.0;
	
	complx anisV1, anisV2;
	TinyVector<complx,3>  K;
	TinyVector<complx,3> V, VcrossK;
	
	DP Kmag, theta, Kperp;
	DP V2sqr;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
									
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					
					V = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
					
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					Wavenumber(lx, ly, lz, K);
					
					Kperp = AnisKperp(lx, ly, lz);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					// k on the anisotropy axis
					if (Kperp < MYEPS)	{
					
						if (global.field.anisotropy_dirn == 1) {
							anisV1 = V(1); 	anisV2 = V(2);
						}

						else if (global.field.anisotropy_dirn == 2) {
							anisV1 = V(2);	anisV2 = V(0);			
						}	
						
						else if (global.field.anisotropy_dirn == 3) {
							anisV1 = V(0);	anisV2 = V(1);
						}	

						local_S1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(anisV1);
						local_S2k(shell_index, sector_index) += factor*my_pow(Kmag,n)* Vsqr(anisV2);					
					}

					else {
						if (Ny > 1) {
							VcrossK = cross(V, K);
					
							if (global.field.anisotropy_dirn == 1)
								anisV1 = VcrossK(0)/Kperp;

							else if (global.field.anisotropy_dirn == 2)
								anisV1 = VcrossK(1)/Kperp;

							else if (global.field.anisotropy_dirn == 3)
								anisV1 = VcrossK(2)/Kperp;


							local_S1k(shell_index, sector_index) += factor*my_pow(Kmag,n)*Vsqr(anisV1);
							
							V2sqr = Vsqr(V) -Vsqr(anisV1);
							
							local_S2k(shell_index, sector_index) +=  factor*my_pow(Kmag,n)*V2sqr;
						}
						else {
							local_S1k(shell_index, sector_index) +=  factor*my_pow(Kmag,n)*Vsqr(V(1));
							// perp to the plane
							
							V2sqr = Vsqr(V(0)) + Vsqr(V(2));
							local_S2k(shell_index, sector_index) +=  factor*my_pow(Kmag,n)*V2sqr;
						}
						
					}	// of inner else
				}		// of if
			}			// of for loop
*/
}


//*********************************************************************************************

void CFFF_SLAB::Compute_ring_spectrum
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	int n, 
	Array<DP,2> S1k, Array<DP,2> S2k
)
{

/*	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_ring_spectrum(Ax, Ay, Az, n,local_S1k, local_S2k);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S1k.data()), reinterpret_cast<DP*>(S1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<DP*>(local_S2k.data()), reinterpret_cast<DP*>(S2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	*/				 

}



//*********************************************************************************************
// 
//  A.B
//


void CFFF_SLAB::Compute_local_ring_spectrum
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	int n, Array<DP,2> local_S1k, Array<DP,2> local_S2k
)											
{

/*	local_S1k = 0.0;
	local_S2k = 0.0;
	
	
	complx anisV1, anisV2, anisW1, anisW2;
	TinyVector<complx,3>  K, e1K;
	TinyVector<complx,3> V, VcrossK, W, WcrossK;
	
	DP Kmag, theta, Kh1, Kh2, Kperp;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					V = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
					W = Bx(lx, ly, lz), By(lx, ly, lz), Bz(lx, ly, lz);
					
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					Wavenumber(lx, ly, lz, K);
				
					Kperp = AnisKperp(lx, ly, lz);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					// k on the anisotropy axis.
					if (Kperp < MYEPS)	{
					
						if (global.field.anisotropy_dirn == 1)	{
							anisV1 = V(1); 	anisV2 = V(2);
							anisW1 = W(1); 	anisW2 = W(2);
						}	
						
						else if (global.field.anisotropy_dirn == 2) {
							anisV1 = V(2);	anisV2 = V(0);	
							anisW1 = W(2); 	anisW2 = W(0);
						}	
						
						else if (global.field.anisotropy_dirn == 3) {
							anisV1 = V(0);	anisV2 = V(1);
							anisW1 = W(0); 	anisW2 = W(1);
						}	

						local_S1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(anisV1, anisW1);
						local_S2k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(anisV2, anisW2);				
					}
				
					else {
						
						if (Ny > 1){
							VcrossK = cross(V, K);
							WcrossK = cross(W, K);
							
							Kh1 = AnisKh1(lx, ly, lz);
							Kh2 = AnisKh2(lx, ly, lz);
						
							if (global.field.anisotropy_dirn == 1) {
								anisV1 = VcrossK(0)/Kperp;
								anisW1 = WcrossK(0)/Kperp;
								
								e1K = 0.0, Kh2/Kperp, -Kh1/Kperp;
								
								anisV2 = dot(VcrossK, e1K);
								anisW2 = dot(WcrossK, e1K);
							}	
						
							else if (global.field.anisotropy_dirn == 2) {
								anisV1 = VcrossK(1)/Kperp;
								anisW1 = WcrossK(1)/Kperp;
								
								e1K = -Kh1/Kperp, 0.0, Kh2/Kperp;
								
								anisV2 = dot(VcrossK, e1K);
								anisW2 = dot(WcrossK, e1K);
							}	
							
							else if (global.field.anisotropy_dirn == 3) {
								anisV1 = VcrossK(2)/Kperp;
								anisW1 = WcrossK(2)/Kperp;
								
								e1K = Kh2/Kperp, -Kh1/Kperp, 0.0;
								
								anisV2 = dot(VcrossK, e1K);
								anisW2 = dot(WcrossK, e1K);
							}	
						
							local_S1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(anisV1, anisW1);
							local_S2k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(anisV2, anisW2);
						}
						// 2D
						else {
							local_S1k(shell_index, sector_index) += factor*my_pow(Kmag,n)* real_cprod(V(1),W(1));
                            real(V(1)*conj(W(1)));
						
							DP temp = real_cprod(V(0),W(0)) + real_cprod(V(2),W(2)); 
							local_S1k(shell_index, sector_index) = factor * my_pow(Kmag,n)*temp; 
						}
				
					}	// of inner else
				}		// of if (shell_index <= Kmax)
			}			// of for loop
	*/		
}


//*********************************************************************************************

void CFFF_SLAB::Compute_ring_spectrum
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	int n, Array<DP,2> S1k, Array<DP,2> S2k
)
{
/*
	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_ring_spectrum(Ax, Ay, Az, Bx, By, Bz, n, local_S1k, local_S2k);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S1k.data()), reinterpret_cast<DP*>(S1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<DP*>(local_S2k.data()), reinterpret_cast<DP*>(S2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);					 
*/
}


//*********************************************************************************************
// scalar
//

void CFFF_SLAB::Compute_local_ring_spectrum
(
	Array<complx,3> F, 
	int n, 
	Array<DP,2> local_Sk
)										
{

/*	local_Sk = 0.0;
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
			
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
				
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
				
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)* Vsqr(F(lx,ly,lz));
				}	
			} */
}
					


void CFFF_SLAB::Compute_ring_spectrum
( 
	Array<complx,3> F, 
	int n, 
	Array<DP,2> Sk
)
{

/*	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_ring_spectrum(F, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()),reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
 */

}


//
// F*G
//

void CFFF_SLAB::Compute_local_ring_spectrum
(
	Array<complx,3> F, 
	Array<complx,3> G, 
	int n, 
	Array<DP,2> local_Sk
)										
{
/*
	local_Sk = 0.0;	
	
	DP Kmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	Kmax = Max_radius_inside();
			
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
				
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = Multiplicity_factor(lx, ly, lz);
				
					local_Sk(shell_index, sector_index) +=  factor*my_pow(Kmag,n)*real_cprod(F(lx,ly,lz),G(lx,ly,lz));
				}						
			}	
*/	

}


void CFFF_SLAB::Compute_ring_spectrum
(
	Array<complx,3> F,
	Array<complx,3> G,  
	int n, 
	Array<DP,2> Sk
)
{
/*
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
			
	Compute_local_ring_spectrum(F, G, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
*/
}



//*********************************************************************************************

// Not for 2D
void CFFF_SLAB::Compute_local_ring_spectrum_helicity
( 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,2> local_H1k
)												
{
/*
	local_H1k = 0.0;	

	DP Kmag, theta;														
	DP modal_helicity;
	DP factor;
	int shell_index, sector_index;
	

	int	Kmax = Max_radius_inside();
		
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz); 
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = 2*Multiplicity_factor(lx, ly, lz); 
					
					modal_helicity = factor* Get_Modal_helicity(lx, ly, lz, Ax, Ay, Az);
									
					local_H1k(shell_index, sector_index) +=  modal_helicity;
				}	
			} 
*/
}

//
//

void CFFF_SLAB::Compute_ring_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,2> H1k
)	
{
/*
	static Array<DP,2> local_H1k(H1k.length()(0), (H1k.length())(1));
	
	local_H1k = 0.0;	
	
	Compute_local_ring_spectrum_helicity(Ax, Ay, Az, local_H1k);
	
	int data_size = (H1k.length())(0) * (H1k.length())(1);
				
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k.data()), reinterpret_cast<DP*>(H1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
	*/
}


/**********************************************************************************************
	
								CYLINDERICAL RING SPECTRUM
 
									Not for 2D
				
***********************************************************************************************/


void CFFF_SLAB::Compute_local_cylindrical_ring_spectrum
( 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	int n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k
)						
{
/*
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;

	complx anisV1;
	
	DP Kmag, Kpll, Kperp;										
	DP V2sqr;
	DP factor;
	int shell_index, slab_index;
	TinyVector<complx,3> V;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				
				V = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {					
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
				
					if (global.field.anisotropy_dirn == 1)
						anisV1 = V(0);
					
					else if (global.field.anisotropy_dirn == 2)
						anisV1 = V(1);

					else if (global.field.anisotropy_dirn == 2)
						anisV1 = V(2);
					
					
					local_S1k(shell_index, slab_index) += factor*my_pow(Kmag,n)* Vsqr(anisV1);
					
					V2sqr = Vsqr(V)- Vsqr(anisV1);
					local_S2k(shell_index, slab_index) += factor*my_pow(Kmag,n)*V2sqr;
				}	
			}
	*/
}

//
//
void CFFF_SLAB::Compute_cylindrical_ring_spectrum
( 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	int n, 
	Array<DP,2> S1k, Array<DP,2> S2k
)
{
/*
	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, n, local_S1k, local_S2k);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S1k.data()), reinterpret_cast<DP*>(S1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<DP*>(local_S2k.data()), reinterpret_cast<DP*>(S2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); */

}

//*********************************************************************************************

void CFFF_SLAB::Compute_local_cylindrical_ring_spectrum
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	int n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k
)
{
/*	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;

	complx anisV1, anisW1;
	TinyVector<complx,3> V, Vrho, W, Wrho;
	
	DP Kmag, Kpll, Kperp;
	DP factor;
	int shell_index, slab_index;

	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				
				Kmag = Kmagnitude(lx, ly, lz);
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					V = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
					W = Bx(lx, ly, lz), By(lx, ly, lz), Bz(lx, ly, lz);
					
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
			} */

}	



//
//
void CFFF_SLAB::Compute_cylindrical_ring_spectrum
( 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	int n, 
	Array<DP,2> S1k, Array<DP,2> S2k
)
{

/*	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(Ax, Ay, Az, Bx, By, Bz, n, local_S1k, local_S2k);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_S1k.data()), reinterpret_cast<DP*>(S1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<DP*>(local_S2k.data()), reinterpret_cast<DP*>(S2k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);		*/			 

}


//*********************************************************************************************

void CFFF_SLAB::Compute_local_cylindrical_ring_spectrum
(
	Array<complx,3> F, 
	int n,   
	Array<DP,2> local_Sk
)
{
/*
	local_Sk = 0.0;
	
	DP Kmag, Kperp, Kpll;
	DP factor;
	int shell_index, slab_index;

	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)* Vsqr(F(lx,ly,lz));
				}	
			}	
*/
}

//
//
void CFFF_SLAB::Compute_cylindrical_ring_spectrum
(
	Array<complx,3> F, 
	int n,   
	Array<DP,2> Sk
)
{
/*
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(F, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
*/
}


//*********************************************************************************************

void CFFF_SLAB::Compute_local_cylindrical_ring_spectrum
(
	Array<complx,3> F, 
	Array<complx,3> G,
	int n,   
	Array<DP,2> local_Sk
)
{
/*
	local_Sk = 0.0;
	
	DP Kmag, Kperp, Kpll;
	DP factor;
	int shell_index, slab_index;

	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = Multiplicity_factor(lx, ly, lz);
					
					local_Sk(shell_index, slab_index) += factor*my_pow(Kmag,n)*real_cprod(F(lx,ly,lz),G(lx,ly,lz));
				}	
			}	
*/
}

//
//
void CFFF_SLAB::Compute_cylindrical_ring_spectrum
(
	Array<complx,3> F, 
	Array<complx,3> G,
	int n,   
	Array<DP,2> Sk
)
{
/*	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_cylindrical_ring_spectrum(F, G, n, local_Sk);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
*/
}


//*********************************************************************************************

void CFFF_SLAB::Compute_local_cylindrical_ring_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,2> local_H1k
)
{
/*
	local_H1k = 0.0;  											

	DP Kmag, Kpll, Kperp;															
	DP modal_helicity;
	DP factor;
	
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	
	for (int lx=0; lx<local_Nx; lx++)										
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {			
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = 2*Multiplicity_factor(lx, ly, lz); 
					
					modal_helicity = factor* Get_Modal_helicity(lx, ly, lz, Ax, Ay, Az);				
					local_H1k(shell_index, slab_index) +=  modal_helicity;
				}	
			}	
*/
}

//
//
void CFFF_SLAB::Compute_cylindrical_ring_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,2> H1k
)
{
/*	static Array<DP,2> local_H1k(H1k.length()(0), (H1k.length())(1));
	
	local_H1k = 0.0;	
	
	Compute_local_ring_spectrum_helicity(Ax, Ay, Az, local_H1k);
	
	int data_size = (H1k.length())(0) * (H1k.length())(1);
				
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k.data()), reinterpret_cast<DP*>(H1k.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
	*/
}

/**********************************************************************************************

			(B0.k) Im[vectA. conj(vectB)]

***********************************************************************************************/

void CFFF_SLAB::Compute_local_imag_shell_spectrum_B0
( 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> local_Sk
)
{
/*	DP Kmag;
	int shell_index;
	TinyVector<DP,3> K;
	TinyVector<complx,3> V, W;
	
	
	local_Sk = 0.0;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				
				shell_index = (int) ceil(Kmag);
				
				if (shell_index <= Kmax) {
					Wavenumber(lx, ly, lz, K);
					V = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
					W = Bx(lx, ly, lz), By(lx, ly, lz), Bz(lx, ly, lz);
					
					local_Sk(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);
													
					// factor multiplied by 2. recall the defn of energy spectrum and factor	
				}									
			}
	*/				
}


//
//
void CFFF_SLAB::Compute_imag_shell_spectrum_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> Sk
)
{

/*	static Array<DP,1> local_Sk(Sk.length());
	
	local_Sk = 0.0;	
	
	Compute_local_imag_shell_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz,local_Sk);

	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
 */
					
}	
				
//******************************** for rings **************************************************

void CFFF_SLAB::Compute_local_imag_ring_spectrum_B0
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
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
				
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);	
					
					Wavenumber(lx, ly, lz, K);
					V = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
					W = Bx(lx, ly, lz), By(lx, ly, lz), Bz(lx, ly, lz);
					
					local_Sk(shell_index, sector_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* mydot_imag(V,W);		
				}									
			}
 */
					
}

//
//
void CFFF_SLAB::Compute_imag_ring_spectrum_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> Sk
)
{

/*	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	local_Sk = 0.0;	
	
	Compute_local_imag_ring_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz,local_Sk);

	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
	*/				
}	


//******************************** for cylinder rings *****************************************

// Not for 2D
void CFFF_SLAB::Compute_local_imag_cylindrical_ring_spectrum_B0
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
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<Ny; ly++) 
			for (int lz=0; lz<=Nz/2; lz++)  {
				Kpll = AnisKpll(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					Wavenumber(lx, ly, lz, K);
					V = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
					W = Bx(lx, ly, lz), By(lx, ly, lz), Bz(lx, ly, lz);
					
					local_Sk(shell_index, slab_index) +=  2* Multiplicity_factor(lx, ly, lz)* dot(B0,K)* mydot_imag(V,W);
														
				}									
			}
*/
}

//
//
void CFFF_SLAB::Compute_imag_cylindrical_ring_spectrum_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP,2> Sk
)
{
/*
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	local_Sk = 0.0;	
	
	Compute_local_imag_cylindrical_ring_spectrum_B0(B0, Ax, Ay, Az, Bx, By, Bz, local_Sk);

	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_Sk.data()), reinterpret_cast<DP*>(Sk.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
*/					
}	

//*****************************  End of four_energy.cc ****************************************	


