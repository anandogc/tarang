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

#include "SFF_slab.h"
#include "SFF_slab_inline.h"
#include "shell_etc_indices.h"

/**********************************************************************************************

 	Computes total energy |f(m,k)|^2  without k=0.
        Factor 2 because of negative k modes are not stored.
		
	for ky plane Eky = 2*sum(A^2) - sum(A(m,ky,0)^2) - sum(A(0,ky,kz)^2) + A(0,ky,0)^2
	
	for ky = 0 plane: subtract (1/2)A(0)^2 (origin) from the above
		
***********************************************************************************************/



DP SFF_SLAB::Get_local_energy_real_space(Array<DP,3> Ar)
{

	return Array_sqr(Ar(Range::all(),Range::all(),Range(0,Nz-1))) / (DP(Nx)*DP(Ny)*DP(Nz));
}



// 2D included in this
DP SFF_SLAB::Get_local_energy(Array<complx,3> A)  
{
	DP total = 2*Array_sqr(A);
	
	// Subtract kz=0 plane
	total -= Array_sqr(A(Range::all(),Range::all(),0));
	
	// kx=0
	if (my_id == 0) {
		total -= Array_sqr(A(0,Range::all(),Range::all()));
		
			// kz=0, kx=0 (ADD since it has been subtracted twice)
		total += Array_sqr(A(0,Range::all(),0))/2;
	} 

    return total;
}

/*
// Exclude |A(kx,0,0)|^2
DP SFF_SLAB::Get_total_energy_residual(Array<complx,3> A)
{
	
	DP local_total = Get_local_energy(A);
	local_total = local_total - sum(sqr(abs((A)(Range::all(),0,0))));

	DP total;
	
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	
	
	if (my_id == master_id) 
		return total - pow2(abs( A(0,0,0) ))/2;			// Subtract the origin contributions
	
	else
		return 0;	
}
 */

/**********************************************************************************************

		Computes total A(k)*conj(B(k)) without k=0.
	
		
***********************************************************************************************/

DP SFF_SLAB::Get_local_energy_real_space(Array<DP,3> Ar, Array<DP,3> Br)
{
	DP ans= mydot(Ar(Range::all(),Range::all(),Range(0,Nz-1)), Br(Range::all(),Range::all(),Range(0,Nz-1)));
	
	return ans/ (DP(Nx)*DP(Ny)*DP(Nz));
}

// 2D included in this
DP SFF_SLAB::Get_local_energy(Array<complx,3> A, Array<complx,3> B)
{
	
	DP total = 2*mydot(A, B);
	
		// Subtract kz=0 plane
	total -= mydot(A(Range::all(),Range::all(),0), B(Range::all(),Range::all(),0));
	
		// kx=0
	if (my_id == 0) {
		total -= mydot(A(0,Range::all(),Range::all()), B(0,Range::all(),Range::all()));
		
		// kz=0, kx=0 (ADD since it has been subtracted twice)
		total += mydot(A(0,Range::all(),0), B(0,Range::all(),0))/2;
	}
    return total;
}


/**********************************************************************************************
 
 Computes sum[ Kperp^2/K^n |A(k)|^2];   k=0 excluded
 
 ***********************************************************************************************/
/*
DP SFF_SLAB::Get_local_Sn_anis(Array<complx,3> A, DP n)
{
	DP Sn = 0.0;
	DP Kmag, Kperpsqr;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx< local_Nx; lx++) 
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++)	{
				Kmag = Kmagnitude(lx, ly, lz);
				Kperpsqr = pow2(Get_ky(ly)*kfactor[2]) + pow2(lz*kfactor[3]);
				
				if ((Kmag <= Kmax) && (Kmag > MYEPS))
					Sn +=  Multiplicity_factor(lx, ly, lz) *(Kperpsqr/my_pow(Kmag,n))*Vsqr(A(lx,ly,lz));
			} 
	
	return Sn; 	
}

	//

DP SFF_SLAB::Get_total_Sn_anis(Array<complx,3> A, DP n)
{
	
	DP local_total = Get_local_Sn_anis(A, n);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);	
	
	if (my_id == master_id) 
		return total;
	
	else
		return 0;
}

*/
/**********************************************************************************************

	Compute Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2
	Multiplication factor = 1 for kz>0 because of energy spectrum details.
 
	not for 2D

***********************************************************************************************/


void SFF_SLAB::Compute_local_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &local_helicity1, DP &local_helicity2, 
	DP &local_k2H1, DP &local_k2H2
)
{
	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, K;
	DP modal_helicity, Kmag, Ksqr;
	
	
	local_helicity1 = local_helicity2 = 0.0;
	local_k2H1 = 0.0;
	
	int	Kmax = Min_radius_outside();
	
	for (int lx=0; lx<Ax.extent(0); lx++)
        for (int ly=0; ly<Ax.extent(1); ly++)
            for (int lz=0; lz<Ax.extent(2); lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax) {
					modal_helicity = 2*Multiplicity_factor(lx, ly, lz)* Get_Modal_helicity(lx, ly, lz, Ax, Ay, Az);
												
					local_helicity1 += modal_helicity;
					
					Ksqr = pow2(Kmag);
					if (Ksqr > MYEPS)
						local_helicity2 += modal_helicity / Ksqr;
						
					local_k2H1 += Ksqr * modal_helicity;
				}	
			}
			
	local_k2H1 *= 2.0;
	local_k2H2 = 2.0 * local_helicity1;				
}

//

void SFF_SLAB::Compute_total_helicity
( 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &total_k2H1, DP &total_k2H2
)
{
	DP local_helicity1, local_helicity2;
	DP local_k2H1, local_k2H2;
	
	Compute_local_helicity(Ax, Ay, Az, local_helicity1, local_helicity2, local_k2H1, local_k2H2);
	
	MPI_Reduce(&local_helicity1, &total_helicity1, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
								
	MPI_Reduce(&local_helicity2, &total_helicity2, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
								
	MPI_Reduce(&local_k2H1, &total_k2H1, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
									
	MPI_Reduce(&local_k2H2, &total_k2H2, 1, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
}



/**********************************************************************************************

	Compute helicity spectrum
	Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2

	Not for 2D
***********************************************************************************************/

void SFF_SLAB::Compute_local_shell_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3, 
	Array<DP,1> local_H1k_count
)
{
	local_H1k_count = 0.0;
	local_H1k1 = 0.0;	
	local_H1k2 = 0.0;
	local_H1k3 = 0.0;

	TinyVector<complx,3> V, VFour;
	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, K;
	DP Kmag;													// kmag = sqrt(Kx^2+Ky^2+Kz^2)
	int index;
	DP factor;
	
	int	Kmax = Min_radius_outside();

	
	for (int lx=0; lx<Ax.extent(0); lx++)
        for (int ly=0; ly<Ax.extent(1); ly++)
            for (int lz=0; lz<Ax.extent(2); lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				index = (int) ceil(Kmag);
				
				if (index <= Kmax)  {
					factor = 2*Multiplicity_factor(lx, ly, lz);
					
					V = Ax(lx,ly,lz), Ay(lx,ly,lz), Az(lx,ly,lz);
					Convert_to_Fourier_space(V, VFour);
					
					Vreal = real(VFour(0)), real(VFour(1)), real(VFour(2));
					Vimag = imag(VFour(0)), imag(VFour(1)), imag(VFour(2));
					
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber(lx, ly, lz, K);
					
					// modal_helicity = factor * dot(K, VrcrossVi);	
					local_H1k1(index) += factor* (K(0)*VrcrossVi(0));
					local_H1k2(index) += factor* (K(1)*VrcrossVi(1));
					local_H1k3(index) += factor* (K(2)*VrcrossVi(2));
					
					local_H1k_count(index) += 2*factor;
				}	
			}

}

//

void SFF_SLAB::Compute_shell_spectrum_helicity
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3
)	
{

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
				H1k1(index) = H1k1(index) *Approx_number_modes_in_shell(index)/ H1k1_count(index); 
										
				H1k2(index) = H1k2(index) *Approx_number_modes_in_shell(index)/ H1k1_count(index); 
										
				H1k3(index) = H1k3(index) *Approx_number_modes_in_shell(index)/ H1k1_count(index); 
			}	
			
	}
}

//*********************************************************************************************
// Not for 2D
void SFF_SLAB::Compute_local_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,2> local_H1k1, Array<DP,2> local_H1k2, Array<DP,2> local_H1k3
 )
{
	
	local_H1k1 = 0.0;
	local_H1k2 = 0.0;
	local_H1k3 = 0.0;
	
	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, K;
	
	DP Kmag, theta;
	DP modal_helicity;
	DP factor;
	int shell_index, sector_index;
	
	
	int	Kmax = Max_radius_inside();
	
	for (int lx=0; lx<Ax.extent(0); lx++)
        for (int ly=0; ly<Ax.extent(1); ly++)
            for (int lz=0; lz<Ax.extent(2); lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				shell_index = (int) ceil(Kmag);
				
				if ((Kmag > MYEPS) && (shell_index <= Kmax)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					sector_index = Get_sector_index(theta, global.spectrum.ring.sector_angles);
					
					factor = 2*Multiplicity_factor(lx, ly, lz);
					
					Vreal = real(Ax(lx,ly,lz)), real(Ay(lx,ly,lz)), real(Az(lx,ly,lz));
					Vimag = imag(Ax(lx,ly,lz)), imag(Ay(lx,ly,lz)), imag(Az(lx,ly,lz));
					
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber(lx, ly, lz, K);
					
					// modal_helicity = factor * dot(K, VrcrossVi);
					local_H1k1(shell_index, sector_index) += factor* (K(0)*VrcrossVi(0));
					local_H1k2(shell_index, sector_index) += factor* (K(1)*VrcrossVi(1));
					local_H1k3(shell_index, sector_index) += factor* (K(2)*VrcrossVi(2));
				}
			}
	
}

//
//

void SFF_SLAB::Compute_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,2> H1k1, Array<DP,2> H1k2, Array<DP,2> H1k3
 )
{
	static Array<DP,2> local_H1k1(H1k1.shape());
	static Array<DP,2> local_H1k2(H1k1.shape());
	static Array<DP,2> local_H1k3(H1k1.shape());
	
	Compute_local_ring_spectrum_helicity(Ax, Ay, Az, local_H1k1, local_H1k2, local_H1k3);
	
	int data_size = H1k1.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k1.data()), reinterpret_cast<DP*>(H1k1.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k2.data()), reinterpret_cast<DP*>(H1k2.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k3.data()), reinterpret_cast<DP*>(H1k3.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	
}


//*********************************************************************************************

void SFF_SLAB::Compute_local_cylindrical_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,2> local_H1k1,  Array<DP,2> local_H1k2
 )
{
	
	local_H1k1 = 0.0;
	local_H1k2 = 0.0;
	
	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, K;
	
	DP Kmag, Kpll, Kperp;
	DP modal_helicity;
	DP factor;
	
	int shell_index, slab_index;
	
	int	Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<Ax.extent(0); lx++)
        for (int ly=0; ly<Ax.extent(1); ly++)
            for (int lz=0; lz<Ax.extent(2); lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				Kperp = AnisKperp(lx, ly, lz);
				
				shell_index = (int) ceil(Kperp);
				
				if (shell_index <= Kperp_max) {
					Kpll = AnisKpll(lx, ly, lz);
					
					slab_index = Get_slab_index(Kpll, Kperp, global.spectrum.cylindrical_ring.kpll_array);
					
					factor = 2*Multiplicity_factor(lx, ly, lz);
					
					Vreal = real(Ax(lx,ly,lz)), real(Ay(lx,ly,lz)), real(Az(lx,ly,lz));
					Vimag = imag(Ax(lx,ly,lz)), imag(Ay(lx,ly,lz)), imag(Az(lx,ly,lz));
					
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber(lx, ly, lz, K);
					
					// modal_helicity = factor * dot(K, VrcrossVi);
					if (global.field.anisotropy_dirn == 1) {
						local_H1k1(shell_index, slab_index) += factor* (K(0)*VrcrossVi(0));
						local_H1k2(shell_index, slab_index) += factor* (K(1)*VrcrossVi(1) + K(2)*VrcrossVi(2));
					}
					
					else if (global.field.anisotropy_dirn == 2) {
						local_H1k1(shell_index, slab_index) += factor* (K(1)*VrcrossVi(1));
						local_H1k2(shell_index, slab_index) += factor* (K(0)*VrcrossVi(0) + K(2)*VrcrossVi(2));
					}
					
					else if (global.field.anisotropy_dirn == 3) {
						local_H1k1(shell_index, slab_index) += factor* (K(2)*VrcrossVi(2));
						local_H1k2(shell_index, slab_index) += factor* (K(0)*VrcrossVi(0) + K(1)*VrcrossVi(1));
					}
				}
			}
	
}

//
//
void SFF_SLAB::Compute_cylindrical_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,2> H1k1,  Array<DP,2> H1k2
 )
{
	static Array<DP,2> local_H1k1(H1k1.shape());
	static Array<DP,2> local_H1k2(H1k1.shape());
	
	Compute_local_cylindrical_ring_spectrum_helicity(Ax, Ay, Az, local_H1k1, local_H1k2);
	
	int data_size = H1k1.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k1.data()), reinterpret_cast<DP*>(H1k1.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
	MPI_Reduce(reinterpret_cast<DP*>(local_H1k2.data()), reinterpret_cast<DP*>(H1k2.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
}




//*****************************  End of sff_slab_energy.cc ****************************


