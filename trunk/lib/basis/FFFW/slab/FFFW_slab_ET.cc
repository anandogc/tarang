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


/*! \file four_ET.cc
 * 
 * @sa four_ET.h
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 22/08/2008
 * @bug no known bug
 */ 

#include "FFFW_slab.h"
#include "FFFW_slab_inline.h"
#include "shell_etc_indices.h"

/**********************************************************************************************

					SHELL MULT  Real(A.B*) 

***********************************************************************************************/



DP FFFW_SLAB::Local_shell_mult_single
( 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius
)
{
	DP Kmag;
	DP local_result = 0.0;
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) ) 
					local_result += 2 * Multiplicity_factor(lx, ly, lz)*real( A(lx,ly,lz)*conj(B(lx,ly,lz)) );   
			}
		
	return local_result;
}


//
//
DP FFFW_SLAB::Shell_mult_single
( 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius
)
{
	DP local_total;
	local_total = Local_shell_mult_single(A, B, inner_radius, outer_radius);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR);	

	if (my_id == master_id) 
		return total;
		
	else
		return 0;	
}				


//*********************************************************************************************


void FFFW_SLAB::Local_shell_mult_all
( 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array,
	Array<DP,1> local_result
)
{
	DP Kmag;
	int shell_index;
	
	local_result = 0.0;
	 
	 int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++)	{
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
						shell_index = Get_shell_index(Kmag, shell_radius_array);
					
						local_result(shell_index) += 2*Multiplicity_factor(lx,ly,lz)* real_cprod(A(lx,ly,lz),B(lx,ly,lz));
				}									   
			}
				
}


//
//
void FFFW_SLAB::Shell_mult_all
( 
	Array<complx,3> A, Array<complx,3> B,
	Array<DP, 1> shell_radius_array,
	Array<DP,1> result
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_all(A, B, shell_radius_array, local_result);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}				


//*********************************************************************************************


void FFFW_SLAB::Local_shell_mult_all
(	 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> shell_radius_array,
	Array<DP,1> local_result
)
{
	DP Kmag;
	int shell_index;
	DP AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
                    shell_index = Get_shell_index(Kmag, shell_radius_array);
					
                    AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));
					
                    local_result(shell_index) += 2*Multiplicity_factor(lx,ly,lz)* AdotB;  
				}										
			} 
				
}


//
//
void FFFW_SLAB::Shell_mult_all
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP, 1> shell_radius_array,
	Array<DP,1> result
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_all(Ax, Ay, Az, Bx, By, Bz, shell_radius_array, local_result);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}	


/**********************************************************************************************

			RING MULT

***********************************************************************************************/


void FFFW_SLAB::Local_ring_mult_all
( 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP,2> local_result
)
{
	DP Kmag, theta;
	int shell_index, sector_index;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					local_result(shell_index, sector_index) += 2 *Multiplicity_factor(lx,ly,lz)* real_cprod(A(lx,ly,lz),B(lx,ly,lz));
				}									   
			}

}


//
//
void FFFW_SLAB::Ring_mult_all
( 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP,2> result
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_ring_mult_all(A, B, local_result);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}				


//*********************************************************************************************


void FFFW_SLAB::Local_ring_mult_all
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> local_result
)
{
	DP Kmag, theta;
	int shell_index, sector_index;
	DP AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));
					
					local_result(shell_index, sector_index) +=  2* Multiplicity_factor(lx, ly, lz)* AdotB;  
				}											
			}
			
}


//
//
void FFFW_SLAB::Ring_mult_all
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP,2> result
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_ring_mult_all(Ax, Ay, Az, Bx, By, Bz, local_result);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}				
			

/**********************************************************************************************
	
						CYLINDRICAL RING MULTIPLICATION OD ARRAYS
 
								Not for 2D
				
***********************************************************************************************/

void FFFW_SLAB::Local_cyl_ring_mult_all
(
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP,2> local_result
)
{

	DP  Kpll,  Kperp;
	int shell_index, slab_index;
	
	local_result = 0.0;
		
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index( Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);
					
					local_result(shell_index, slab_index) +=  2* Multiplicity_factor(lx,ly,lz)* real_cprod(A(lx,ly,lz),B(lx,ly,lz));
				}										
			}
		
}

//
//

void FFFW_SLAB::Cyl_ring_mult_all
(
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP,2> result
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_cyl_ring_mult_all(A, B, local_result);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}



//*********************************************************************************************

void FFFW_SLAB::Local_cyl_ring_mult_all
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP,2> local_result
)
{

	DP  Kpll,  Kperp, AdotB;
	int shell_index, slab_index;
	
	local_result = 0.0;
	
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index( Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);
					
					AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));
					
					local_result(shell_index, slab_index) +=  2 *Multiplicity_factor(lx,ly,lz)* AdotB;
				}									
			}

}


void FFFW_SLAB::Cyl_ring_mult_all
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> result
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_cyl_ring_mult_all(Ax, Ay, Az, Bx, By, Bz, local_result);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}

/**********************************************************************************************

			(B0.k) Im[vectA. conj(vectB)]

***********************************************************************************************/


void FFFW_SLAB::Local_shell_mult_all_imagVW_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> local_result
)
{
	DP Kmag;
	int shell_index;
	TinyVector<DP,3>  K;
	DP AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					Wavenumber(lx, ly, lz, K);	
					
					AdotB = mydot_imag(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));
					
					local_result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)*  dot(B0,K)* AdotB;
				}									
			}
					
}

//
//

void FFFW_SLAB::Shell_mult_all_imagVW_B0
( 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP,1> result
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_all_imagVW_B0(B0, Ax, Ay, Az, Bx, By, Bz, local_result);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}


//*********************************************************************************************

void FFFW_SLAB::Local_ring_mult_all_imagVW_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> local_result
)
{
	DP Kmag, theta, AdotB;
	int shell_index, sector_index;
	TinyVector<DP,3>  K;
	
	local_result = 0.0;
		
	int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					Wavenumber(lx, ly, lz, K);
					AdotB = mydot_imag(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));
					
					local_result(shell_index, sector_index)  +=	 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* AdotB;
				}									
			}
					
}

//
//
//

void FFFW_SLAB::Ring_mult_all_imagVW_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP,2> result
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_ring_mult_all_imagVW_B0(B0, Ax, Ay, Az, Bx, By, Bz, local_result);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD);
}


//*********************************************************************************************

void FFFW_SLAB::Local_cyl_ring_mult_all_imagVW_B0
( 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> local_result
)
{

	DP  Kpll,  Kperp, AdotB;
	int shell_index, slab_index;
	TinyVector<DP,3>  K;
	
	local_result = 0.0;
	
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++)  {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index( Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);

					Wavenumber(lx, ly, lz, K);
					AdotB = mydot_imag(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));

					local_result(shell_index, slab_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* AdotB;
				}									
			}

}

//
//

void FFFW_SLAB::Cyl_ring_mult_all_imagVW_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> result
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_cyl_ring_mult_all_imagVW_B0(B0, Ax, Ay, Az, Bx, By, Bz, local_result);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()),data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}



//*********************************************************************************************

//
//  sum(A.curl(B)) for a given shell
// 
DP FFFW_SLAB::Local_shell_mult_vorticity
(
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
    DP inner_radius, DP outer_radius
)
{
	
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP Kmag, AdotB;
	DP result = 0.0;
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
					if (N[2] > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 below to account for the absent complex conj modes. 
						result += 2* Multiplicity_factor(lx,ly,lz)* AdotB;
										
					}	
					
					else if (N[2] ==1) {
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						// nlin2 contains U.grad omega
						// factor 2 to account for the absent complex conj modes.
						result += 2* Multiplicity_factor(lx,0,lz)*real_cprod(Ay(lx,ly,lz), vort_y);
					}	
				}	
			}
	
	return result;
}

//
//


DP FFFW_SLAB::Shell_mult_vorticity
(
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
    DP inner_radius, DP outer_radius
)
{
	DP local_total;
	local_total = Local_shell_mult_vorticity(Ax, Ay, Az, Bx, By, Bz, inner_radius, outer_radius);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR);	
	
	if (my_id == master_id) 
		return total;
	
	else
		return 0;	
}	

//
//	sum(A.curl(B)/k^2) for a given shell
//
DP FFFW_SLAB::Local_shell_mult_vector_potential
(
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
    DP inner_radius, DP outer_radius
)
{
	TinyVector<complx,3> vorticity;
	complx vort_y;
	DP AdotB;
	
	DP Kmag;
	DP result = 0.0;
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++)  {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) )  {
					if (N[2] > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						result += 2* Multiplicity_factor(lx,ly,lz)/pow2(Kmag)* AdotB;
							// factor 2 to account for the absent complex conj modes. 
					}	
					
					else if (N[2] ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						result += 2 * Multiplicity_factor(lx,0,lz)*real_cprod(Ay(lx,ly,lz), vort_y)/pow2(Kmag);
							// nlin2 contains U.grad a; a(k) = curl(b)/k^2
							// factor 2 to account for the absent complex conj modes. 
					}	
					
				}	
			}
	
	return result;
}


//
//
//
DP FFFW_SLAB::Shell_mult_vector_potential
( 
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
    DP inner_radius, DP outer_radius
)
{
	DP local_total;
	local_total = Local_shell_mult_vector_potential(Ax, Ay, Az, Bx, By, Bz, inner_radius, outer_radius);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR);	
	
	if (my_id == master_id) 
		return total;
	
	else
		return 0;	
}	


//
//	sum(A.curl(B)) for all shell
//
void FFFW_SLAB::Local_shell_mult_vorticity_all
( 
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
    Array<DP,1> result
)
{
	
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP Kmag, AdotB;
	int shell_index;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					if (N[2] > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* AdotB;
						
					}	
					
					else if (N[2] ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						// nlin2 contains U.grad omega
						// factor 2 to account for the absent complex conj modes. 
						result(shell_index) += 2* Multiplicity_factor(lx,0,lz)*real_cprod(Ay(lx,ly,lz), vort_y);
					}
				}							
			}
}

//
//
void FFFW_SLAB::Shell_mult_vorticity_all
(
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
    Array<DP,1> result
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_vorticity_all(Ax, Ay, Az, Bx, By, Bz, local_result);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}	


//
//	sum(A.curl(B)/k^2) for all shell
//
void FFFW_SLAB::Local_shell_mult_vector_potential_all
(
     Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
     Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
     Array<DP,1> result
)
{
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP Kmag, AdotB;
	int shell_index;
	
	result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int lx=0; lx<local_Nx; lx++)				
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					if (N[2] > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* AdotB/pow2(Kmag);
					}	
					
					else if (N[2] ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						result(shell_index) += 2* Multiplicity_factor(lx,0,lz)* real_cprod(Ay(lx,ly,lz),vort_y)/ pow2(Kmag);
						
					}
				}							
			}
	
}


//
//
void FFFW_SLAB::Shell_mult_vector_potential_all
(
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
    Array<DP,1> result
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_vector_potential_all(Ax, Ay, Az, Bx, By, Bz, local_result);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<DP*>(local_result.data()), reinterpret_cast<DP*>(result.data()), data_size, MPI_DP, MPI_SUM, master_id, MPI_COMMUNICATOR); 						
}	


//*********************************************************************************************

void FFFW_SLAB::Fill_array_shell(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius)					
{	
	DP Kmag;
	
	B = 0.0;
	
	for (int lx=0; lx<local_Nx; lx++)			
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) ) 
					B(lx,ly,lz) = A(lx,ly,lz); 
				}
}


void FFFW_SLAB::Fill_array_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP left_angle, DP right_angle)		
{	
	DP Kmag, theta;
	
	B = 0.0;
	
	for (int lx=0; lx<local_Nx; lx++)			
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				theta = AnisKvect_polar_angle(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) 
					if ( ((theta > left_angle) && (theta <= right_angle)) || ((abs(theta) < MYEPS) && (left_angle < MYEPS)) )
						B(lx,ly,lz) = A(lx,ly,lz); 
			}
}


	
void FFFW_SLAB::Fill_array_cylindrical_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP h_lower, DP h_upper)		
{	
	DP Kmag, Kpll;
	
	B = 0.0;
	
	for (int lx=0; lx<local_Nx; lx++)			
		for (int ly=0; ly<N[2]; ly++) 
			for (int lz=0; lz<=N[3]/2; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				Kpll = AnisKpll(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) 
					if ( ((Kpll > h_lower) && (Kpll <= h_upper)) || ((abs(Kpll) < MYEPS) && (h_lower < MYEPS)) )
						B(lx,ly,lz) = A(lx,ly,lz); 
			}
}








//****************************  End of four_ET.cc   *******************************************





