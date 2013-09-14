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

#include "SSS_pencil.h"
#include "SSS_pencil_inline.h"
#include "shell_etc_indices.h"

/**********************************************************************************************

					SHELL MULT  Real(A.B*) 

***********************************************************************************************/



DP SSS_PENCIL::Local_shell_mult_single
( 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius
)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP Kmag;
	DP local_result = 0.0;
	
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) ) 
					local_result += 2 * Multiplicity_factor(lx, ly, lz)* Ar(ly,lz,lx)*Br(ly,lz,lx);
			}
		
	return local_result;
}
			


//*********************************************************************************************


void SSS_PENCIL::Local_shell_mult_all
( 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array,
	Array<DP,1> local_result
)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP Kmag;
	int shell_index;
	
	local_result = 0.0;
	 
	 int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
						shell_index = Get_shell_index(Kmag, shell_radius_array);
					
						local_result(shell_index) += 2*Multiplicity_factor(lx,ly,lz)* Ar(ly,lz,lx)*Br(ly,lz,lx);
				}									   
			}
				
}			


//*********************************************************************************************


void SSS_PENCIL::Local_shell_mult_all
(	 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> shell_radius_array,
	Array<DP,1> local_result
)
{
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	Array<DP,3> Bxr=Array<DP,3>(reinterpret_cast<DP*>(Bx.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Byr=Array<DP,3>(reinterpret_cast<DP*>(By.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	Array<DP,3> Bzr=Array<DP,3>(reinterpret_cast<DP*>(Bz.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP Kmag;
	int shell_index;
	DP AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Axr.extent(0); ly++)
        for (int lz=0; lz<Axr.extent(1); lz++)
            for (int lx=0; lx<Axr.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
                    shell_index = Get_shell_index(Kmag, shell_radius_array);
					
                    AdotB = Axr(ly,lz,lx)*Bxr(ly,lz,lx) + Ayr(ly,lz,lx)*Byr(ly,lz,lx)  + Azr(ly,lz,lx)*Bzr(ly,lz,lx);
					
                    local_result(shell_index) += 2*Multiplicity_factor(lx,ly,lz)* AdotB;  
				}										
			} 
				
}



/**********************************************************************************************

			RING MULT

***********************************************************************************************/


void SSS_PENCIL::Local_ring_mult_all
( 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP,2> local_result
)
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP Kmag, theta;
	int shell_index, sector_index;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					local_result(shell_index, sector_index) += 2 *Multiplicity_factor(lx,ly,lz)* Ar(ly,lz,lx)*Br(ly,lz,lx);
				}									   
			}

}


//*********************************************************************************************


void SSS_PENCIL::Local_ring_mult_all
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> local_result
)
{
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	Array<DP,3> Bxr=Array<DP,3>(reinterpret_cast<DP*>(Bx.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Byr=Array<DP,3>(reinterpret_cast<DP*>(By.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	Array<DP,3> Bzr=Array<DP,3>(reinterpret_cast<DP*>(Bz.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP Kmag, theta;
	int shell_index, sector_index;
	DP AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Axr.extent(0); ly++)
        for (int lz=0; lz<Axr.extent(1); lz++)
            for (int lx=0; lx<Axr.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					AdotB = Axr(ly,lz,lx)*Bxr(ly,lz,lx) + Ayr(ly,lz,lx)*Byr(ly,lz,lx)  + Azr(ly,lz,lx)*Bzr(ly,lz,lx);
					
					local_result(shell_index, sector_index) +=  2* Multiplicity_factor(lx, ly, lz)* AdotB;  
				}											
			}
			
}
		
			

/**********************************************************************************************
	
						CYLINDRICAL RING MULTIPLICATION OD ARRAYS
 
								Not for 2D
				
***********************************************************************************************/

void SSS_PENCIL::Local_cyl_ring_mult_all
(
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP,2> local_result
)
{

	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP  Kpll,  Kperp;
	int shell_index, slab_index;
	
	local_result = 0.0;
		
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index(Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);
					
					local_result(shell_index, slab_index) +=  2* Multiplicity_factor(lx,ly,lz)* Ar(ly,lz,lx)*Br(ly,lz,lx);
				}										
			}
		
}


//*********************************************************************************************

void SSS_PENCIL::Local_cyl_ring_mult_all
(
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP,2> local_result
)
{
	Array<DP,3> Axr=Array<DP,3>(reinterpret_cast<DP*>(Ax.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Ayr=Array<DP,3>(reinterpret_cast<DP*>(Ay.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	Array<DP,3> Azr=Array<DP,3>(reinterpret_cast<DP*>(Az.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	Array<DP,3> Bxr=Array<DP,3>(reinterpret_cast<DP*>(Bx.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Byr=Array<DP,3>(reinterpret_cast<DP*>(By.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	Array<DP,3> Bzr=Array<DP,3>(reinterpret_cast<DP*>(Bz.data()), shape_complex_array*shape(1,2,1), neverDeleteData);

	DP  Kpll,  Kperp, AdotB;
	int shell_index, slab_index;
	
	local_result = 0.0;
	
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<Axr.extent(0); ly++)
        for (int lz=0; lz<Axr.extent(1); lz++)
            for (int lx=0; lx<Axr.extent(2); lx++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index( Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);
					
					AdotB = Axr(ly,lz,lx)*Bxr(ly,lz,lx) + Ayr(ly,lz,lx)*Byr(ly,lz,lx)  + Azr(ly,lz,lx)*Bzr(ly,lz,lx);
					
					local_result(shell_index, slab_index) +=  2 *Multiplicity_factor(lx,ly,lz)* AdotB;
				}									
			}

}

/**********************************************************************************************

			(B0.k) Im[vectA. conj(vectB)]

***********************************************************************************************/


void SSS_PENCIL::Local_shell_mult_all_imagVW_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> local_result
)
{
/*	DP Kmag;
	int shell_index;
	TinyVector<DP,3>  K;
	DP AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					Wavenumber(lx, ly, lz, K);	
					
					AdotB = mydot_imag(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx),Bx(ly,lz,lx),By(ly,lz,lx),Bz(ly,lz,lx));
					
					local_result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)*  dot(B0,K)* AdotB;
				}									
			}
 */
					
}



//*********************************************************************************************

void SSS_PENCIL::Local_ring_mult_all_imagVW_B0
(
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> local_result
)
{
/*	DP Kmag, theta, AdotB;
	int shell_index, sector_index;
	TinyVector<DP,3>  K;
	
	local_result = 0.0;
		
	int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					Wavenumber(lx, ly, lz, K);
					AdotB = mydot_imag(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx),Bx(ly,lz,lx),By(ly,lz,lx),Bz(ly,lz,lx));
					
					local_result(shell_index, sector_index)  +=	 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* AdotB;
				}									
			}
 */
					
}

//
//
//


//*********************************************************************************************

void SSS_PENCIL::Local_cyl_ring_mult_all_imagVW_B0
( 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,2> local_result
)
{

/*	DP  Kpll,  Kperp, AdotB;
	int shell_index, slab_index;
	TinyVector<DP,3>  K;
	
	local_result = 0.0;
	
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index( Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);

					Wavenumber(lx, ly, lz, K);
					AdotB = mydot_imag(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx),Bx(ly,lz,lx),By(ly,lz,lx),Bz(ly,lz,lx));

					local_result(shell_index, slab_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* AdotB;
				}									
			}
*/
}




//*********************************************************************************************

//
//  sum(A.curl(B)) for a given shell
// 
DP SSS_PENCIL::Local_shell_mult_vorticity
(
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
    DP inner_radius, DP outer_radius
)
{
	
/*	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP Kmag, AdotB;
	DP result = 0.0;
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 below to account for the absent complex conj modes. 
						result += 2* Multiplicity_factor(lx,ly,lz)* AdotB;
										
					}	
					
					else if (Ny ==1) {
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						// nlin2 contains U.grad omega
						// factor 2 to account for the absent complex conj modes.
						result += 2* Multiplicity_factor(lx,0,lz)*real_cprod(Ay(ly,lz,lx), vort_y);
					}	
				}	
			}
	
	return result;
 */

	return 0;
}

//
//	sum(A.curl(B)/k^2) for a given shell
//
DP SSS_PENCIL::Local_shell_mult_vector_potential
(
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
    DP inner_radius, DP outer_radius
)
{
/*	TinyVector<complx,3> vorticity;
	complx vort_y;
	DP AdotB;
	
	DP Kmag;
	DP result = 0.0;
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) )  {
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx), vorticity(0),vorticity(1),vorticity(2));
						
						result += 2* Multiplicity_factor(lx,ly,lz)/pow2(Kmag)* AdotB;
							// factor 2 to account for the absent complex conj modes. 
					}	
					
					else if (Ny ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						result += 2 * Multiplicity_factor(lx,0,lz)*real_cprod(Ay(ly,lz,lx), vort_y)/pow2(Kmag);
							// nlin2 contains U.grad a; a(k) = curl(b)/k^2
							// factor 2 to account for the absent complex conj modes. 
					}	
					
				}	
			}
	
	return result;
 */

	return 0;
}



//
//	sum(A.curl(B)) for all shell
//
void SSS_PENCIL::Local_shell_mult_vorticity_all
( 
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
    Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
    Array<DP,1> result
)
{
	
/*	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP Kmag, AdotB;
	int shell_index;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* AdotB;
						
					}	
					
					else if (Ny ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						// nlin2 contains U.grad omega
						// factor 2 to account for the absent complex conj modes. 
						result(shell_index) += 2* Multiplicity_factor(lx,0,lz)*real_cprod(Ay(ly,lz,lx), vort_y);
					}
				}							
			}
 */

	result =0;
}


//
//	sum(A.curl(B)/k^2) for all shell
//
void SSS_PENCIL::Local_shell_mult_vector_potential_all
(
     Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
     Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
     Array<DP,1> result
)
{
/*	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP Kmag, AdotB;
	int shell_index;
	
	result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
	for (int ly=0; ly<Ax.extent(0); ly++)
        for (int lz=0; lz<Ax.extent(1); lz++)
            for (int lx=0; lx<Ax.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(ly,lz,lx),Ay(ly,lz,lx),Az(ly,lz,lx), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* AdotB/pow2(Kmag);
					}	
					
					else if (Ny ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						result(shell_index) += 2* Multiplicity_factor(lx,0,lz)* real_cprod(Ay(ly,lz,lx),vort_y)/ pow2(Kmag);
						
					}
				}							
			}
 */
	result =0;
	
}


//*********************************************************************************************

void SSS_PENCIL::Fill_array_shell(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius)					
{
	
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP Kmag;
	
	B = 0.0;
	
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) ) 
					Br(ly,lz,lx) = Ar(ly,lz,lx);
				}
}


void SSS_PENCIL::Fill_array_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP left_angle, DP right_angle)		
{
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	DP Kmag, theta;
	
	B = 0.0;
	
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				theta = AnisKvect_polar_angle(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) 
					if ( ((theta > left_angle) && (theta <= right_angle)) || ((abs(theta) < MYEPS) && (left_angle < MYEPS)) )
						Br(ly,lz,lx) = Ar(ly,lz,lx);
			}
}


	
void SSS_PENCIL::Fill_array_cylindrical_ring(Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, DP h_lower, DP h_upper)		
{
	
	Array<DP,3> Ar=Array<DP,3>(reinterpret_cast<DP*>(A.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
    Array<DP,3> Br=Array<DP,3>(reinterpret_cast<DP*>(B.data()), shape_complex_array*shape(1,2,1), neverDeleteData);
	
	
	DP Kmag, Kpll;
	
	B = 0.0;
	
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
				Kmag = Kmagnitude(lx, ly, lz);
				Kpll = AnisKpll(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) 
					if ( ((Kpll > h_lower) && (Kpll <= h_upper)) || ((abs(Kpll) < MYEPS) && (h_lower < MYEPS)) )
						Br(ly,lz,lx) = Ar(ly,lz,lx);
			}
}








//****************************  End of four_ET.cc   *******************************************





