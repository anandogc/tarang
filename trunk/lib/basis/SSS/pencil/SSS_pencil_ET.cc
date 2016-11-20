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



Real SSS_PENCIL::Local_shell_mult_single
( 
	Array<Complex,3> A, Array<Complex,3> B, 
	Real inner_radius, Real outer_radius
)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Br=Array<Real,3>(reinterpret_cast<Real*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kmag;
	Real local_result = 0.0;
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) ) 
					local_result += 2 * Multiplicity_factor(lx, ly, lz)* Ar(lx,ly,lz)*Br(lx,ly,lz);
			}
		
	return local_result;
}
			


//*********************************************************************************************


void SSS_PENCIL::Local_shell_mult_all
( 
	Array<Complex,3> A, Array<Complex,3> B, 
	Array<Real, 1> shell_radius_array,
	Array<Real,1> local_result
)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Br=Array<Real,3>(reinterpret_cast<Real*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kmag;
	int shell_index;
	
	local_result = 0.0;
	 
	 int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
						shell_index = Get_shell_index(Kmag, shell_radius_array);
					
						local_result(shell_index) += 2*Multiplicity_factor(lx,ly,lz)* Ar(lx,ly,lz)*Br(lx,ly,lz);
				}									   
			}
				
}			


//*********************************************************************************************


void SSS_PENCIL::Local_shell_mult_all
(	 
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
	Array<Real, 1> shell_radius_array,
	Array<Real,1> local_result
)
{
	Array<Real,3> Axr=Array<Real,3>(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Ayr=Array<Real,3>(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Azr=Array<Real,3>(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Array<Real,3> Bxr=Array<Real,3>(reinterpret_cast<Real*>(Bx.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Byr=Array<Real,3>(reinterpret_cast<Real*>(By.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Bzr=Array<Real,3>(reinterpret_cast<Real*>(Bz.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kmag;
	int shell_index;
	Real AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
                    shell_index = Get_shell_index(Kmag, shell_radius_array);
					
                    AdotB = Axr(lx,ly,lz)*Bxr(lx,ly,lz) + Ayr(lx,ly,lz)*Byr(lx,ly,lz)  + Azr(lx,ly,lz)*Bzr(lx,ly,lz);
					
                    local_result(shell_index) += 2*Multiplicity_factor(lx,ly,lz)* AdotB;  
				}										
			} 
				
}



/**********************************************************************************************

			RING MULT

***********************************************************************************************/


void SSS_PENCIL::Local_ring_mult_all
( 
	Array<Complex,3> A, Array<Complex,3> B, 
	Array<Real,2> local_result
)
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Br=Array<Real,3>(reinterpret_cast<Real*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kmag, theta;
	int shell_index, sector_index;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					local_result(shell_index, sector_index) += 2 *Multiplicity_factor(lx,ly,lz)* Ar(lx,ly,lz)*Br(lx,ly,lz);
				}									   
			}

}


//*********************************************************************************************


void SSS_PENCIL::Local_ring_mult_all
(
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
	Array<Real,2> local_result
)
{
	Array<Real,3> Axr=Array<Real,3>(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Ayr=Array<Real,3>(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Azr=Array<Real,3>(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Array<Real,3> Bxr=Array<Real,3>(reinterpret_cast<Real*>(Bx.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Byr=Array<Real,3>(reinterpret_cast<Real*>(By.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Bzr=Array<Real,3>(reinterpret_cast<Real*>(Bz.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kmag, theta;
	int shell_index, sector_index;
	Real AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					AdotB = Axr(lx,ly,lz)*Bxr(lx,ly,lz) + Ayr(lx,ly,lz)*Byr(lx,ly,lz)  + Azr(lx,ly,lz)*Bzr(lx,ly,lz);
					
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
	Array<Complex,3> A, Array<Complex,3> B, 
	Array<Real,2> local_result
)
{

	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Br=Array<Real,3>(reinterpret_cast<Real*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real  Kpll,  Kperp;
	int shell_index, slab_index;
	
	local_result = 0.0;
		
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index(Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);
					
					local_result(shell_index, slab_index) +=  2* Multiplicity_factor(lx,ly,lz)* Ar(lx,ly,lz)*Br(lx,ly,lz);
				}										
			}
		
}


//*********************************************************************************************

void SSS_PENCIL::Local_cyl_ring_mult_all
(
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
	Array<Real,2> local_result
)
{
	Array<Real,3> Axr=Array<Real,3>(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Ayr=Array<Real,3>(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Azr=Array<Real,3>(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Array<Real,3> Bxr=Array<Real,3>(reinterpret_cast<Real*>(Bx.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Byr=Array<Real,3>(reinterpret_cast<Real*>(By.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Bzr=Array<Real,3>(reinterpret_cast<Real*>(Bz.data()), shape_complex_array*shape(1,1,2), neverDeleteData);

	Real  Kpll,  Kperp, AdotB;
	int shell_index, slab_index;
	
	local_result = 0.0;
	
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index( Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);
					
					AdotB = Axr(lx,ly,lz)*Bxr(lx,ly,lz) + Ayr(lx,ly,lz)*Byr(lx,ly,lz)  + Azr(lx,ly,lz)*Bzr(lx,ly,lz);
					
					local_result(shell_index, slab_index) +=  2 *Multiplicity_factor(lx,ly,lz)* AdotB;
				}									
			}

}

/**********************************************************************************************

			(B0.k) Im[vectA. conj(vectB)]

***********************************************************************************************/


void SSS_PENCIL::Local_shell_mult_all_imagVW_B0
(
	TinyVector<Real,3> B0,
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
	Array<Real,1> local_result
)
{
/*	Real Kmag;
	int shell_index;
	TinyVector<Real,3>  K;
	Real AdotB;
	
	local_result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					Wavenumber(lx, ly, lz, K);	
					
					AdotB = mydot_imag(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));
					
					local_result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)*  dot(B0,K)* AdotB;
				}									
			}
 */
					
}



//*********************************************************************************************

void SSS_PENCIL::Local_ring_mult_all_imagVW_B0
(
	TinyVector<Real,3> B0,
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
	Array<Real,2> local_result
)
{
/*	Real Kmag, theta, AdotB;
	int shell_index, sector_index;
	TinyVector<Real,3>  K;
	
	local_result = 0.0;
		
	int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag <= Kmax_inside) && (Kmag > MYEPS)) {
					theta = AnisKvect_polar_angle(lx, ly, lz);
					
					Compute_ring_index(Kmag, theta, global.energy_transfer.ring_to_ring.radii, global.energy_transfer.ring_to_ring.sector_angles, shell_index, sector_index);
					
					Wavenumber(lx, ly, lz, K);
					AdotB = mydot_imag(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));
					
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
	TinyVector<Real,3> B0,
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
	Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz, 
	Array<Real,2> local_result
)
{

/*	Real  Kpll,  Kperp, AdotB;
	int shell_index, slab_index;
	TinyVector<Real,3>  K;
	
	local_result = 0.0;
	
	int	 Kperp_max = Anis_max_Krho_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				 Kpll = AnisKpll(lx, ly, lz);
				
				 Kperp = AnisKperp(lx, ly, lz);
				
				if ( Kperp <=  Kperp_max) {
					Compute_cylindrical_ring_index( Kpll,  Kperp, global.energy_transfer.cylindrical_ring_to_ring.radii, global.energy_transfer.cylindrical_ring_to_ring.kpll_array, shell_index, slab_index);

					Wavenumber(lx, ly, lz, K);
					AdotB = mydot_imag(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz),Bx(lx,ly,lz),By(lx,ly,lz),Bz(lx,ly,lz));

					local_result(shell_index, slab_index) += 2* Multiplicity_factor(lx,ly,lz)* dot(B0,K)* AdotB;
				}									
			}
*/
}




//*********************************************************************************************

//
//  sum(A.curl(B)) for a given shell
// 
Real SSS_PENCIL::Local_shell_mult_vorticity
(
    Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
    Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
    Real inner_radius, Real outer_radius
)
{
	
/*	TinyVector<Complex,3> vorticity;
	Complex vort_y;
	
	Real Kmag, AdotB;
	Real result = 0.0;
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 below to account for the absent complex conj modes. 
						result += 2* Multiplicity_factor(lx,ly,lz)* AdotB;
										
					}	
					
					else if (Ny ==1) {
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						// nlin2 contains U.grad omega
						// factor 2 to account for the absent complex conj modes.
						result += 2* Multiplicity_factor(lx,0,lz)*real_cprod(Ay(lx,ly,lz), vort_y);
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
Real SSS_PENCIL::Local_shell_mult_vector_potential
(
    Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
    Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
    Real inner_radius, Real outer_radius
)
{
/*	TinyVector<Complex,3> vorticity;
	Complex vort_y;
	Real AdotB;
	
	Real Kmag;
	Real result = 0.0;
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) )  {
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						result += 2* Multiplicity_factor(lx,ly,lz)/pow2(Kmag)* AdotB;
							// factor 2 to account for the absent complex conj modes. 
					}	
					
					else if (Ny ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						result += 2 * Multiplicity_factor(lx,0,lz)*real_cprod(Ay(lx,ly,lz), vort_y)/pow2(Kmag);
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
    Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
    Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,
    Array<Real,1> result
)
{
	
/*	TinyVector<Complex,3> vorticity;
	Complex vort_y;
	
	Real Kmag, AdotB;
	int shell_index;
	
	int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* AdotB;
						
					}	
					
					else if (Ny ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						// nlin2 contains U.grad omega
						// factor 2 to account for the absent complex conj modes. 
						result(shell_index) += 2* Multiplicity_factor(lx,0,lz)*real_cprod(Ay(lx,ly,lz), vort_y);
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
     Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
     Array<Complex,3> Bx, Array<Complex,3> By, Array<Complex,3> Bz,  
     Array<Real,1> result
)
{
/*	TinyVector<Complex,3> vorticity;
	Complex vort_y;
	
	Real Kmag, AdotB;
	int shell_index;
	
	result = 0.0;
	
	int	Kmax_inside = Max_radius_inside();
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if (Kmag <= Kmax_inside) {
					shell_index = Get_shell_index(Kmag, global.energy_transfer.shell_to_shell.radii);
					
					if (Ny > 1) {	
						Compute_Modal_vorticity(lx, ly, lz, Bx, By, Bz, vorticity);
						
						AdotB = mydot(Ax(lx,ly,lz),Ay(lx,ly,lz),Az(lx,ly,lz), vorticity(0),vorticity(1),vorticity(2));
						
						// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2* Multiplicity_factor(lx,ly,lz)* AdotB/pow2(Kmag);
					}	
					
					else if (Ny ==1)	{
						Compute_Modal_vorticity_y_component(lx, 0, lz, Bx, By, Bz, vort_y);
						
						result(shell_index) += 2* Multiplicity_factor(lx,0,lz)* real_cprod(Ay(lx,ly,lz),vort_y)/ pow2(Kmag);
						
					}
				}							
			}
 */
	result =0;
	
}


//*********************************************************************************************

void SSS_PENCIL::Fill_array_shell(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius)					
{
	
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Br=Array<Real,3>(reinterpret_cast<Real*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kmag;
	
	B = 0.0;
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				
				if ( (Kmag > inner_radius) && (Kmag <= outer_radius) ) 
					Br(lx,ly,lz) = Ar(lx,ly,lz);
				}
}


void SSS_PENCIL::Fill_array_ring(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius, Real left_angle, Real right_angle)		
{
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Br=Array<Real,3>(reinterpret_cast<Real*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	Real Kmag, theta;
	
	B = 0.0;
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				theta = AnisKvect_polar_angle(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) 
					if ( ((theta > left_angle) && (theta <= right_angle)) || ((abs(theta) < MYEPS) && (left_angle < MYEPS)) )
						Br(lx,ly,lz) = Ar(lx,ly,lz);
			}
}


	
void SSS_PENCIL::Fill_array_cylindrical_ring(Array<Complex,3> A, Array<Complex,3> B, Real inner_radius, Real outer_radius, Real h_lower, Real h_upper)		
{
	
	Array<Real,3> Ar=Array<Real,3>(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
    Array<Real,3> Br=Array<Real,3>(reinterpret_cast<Real*>(B.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	
	Real Kmag, Kpll;
	
	B = 0.0;
	
    for (int lx=0; lx<maxlx; lx++)
		for (int ly=0; ly<maxly; ly++)
	        for (int lz=0; lz<2*maxlz; lz++) {
				Kmag = Kmagnitude(lx, ly, lz);
				Kpll = AnisKpll(lx, ly, lz);
				
				if ((Kmag > inner_radius) && (Kmag <= outer_radius)) 
					if ( ((Kpll > h_lower) && (Kpll <= h_upper)) || ((abs(Kpll) < MYEPS) && (h_lower < MYEPS)) )
						Br(lx,ly,lz) = Ar(lx,ly,lz);
			}
}








//****************************  End of four_ET.cc   *******************************************





