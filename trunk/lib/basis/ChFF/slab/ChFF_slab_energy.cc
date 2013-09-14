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

#include "ChFF_slab.h"
#include "ChFF_slab_inline.h"
#include "shell_etc_indices.h"



/**********************************************************************************************
 
 Computes real-space total energy int |f(x)|^2/2 dx;
 Note dx=sin(theta) dtheta; dtheta = pi/N with equispaced theta
 
 ***********************************************************************************************/

DP ChFF_SLAB::Get_local_energy_real_space(Array<DP,3> Ar)
{
    int rx;
    DP factor,x;
    DP sum = 0.0;

	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
                rx = Get_rx_real_space(lx);
                factor =  sin((rx)*M_PI/Nx);
                sum += factor*my_pow(Ar(ly,lz,lx),2);
            }
	
	return sum*(M_PI/(2* DP(Nx) * DP(Ny) * DP(Nz)));
}


DP ChFF_SLAB::Get_local_energy(Array<complx,3> A)  
{
	DP total = 2*Array_sqr(A);
	
	// Subtract kz=0 plane
	total -= Array_sqr(A(Range::all(), 0, Range::all()));
	
	// kx=0
	if (my_id == 0) {
		total -= Array_sqr(A(Range::all(), Range::all(), 0));
		
		// kz=0, kx=0 (ADD since it has been subtracted twice)
		total += Array_sqr(A(Range::all(), 0, 0))/2;
	}
    
    return total;
}



/**********************************************************************************************

		Computes total A(k)*conj(B(k)) without k=0.
	
		
***********************************************************************************************/




DP ChFF_SLAB::Get_local_energy_real_space(Array<DP,3> Ar, Array<DP,3> Br)
{
    int rx;
    DP factor,x;
    DP sum = 0.0;
    
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
                rx = Get_rx_real_space(lx);
                factor =  sin((rx)*M_PI/Nx);
                sum += factor*Ar(ly,lz,lx)*Br(ly,lz,lx);
            }
	
	return sum*(M_PI/(2* DP(Nx) * DP(Ny) * DP(Nz)));
}

DP ChFF_SLAB::Get_local_energy(Array<complx,3> A, Array<complx,3> B)
{
	DP total = 2*mydot(A, B);
	
	// Subtract kz=0 plane
	total -= mydot(A(Range::all(),0,Range::all()), B(Range::all(),0,Range::all()));
	
	// kx=0
	if (my_id == 0) {
		total -= mydot(A(Range::all(),Range::all(),0), B(Range::all(),Range::all(),0));
		
		// kz=0, kx=0 (ADD since it has been subtracted twice)
		total += mydot(A(Range::all(),0,0), B(Range::all(),0,0))/2;
	}
    return total;
}

//******************************************************************************

void ChFF_SLAB::Compute_local_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 DP &local_helicity1, DP &local_helicity2,
 DP &local_H1k1, DP &local_H1k2
)
{}

void ChFF_SLAB::Compute_total_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 DP &total_helicity1, DP &total_helicity2,
 DP &total_H1k1, DP &total_H1k2
 )
{}

void ChFF_SLAB::Compute_local_shell_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3,
 Array<DP,1> local_H1k_count
 )
{}


void ChFF_SLAB::Compute_shell_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3
 )
{}

//*********************************************************************************************

// Not for 2D
void ChFF_SLAB::Compute_local_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,2> local_H1k1, Array<DP,2> local_H1k2, Array<DP,2> local_H1k3
 )
{}

//
//

void ChFF_SLAB::Compute_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
Array<DP,2> H1k1, Array<DP,2> H1k2, Array<DP,2> H1k3
 )
{}


//*********************************************************************************************

void ChFF_SLAB::Compute_local_cylindrical_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
  Array<DP,2> local_H1k1,  Array<DP,2> local_H1k2
 )
{}

//
//
void ChFF_SLAB::Compute_cylindrical_ring_spectrum_helicity
(
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<DP,2> H1k1,  Array<DP,2> H1k2
 )
{}



//*****************************  End of ChFF_slab_energy.cc ****************************







//*****************************  End of scft_energy.cc ****************************************






