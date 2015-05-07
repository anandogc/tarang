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

#include "ChSF_slab.h"
#include "ChSF_slab_inline.h"
#include "shell_etc_indices.h"



/**********************************************************************************************
 
 Computes real-space total energy int |f(x)|^2/2 dx;
 Note dx=sin(theta) dtheta; dtheta = pi/N with equispaced theta
 
 ***********************************************************************************************/

Real ChSF_SLAB::Get_local_energy_real_space(Array<Real,3> Ar)
{
    int rx;
    Real factor,x;
    Real sum = 0.0;

	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
                rx = Get_rx_real_space(lx);
                factor = sin((rx)*M_PI/Nx);
                sum += my_pow(Ar(ly,lx,lz),2)*factor;
            }
	
	return sum*(M_PI/(2* Real(Nx) * Real(Ny) * Real(Nz)));
}


Real ChSF_SLAB::Get_local_energy(Array<Complex,3> A)  
{
	
	Real total = 4*Array_sqr(A);
	
	// subtractions |A(ky=0,:,:)|^2
	total -= 2*Array_sqr(A(0, Range::all(), Range::all()));
	
	// subtractions |A(:, kz=0, :)|^2
	total -= 2*Array_sqr(A(Range::all(), 0, Range::all()));
	
	// But the edge |A(:,0,0)|^2 should be added with a factor 2x0.5=1
	total += Array_sqr(A(0, 0, Range::all()));
	
	// subtractions |A(0,:,:)|^2
	if (my_id == 0) {
		total -= 2*Array_sqr(A(Range::all(), Range::all(), 0));
		
		// But the edges |A(0,0,:)|^2, |A(0,:,0)|^2 should be added with a factor 2x0.5=1
		total += Array_sqr(A(0, Range::all(), 0));
		total += Array_sqr(A(Range::all(), 0, 0));
		total -= my_pow(real(A(0,0,0)),2)/2;
	}
	
	return total;
}



/**********************************************************************************************

		Computes total A(k)*conj(B(k)) without k=0.
	
		
***********************************************************************************************/


Real ChSF_SLAB::Get_local_energy_real_space(Array<Real,3> Ar, Array<Real,3> Br)
{
    int rx;
    Real factor,x;
    Real sum = 0.0;
    
	for (int ly=0; ly<Ar.extent(0); ly++)
        for (int lz=0; lz<Ar.extent(1); lz++)
            for (int lx=0; lx<Ar.extent(2); lx++) {
                rx = Get_rx_real_space(lx);
                factor =  sin((rx)*M_PI/Nx);
                sum += factor*Ar(ly,lz,lx)*Br(ly,lz,lx);
            }
	
	return sum*(M_PI/(2* Real(Nx) * Real(Ny) * Real(Nz)));
}


Real ChSF_SLAB::Get_local_energy(Array<Complex,3> A, Array<Complex,3> B)
{
	Real  total = 4*mydot(A,B);
	
	// subtractions |A(:, :, 0)|^2
	total -= 2*mydot(A(Range::all(), 0, Range::all()), B(Range::all(), 0, Range::all()));
	
	// subtractions |A(:,0,:)|^2
	total -= 2*mydot(A(0, Range::all(), Range::all()), B(0, Range::all(), Range::all()));
	
	total += mydot(A(0, 0, Range::all()), B(0, 0, Range::all()));
	
	// subtractions |A(0,:,:)|^2
	if (my_id == 0) {
		total -= 2*mydot(A(Range::all(), Range::all(), 0), B(Range::all(), Range::all(), 0));
		
		total += mydot(A(0, Range::all(), 0), B(0, Range::all(), 0));
		total += mydot(A(Range::all(), 0, 0), B(Range::all(), 0, 0));
		total -= real(A(0,0,0))*real(B(0,0,0))/2;
	}
	
	return total;
}


void ChSF_SLAB::Compute_local_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Real &local_helicity1, Real &local_helicity2,
 Real &local_H1k1, Real &local_H1k2
)
{}

void ChSF_SLAB::Compute_total_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Real &total_helicity1, Real &total_helicity2,
 Real &total_H1k1, Real &total_H1k2
 )
{}

void ChSF_SLAB::Compute_local_shell_spectrum_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Real,1> local_H1k1, Array<Real,1> local_H1k2, Array<Real,1> local_H1k3,
 Array<Real,1> local_H1k_count
 )
{}


void ChSF_SLAB::Compute_shell_spectrum_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Real,1> H1k1, Array<Real,1> H1k2, Array<Real,1> H1k3
 )
{}

//*********************************************************************************************

// Not for 2D
void ChSF_SLAB::Compute_local_ring_spectrum_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Real,2> local_H1k1, Array<Real,2> local_H1k2, Array<Real,2> local_H1k3
 )
{}

//
//

void ChSF_SLAB::Compute_ring_spectrum_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
Array<Real,2> H1k1, Array<Real,2> H1k2, Array<Real,2> H1k3
 )
{}


//*********************************************************************************************

void ChSF_SLAB::Compute_local_cylindrical_ring_spectrum_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
  Array<Real,2> local_H1k1,  Array<Real,2> local_H1k2
 )
{}

//
//
void ChSF_SLAB::Compute_cylindrical_ring_spectrum_helicity
(
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az,
 Array<Real,2> H1k1,  Array<Real,2> H1k2
 )
{}



//*****************************  End of ChFF_slab_energy.cc ****************************







//*****************************  End of scft_energy.cc ****************************************






