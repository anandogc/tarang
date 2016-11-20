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


/*! \file array_basic.cc
 *
 * @sa field_basic.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */


#include "ArrayOps.h"


/**********************************************************************************************
 
 C = A.B (assuming A and B are real)
 Term by term multiplication.
 
 ***********************************************************************************************/

void ArrayOps::Real_space_multiply(Array<Real,3> A, Array<Real,3> B, Array<Real,3> C)
{
	C = A * B;
	
 /*   #pragma omp parallel for
    for (int ly=0; ly<A.extent(0); ly++)
        for (int lz=0; lz<A.extent(1); lz++)
            for (int lx=0; lx<A.extent(2); lx++){
                real(C(ly, lz, lx)) = real(A(ly, lz, lx)) * real(B(ly, lz, lx));
                imag(C(ly, lz, lx)) = imag(A(ly, lz, lx)) * imag(B(ly, lz, lx));
            } */
//	real(C) = real(A) * real(B);
//	imag(C) = imag(A) * imag(B);
  //  cout << "A,B,C = " << A.data() << " " << B.data() << " " << C.data() << endl;
}

void ArrayOps::Real_space_divide(Array<Real,3> A, Array<Real,3> B, Array<Real,3> C)
{

	C = A / B;
	
   /* #pragma omp parallel for
    for (int ly=0; ly<A.extent(0); ly++)
        for (int lz=0; lz<A.extent(1); lz++)
            for (int lx=0; lx<A.extent(2); lx++){
                real(C(ly, lz, lx)) = real(A(ly, lz, lx)) / real(B(ly, lz, lx));
                imag(C(ly, lz, lx)) = imag(A(ly, lz, lx)) / imag(B(ly, lz, lx));
            } */
	// real(C) = real(A) / real(B);
	// imag(C) = imag(A) / imag(B);
}
