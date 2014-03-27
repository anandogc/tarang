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

void ArrayOps::Real_space_multiply(Array<DP,3> A, Array<DP,3> B, Array<DP,3> C)
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

void ArrayOps::Real_space_divide(Array<DP,3> A, Array<DP,3> B, Array<DP,3> C)
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


/**********************************************************************************************

	Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
	
	// For kz =0 and kz = N[3]/2 planes
	// f(-kx, -ky, 0) = conj(f(kx, ky, 0))- do for kz=N[3]/2

***********************************************************************************************/ 


/**********************************************************************************************
 
 Collect pieced of kz planes from all the procs and put them in plane
 
 ***********************************************************************************************/


void ArrayOps::Get_XY_plane(Array<complx,3> A, Array<complx,2> plane_xy, int kz)
{

    int lz;
	int data_size, full_data_size;
	
	if (global.program.decomposition == "SLAB") {
		if (numprocs == 1) { // for SLAB
			plane_xy= A(Range::all(),Range::all(),kz);
			return;
		}
		
		else {
			int num_blocks;
			int extent;
			int stride;
			int block_lengths[2];
			
			MPI_Aint displacements[2];
			MPI_Datatype types[2];
			
			MPI_Datatype MPI_Vector_block;
			MPI_Datatype MPI_Struct_block;
			
			
			num_blocks = Nx;
			extent = 2*local_Ny;
			stride = 2*Ny;
			
			MPI_Type_vector(num_blocks,extent,stride,MPI_DP,&MPI_Vector_block);
			MPI_Type_commit(&MPI_Vector_block);
			
			
			block_lengths[0]=1;
			block_lengths[1]=1;
			
			displacements[0]=0;
			displacements[1]=local_Ny*sizeof(complx);
			
			types[0]=MPI_Vector_block;
			types[1]=MPI_UB;
			
			MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_block);
			MPI_Type_commit(&MPI_Struct_block);
				
			
			lz = universal->Get_lz(kz);
			
			global.temp_array.plane_xy_inproc = A(Range::all(),Range::all(),lz);
		
			data_size = 2*local_Ny*Nx;
			
			full_data_size = 2*Nx*Ny;
			
			MPI_Gather(reinterpret_cast<DP*>(global.temp_array.plane_xy_inproc.data()), data_size, MPI_DP, reinterpret_cast<DP*>(plane_xy.data()), 1, MPI_Struct_block, 0, MPI_COMM_WORLD);
			
			MPI_Bcast(reinterpret_cast<DP*>(global.temp_array.plane_xy.data()), full_data_size, MPI_DP, 0, MPI_COMM_WORLD);
		}
	}
	
	else if (global.program.decomposition == "PENCIL") {
		if ((num_p_vert > 1)) {
			lz = universal->Get_lz(kz);
			
			global.temp_array.plane_xy_inproc = A(Range::all(), lz, Range::all());
			
			data_size = 2*local_Ny_vert*Nx;
			full_data_size = 2*Nx*Ny;
			
			MPI_Gather(reinterpret_cast<DP*>(global.temp_array.plane_xy_inproc.data()), data_size, MPI_DP, reinterpret_cast<DP*>(plane_xy.data()), data_size, MPI_DP, 0, global.mpi.MPI_COMM_VERT_SEGMENT);
			
			MPI_Bcast(reinterpret_cast<DP*>(global.temp_array.plane_xy.data()), full_data_size, MPI_DP, 0, global.mpi.MPI_COMM_VERT_SEGMENT);
		}
		
		else {
			plane_xy = A(Range::all(),kz,Range::all());
		}
		
	}
}


