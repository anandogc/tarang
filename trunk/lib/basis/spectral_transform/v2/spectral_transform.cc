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
/*! \file fourier.cc 
 * 
 * @sa fourier.hFTr2c_yz(plane_yz);
 * 
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		ArrayIFFT
 */ 


#include "spectral_transform.h"


//*********************************************************************************************

void SpectralTransform::Init(string basis, string decomposition, int plan_id, int Nx, int Ny, int Nz, int num_p_hor)
{

	int my_id, numprocs;	

	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


	if (basis=="FFF") {
		switch (plan_id) {
			case 1:
				plan = new FFF_slab_Isend_Recv_overlap_Isend_forward_3D(my_id, numprocs, Nx, Ny, Nz);
				break;

			case 2:
				plan = new FFF_slab_Isend_Recv_overlap_Isend_both_3D(my_id, numprocs, Nx, Ny, Nz);
			
			case 3:
				plan = new FFF_slab_Alltoall_3D(my_id, numprocs, Nx, Ny, Nz);
				break;
		}
	}		
	else if (basis=="SFF") { 
		switch (plan_id) {
			case 1:
				plan = new SFF_slab_Isend_Recv_3D(my_id, numprocs, Nx, Ny, Nz);
				break;
			
			case 2:
				plan = new SFF_slab_Alltoall_3D(my_id, numprocs, Nx, Ny, Nz);
				break;
		}
	}
	
	Nx=plan->Nx;
	Ny=plan->Ny;
	Nz=plan->Nz;
	
	local_Nx=plan->local_Nx;
	local_Ny=plan->local_Ny;
	local_Nz=plan->local_Nz;
	
	local_Nx_start=plan->local_Nx_start;
	local_Ny_start=plan->local_Ny_start;
	local_Nz_start=plan->local_Nz_start;
}

//********************************	End of four_tr.cc *****************************************


//3D
void SpectralTransform::Forward_transform(Array<DP,3> Ar, Array<complx,3> A){
	plan->Forward_transform(Ar, A);
}
void SpectralTransform::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar){
	plan->Inverse_transform(A, Ar);
}

void SpectralTransform::Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A){
	plan->Forward_transform(sincostr_option, Ar, A);
}
void SpectralTransform::Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar){
	plan->Inverse_transform(sincostr_option, A, Ar);
}

void SpectralTransform::Transpose(Array<DP,3> Ar, Array<complx,3> A){
	plan->Transpose(Ar, A);
}
void SpectralTransform::Transpose(Array<complx,3> A, Array<DP,3> Ar){
	plan->Transpose(A, Ar);
}

//2D
void SpectralTransform::Forward_transform(Array<DP,2> Ar, Array<complx,2> A){
	plan->Forward_transform(Ar, A);
}
void SpectralTransform::Inverse_transform(Array<complx,2> A, Array<DP,2> Ar){
	plan->Inverse_transform(A, Ar);
}

void SpectralTransform::Forward_transform(string sincostr_option, Array<DP,2> Ar, Array<complx,2> A){
	plan->Forward_transform(sincostr_option, Ar, A);
}
void SpectralTransform::Inverse_transform(string sincostr_option, Array<complx,2> A, Array<DP,2> Ar){
	plan->Inverse_transform(sincostr_option, A, Ar);
}

void SpectralTransform::Transpose(Array<DP,2> Ar, Array<complx,2> A){
	plan->Transpose(Ar, A);
}
void SpectralTransform::Transpose(Array<complx,2> A, Array<DP,2> Ar){
	plan->Transpose(A, Ar);
}

