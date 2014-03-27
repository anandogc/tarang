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


/*! \file  spectral_plan.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  A. G. Chatterjee, M. K. Verma
 * @version 2
 * @date March 2014
 * @bug	No known bugs
 */

#include "spectral_plan.h"

SpectralPlan::SpectralPlan(int my_id, int numprocs, size_t Nx, size_t Ny, size_t Nz): 
		my_id(my_id), numprocs(numprocs), Nx(Nx), Ny(Ny), Nz(Nz)
{
	local_Nx = Nx/numprocs;
	local_Ny = Ny/numprocs;
	local_Nz = (Nz/2+1);

	local_Nx_start = my_id * local_Nx;
	local_Ny_start = my_id * local_Ny;
	local_Nz_start = 0;

	
	request = new MPI_Request[max(Nx,Ny)];
	status = new MPI_Status[max(Nx,Ny)];

	//Vector types required during Isend-Recv
	MPI_Type_vector(local_Nx,Nz+2,local_Ny*(Nz+2),MPI_DP,&MPI_Vector_z_strip_send);
	MPI_Type_commit(&MPI_Vector_z_strip_send);

	MPI_Type_vector(local_Nx,Nz+2,Ny*(Nz+2),MPI_DP,&MPI_Vector_z_strip_recv);
	MPI_Type_commit(&MPI_Vector_z_strip_recv);


	//Required for Alltoall transpose
	int block_lengths[2];
	
	MPI_Aint displacements[2];
	MPI_Datatype types[2];

	block_lengths[0]=(Nz+2)*local_Ny;
	block_lengths[1]=1;
	
	displacements[0]=0;
	displacements[1]=(Nz+2)*local_Ny*local_Nx*sizeof(DP);

	types[0]=MPI_DP;
	types[1]=MPI_UB;

	MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_yz_plane_block);
	MPI_Type_commit(&MPI_Struct_yz_plane_block);

}

void SpectralPlan::Zero_pad_last_plane(Array<DP,3> Ar)
{
	Ar(Range::all(),Range::all(),Range(Nz,Nz+1)) = DP(0.0);
}

//After SinCos transform, shift right
void SpectralPlan::ArrayShiftRight(Array<DP,3> Ar, int axis)
{	
	Ar(Range(Nx-1,1,-1),Range::all(),Range::all()) = Ar(Range(Nx-2,0,-1),Range::all(),Range::all());
	Ar(0,Range::all(),Range::all()) = 0.0;
}


//Before inverse SinCos transform, shift left
void SpectralPlan::ArrayShiftLeft(Array<DP,3> Ar, int axis)
{	
	Ar(Range(0,Nx-2),Range::all(),Range::all()) = Ar(Range(1,Nx-1),Range::all(),Range::all());
	Ar(Nx-1,Range::all(),Range::all()) = 0.0;
}

//3D FFF
void SpectralPlan::Forward_transform(Array<DP,3> Ar, Array<complx,3> A){}
void SpectralPlan::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar){}

//3D SFF, SSF, SSS
void SpectralPlan::Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A){}
void SpectralPlan::Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar){}

//3D Transpose
void SpectralPlan::Transpose(Array<DP,3> Ar, Array<complx,3> A){}
void SpectralPlan::Transpose(Array<complx,3> A, Array<DP,3> Ar){}


//2D FFF
void SpectralPlan::Forward_transform(Array<DP,2> Ar, Array<complx,2> A){}
void SpectralPlan::Inverse_transform(Array<complx,2> A, Array<DP,2> Ar){}

//2D SFF, SSF, SSS
void SpectralPlan::Forward_transform(string sincostr_option, Array<DP,2> Ar, Array<complx,2> A){}
void SpectralPlan::Inverse_transform(string sincostr_option, Array<complx,2> A, Array<DP,2> Ar){}

//2D Transpose
void SpectralPlan::Transpose(Array<DP,2> Ar, Array<complx,2> A){}
void SpectralPlan::Transpose(Array<complx,2> A, Array<DP,2> Ar){}



//******************************** End of field_basic.h  **************************************



