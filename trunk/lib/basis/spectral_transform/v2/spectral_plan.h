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


#ifndef _SPECTRAL_PLAN_H
#define _SPECTRAL_PLAN_H

#include "spectral_def_vars.h"


//*********************************************************************************************	

class SpectralPlan
{
protected:

	enum orientation{
		X,
		Y,
		Z,
		XZ,
		YZ
	};

	Array<complx,3> X_3d;

	MPI_Request *request;
	MPI_Status *status;


	MPI_Datatype MPI_Vector_z_strip_send;
	MPI_Datatype MPI_Vector_z_strip_recv;

	MPI_Datatype MPI_Struct_yz_plane_block;


	template<typename T>
	void Isend_yz_x(Array<T,3> A);

	template<typename T>
	void Isend_yz_y(Array<T,3> Ar);

	template<typename T>

	void Isend_xz_y(Array<T,3> A);

	template<typename T>
	void Recv_yz_y(Array<T,3> A);

	template<typename T>
	void Recv_xz_x(Array<T,3> A);
	
	template<typename T>
	void Recv_yz_x(Array<T,3> A);

	template<typename T1, typename T2>
	void Alltoall(Array<T1,3> A1, Array<T2,3> A2, int orientation);

	void Zero_pad_last_plane(Array<DP,3> Ar);

	void SinCostr_x(char sincostr_switch, Array<DP,3> Ar);
	void ISinCostr_x(char sincostr_switch, Array<DP,3> Ar);

	void ArrayShiftRight(Array<DP,3> Ar, int axis);
	void ArrayShiftLeft(Array<DP,3> Ar, int axis);

public:
	int my_id;
	int numprocs;

	size_t Nx;
	size_t Ny;
	size_t Nz;

	size_t local_Nx;
	size_t local_Ny;
	size_t local_Nz;

	size_t local_Nx_start;
	size_t local_Ny_start;
	size_t local_Nz_start;



	SpectralPlan(int my_id, int numprocs, size_t Nx, size_t Ny, size_t Nz);

	//3D FFF
	virtual void Forward_transform(Array<DP,3> Ar, Array<complx,3> A);
	virtual void Inverse_transform(Array<complx,3> A, Array<DP,3> Ar);

	//3D SFF, SSF, SSS
	virtual void Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A);
	virtual void Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar);

	//3D Transpose
	virtual void Transpose(Array<DP,3> Ar, Array<complx,3> A);
	virtual void Transpose(Array<complx,3> A, Array<DP,3> Ar);

	//2D FFF
	virtual void Forward_transform(Array<DP,2> Ar, Array<complx,2> A);
	virtual void Inverse_transform(Array<complx,2> A, Array<DP,2> Ar);

	//2D SFF, SSF, SSS
	virtual void Forward_transform(string sincostr_option, Array<DP,2> Ar, Array<complx,2> A);
	virtual void Inverse_transform(string sincostr_option, Array<complx,2> A, Array<DP,2> Ar);
	
	//2D Transpose
	virtual void Transpose(Array<DP,2> Ar, Array<complx,2> A);
	virtual void Transpose(Array<complx,2> A, Array<DP,2> Ar);

};

//Isend
template<typename T>
void SpectralPlan::Isend_yz_x(Array<T,3> A)
{
	for (int lx=0; lx<local_Nx; lx++)
		for (int p=0; p<numprocs; p++)
			MPI_Isend(A(lx, Range(p*local_Ny, (p+1)*local_Ny-1), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, &request[p*local_Nx + lx]);
}

template<typename T>
void SpectralPlan::Isend_yz_y(Array<T,3> Ar)
{
	for (int lx=0; lx<local_Nx; lx++)
		for (int p=0; p<numprocs; p++)
			MPI_Isend(Ar(p*local_Nx + lx, Range::all(), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, &request[p*local_Nx + lx]);
}

template<typename T>
void SpectralPlan::Isend_xz_y(Array<T,3> A)
{
 	for (int ly=0; ly<local_Ny; ly++)
		for (int p=0; p<numprocs; p++)
			MPI_Isend(A(Range(p*local_Nx,(p+1)*local_Nx-1),ly,Range::all()).data(), 1, MPI_Vector_z_strip_send, p, ly, MPI_COMM_WORLD, &request[p*local_Ny + ly]);
}


//Recv
template<typename T>
void SpectralPlan::Recv_xz_x(Array<T,3> A)
{
	for (int ly=0; ly<local_Ny; ly++)
		for (int p=0; p<numprocs; p++)
			MPI_Recv(A(Range::all(),p*local_Ny+ly,Range::all()).data(), 1, MPI_Vector_z_strip_recv, p, ly, MPI_COMM_WORLD, status);

	MPI_Waitall(Ny, request, status);
}

template<typename T>
void SpectralPlan::Recv_yz_x(Array<T,3> A)
{
	for (int lx=0; lx<local_Nx; lx++)
	{
		for (int p=0; p<numprocs; p++)
		{
			MPI_Recv(A(lx, Range(p*local_Ny, (p+1)*local_Ny-1), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, status);

		}
	}
	MPI_Waitall(Nx, request, status);
}

template<typename T>
void SpectralPlan::Recv_yz_y(Array<T,3> A)
{
	for (int lx=0; lx<local_Nx; lx++)
		for (int p=0; p<numprocs; p++) {
			MPI_Recv(A(p*local_Nx + lx, Range::all(), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, status);
		}
	MPI_Waitall(Nx, request, status);
}




//Alltoall
template<typename T1, typename T2>
void SpectralPlan::Alltoall(Array<T1,3> A1, Array<T2,3> A2, int orientation)
{
	switch (orientation)
	{
		case YZ:
			for (int lx=0; lx<local_Nx; lx++)
				MPI_Alltoall(A1(lx,Range::all(),Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, A2(lx,Range::all(),Range::all()).data(), 1, MPI_Struct_yz_plane_block, MPI_COMM_WORLD);
			break;

		case XZ:
			for (int lx=0; lx<local_Nx; lx++)
				MPI_Alltoall(A1(lx,Range::all(),Range::all()).data(), 1, MPI_Struct_yz_plane_block, A2(lx,Range::all(),Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, MPI_COMM_WORLD);
			break;
	}
}

#endif

//******************************** End of field_basic.h  **************************************



