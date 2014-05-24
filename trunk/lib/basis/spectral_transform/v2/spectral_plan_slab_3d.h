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


#ifndef _SPECTRAL_PLAN_SLAB_3D_H
#define _SPECTRAL_PLAN_SLAB_3D_H

#include "spectral_plan.h"


//*********************************************************************************************	

class SpectralPlan_Slab_3D : public SpectralPlan
{
protected:

	enum orientation{
		XZ,
		YZ
	};


	Array<complx,3> X_3d;
	Array<DP,3> Xr_3d;

	MPI_Datatype MPI_Vector_z_strip_send;
	MPI_Datatype MPI_Vector_z_strip_recv;

	MPI_Datatype MPI_Vector_yz_plane_block;
	MPI_Datatype MPI_Vector_resized_yz_plane_block;


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

	void SinCostr_x(char sincostr_option, Array<DP,3> Ar);
	void ISinCostr_x(char sincostr_option, Array<DP,3> Ar);

	void ArrayShiftRight_basic(void *data, TinyVector<int,3> shape, int axis);
	void ArrayShiftRight(Array<DP,3> Ar, int axis);
	void ArrayShiftRight(Array<complx,3> Ar, int axis);

	void ArrayShiftLeft_basic(void *data, TinyVector<int,3> shape, int axis);
	void ArrayShiftLeft(Array<DP,3> Ar, int axis);
	void ArrayShiftLeft(Array<complx,3> Ar, int axis);

public:
	SpectralPlan_Slab_3D(string basis, int my_id, int numprocs, size_t Nx, size_t Ny, size_t Nz);

	void Evaluate_time_per_step(double num_iter);
	void Evaluate_time_per_step(string sincostr_option, double num_iter);
};

//Isend
template<typename T>
void SpectralPlan_Slab_3D::Isend_yz_x(Array<T,3> A)
{
	for (int lx=0; lx<local_Nx; lx++)
		for (int p=0; p<numprocs; p++)
			MPI_Isend(A(lx, Range(p*local_Ny, (p+1)*local_Ny-1), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, &request[p*local_Nx + lx]);
}

template<typename T>
void SpectralPlan_Slab_3D::Isend_yz_y(Array<T,3> Ar)
{
	for (int lx=0; lx<local_Nx; lx++)
		for (int p=0; p<numprocs; p++)
			MPI_Isend(Ar(p*local_Nx + lx, Range::all(), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, &request[p*local_Nx + lx]);
}

template<typename T>
void SpectralPlan_Slab_3D::Isend_xz_y(Array<T,3> A)
{
 	for (int ly=0; ly<local_Ny; ly++)
		for (int p=0; p<numprocs; p++)
			MPI_Isend(A(Range(p*local_Nx,(p+1)*local_Nx-1),ly,Range::all()).data(), 1, MPI_Vector_z_strip_send, p, ly, MPI_COMM_WORLD, &request[p*local_Ny + ly]);
}


//Recv
template<typename T>
void SpectralPlan_Slab_3D::Recv_xz_x(Array<T,3> A)
{
	for (int ly=0; ly<local_Ny; ly++)
		for (int p=0; p<numprocs; p++)
			MPI_Recv(A(Range::all(),p*local_Ny+ly,Range::all()).data(), 1, MPI_Vector_z_strip_recv, p, ly, MPI_COMM_WORLD, status);

	MPI_Waitall(Ny, request, status);
}

template<typename T>
void SpectralPlan_Slab_3D::Recv_yz_x(Array<T,3> A)
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
void SpectralPlan_Slab_3D::Recv_yz_y(Array<T,3> A)
{
	for (int lx=0; lx<local_Nx; lx++)
		for (int p=0; p<numprocs; p++) {
			MPI_Recv(A(p*local_Nx + lx, Range::all(), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, status);
		}
	MPI_Waitall(Nx, request, status);
}




//Alltoall
/*template<typename T1, typename T2>
void SpectralPlan_Slab_3D::Alltoall(Array<T1,3> A1, Array<T2,3> A2, int orientation)
{
}*/

template<typename T1, typename T2>
void SpectralPlan_Slab_3D::Alltoall(Array<T1,3> A1, Array<T2,3> A2, int orientation)
{
	if (basis != "SSS")
	{
		switch (orientation)
		{
			case YZ:
				for (int lx=0; lx<local_Nx; lx++)
					MPI_Alltoall(A1(lx,Range::all(),Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, A2(lx,Range::all(),Range::all()).data(), 1, MPI_Vector_resized_yz_plane_block, MPI_COMM_WORLD);
				break;

			case XZ:
				for (int lx=0; lx<local_Nx; lx++)
					MPI_Alltoall(A1(lx,Range::all(),Range::all()).data(), 1, MPI_Vector_resized_yz_plane_block, A2(lx,Range::all(),Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, MPI_COMM_WORLD);
				break;
		}
	}
	else
	{
		switch (orientation)
		{
			case YZ:
				for (int lx=0; lx<local_Nx; lx++)
					MPI_Alltoall(A1(lx,Range::all(),Range::all()).data(), Nz*local_Ny, MPI_DP, A2(lx,Range::all(),Range::all()).data(), 1, MPI_Vector_resized_yz_plane_block, MPI_COMM_WORLD);
				break;

			case XZ:
				for (int lx=0; lx<local_Nx; lx++)
					MPI_Alltoall(A1(lx,Range::all(),Range::all()).data(), 1, MPI_Vector_resized_yz_plane_block, A2(lx,Range::all(),Range::all()).data(), Nz*local_Ny, MPI_DP, MPI_COMM_WORLD);
				break;
		}
	}
}

#endif

//******************************** End of field_basic.h  **************************************



