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


#ifndef _SPECTRAL_PLAN_SLAB_2D_H
#define _SPECTRAL_PLAN_SLAB_2D_H

#include "spectral_plan.h"


//*********************************************************************************************	

class SpectralPlan_Slab_2D : public SpectralPlan
{
protected:

	enum orientation{
		Z,
		X
	};


	Array<complx,2> X_2d;
	Array<DP,2> Xr_2d;

	MPI_Datatype MPI_Vector_z_strip_send;
	MPI_Datatype MPI_Vector_z_strip_recv;

	MPI_Datatype MPI_Vector_xz_plane_block;
	MPI_Datatype MPI_Struct_xz_plane_block;


	void Isend_x(Array<DP,2> A);
	void Isend_x(Array<complx,2> A);

	void Isend_z(Array<DP,2> A);
	void Isend_z(Array<complx,2> A);

	void Recv_x(Array<DP,2> A);
	void Recv_x(Array<complx,2> A);

	void Recv_z(Array<DP,2> A);
	void Recv_z(Array<complx,2> A);

	template<typename T1, typename T2>
	void Alltoall(Array<T1,2> A1, Array<T2,2> A2, int orientation);

	void Zero_pad_last_plane(Array<DP,2> Ar);
	void Zero_pad_last_plane(Array<complx,2> Ar);

	void SinCostr_x(char sincostr_option, Array<DP,2> Ar);
	void ISinCostr_x(char sincostr_option, Array<DP,2> Ar);

	void ArrayShiftRight_basic(void *data, TinyVector<int,2> shape, int axis);
	void ArrayShiftLeft_basic(void *data, TinyVector<int,2> shape, int axis);

	void ArrayShiftRight(Array<DP,2> Ar, int axis);
	void ArrayShiftLeft(Array<DP,2> Ar, int axis);

	void ArrayShiftRight(Array<complx,2> Ar, int axis);
	void ArrayShiftLeft(Array<complx,2> Ar, int axis);

public:
	SpectralPlan_Slab_2D(string basis, int my_id, int numprocs, int Nx, int Nz);

	void Evaluate_time_per_step(double num_iter);
	void Evaluate_time_per_step(string sincostr_option, double num_iter);
};


//Alltoall
template<typename T1, typename T2>
void SpectralPlan_Slab_2D::Alltoall(Array<T1,2> A1, Array<T2,2> A2, int orientation)
{
	switch (orientation)
	{
		case Z:
			MPI_Alltoall(A1.data(), 1, MPI_Struct_xz_plane_block, A2.data(), 2*local_Nx*local_Nz, MPI_DP, MPI_COMM_WORLD);
			break;

		case X:
			MPI_Alltoall(A1.data(), 2*local_Nx*local_Nz, MPI_DP, A2.data(), 1, MPI_Struct_xz_plane_block, MPI_COMM_WORLD);
			break;
	}
}

#endif

//******************************** End of field_basic.h  **************************************



