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

#include "spectral_plan_slab_2d.h"

SpectralPlan_Slab_2D::SpectralPlan_Slab_2D(string basis, int my_id, int numprocs, int Nx, int Nz): SpectralPlan(basis, my_id, numprocs, Nx, Nz)
{
	if (basis!="SS")
	{
		if ((Nz/2+1)/numprocs <= 0)
		{
			if (my_id==0) cerr << "Unable to divide data along last axis.\n\n"<<endl;
			exit(1);
		}

		local_Nx = Nx/numprocs;
		local_Nz = (Nz/2+1)/numprocs;

		local_Nx_start = my_id * local_Nx;
		local_Nz_start = my_id * local_Nz;

		request = new MPI_Request[numprocs];
		status = new MPI_Status[numprocs];

		//Required for Alltoall transpose
		int count=local_Nx;
		int blocklength=2*local_Nz;
		int stride=Nz+2;

		MPI_Type_vector(count, blocklength, stride, MPI_DP, &MPI_Vector_xz_plane_block);
		MPI_Type_commit(&MPI_Vector_xz_plane_block);

		int block_lengths[2];
		MPI_Aint displacements[2];
		MPI_Datatype types[2];

		block_lengths[0]=1;
		block_lengths[1]=1;
		
		displacements[0]=0;
		displacements[1]=local_Nz*sizeof(complx);

		types[0]=MPI_Vector_xz_plane_block;
		types[1]=MPI_UB;

		MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_xz_plane_block);
		MPI_Type_commit(&MPI_Struct_xz_plane_block);
	}
	else if (basis=="SS")
	{
		if ((Nz/2)/numprocs <= 0)
		{
			if (my_id==0) cerr << "Unable to divide data along last axis.\n\n"<<endl;
			exit(1);
		}

		local_Nx = Nx/numprocs;
		local_Nz = (Nz/2)/numprocs;

		local_Nx_start = my_id * local_Nx;
		local_Nz_start = my_id * local_Nz;

		request = new MPI_Request[numprocs];
		status = new MPI_Status[numprocs];

		//Required for Alltoall transpose
		int count=local_Nx;
		int blocklength=2*local_Nz;
		int stride=Nz;

		MPI_Type_vector(count, blocklength, stride, MPI_DP, &MPI_Vector_xz_plane_block);
		MPI_Type_commit(&MPI_Vector_xz_plane_block);

		int block_lengths[2];
		MPI_Aint displacements[2];
		MPI_Datatype types[2];

		block_lengths[0]=1;
		block_lengths[1]=1;
		
		displacements[0]=0;
		displacements[1]=local_Nz*sizeof(complx);

		types[0]=MPI_Vector_xz_plane_block;
		types[1]=MPI_UB;

		MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_xz_plane_block);
		MPI_Type_commit(&MPI_Struct_xz_plane_block);
	}
}


void SpectralPlan_Slab_2D::Evaluate_time_per_step(double num_iter)
{

	DP local_time_per_step=-MPI_Wtime();
	for (double n=0; n<num_iter; n++)
	{
		Forward_transform(Xr_2d, X_2d);
		Inverse_transform(X_2d, Xr_2d);
	}
	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	time_per_step/=num_iter;


//For debugging uncomment lower section and comment upper section.
/*
	DP local_energy=(sum(sqr(Xr_2d))/(Nx*Nz));
	DP total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id==0)
		cout << "Initial Energy = " << total_energy << endl;

	Array<DP,2> Xr_2d_init(Xr_2d.shape());
	Xr_2d_init = Xr_2d;
	DP local_time_per_step=-MPI_Wtime();

	//

	// MPI_Barrier(MPI_COMM_WORLD);	
	// for (int i=0; i<numprocs; i++)
	// {
	// 	if (my_id==i)
	// 		cout << fixed << Xr_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }

	X_2d=0;

	Forward_transform(Xr_2d, X_2d);

	// for (int i=0; i<numprocs; i++)
	// {
	// 	if (my_id==i)
	// 		cout << fixed << X_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }


	for (int lx=0; lx<X_2d.extent(0); lx++)
		for (int lz=0; lz<X_2d.extent(1); lz++)
			if (abs(X_2d(lx,lz))>1e-6)
				cout << "my_id = " << my_id << " X(" << lx << "," << lz << ") = " << fixed << X_2d(lx,lz) << endl;
	Xr_2d=0;

	Inverse_transform(X_2d, Xr_2d);


	//MPI_Barrier(MPI_COMM_WORLD);	

	// for (int i=0; i<numprocs; i++)
	// {
	// 	if (my_id==i)
	// 		cout << fixed << Xr_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }

	//

	local_energy=(sum(sqr(Xr_2d))/(Nx*Nz));
	total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	DP local_sum_sqr_diff=(sum(sqr(Xr_2d_init - Xr_2d)));
	DP total_sum_sqr_diff;
	MPI_Reduce(&local_sum_sqr_diff, &total_sum_sqr_diff, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id==0)
	{
		cout << "Final Energy = " << total_energy << endl;
		cout << "(sum(sqr(Xr_2d_init-Xr_2d))) = " << scientific << total_sum_sqr_diff << endl;

		// for (int rx=0; rx<local_Nx; rx++)
		// 	for (int rz=0; rz<Nz+2; rz++)
		// 		if (abs(Xr_2d_init(rx,rz) - Xr_2d(rx,rz))>1e-12)
		// 			cout << fixed << "Xr_2d_init("<<rx<<","<<rz<<") , Xr_2d("<<rx<<","<<rz<<") = " << Xr_2d_init(rx,rz) << ", " << Xr_2d(rx,rz) << endl;
	}

	//MPI_Barrier(MPI_COMM_WORLD);

	//if (local_sum_sqr_diff>1e-12)
	//	cout << "sum not matching in " << my_id << endl;

	//MPI_Barrier(MPI_COMM_WORLD);

	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	// time_per_step/=num_iter;

*/
}


void SpectralPlan_Slab_2D::Evaluate_time_per_step(string sincostr_option, double num_iter)
{
	DP local_time_per_step=-MPI_Wtime();
	for (double n=0; n<num_iter; n++) {
		Forward_transform(sincostr_option, Xr_2d, X_2d);
		Inverse_transform(sincostr_option, X_2d, Xr_2d);
	}
	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	time_per_step/=num_iter;

//For debugging uncomment lower section and comment upper section.
/*
	DP local_energy=sum(sqr(Xr_2d));
	DP total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id==0)
		cout << "Initial Energy = " << total_energy << endl;

	Array<DP,2> Xr_2d_init(Xr_2d.shape());
	Xr_2d_init = Xr_2d;
	DP local_time_per_step=-MPI_Wtime();

	//

	// MPI_Barrier(MPI_COMM_WORLD);	
	// for (int i=0; i<numprocs; i++)
	// {
	// 	if (my_id==i)
	// 		cout << fixed << Xr_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }

	X_2d=0;

	Forward_transform(sincostr_option, Xr_2d, X_2d);

	// for (int i=0; i<numprocs; i++)
	// {
	// 	if (my_id==i)
	// 		cout << fixed << X_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }


	for (int lx=0; lx<X_2d.extent(0); lx++)
		for (int lz=0; lz<X_2d.extent(1); lz++)
			if (abs(X_2d(lx,lz))>1e-6)
				cout << "my_id = " << my_id << " X(" << lx << "," << lz << ") = " << fixed << X_2d(lx,lz) << endl;
	Xr_2d=0;

	Inverse_transform(sincostr_option, X_2d, Xr_2d);


	//MPI_Barrier(MPI_COMM_WORLD);	

	// for (int i=0; i<numprocs; i++){
	// 	if (my_id==i)
	// 		cout << fixed << Xr_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }

	//

	local_energy=sum(sqr(Xr_2d));
	total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	DP local_sum_sqr_diff=(sum(sqr(Xr_2d_init - Xr_2d)));
	DP total_sum_sqr_diff;
	MPI_Reduce(&local_sum_sqr_diff, &total_sum_sqr_diff, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id==0)
	{
		cout << "Final Energy = " << total_energy << endl;
		cout << "(sum(sqr(Xr_2d_init-Xr_2d))) = " << scientific << total_sum_sqr_diff << endl;

		// for (int rx=0; rx<local_Nx; rx++)
		// 	for (int rz=0; rz<Nz+2; rz++)
		// 		if (abs(Xr_2d_init(rx,rz) - Xr_2d(rx,rz))>1e-12)
		// 			cout << fixed << "Xr_2d_init("<<rx<<","<<rz<<") , Xr_2d("<<rx<<","<<rz<<") = " << Xr_2d_init(rx,rz) << ", " << Xr_2d(rx,rz) << endl;
	}

	//MPI_Barrier(MPI_COMM_WORLD);

	//if (local_sum_sqr_diff>1e-12)
	//	cout << "sum not matching in " << my_id << endl;

	//MPI_Barrier(MPI_COMM_WORLD);

	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
*/
}

//Isend
void SpectralPlan_Slab_2D::Isend_x(Array<DP,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Isend(A(Range(p*local_Nx, (p+1)*local_Nx-1), Range::all()).data(), 2*local_Nx*local_Nz, MPI_DP, p, 0, MPI_COMM_WORLD, &request[p]);
	
}

void SpectralPlan_Slab_2D::Isend_x(Array<complx,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Isend(A(Range(p*local_Nx, (p+1)*local_Nx-1), Range::all()).data(), 2*local_Nx*local_Nz, MPI_DP, p, 0, MPI_COMM_WORLD, &request[p]);
	
}

void SpectralPlan_Slab_2D::Isend_z(Array<DP,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Isend(A(Range::all(), Range(2*p*local_Nz, 2*(p+1)*local_Nz-1)).data(), 1, MPI_Vector_xz_plane_block, p, 0, MPI_COMM_WORLD, &request[p]);
	
}

void SpectralPlan_Slab_2D::Isend_z(Array<complx,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Isend(A(Range::all(), Range(p*local_Nz, (p+1)*local_Nz-1)).data(), 1, MPI_Vector_xz_plane_block, p, 0, MPI_COMM_WORLD, &request[p]);
	
}

//Recv
void SpectralPlan_Slab_2D::Recv_x(Array<DP,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Recv(A(Range(p*local_Nx, (p+1)*local_Nx-1),Range::all()).data(), 2*local_Nx*local_Nz, MPI_DP, p, 0, MPI_COMM_WORLD, status);

	MPI_Waitall(numprocs, request, status);
}

void SpectralPlan_Slab_2D::Recv_x(Array<complx,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Recv(A(Range(p*local_Nx, (p+1)*local_Nx-1),Range::all()).data(), 2*local_Nx*local_Nz, MPI_DP, p, 0, MPI_COMM_WORLD, status);

	MPI_Waitall(numprocs, request, status);
}


void SpectralPlan_Slab_2D::Recv_z(Array<DP,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Recv(A(Range::all(), Range(2*p*local_Nz, 2*(p+1)*local_Nz-1)).data(), 1, MPI_Vector_xz_plane_block, p, 0, MPI_COMM_WORLD, status);

	MPI_Waitall(numprocs, request, status);
}

void SpectralPlan_Slab_2D::Recv_z(Array<complx,2> A)
{
	for (int p=0; p<numprocs; p++)
		MPI_Recv(A(Range::all(), Range(p*local_Nz, (p+1)*local_Nz-1)).data(), 1, MPI_Vector_xz_plane_block, p, 0, MPI_COMM_WORLD, status);

	MPI_Waitall(numprocs, request, status);
}

void SpectralPlan_Slab_2D::Zero_pad_last_plane(Array<DP,2> Ar)
{
	Ar(Range::all(),Range(Nz,Nz+1)) = DP(0.0);
}

void SpectralPlan_Slab_2D::Zero_pad_last_plane(Array<complx,2> A)
{
	A(Range::all(),Nz/2) = complx(0.0,0.0);
}

//After SinCos transform, shift right
void SpectralPlan_Slab_2D::ArrayShiftRight_basic(void *data, TinyVector<int,2> shape, int axis)
{
	Array<DP,2> Ar((DP*)(data), shape, neverDeleteData);

	if (axis=='X') {
		Ar(Range(Nx-1,1,-1),Range::all()) = Ar(Range(Nx-2,0,-1),Range::all());
		Ar(0,Range::all()) = 0.0;
	}
	else if (axis=='Z') {
		Ar(Range::all(),Range(Nz-1,1,-1)) = Ar(Range::all(),Range(Nx-2,0,-1));
		Ar(Range::all(),0) = 0.0;
	}
}

void SpectralPlan_Slab_2D::ArrayShiftRight(Array<DP,2> Ar, int axis)
{
	ArrayShiftRight_basic(Ar.data(), Ar.shape(), axis);
}

void SpectralPlan_Slab_2D::ArrayShiftRight(Array<complx,2> A, int axis)
{
	ArrayShiftRight_basic(A.data(), A.shape()*shape(1,1,2), axis);
}


//Before inverse SinCos transform, shift left
void SpectralPlan_Slab_2D::ArrayShiftLeft_basic(void *data, TinyVector<int,2> shape, int axis)
{
	Array<DP,2> Ar((DP*)(data), shape, neverDeleteData);

	if (axis=='X')	{
		Ar(Range(0,Nx-2),Range::all()) = Ar(Range(1,Nx-1),Range::all());
		Ar(Nx-1,Range::all()) = 0.0;
	}
	else if (axis=='Z') {
		Ar(Range::all(),Range(0,Nz-2)) = Ar(Range::all(),Range(1,Nz-1));
		Ar(Range::all(),Nz-1) = 0.0;
	}
}

void SpectralPlan_Slab_2D::ArrayShiftLeft(Array<DP,2> Ar, int axis)
{	
	ArrayShiftLeft_basic(Ar.data(), Ar.shape(), axis);
}

void SpectralPlan_Slab_2D::ArrayShiftLeft(Array<complx,2> A, int axis)
{	
	ArrayShiftLeft_basic(A.data(), A.shape()*shape(1,1,2), axis);
}



//******************************** End of field_basic.h  **************************************



