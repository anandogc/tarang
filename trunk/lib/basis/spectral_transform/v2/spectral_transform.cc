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
void SpectralTransform::Init(string basis, int Nx, int Nz)
{

	int my_id, numprocs;	

	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	double num_iter=1;

	plan=NULL;

	vector<SpectralPlan*> test_plan;


	if (basis=="FFW") {
		test_plan.resize(1);
		test_plan[0] = new FFFW_slab_transposed_order_2D(my_id, numprocs, num_iter, Nx, Nz);
	}
	else if (basis=="FF") {
		test_plan.resize(1);
		test_plan[0] = new FFF_slab_Alltoall_2D(my_id, numprocs, num_iter, Nx, Nz);
	}
	else if (basis=="SF") {
		test_plan.resize(1);
		test_plan[0] = new SF_slab_Alltoall_2D(my_id, numprocs, num_iter, Nx, Nz);
	}
	else if (basis=="SS") {
		test_plan.resize(1);
		test_plan[0] = new SS_slab_Alltoall_2D(my_id, numprocs, num_iter, Nx, Nz);
	}

	int min_time_plan_index = -1;
	double min_time = numeric_limits<double>::max();

	for (int i=0; i<test_plan.size(); i++)
		if (test_plan[i]->time_per_step < min_time)
			min_time_plan_index=i;

	for (int i=0; i<test_plan.size(); i++)
	{
		if (my_id==0) cout << "Algo " << i << " time = " << test_plan[i]->time_per_step << '\n';
		if (i!=min_time_plan_index)
		{
			// if (my_id==0) cout << "deleting plan " << i << endl;
			delete test_plan[i];
		}
	}

	// if (my_id==0) cout << "Selected plan: " << basis << ", algo " << min_time_plan_index << endl;


	if (test_plan.size()>0) // A plan has been selected
	{
		plan=test_plan[min_time_plan_index];

		Nx=plan->Nx;
		Ny=1;
		Nz=plan->Nz;
		
		local_Nx=plan->local_Nx;
		local_Ny=1;
		local_Nz=plan->local_Nz;
		
		local_Nx_start=plan->local_Nx_start;
		local_Ny_start=0;
		local_Nz_start=plan->local_Nz_start;
	}
	else 
	{
		if (my_id==0)
			cerr << "Invalid parameters. No plan Initialized. " << endl;
		exit(1);
	}
}

void SpectralTransform::Init(string basis, int Nx, int Ny, int Nz)
{

	int my_id, numprocs;	

	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int num_iter=1;

	// if (my_id==0) cout << "num_iter = " << num_iter << '\n';

	vector<SpectralPlan*> test_plan;


	if (basis=="FFFW") {
		test_plan.resize(1);
		test_plan[0] = new FFFW_slab_transposed_order_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
	}
	else if (basis=="FFF") {
		test_plan.resize(3);
		test_plan[0] = new FFF_slab_Isend_Recv_overlap_Isend_forward_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
		test_plan[1] = new FFF_slab_Isend_Recv_overlap_Isend_both_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
		test_plan[2] = new FFF_slab_Alltoall_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
	}
	else if (basis=="SFF") {
		test_plan.resize(2);
		test_plan[0] = new SFF_slab_Isend_Recv_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
		test_plan[1] = new SFF_slab_Alltoall_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
	}
	else if (basis=="SSF") {
		test_plan.resize(1);
		test_plan[0] = new SSF_slab_Alltoall_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
	}
	else if (basis=="SSS") {
		test_plan.resize(1);
		test_plan[0] = new SSS_slab_Alltoall_3D(my_id, numprocs, num_iter, Nx, Ny, Nz);
	}
	int min_time_plan_index = -1;
	double min_time = numeric_limits<double>::max();

	for (int i=0; i<test_plan.size(); i++)
		if (test_plan[i]->time_per_step < min_time)
			min_time_plan_index=i;

	for (int i=0; i<test_plan.size(); i++)
	{
		// if (my_id==0) cout << "Algo " << i << " time = " << test_plan[i]->time_per_step << '\n';
		if (i!=min_time_plan_index)
		{
			// if (my_id==0) cout << "deleting plan " << i << endl;
			delete test_plan[i];
		}
	}

	// if (my_id==0) cout << "Selected plan: " << basis << ", algo " << min_time_plan_index << endl;


	if (test_plan.size()>0) // A plan has been selected
	{
		plan=test_plan[min_time_plan_index];

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
	else {
		if (my_id==0)
			cerr << "Invalid parameters. No plan Initialized. " << endl;
		exit(1);
	}
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

