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

#include "spectral_plan_slab_3d.h"

SpectralPlan_Slab_3D::SpectralPlan_Slab_3D(string basis, int my_id, int numprocs, size_t Nx, size_t Ny, size_t Nz): SpectralPlan(basis, my_id, numprocs, Nx, Ny, Nz)
{
	if (basis != "SSS")
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
		int count=local_Nx;
		int blocklength=Nz+2;
		int stride=local_Ny*(Nz+2);
		MPI_Type_vector(count,blocklength,stride,MPI_DP,&MPI_Vector_z_strip_send);
		MPI_Type_commit(&MPI_Vector_z_strip_send);

		count=local_Nx;
		blocklength=Nz+2;
		stride=Ny*(Nz+2);
		MPI_Type_vector(count,blocklength,stride,MPI_DP,&MPI_Vector_z_strip_recv);
		MPI_Type_commit(&MPI_Vector_z_strip_recv);


		//Required for Alltoall transpose
		count=1;
		blocklength=(Nz+2)*local_Ny;
		stride=1;
		MPI_Type_vector(1,(Nz+2)*local_Ny,1,MPI_DP,&MPI_Vector_yz_plane_block);
		MPI_Type_commit(&MPI_Vector_yz_plane_block);

		int lower_bound=0;
		int extent=(Nz+2)*local_Ny*local_Nx*sizeof(DP);
		MPI_Type_create_resized(MPI_Vector_yz_plane_block, 0, (Nz+2)*local_Ny*local_Nx*sizeof(DP), &MPI_Vector_resized_yz_plane_block);
		MPI_Type_commit(&MPI_Vector_resized_yz_plane_block);

	}
	else
	{
		local_Nx = Nx/numprocs;
		local_Ny = Ny/numprocs;
		local_Nz = Nz/2;

		local_Nx_start = my_id * local_Nx;
		local_Ny_start = my_id * local_Ny;
		local_Nz_start = 0;

		
		/*request = new MPI_Request[max(Nx,Ny)];
		status = new MPI_Status[max(Nx,Ny)];

		//Vector types required during Isend-Recv
		MPI_Type_vector(local_Nx,Nz+2,local_Ny*(Nz+2),MPI_DP,&MPI_Vector_z_strip_send);
		MPI_Type_commit(&MPI_Vector_z_strip_send);

		MPI_Type_vector(local_Nx,Nz+2,Ny*(Nz+2),MPI_DP,&MPI_Vector_z_strip_recv);
		MPI_Type_commit(&MPI_Vector_z_strip_recv);*/


		//Required for Alltoall transpose
		MPI_Type_vector(1,Nz*local_Ny,1,MPI_DP,&MPI_Vector_yz_plane_block);
		MPI_Type_commit(&MPI_Vector_yz_plane_block);

		MPI_Type_create_resized(MPI_Vector_yz_plane_block, 0, Nz*local_Ny*local_Nx*sizeof(DP), &MPI_Vector_resized_yz_plane_block);
		MPI_Type_commit(&MPI_Vector_resized_yz_plane_block);
	}
}


void SpectralPlan_Slab_3D::Evaluate_time_per_step(double num_iter){
	DP local_time_per_step=-MPI_Wtime();
	for (double n=0; n<num_iter; n++) {
		Forward_transform(Xr_3d, X_3d);
		Inverse_transform(X_3d, Xr_3d);
	}
	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	time_per_step/=num_iter;

//For debugging uncomment lower section and comment upper section.
/*
	DP local_energy=sum(sqr(Xr_3d));
	DP total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id==0)
		cout << "Initial Energy = " << fixed << total_energy << endl;

	Array<DP,3> Xr_3d_init(Xr_3d.shape());
	Xr_3d_init = Xr_3d;

	// MPI_Barrier(MPI_COMM_WORLD);	
	// for (int i=0; i<numprocs; i++){
	// 	if (my_id==i)
	// 		cout << fixed << Xr_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }


	DP local_time_per_step=-MPI_Wtime();
	for (double n=0; n<num_iter; n++) {
		X_3d=0;
		Forward_transform(Xr_3d, X_3d);

		// for (int i=0; i<numprocs; i++){
		// 	if (my_id==i)
		// 		cout << fixed << X_2d << endl;
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// }


		for (int lx=0; lx<X_3d.extent(0); lx++)
			for (int ly=0; ly<X_3d.extent(1); ly++)
				for (int lz=0; lz<X_3d.extent(2); lz++)
					if (abs(X_3d(lx,ly,lz))>1e-6)
						cout << "my_id = " << my_id << " X(" << lx << "," << ly << "," << lz << ") = " << fixed << X_3d(lx,ly,lz) << endl;


		Xr_3d=0;
		Inverse_transform(X_3d, Xr_3d);
	}
	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	time_per_step/=num_iter;


	local_energy=sum(sqr(Xr_3d));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	DP local_sum_sqr_diff=(sum(sqr(Xr_3d_init - Xr_3d)));
	DP total_sum_sqr_diff;
	MPI_Reduce(&local_sum_sqr_diff, &total_sum_sqr_diff, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id==0) {
		cout << "Final Energy = " << fixed << total_energy << endl;
		cout << "(sum(sqr(Xr_3d_init-Xr_3d))) = " << scientific << total_sum_sqr_diff << endl;

		// for (int rx=0; rx<Nx; rx++)
		// 	for (int ry=0; ry<local_Ny; ry++)
		// 		for (int rz=0; rz<Nz+2; rz++)
		// 			if (abs(Xr_3d_init(rx,ry,rz) - Xr_3d(rx,ry,rz))>1e-12)
		// 				cout << fixed << "Xr_3d_init("<<rx<<","<<ry<<","<<rz<<") , Xr_3d("<<rx<<","<<ry<<","<<rz<<") = " << Xr_3d_init(rx,ry,rz) << ", " << Xr_3d(rx,ry,rz) << endl;
	}

	//MPI_Barrier(MPI_COMM_WORLD);

	//if (local_sum_sqr_diff>1e-12)
	//	cout << "sum not matching in " << my_id << endl;

	//MPI_Barrier(MPI_COMM_WORLD);
*/
}


void SpectralPlan_Slab_3D::Evaluate_time_per_step(string sincostr_option, double num_iter){
	DP local_time_per_step=-MPI_Wtime();
	for (double n=0; n<num_iter; n++) {
		Forward_transform(sincostr_option, Xr_3d, X_3d);
		Inverse_transform(sincostr_option, X_3d, Xr_3d);
	}
	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	time_per_step/=num_iter;

//For debugging uncomment lower section and comment upper section.
/*
	DP local_energy=sum(sqr(Xr_3d));//(8.0 * DP(Nx) * DP(Ny) * DP(Nz));
	DP total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id==0)
		cout << "Initial Energy = " << fixed << total_energy << endl;

	Array<DP,3> Xr_3d_init(Xr_3d.shape());
	Xr_3d_init = Xr_3d;

	// MPI_Barrier(MPI_COMM_WORLD);	
	// for (int i=0; i<numprocs; i++){
	// 	if (my_id==i)
	// 		cout << fixed << Xr_2d << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }


	DP local_time_per_step=-MPI_Wtime();
	for (double n=0; n<num_iter; n++) {
		X_3d=0;
		Forward_transform(sincostr_option, Xr_3d, X_3d);

		// for (int i=0; i<numprocs; i++){
		// 	if (my_id==i)
		// 		cout << fixed << X_2d << endl;
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// }

		for (int lx=0; lx<local_Nx; lx++)
			for (int ly=0; ly<Ny; ly++)
				for (int lz=0; lz<Nz/2; lz++)
					if (abs(X_3d(lx,ly,lz))>1e-6)
						cout << "my_id = " << my_id << " X(" << lx << "," << ly << "," << lz << ") = " << fixed << X_3d(lx,ly,lz) << endl;


		Xr_3d=0;
		Inverse_transform(sincostr_option, X_3d, Xr_3d);
	}
	local_time_per_step += MPI_Wtime();
	MPI_Allreduce(&local_time_per_step, &time_per_step, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	time_per_step/=num_iter;


	local_energy=sum(sqr(Xr_3d));//(8.0 * DP(Nx) * DP(Ny) * DP(Nz));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	DP local_sum_sqr_diff=(sum(sqr(Xr_3d_init - Xr_3d)));
	DP total_sum_sqr_diff;
	MPI_Reduce(&local_sum_sqr_diff, &total_sum_sqr_diff, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id==0) {
		cout << "Final Energy = " << fixed << total_energy << endl;
		cout << "(sum(sqr(Xr_3d_init-Xr_3d))) = " << scientific << total_sum_sqr_diff << endl;

		// for (int rx=0; rx<Nx; rx++)
		// 	for (int ry=0; ry<local_Ny; ry++)
		// 		for (int rz=0; rz<Nz+2; rz++)
		// 			if (abs(Xr_3d_init(rx,ry,rz) - Xr_3d(rx,ry,rz))>1e-12)
		// 				cout << fixed << "Xr_3d_init("<<rx<<","<<ry<<","<<rz<<") , Xr_3d("<<rx<<","<<ry<<","<<rz<<") = " << Xr_3d_init(rx,ry,rz) << ", " << Xr_3d(rx,ry,rz) << endl;
	}

	//MPI_Barrier(MPI_COMM_WORLD);

	//if (local_sum_sqr_diff>1e-12)
	//	cout << "sum not matching in " << my_id << endl;

	//MPI_Barrier(MPI_COMM_WORLD);
*/
}


void SpectralPlan_Slab_3D::Zero_pad_last_plane(Array<DP,3> Ar)
{
	Ar(Range::all(),Range::all(),Range(Nz,Nz+1)) = DP(0.0);
}

//After SinCos transform, shift right
void SpectralPlan_Slab_3D::ArrayShiftRight_basic(void *data, TinyVector<int,3> shape, int axis)
{	
	Array<DP,3> Ar((DP*)(data), shape, neverDeleteData);

	if (axis == 'X') {
		Ar(Range(Nx-1,1,-1),Range::all(),Range::all()) = Ar(Range(Nx-2,0,-1),Range::all(),Range::all());
		Ar(0,Range::all(),Range::all()) = 0.0;
	}
	else if (axis == 'Y') {
		Ar(Range::all(),Range(Ny-1,1,-1),Range::all()) = Ar(Range::all(),Range(Ny-2,0,-1),Range::all());
		Ar(Range::all(),0,Range::all()) = 0.0;	
	}
	else if (axis == 'Z') {
		Ar(Range::all(),Range::all(),Range(Nz-1,1,-1)) = Ar(Range::all(),Range::all(),Range(Nz-2,0,-1));
		Ar(Range::all(),Range::all(),0) = 0.0;	
	}
}

void SpectralPlan_Slab_3D::ArrayShiftRight(Array<DP,3> Ar, int axis)
{	
	ArrayShiftRight_basic(Ar.data(), Ar.shape(), axis);
}

void SpectralPlan_Slab_3D::ArrayShiftRight(Array<complx,3> Ar, int axis)
{	
	ArrayShiftRight_basic(Ar.data(), Ar.shape()*shape(1,1,2), axis);
}

//Before inverse SinCos transform, shift left
void SpectralPlan_Slab_3D::ArrayShiftLeft_basic(void *data, TinyVector<int,3> shape, int axis)
{
	Array<DP,3> Ar((DP*)(data), shape, neverDeleteData);

	if (axis == 'X') {
		Ar(Range(0,Nx-2),Range::all(),Range::all()) = Ar(Range(1,Nx-1),Range::all(),Range::all());
		Ar(Nx-1,Range::all(),Range::all()) = 0.0;	
	}
	else if (axis == 'Y') {
		Ar(Range::all(),Range(0,Ny-2),Range::all()) = Ar(Range::all(),Range(1,Ny-1),Range::all());
		Ar(Range::all(),Ny-1,Range::all()) = 0.0;

	}
	else if (axis == 'Z') {
		Ar(Range::all(),Range::all(),Range(0,Nz-2)) = Ar(Range::all(),Range::all(),Range(1,Nz-1));
		Ar(Range::all(),Range::all(),Nz-1) = 0.0;	
	}
}

void SpectralPlan_Slab_3D::ArrayShiftLeft(Array<DP,3> Ar, int axis)
{	
	ArrayShiftLeft_basic(Ar.data(), Ar.shape(), axis);
}

void SpectralPlan_Slab_3D::ArrayShiftLeft(Array<complx,3> Ar, int axis)
{	
	ArrayShiftLeft_basic(Ar.data(), Ar.shape()*shape(1,1,2), axis);
}


//******************************** End of field_basic.h  **************************************
