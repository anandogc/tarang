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
 * along with Tarang-2; if 2D Not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */
/*! \file fourier.cc 
 * 
 * @sa fourier.h
 * 
 * @author  A. G. Chatterjee
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		ArrayIFFT
 */ 

#include "spectral_transform.h"
#include <iomanip>

SpectralTransform spectralTransform;

int my_id,numprocs;
int N_iter;
int Nx,Ny,Nz;

string basis_type, sincostr_option;

void test_FFFW_2d(char** argv);
void test_FFFW_3d(char** argv);

void test_FFF_slab_2d(char** argv);
void test_FFF_slab_3d(char** argv);
void test_FFF_pencil_3d(char** argv);

void test_SFF_2d(char** argv);
void test_SFF_3d(char** argv);
void test_SFF_pencil_3d(char** argv);

void test_SSF_3d(char** argv);
void test_SSF_pencil_3d(char** argv);

void test_SS_2d(char** argv);
void test_SSS_3d(char** argv);
void test_SSS_pencil_3d(char** argv);

DP f(int rx, int ry, int rz);


int main(int argc, char** argv)
{

	MPI_Init(&argc, &argv);



	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	basis_type=argv[1];
	
	if (basis_type=="FFW")
		test_FFFW_2d(argv);
	else if (basis_type=="FF")
		test_FFF_slab_2d(argv);
	else if (basis_type=="SF")
		test_SFF_2d(argv);
	else if (basis_type=="SS")
		test_SS_2d(argv);

	else if (basis_type=="FFFW")
		test_FFFW_3d(argv);
	else if (basis_type=="FFF")
		// test_FFF_slab_3d(argv);
		test_FFF_pencil_3d(argv);
	else if (basis_type=="SFF")
		//test_SFF_3d(argv);
		test_SFF_pencil_3d(argv);
	else if (basis_type=="SSF")
		//test_SSF_3d(argv);
		test_SSF_pencil_3d(argv);
	else if (basis_type=="SSS")
		//test_SSS_3d(argv);
		test_SSS_pencil_3d(argv);

	MPI_Finalize();
	return 0;
}

void test_FFFW_2d(char** argv)
{
	Nx=atoi(argv[2]);
	Nz=atoi(argv[3]);


	if (my_id==0) cout << "numprocs, Nx, Nz = " << numprocs << " " << Nx << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Nz);
}

void test_FFFW_3d(char** argv)
{
	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);

	//int N_iter=atoi(argv[3]);


	if (my_id==0) cout << "numprocs, Nx, Ny, Nz = " << numprocs << " " << Nx << " " << Ny << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Ny, Nz);


}

void test_FFF_slab_2d(char** argv)
{
	Nx=atoi(argv[2]);
	Nz=atoi(argv[3]);


	if (my_id==0) cout << "numprocs, Nx, Nz = " << numprocs << " " << Nx << " " << Nz << endl;

	spectralTransform.Init("FF", Nx, Nz);

	return;

	int local_Nx = spectralTransform.local_Nx;
	int local_Nz = spectralTransform.local_Nz;

	int local_Nx_start = spectralTransform.local_Nx_start;
	int local_Nz_start = spectralTransform.local_Nz_start;


	Array<complx, 2> A(Nx, Nz/2+1);

	Array<DP, 2> Ar_init(local_Nx, Nz+2);
	Array<DP, 2> Ar(local_Nx, Nz+2);

	/*Ar_init=0;

	for (int rx=0; rx<local_Nx; rx++)
			for (int rz=0; rz<Nz; rz++){
				Ar_init(rx, rz) = f(local_Nx_start + rx,rz);
			}



*/
	DP local_energy = sum(sqr(Ar_init));
	DP total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;

    DP start_time;
    DP end_time;

	Ar = Ar_init;

    DP local_time=0, local_time_max;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

	local_time -= MPI_Wtime();
	for (int iter=0; iter<N_iter; iter++) {

		spectralTransform.Forward_transform(Ar, A);

		local_energy = sum(sqr(abs(A)));
		MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

		if (my_id==0)	cout << "After F transform energy = " << total_energy << endl;
		
		//FFF	
		for (int lx=0; lx<Nx; lx++)
				for (int lz=0; lz<local_Nz; lz++)
					if (abs(A(lx,lz)) > 1E-3)
						cout << "Large modes: " << my_id << " " << local_Nx_start + lx << " " << lz << " " << A(lx, lz) << endl;
		Ar=0;

		spectralTransform.Inverse_transform(A, Ar);
	}
	local_time += MPI_Wtime();
	MPI_Reduce(&local_time, &local_time_max, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);
		
	local_energy = sum(sqr(Ar));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "After I transform energy = " << total_energy << endl;

	DP local_diff_energy = sum(sqr(Ar-Ar_init));
	DP total_diff_energy = sum(sqr(Ar-Ar_init));
	MPI_Reduce(&local_diff_energy, &total_diff_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "sum(sqr(Ar-Ar_init)) = " << total_diff_energy << endl;

	if (my_id == 0){
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << local_time_max << endl;
	}
}

void test_FFF_slab_3d(char** argv)
{
	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);

	int N_iter=atoi(argv[5]);


	if (my_id==0) cout << "numprocs, Nx, Ny, Nz = " << numprocs << " " << Nx << " " << Ny << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Ny, Nz);


	int local_Nx = spectralTransform.local_Nx;
	int local_Ny = spectralTransform.local_Ny;
	int local_Nz = spectralTransform.local_Nz;

	int local_Nx_start = spectralTransform.local_Nx_start;
	int local_Ny_start = spectralTransform.local_Ny_start;
	int local_Nz_start = spectralTransform.local_Nz_start;


	Array<complx, 3> A(Nx, local_Ny, Nz/2+1);

	Array<DP, 3> Ar_init(local_Nx, Ny, Nz+2);
	Array<DP, 3> Ar(local_Nx, Ny, Nz+2);

	Ar_init=0;

	for (int rx=0; rx<local_Nx; rx++)
		for (int ry=0; ry<Ny; ry++)
			for (int rz=0; rz<Nz; rz++){
				Ar_init(rx, ry, rz) = f(local_Nx_start + rx, ry,rz);
			}



	DP local_energy = sum(sqr(Ar_init));
	DP total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;

    DP start_time;
    DP end_time;

	Ar = Ar_init;

    DP local_time=0, local_time_max;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

	local_time -= MPI_Wtime();
	for (int iter=0; iter<N_iter; iter++) {

		spectralTransform.Forward_transform(Ar, A);

/*		local_energy = sum(sqr(abs(A)));
		MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

		if (my_id==0)	cout << "After F transform energy = " << total_energy << endl;
		
		//FFF	
		for (int lx=0; lx<Nx; lx++)
			for (int ly=0; ly<local_Ny; ly++)
				for (int lz=0; lz<Nz/2+1; lz++)
					if (abs(A(lx,ly,lz)) > 1E-3)
						cout << "Large modes: " << my_id << " " << lx << " " << local_Ny_start + ly << " " << lz << " " << A(lx, ly, lz) << endl;
		Ar=0;
*/
		spectralTransform.Inverse_transform(A, Ar);
	}
	local_time += MPI_Wtime();
	MPI_Reduce(&local_time, &local_time_max, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);
		
	local_energy = sum(sqr(Ar));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "After I transform energy = " << total_energy << endl;

	DP local_diff_energy = sum(sqr(Ar-Ar_init));
	DP total_diff_energy = sum(sqr(Ar-Ar_init));
	MPI_Reduce(&local_diff_energy, &total_diff_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "sum(sqr(Ar-Ar_init)) = " << total_diff_energy << endl;

	if (my_id == 0){
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << local_time_max << endl;
	}
}

void test_FFF_pencil_3d(char** argv)
{
	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);
	int num_p_cols=atoi(argv[5]);

	int N_iter=atoi(argv[6]);



	// spectralTransform.select_plan_id = atoi(argv[7]);

	if (my_id==0) cout << "FFF pencil: numprocs, Nx, Ny, Nz, iter = " << numprocs << " " << Nx << " " << Ny << " " << Nz << " " << N_iter << " " << spectralTransform.select_plan_id << endl;

	spectralTransform.Init(basis_type, Nx, Ny, Nz, num_p_cols);

	int num_p_rows = spectralTransform.num_p_rows;
	if (my_id==0) cout << "FFF pencil: Initialization done with proc_grid = " << num_p_rows << "x" << num_p_cols << endl;

	int local_Nx_col = spectralTransform.local_Nx_col;
	int local_Ny_col = spectralTransform.local_Ny_col;
	int local_Nz_col = spectralTransform.local_Nz_col;

	int local_Nx_row = spectralTransform.local_Nx_row;
	int local_Ny_row = spectralTransform.local_Ny_row;
	int local_Nz_row = spectralTransform.local_Nz_row;

	int local_Nx_start_col = spectralTransform.local_Nx_start_col;
	int local_Ny_start_col = spectralTransform.local_Ny_start_col;
	int local_Nz_start_col = spectralTransform.local_Nz_start_col;

	int local_Nx_start_row = spectralTransform.local_Nx_start_row;
	int local_Ny_start_row = spectralTransform.local_Ny_start_row;
	int local_Nz_start_row = spectralTransform.local_Nz_start_row;


	Array<complx, 3> A(Nx, local_Ny_row, local_Nz_col);
	Array<DP, 3> Ar_init(local_Nx_row, local_Ny_col, Nz+2);
	Array<DP, 3> Ar(local_Nx_row, local_Ny_col, Nz+2);

	Ar_init=0;


	//return;

	for (int rx=0; rx<local_Nx_row; rx++)
		for (int ry=0; ry<local_Ny_col; ry++)
			for (int rz=0; rz<2*local_Nz_col*num_p_rows; rz++){
				Ar_init(rx, ry, rz) = f(local_Nx_start_row + rx, local_Ny_start_col + ry, rz);//my_id*Nx*Ny*(Nz+2) + (rx*Ar_init.extent(1)*Ar_init.extent(2)) + ry*Ar_init.extent(2) + rz;
			}

	//Zero pad
	Ar_init(Range::all(),Range::all(),Range(Nz,Nz+1)) = 0;

/*
	for (int rx=0; rx<Ar_init.extent(0); rx++)
		for (int ry=0; ry<Ar_init.extent(1); ry++)
			for (int rz=0; rz<Ar_init.extent(1); rz++){
				// if (abs(Ar_init(rx,ry,rz)) > 1E-3)
					cout << "Init: " << my_id << " " << rx << " " << local_Ny_start + ry << " " << 2*local_Nz_start + rz << " " << Ar_init(rx, ry, rz) << " " << endl;
				}
*/

	/*for (int i=0; i<numprocs; i++) {
		if (i==my_id)
			cout << Ar_init(Range::all(), Range::all(), 0);
		MPI_Barrier(MPI_COMM_WORLD);
	}
*/
	double local_energy = sum(sqr(Ar_init));
	double total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;

    DP start_time;
    DP end_time;

	Ar = Ar_init;

    DP local_time=0, local_time_max;
    MPI_Barrier(MPI_COMM_WORLD);

	spectralTransform.Forward_transform(Ar, A);
	spectralTransform.Inverse_transform(A, Ar);


	double *timers=spectralTransform.plan->timers;
	int num_timers=spectralTransform.plan->num_timers;

	for (int i=0; i<22; i++)
		timers[i]=0;


	for (int iter=0; iter<N_iter; iter++) {

		// A=0;
		local_time -= MPI_Wtime();
		spectralTransform.Forward_transform(Ar, A);
		local_time += MPI_Wtime();
		
		/*for (int lx=0; lx<A.extent(0); lx++)
			for (int ly=0; ly<A.extent(1); ly++)
				for (int lz=0; lz<A.extent(2); lz++)
					if (abs(A(lx,ly,lz))>1e-3)
						cout << "Large modes: " << my_id << " " << lx << " " << local_Ny_start_col + ly << " " << local_Nz_start_row + lz << " " << A(lx, ly, lz) << endl;

		local_energy=sum(sqr(abs(A)));
		MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);*/
		// if (my_id==0)	cout << "After F transform energy = " << total_energy << endl;
		// Ar=0;
		local_time -= MPI_Wtime();
		spectralTransform.Inverse_transform(A, Ar);
		local_time += MPI_Wtime();

	}
	MPI_Reduce(&local_time, &local_time_max, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);
		
	local_energy = sum(sqr(Ar));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "After I transform energy = " << total_energy << endl;

	DP local_diff_energy = sum(sqr(Ar-Ar_init));
	DP total_diff_energy = sum(sqr(Ar-Ar_init));
	MPI_Reduce(&local_diff_energy, &total_diff_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id == 0)
		cout << "sum(sqr(Ar-Ar_init)) = " << total_diff_energy << endl;

	double *timer_sum, *timer_max, *timer_min;
	timer_sum = new double[num_timers];
	timer_max = new double[num_timers];
	timer_min = new double[num_timers];


 	for(int i=0;i < num_timers;i++) {
 		timers[i]/=N_iter;
 	}
    

	MPI_Reduce(timers,timer_sum,num_timers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(timers,timer_max,num_timers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Reduce(timers,timer_min,num_timers,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

 	for(int i=0;i < num_timers;i++)
 		timer_sum[i]/=numprocs;

	if (my_id == 0){
		cout << fixed << setprecision(3) << "Time taken per loop = "<< local_time_max/N_iter << endl;
     	for(int i=0;i < num_timers;i++)
			cout << "timer["<<i<<"] (avg/max/min): " << timer_sum[i] << " " << timer_max[i] << " " << timer_min[i] << endl;

	}
}

void test_SFF_pencil_3d(char** argv)
{
	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);
	int num_p_cols=atoi(argv[5]);

	int N_iter=atoi(argv[6]);

	sincostr_option = "SFF";

	if (my_id==0) cout << "SFF pencil: numprocs, Nx, Ny, Nz = " << numprocs << " " << Nx << " " << Ny << " " << Nz << endl;


	spectralTransform.Init(basis_type, Nx, Ny, Nz, num_p_cols);

	int num_p_rows = spectralTransform.num_p_rows;

	int my_id_row = spectralTransform.my_id_row;
	int my_id_col = spectralTransform.my_id_col;

	int local_Nx_col = spectralTransform.local_Nx_col;
	int local_Ny_col = spectralTransform.local_Ny_col;
	int local_Nz_col = spectralTransform.local_Nz_col;

	int local_Nx_row = spectralTransform.local_Nx_row;
	int local_Ny_row = spectralTransform.local_Ny_row;
	int local_Nz_row = spectralTransform.local_Nz_row;

	int local_Nx_start_col = spectralTransform.local_Nx_start_col;
	int local_Ny_start_col = spectralTransform.local_Ny_start_col;
	int local_Nz_start_col = spectralTransform.local_Nz_start_col;

	int local_Nx_start_row = spectralTransform.local_Nx_start_row;
	int local_Ny_start_row = spectralTransform.local_Ny_start_row;
	int local_Nz_start_row = spectralTransform.local_Nz_start_row;


	Array<complx, 3> A(local_Nx_row, Ny, local_Nz_col);
	Array<DP, 3> Ar_init(Nx, local_Ny_col, 2*local_Nz_row);
	Array<DP, 3> Ar(Nx, local_Ny_col, 2*local_Nz_row);

	Ar_init=0;


	for (int rx=0; rx<Nx; rx++)
		for (int ry=0; ry<local_Ny_col; ry++)
			for (int rz=0; rz<2*local_Nz_row; rz++){
				Ar_init(rx, ry, rz) = f(rx, local_Ny_start_col + ry, 2*local_Nz_start_row + rz); //my_id*Nx*Ny*(Nz+2) + (rx*Ar_init.extent(1)*Ar_init.extent(2)) + ry*Ar_init.extent(2) + rz;
			}

	//ZeroPad only when Nz is divisible
	if ( (my_id_row == num_p_cols-1) && (2*local_Nz_row*num_p_cols == Nz+2) ) {
		Ar_init(Range::all(),Range::all(),Range(2*local_Nz_row-2,2*local_Nz_row-1)) = 0;
	}

	//Zero pad
	//Ar_init(Range::all(),Range::all(),Range(Nz,Nz+1)) = 0;

/*
	for (int rx=0; rx<Ar_init.extent(0); rx++)
		for (int ry=0; ry<Ar_init.extent(1); ry++)
			for (int rz=0; rz<Ar_init.extent(1); rz++){
				// if (abs(Ar_init(rx,ry,rz)) > 1E-3)
					cout << "Init: " << my_id << " " << rx << " " << local_Ny_start + ry << " " << 2*local_Nz_start + rz << " " << Ar_init(rx, ry, rz) << " " << endl;
				}
*/

	/*for (int i=0; i<numprocs; i++) {
		if (i==my_id)
			cout << Ar_init(Range::all(), Range::all(), 0);
		MPI_Barrier(MPI_COMM_WORLD);
	}
*/
	double local_energy = sum(sqr(Ar_init));
	double total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;

    DP start_time;
    DP end_time;

	Ar = Ar_init;

    DP local_time=0, local_time_max;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

	local_time -= MPI_Wtime();
	for (int iter=0; iter<N_iter; iter++) {

		A=0;
		spectralTransform.Forward_transform(sincostr_option, Ar, A);

		for (int lx=0; lx<A.extent(0); lx++)
			for (int ly=0; ly<A.extent(1); ly++)
				for (int lz=0; lz<A.extent(2); lz++)
					if (abs(A(lx,ly,lz))>1e-3)
						cout << "Large modes: " << my_id << " " << lx << " " << local_Ny_start_col + ly << " " << local_Nz_start_row + lz << " " << A(lx, ly, lz) << endl;

		local_energy=sum(sqr(abs(A)));
		MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
		if (my_id==0)	cout << "After F transform energy = " << total_energy << endl;
		Ar=0;

		spectralTransform.Inverse_transform(sincostr_option, A, Ar);

	}
	local_time += MPI_Wtime();
	MPI_Reduce(&local_time, &local_time_max, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);

	local_energy = sum(sqr(Ar));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "After I transform energy = " << total_energy << endl;

	DP local_diff_energy = sum(sqr(Ar-Ar_init));
	DP total_diff_energy = sum(sqr(Ar-Ar_init));
	MPI_Reduce(&local_diff_energy, &total_diff_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "sum(sqr(Ar-Ar_init)) = " << total_diff_energy << endl;

	if (my_id == 0){
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << local_time_max << endl;
	}
}

void test_SSF_pencil_3d(char** argv)
{
	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);
	int num_p_cols=atoi(argv[5]);

	int N_iter=atoi(argv[6]);

	sincostr_option = "SSF";

	if (my_id==0) cout << "SSF pencil: numprocs, Nx, Ny, Nz = " << numprocs << " " << Nx << " " << Ny << " " << Nz << endl;


	spectralTransform.Init(basis_type, Nx, Ny, Nz, num_p_cols);

	int num_p_rows = spectralTransform.num_p_rows;

	int my_id_row = spectralTransform.my_id_row;
	int my_id_col = spectralTransform.my_id_col;

	int local_Nx_col = spectralTransform.local_Nx_col;
	int local_Ny_col = spectralTransform.local_Ny_col;
	int local_Nz_col = spectralTransform.local_Nz_col;

	int local_Nx_row = spectralTransform.local_Nx_row;
	int local_Ny_row = spectralTransform.local_Ny_row;
	int local_Nz_row = spectralTransform.local_Nz_row;

	int local_Nx_start_col = spectralTransform.local_Nx_start_col;
	int local_Ny_start_col = spectralTransform.local_Ny_start_col;
	int local_Nz_start_col = spectralTransform.local_Nz_start_col;

	int local_Nx_start_row = spectralTransform.local_Nx_start_row;
	int local_Ny_start_row = spectralTransform.local_Ny_start_row;
	int local_Nz_start_row = spectralTransform.local_Nz_start_row;


	Array<complx, 3> A(local_Nx_row, local_Ny_col, Nz/2+1);
	Array<DP, 3> Ar_init(Nx, local_Ny_row, 2*local_Nz_col);
	Array<DP, 3> Ar(Nx, local_Ny_row, 2*local_Nz_col);

	Ar_init=0;


	for (int rx=0; rx<Nx; rx++)
		for (int ry=0; ry<local_Ny_row; ry++)
			for (int rz=0; rz<2*local_Nz_col; rz++){
				Ar_init(rx, ry, rz) = f(rx, local_Ny_start_row + ry, 2*local_Nz_start_col + rz);//my_id*Nx*Ny*(Nz+2) + (rx*Ar_init.extent(1)*Ar_init.extent(2)) + ry*Ar_init.extent(2) + rz;
			}

	//ZeroPad only when Nz+2 is divisible
	if ( (my_id_col == num_p_rows-1) && (2*local_Nz_col*num_p_rows == Nz+2) ) {
		Ar_init(Range::all(),Range::all(),Range(2*local_Nz_col-2,2*local_Nz_col-1)) = 0;
	}

/*
	for (int rx=0; rx<Ar_init.extent(0); rx++)
		for (int ry=0; ry<Ar_init.extent(1); ry++)
			for (int rz=0; rz<Ar_init.extent(1); rz++){
				// if (abs(Ar_init(rx,ry,rz)) > 1E-3)
					cout << "Init: " << my_id << " " << rx << " " << local_Ny_start + ry << " " << 2*local_Nz_start + rz << " " << Ar_init(rx, ry, rz) << " " << endl;
				}
*/

	/*for (int i=0; i<numprocs; i++) {
		if (i==my_id)
			cout << Ar_init(Range::all(), Range::all(), 0);
		MPI_Barrier(MPI_COMM_WORLD);
	}
*/
	double local_energy = sum(sqr(Ar_init));
	double total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;

    DP start_time;
    DP end_time;

	Ar = Ar_init;

    DP local_time=0, local_time_max;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

	local_time -= MPI_Wtime();
	for (int iter=0; iter<N_iter; iter++) {

		A=0;
		spectralTransform.Forward_transform(sincostr_option, Ar, A);

		for (int lx=0; lx<A.extent(0); lx++)
			for (int ly=0; ly<A.extent(1); ly++)
				for (int lz=0; lz<A.extent(2); lz++)
					if (abs(A(lx,ly,lz))>1e-3)
						cout << "Large modes: " << my_id << " " << lx << " " << local_Ny_start_col + ly << " " << local_Nz_start_row + lz << " " << A(lx, ly, lz) << endl;

		local_energy=sum(sqr(abs(A)));
		MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
		if (my_id==0)	cout << "After F transform energy = " << total_energy << endl;
		Ar=0;

		spectralTransform.Inverse_transform(sincostr_option, A, Ar);

	}
	local_time += MPI_Wtime();
	MPI_Reduce(&local_time, &local_time_max, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);

	local_energy = sum(sqr(Ar));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "After I transform energy = " << total_energy << endl;

	DP local_diff_energy = sum(sqr(Ar-Ar_init));
	DP total_diff_energy = sum(sqr(Ar-Ar_init));
	MPI_Reduce(&local_diff_energy, &total_diff_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "sum(sqr(Ar-Ar_init)) = " << total_diff_energy << endl;

	if (my_id == 0){
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << local_time_max << endl;
	}
}



void test_SFF_2d(char** argv)
{
	sincostr_option = "CFF";//argv[2];

	Nx=atoi(argv[2]);
	Nz=atoi(argv[3]);

	//N_iter=atoi(argv[5]);


	if (my_id==0) cout << "Nx, Nz = " << Nx << " " << Nz << endl;

	spectralTransform.Init("SF", Nx, Nz);

	return;
}


void test_SFF_3d(char** argv)
{
	sincostr_option = "CFF";//argv[2];

	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);

	//N_iter=atoi(argv[5]);


	if (my_id==0) cout << "Nx, Ny, Nz = " << Nx << " " << Ny << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Ny, Nz);

	return;
	int local_Nx = spectralTransform.local_Nx;
	int local_Ny = spectralTransform.local_Ny;
	int local_Nz = spectralTransform.local_Nz;

	int local_Nx_start = spectralTransform.local_Nx_start;
	int local_Ny_start = spectralTransform.local_Ny_start;
	int local_Nz_start = spectralTransform.local_Nz_start;


	Array<complx, 3> A(local_Nx, Ny, Nz/2+1);

	Array<DP, 3> Ar_init(Nx, local_Ny, Nz+2);
	Array<DP, 3> Ar(Nx, local_Ny, Nz+2);

	Ar_init=0;

	for (int rx=0; rx<Nx; rx++)
		for (int ry=0; ry<local_Ny; ry++)
			for (int rz=0; rz<Nz; rz++){
				Ar_init(rx, ry, rz) = f(rx, local_Ny_start + ry,rz);
			}



	DP local_energy = sum(sqr(Ar_init));
	DP total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;

    DP start_time;
    DP end_time;

	Ar = Ar_init;

    DP local_time_max=0;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

	local_time_max -= MPI_Wtime();
	for (int iter=0; iter<N_iter; iter++) {

		spectralTransform.Forward_transform(sincostr_option, Ar, A);

		local_energy = sum(sqr(abs(A)));
		MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

		if (my_id==0)	cout << "After F transform energy = " << total_energy << endl;

		for (int lx=0; lx<local_Nx; lx++)
			for (int ly=0; ly<Ny; ly++)
				for (int lz=0; lz<Nz/2+1; lz++)
					if (abs(A(lx,ly,lz)) > 1E-3)
						cout << "Large modes: " << my_id << " " << local_Nx_start + lx << " " << ly << " " << lz << " " << A(lx, ly, lz) << endl;
		
		Ar=0;

		spectralTransform.Inverse_transform(sincostr_option, A, Ar);

	}
	local_time_max += MPI_Wtime();
		
	local_energy = sum(sqr(Ar));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "After I transform energy = " << total_energy << endl;

	DP local_diff_energy = sum(sqr(Ar-Ar_init));
	DP total_diff_energy = sum(sqr(Ar-Ar_init));
	MPI_Reduce(&local_diff_energy, &total_diff_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "sum(sqr(Ar-Ar_init)) = " << total_diff_energy << endl;

	if (my_id == 0){
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << local_time_max << endl;
	}
}


void test_SSF_3d(char** argv)
{
	sincostr_option = "CCF";//argv[2];

	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);

	//N_iter=atoi(argv[5]);


	if (my_id==0) cout << "Nx, Ny, Nz = " << Nx << " " << Ny << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Ny, Nz);

	return;
}



void test_SS_2d(char** argv)
{
	sincostr_option = "CCC";// argv[2];

	Nx=atoi(argv[2]);
	Nz=atoi(argv[3]);

	if (my_id==0) cout << "Nx, Nz = " << Nx << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Nz);

	return;
}


void test_SSS_3d(char** argv)
{
	sincostr_option = "CCC";//argv[2];

	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);

	if (my_id==0) cout << "Nx, Ny, Nz = " << Nx << " " << Ny << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Ny, Nz);

	return;
}

void test_SSS_pencil_3d(char** argv)
{
	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);
	int num_p_cols=atoi(argv[5]);

	int N_iter=atoi(argv[6]);

	sincostr_option = "SSS";

	if (my_id==0) cout << "SSS pencil: numprocs, Nx, Ny, Nz = " << numprocs << " " << Nx << " " << Ny << " " << Nz << endl;


	spectralTransform.Init(basis_type, Nx, Ny, Nz, num_p_cols);

	int num_p_rows = spectralTransform.num_p_rows;

	int my_id_row = spectralTransform.my_id_row;
	int my_id_col = spectralTransform.my_id_col;

	int local_Nx_col = spectralTransform.local_Nx_col;
	int local_Ny_col = spectralTransform.local_Ny_col;
	int local_Nz_col = spectralTransform.local_Nz_col;

	int local_Nx_row = spectralTransform.local_Nx_row;
	int local_Ny_row = spectralTransform.local_Ny_row;
	int local_Nz_row = spectralTransform.local_Nz_row;

	int local_Nx_start_col = spectralTransform.local_Nx_start_col;
	int local_Ny_start_col = spectralTransform.local_Ny_start_col;
	int local_Nz_start_col = spectralTransform.local_Nz_start_col;

	int local_Nx_start_row = spectralTransform.local_Nx_start_row;
	int local_Ny_start_row = spectralTransform.local_Ny_start_row;
	int local_Nz_start_row = spectralTransform.local_Nz_start_row;


	Array<complx, 3> A(local_Nx_row, local_Ny_col, Nz/2);
	Array<DP, 3> Ar_init(Nx, local_Ny_row, 2*local_Nz_col);
	Array<DP, 3> Ar(Nx, local_Ny_row, 2*local_Nz_col);

	Ar_init=0;


	for (int rx=0; rx<Nx; rx++)
		for (int ry=0; ry<local_Ny_row; ry++)
			for (int rz=0; rz<2*local_Nz_col; rz++){
				Ar_init(rx, ry, rz) = f(rx, local_Ny_start_row + ry, 2*local_Nz_start_col + rz);//my_id*Nx*Ny*(Nz+2) + (rx*Ar_init.extent(1)*Ar_init.extent(2)) + ry*Ar_init.extent(2) + rz;
			}

/*
	for (int rx=0; rx<Ar_init.extent(0); rx++)
		for (int ry=0; ry<Ar_init.extent(1); ry++)
			for (int rz=0; rz<Ar_init.extent(1); rz++){
				// if (abs(Ar_init(rx,ry,rz)) > 1E-3)
					cout << "Init: " << my_id << " " << rx << " " << local_Ny_start + ry << " " << 2*local_Nz_start + rz << " " << Ar_init(rx, ry, rz) << " " << endl;
				}
*/

	/*for (int i=0; i<numprocs; i++) {
		if (i==my_id)
			cout << Ar_init(Range::all(), Range::all(), 0);
		MPI_Barrier(MPI_COMM_WORLD);
	}
*/
	double local_energy = sum(sqr(Ar_init));
	double total_energy;
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;

    DP start_time;
    DP end_time;

	Ar = Ar_init;

    DP local_time=0, local_time_max;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

	local_time -= MPI_Wtime();
	for (int iter=0; iter<N_iter; iter++) {

		A=0;
		spectralTransform.Forward_transform(sincostr_option, Ar, A);

		for (int lx=0; lx<A.extent(0); lx++)
			for (int ly=0; ly<A.extent(1); ly++)
				for (int lz=0; lz<A.extent(2); lz++)
					if (abs(A(lx,ly,lz))>1e-3)
						cout << "Large modes: " << my_id << " " << lx << " " << local_Ny_start_col + ly << " " << local_Nz_start_row + lz << " " << A(lx, ly, lz) << endl;

		local_energy=sum(sqr(abs(A)));
		MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
		if (my_id==0)	cout << "After F transform energy = " << total_energy << endl;
		Ar=0;

		spectralTransform.Inverse_transform(sincostr_option, A, Ar);

	}
	local_time += MPI_Wtime();
	MPI_Reduce(&local_time, &local_time_max, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);

	local_energy = sum(sqr(Ar));
	MPI_Reduce(&local_energy, &total_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "After I transform energy = " << total_energy << endl;

	DP local_diff_energy = sum(sqr(Ar-Ar_init));
	DP total_diff_energy = sum(sqr(Ar-Ar_init));
	MPI_Reduce(&local_diff_energy, &total_diff_energy, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_id == 0)
		cout << "sum(sqr(Ar-Ar_init)) = " << total_diff_energy << endl;

	if (my_id == 0){
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << local_time_max << endl;
	}
}


DP f(int rx, int ry, int rz)
{
	DP k0 = 1;
	DP x,y,z;
	
	DP Lx=2*M_PI;
	DP Ly=2*M_PI;
	DP Lz=2*M_PI;
	if (basis_type == "FFF") {
		x = rx*Lx/Nx;
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;
		return 8*sin(k0*x)*sin(k0*y)*sin(k0*z);
	}

	else if (basis_type == "FFFW") {
		x = rx*Lx/Nx;
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;
		return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis_type == "SFF") {
		Lx=M_PI;
		x = (rx+0.5)*Lx/Nx;
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;

		if (sincostr_option=="SFF")
			return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
		else if (sincostr_option=="CFF")
			return 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis_type == "SSF") {
		x = (rx+0.5)*Lx/Nx;
		y = (ry+0.5)*Ly/Ny;
		z = rz*Lz/Nz;

		if (sincostr_option=="SSF")
			return 8*sin(k0*x)*sin(k0*y)*cos(k0*z);
		else if (sincostr_option=="CCF")
			return 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis_type == "SSS") {
		x = (rx+0.5)*Lx/Nx;
		y = (ry+0.5)*Ly/Ny;
		z = (rz+0.5)*Lz/Nz;

		if (sincostr_option=="SSS")
			return 8*sin(k0*x)*sin(k0*y)*sin(k0*z);
		else if (sincostr_option=="CCC")
			return 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis_type == "ChFF") {
		x = rx*M_PI/(Nx-1);
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;
	}

	
	
//	return 4*sin(k0*x)*cos(k0*z);
//	return 4*(2*cos(k0*x)*cos(k0*x)-1)*cos(k0*y)*cos(k0*z); // Cheb
//	return 4*(cos(k0*x))*cos(k0*y)*cos(k0*z);  // Cheb
//	return 2*(cos(k0*x))*cos(k0*z);  // Cheb (101)
	// return 8*sin(k0*x)*sin(k0*y)*sin(k0*z);
	// return 4*sin(k0*x)*sin(k0*z);
}
