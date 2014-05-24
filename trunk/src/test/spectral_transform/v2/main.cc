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

#include <iomanip>
#include "spectral_transform.h"

SpectralTransform spectralTransform;

int my_id,numprocs;
int N_iter;
int Nx,Ny,Nz;

string basis_type, sincostr_option;

void test_FFFW_2d(char** argv);
void test_FFFW_3d(char** argv);

void test_FFF_2d(char** argv);
void test_FFF_3d(char** argv);

void test_SFF_2d(char** argv);
void test_SFF_3d(char** argv);

void test_SSF_3d(char** argv);

void test_SS_2d(char** argv);
void test_SSS_3d(char** argv);


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
		test_FFF_2d(argv);
	else if (basis_type=="SF")
		test_SFF_2d(argv);
	else if (basis_type=="SS")
		test_SS_2d(argv);

	else if (basis_type=="FFFW")
		test_FFFW_3d(argv);
	else if (basis_type=="FFF")
		test_FFF_3d(argv);
	else if (basis_type=="SFF")
		test_SFF_3d(argv);
	else if (basis_type=="SSF")
		test_SSF_3d(argv);
	else if (basis_type=="SSS")
		test_SSS_3d(argv);

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

void test_FFF_2d(char** argv)
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

    DP local_time=0, total_time;
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
	MPI_Reduce(&local_time, &total_time, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);
		
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
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << total_time << endl;
	}
}

void test_FFF_3d(char** argv)
{
	Nx=atoi(argv[2]);
	Ny=atoi(argv[3]);
	Nz=atoi(argv[4]);

	//int N_iter=atoi(argv[3]);


	if (my_id==0) cout << "numprocs, Nx, Ny, Nz = " << numprocs << " " << Nx << " " << Ny << " " << Nz << endl;

	spectralTransform.Init(basis_type, Nx, Ny, Nz);


	return;

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

    DP local_time=0, total_time;
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
			for (int ly=0; ly<local_Ny; ly++)
				for (int lz=0; lz<Nz/2+1; lz++)
					if (abs(A(lx,ly,lz)) > 1E-3)
						cout << "Large modes: " << my_id << " " << lx << " " << local_Ny_start + ly << " " << lz << " " << A(lx, ly, lz) << endl;
		Ar=0;

		spectralTransform.Inverse_transform(A, Ar);
	}
	local_time += MPI_Wtime();
	MPI_Reduce(&local_time, &total_time, 1, MPI_DP, MPI_MAX, 0, MPI_COMM_WORLD);
		
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
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << total_time << endl;
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

    DP total_time=0;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

	total_time -= MPI_Wtime();
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
	total_time += MPI_Wtime();
		
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
		cout << fixed << setprecision(3) << "Time taken for "<< N_iter <<" pair of transform = " << total_time << endl;
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
		return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
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
	}
	
	else if (basis_type == "SSS") {
		x = (rx+0.5)*Lx/Nx;
		y = (ry+0.5)*Ly/Ny;
		z = (rz+0.5)*Lz/Nz;
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
