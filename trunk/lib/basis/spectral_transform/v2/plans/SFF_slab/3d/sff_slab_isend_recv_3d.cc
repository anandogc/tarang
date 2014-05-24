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
 * @author  M. K. Verma
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		ArrayIFFT
 */ 


#include "sff_slab_isend_recv_3d.h"


	
//*********************************************************************************************
SFF_slab_Isend_Recv_3D::SFF_slab_Isend_Recv_3D(int my_id, int numprocs, int num_iter, int Nx, int Ny, int Nz): SpectralPlan_Slab_3D("SFF", my_id, numprocs, Nx, Ny, Nz)
{

	X_3d.resize(local_Nx, Ny, Nz/2+1);	   // temp array
	Xr_3d.resize(Nx, local_Ny, Nz+2);

	//Initialize plans
	int Nx_plan[]={Nx};
	int Ny_plan[]={Ny};
	int Nz_plan[]={Nz};

	fftw_r2r_kind kind[1];

	kind[0]=FFTW_RODFT10;
	plan_sintr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT10;
	plan_costr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_RODFT01;
	plan_isintr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);
	
	kind[0]=FFTW_REDFT01;
	plan_icostr_x = FFTW_PLAN_MANY_R2R(1, Nx_plan, (Nz+2)*local_Ny,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    reinterpret_cast<DP*>(X_3d.data()), NULL,
	    (Nz+2)*local_Ny, 1,
	    kind, FFTW_PLAN_FLAG);


	plan_ft_r2c_z = FFTW_PLAN_MANY_DFT_R2C(1, Nz_plan, local_Ny,
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, (Nz+2),
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, (Nz/2+1),
		FFTW_PLAN_FLAG);

	plan_ift_c2r_z = FFTW_PLAN_MANY_DFT_C2R(1, Nz_plan, local_Ny,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		1, (Nz/2+1),
		reinterpret_cast<DP*>(X_3d.data()), NULL,
		1, (Nz+2),
		FFTW_PLAN_FLAG);


	plan_ft_c2c_y = FFTW_PLAN_MANY_DFT(1, Ny_plan, Nz/2+1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		Nz/2+1, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		Nz/2+1, 1,
		FFTW_FORWARD, FFTW_PLAN_FLAG);

	plan_ift_c2c_y = FFTW_PLAN_MANY_DFT(1, Ny_plan, Nz/2+1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		Nz/2+1, 1,
		reinterpret_cast<FFTW_COMPLEX*>(X_3d.data()), NULL,
		Nz/2+1, 1,
		FFTW_BACKWARD, FFTW_PLAN_FLAG);

	Init_array();
	Evaluate_time_per_step("CFF", num_iter);

	X_3d.free();
	Xr_3d.free();
}
//*********************************************************************************************

DP SFF_slab_Isend_Recv_3D::f(string sincostr_option, int rx, int ry, int rz)
{
	DP k0 = 1;
	DP x,y,z;
	
	DP Lx=M_PI;
	DP Ly=2*M_PI;
	DP Lz=2*M_PI;

	x = (rx+0.5)*Lx/Nx;
	y = ry*Ly/Ny;
	z = rz*Lz/Nz;

	if (sincostr_option=="SFF")
		return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
	else if (sincostr_option=="CFF")
		return 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
}
 

void SFF_slab_Isend_Recv_3D::Init_array()
{
	Xr_3d=0;	
	
	for (int rx=0; rx<Nx; rx++)
		for (int ry=0; ry<local_Ny; ry++)
			for (int rz=0; rz<Nz; rz++)
				Xr_3d(rx, ry, rz) = f("CFF", rx,local_Ny_start + ry,rz);
}

//*********************************************************************************************

void SFF_slab_Isend_Recv_3D::Normalize(Array<complx,3> A)
{
	A /= (2.0 * DP(Nx) * DP(Ny) * DP(Nz));
}

//*************************************************************************************


void SFF_slab_Isend_Recv_3D::Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A)
{

	Zero_pad_last_plane(Ar);	// Zero_pad the plane Nz

	SinCostr_x(sincostr_option[0],Ar);

	FT_r2c_z_and_Isend_yz(Ar);

	Recv_yz_and_FT_c2c_y(A);
	
	Normalize(A);
}

//*********************************************************************************************

void SFF_slab_Isend_Recv_3D::Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar)
{
	IFT_c2c_y_and_Isend_yz(A);

	Recv_yz_and_IFT_c2r_z(Ar);

	ISinCostr_x(sincostr_option[0], Ar);

	Zero_pad_last_plane(Ar); 
}

//Transpose
void SFF_slab_Isend_Recv_3D::Transpose(Array<DP,3> Ar, Array<complx,3> A){
	Isend_yz_y(Ar);
	Recv_yz_x(A);
}

void SFF_slab_Isend_Recv_3D::Transpose(Array<complx,3> A, Array<DP,3> Ar){
	Isend_yz_x(A);
	Recv_yz_y(Ar);
}

void SFF_slab_Isend_Recv_3D::SinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
	if (sincostr_option == 'S') {
		FFTW_EXECUTE_R2R(plan_sintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'X');
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_costr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SFF_slab_Isend_Recv_3D::ISinCostr_x(char sincostr_option, Array<DP,3> Ar)
{
    if (sincostr_option == 'S') {
		ArrayShiftLeft(Ar, 'X');
		FFTW_EXECUTE_R2R(plan_isintr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_option == 'C')
		FFTW_EXECUTE_R2R(plan_icostr_x, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SFF_slab_Isend_Recv_3D::FT_r2c_z_and_Isend_yz(Array<DP,3> Ar)
{
	for (int lx=0; lx<local_Nx; lx++)
	{
		for (int p=0; p<numprocs; p++)
		{
			FFTW_EXECUTE_DFT_R2C(plan_ft_r2c_z, reinterpret_cast<DP*>(Ar(p*local_Nx + lx,Range::all(),Range::all()).data()), reinterpret_cast<FFTW_COMPLEX*>(Ar(p*local_Nx + lx,Range::all(),Range::all()).data()));
			
			MPI_Isend(Ar(p*local_Nx + lx, Range::all(), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, &request[p*local_Nx + lx]);
		}
	}
}


void SFF_slab_Isend_Recv_3D::Recv_yz_and_FT_c2c_y(Array<complx,3> A)
{
	for (int lx=0; lx<local_Nx; lx++)
	{
		for (int p=0; p<numprocs; p++)
		{
			MPI_Recv(A(lx, Range(p*local_Ny, (p+1)*local_Ny-1), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, status);

		}

		FFTW_EXECUTE_DFT(plan_ft_c2c_y, reinterpret_cast<FFTW_COMPLEX*>(A(lx,Range::all(),Range::all()).data()), reinterpret_cast<FFTW_COMPLEX*>(A(lx,Range::all(),Range::all()).data()));

			
	}

	MPI_Waitall(Nx, request, status);
}

void SFF_slab_Isend_Recv_3D::IFT_c2c_y_and_Isend_yz(Array<complx,3> A)
{
	for (int lx=0; lx<local_Nx; lx++)
	{
		FFTW_EXECUTE_DFT(plan_ift_c2c_y, reinterpret_cast<FFTW_COMPLEX*>(A(lx,Range::all(),Range::all()).data()), reinterpret_cast<FFTW_COMPLEX*>(A(lx,Range::all(),Range::all()).data()));

		for (int p=0; p<numprocs; p++)
		{
			MPI_Isend(A(lx, Range(p*local_Ny, (p+1)*local_Ny-1), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, &request[p*local_Nx + lx]);

		}			
	}
}


void SFF_slab_Isend_Recv_3D::Recv_yz_and_IFT_c2r_z(Array<DP,3> Ar)
{
	for (int lx=0; lx<local_Nx; lx++)
	{
		for (int p=0; p<numprocs; p++)
		{
			MPI_Recv(Ar(p*local_Nx + lx, Range::all(), Range::all()).data(), (Nz+2)*local_Ny, MPI_DP, p, lx, MPI_COMM_WORLD, status);

			FFTW_EXECUTE_DFT_C2R(plan_ift_c2r_z, reinterpret_cast<FFTW_COMPLEX*>(Ar(p*local_Nx + lx,Range::all(),Range::all()).data()), reinterpret_cast<DP*>(Ar(p*local_Nx + lx,Range::all(),Range::all()).data()));
		}
	}

	MPI_Waitall(Nx, request, status);
}


//********************************	End of four_tr.cc *****************************************


