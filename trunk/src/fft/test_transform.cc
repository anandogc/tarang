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

/*! \file Ifluid_main.cc 
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 


#include "main.h"

//*********************************************************************************************					

DP f(int rx, int ry, int rz)
{
	DP k0 = 1;
	DP x,y,z;
	
	if (global.program.basis_type == "FFF") {
		x = rx*global.field.L[1]/Nx;
		y = ry*global.field.L[2]/Ny;
		z = rz*global.field.L[3]/Nz;
		return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
	}

	else if (global.program.basis_type == "FFFW") {
		x = rx*global.field.L[1]/Nx;
		y = ry*global.field.L[2]/Ny;
		z = rz*global.field.L[3]/Nz;
		return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (global.program.basis_type == "SFF") {
		x = (rx+0.5)*global.field.L[1]/Nx;
		y = ry*global.field.L[2]/Ny;
		z = rz*global.field.L[3]/Nz;
	}
	
	else if (global.program.basis_type == "SSF") {
		x = (rx+0.5)*global.field.L[1]/Nx;
		y = (ry+0.5)*global.field.L[2]/Ny;
		z = rz*global.field.L[3]/Nz;
	}
	
	else if (global.program.basis_type == "SSS") {
		x = (rx+0.5)*global.field.L[1]/Nx;
		y = (ry+0.5)*global.field.L[2]/Ny;
		z = (rz+0.5)*global.field.L[3]/Nz;
	}
	
	else if (global.program.basis_type == "ChFF") {
		x = rx*M_PI/(Nx-1);
		y = ry*global.field.L[2]/Ny;
		z = rz*global.field.L[3]/Nz;
	}

	
//	return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
//	return 4*sin(k0*x)*cos(k0*z);
//	return 4*(2*cos(k0*x)*cos(k0*x)-1)*cos(k0*y)*cos(k0*z); // Cheb
//	return 4*(cos(k0*x))*cos(k0*y)*cos(k0*z);  // Cheb
//	return 2*(cos(k0*x))*cos(k0*z);  // Cheb (101)
	return 8*sin(k0*x)*sin(k0*y)*sin(k0*z);
	// return 4*sin(k0*x)*sin(k0*z);
}
 
//*********************************************************************************************

int test_transform()
{
    DP value, value_find;
	
	// Spectral->real->spectral
	if (master)
		cout << "Testing Spectral->real->spectral : " << endl << endl;
	
	A_old = A;
	if (master)
		cout << "Initial A: " << endl;
	universal->Print_large_Fourier_elements(A);
    
	universal->Inverse_transform(A, Ar);
	
	A=0;

	ArrayOps::Print_array_all_procs(Ar);
    
	universal->Forward_transform(Ar, A);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (master)
		cout << "A After spectral->real->spectral : " << endl;
	universal->Print_large_Fourier_elements(A);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	A = A-A_old;
	DP total_A_minusA_old = universal->Get_total_energy(A);
	if (master)
		cout << "Energy of A-A_old = " << total_A_minusA_old << endl << endl << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	
	return 0;
	// real->spectral->real
	if (master)
		cout << "Testing real->spectral->real : " << endl << endl;
 	// Initialize the array with a function f(x,y,z)
	DP val;
	for (int ry=0; ry<Ny; ry++)
		for (int rz=0; rz<Nz; rz++)
			for (int rx=0; rx<Nx; rx++){
				value = f(rx,ry,rz);
				universal->Assign_real_field(rx, ry, rz, Ar, value);
				//val = universal->Get_real_field(rx,ry,rz,Ar);
                //cout << "my_id, rx, ry, rz, value, value = " << my_id << " " << rx << " " << ry << " " << rz << " " << value  << " " << val << endl;
			}
	
	
	ArrayOps::Print_array_all_procs(Ar);


	Ar_old = Ar;
	
	universal->Forward_transform(Ar, A);
	
	MPI_Barrier(MPI_COMM_WORLD);

    if (master)
       cout << "After forward transform; A = " << endl;
    universal->Print_large_Fourier_elements(A);
	
    Ar = 0; 
	universal->Inverse_transform(A, Ar);
	MPI_Barrier(MPI_COMM_WORLD);

	 if (master)
       cout << "After inverse transform; A = " << endl;
	ArrayOps::Print_array_all_procs(Ar);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	Ar = Ar-Ar_old;
	DP total_Ar_minusAr_old = universal->Get_total_energy_real_space(Ar);
	if (master)
		cout << "Energy of Ar-Ar_old = " << total_Ar_minusAr_old << endl;


	return 0;
	
}


//********************************** End of Ifluid_main.cc ************************************



