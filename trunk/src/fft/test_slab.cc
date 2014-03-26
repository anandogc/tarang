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

void test_CFFF_slab(){
	global.program.basis_type = "CFFF";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "FFF";
	
	
	global.Process_basic_vars();
	universal=new CFFF_SLAB;
	global.Process_advanced_vars();
	
    //global.Print();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	//	result += test_1dtransform();
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,0,1, A, complx(0,-1));
	universal->Assign_spectral_field(-1,0,1, A, complx(0,1));
	
    test_transform();
	
	MPI_Barrier(MPI_COMM_WORLD);
}

void test_FFFW_slab(){
	global.program.basis_type = "FFFW";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "FFF";
	
	
	global.Process_basic_vars();
	universal=new FFFW_SLAB;
	global.Process_advanced_vars();
	
    //global.Print();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	//	result += test_1dtransform();
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,0,1, A, complx(0,-1));
	universal->Assign_spectral_field(-1,0,1, A, complx(0,1));
	
    test_transform();
	
	MPI_Barrier(MPI_COMM_WORLD);
}

void test_FFF_slab(){
	global.program.basis_type = "FFF";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "FFF";
	
	
	global.Process_basic_vars();
	universal=new FFF_SLAB;
	global.Process_advanced_vars();
	
    //global.Print();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	//	result += test_1dtransform();
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,0,1, A, complx(0,-1));
	universal->Assign_spectral_field(-1,0,1, A, complx(0,1));
	
    test_transform();
	
	MPI_Barrier(MPI_COMM_WORLD);
}

void test_SFF_slab(){
	global.program.basis_type = "SFF";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "SFF";
	
	global.Process_basic_vars();
	universal=new SFF_SLAB;
	global.Process_advanced_vars();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,0,1, A, complx(0,-1));
	//	universal->Assign_spectral_field(1,-1,1, A, complx(1,0));
	
    test_transform();
}

void test_SSF_slab(){
	global.program.basis_type = "SSF";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "CSF";
	
	global.Process_basic_vars();
	universal=new SSF_SLAB;
	global.Process_advanced_vars();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,1,1, A, complx(1,0));
	
    test_transform();
}

void test_SSS_slab(){
	global.program.basis_type = "SSS";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "S0S";
	
	global.Process_basic_vars();
	universal=new SSS_SLAB;
	global.Process_advanced_vars();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,0,1, A, 1.0);
	
    test_transform();
}


void test_ChFF_slab(){
	global.program.basis_type = "ChFF";
	global.program.decomposition = "SLAB";
	global.program.sincostr_switch = "C0F";

	global.Process_basic_vars();
	universal=new ChFF_SLAB;
	global.Process_advanced_vars();
	
	if (master) {
		cout <<  "FFT-transform for basis = " << global.program.basis_type << " and decomposition = " << global.program.decomposition <<  endl;
		cout << "L[i] = " << global.field.L[1] << "  " << global.field.L[2] << "  " << global.field.L[3] << endl;
	}
	
	Ar.resize(shape_real_array);
	Ar_old.resize(shape_real_array);
	
	A.resize(shape_complex_array);
	A_old.resize(shape_complex_array);
	
	// Assign nos to kx,ky,kz (first three args)
	universal->Assign_spectral_field(1,1,1, A, complx(1,0));
	universal->Assign_spectral_field(1,-1,1, A, complx(1,0));
	
    test_transform();
}
//********************************** End of Ifluid_main.cc ************************************



