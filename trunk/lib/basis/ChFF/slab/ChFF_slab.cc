/* Tarang-2
 *
 * Copyright_PENCIL( (C) 2008, 2009  Mahendra K. Verma
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

/*! \file ChFF_clab.cc 
 * 
 * SIZE of A(Nx, Ny, Nz/2+1) with Nx=2^n (FFTW convention);
 * FFTW plan uses (Nx+1, col_x, col_x).. Nx+1 points.
 *
 * @sa ChFF_slab.h
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "ChFF_slab.h"
#include "BasicIO.h"

ChFF_SLAB::ChFF_SLAB()
{
		// kfactor, xfactor,  L, Delta_x
	if (global.field.kfactor.size() == 4) {
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		global.field.L[1] = 2.0/global.field.kfactor[1]; //Chebyshev box length = 2
		for (int i=2; i<=3; i++) 
			global.field.L[i] = 2*M_PI/global.field.kfactor[i];
	}
	
	else if (global.field.L.size() == 4) {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		global.field.kfactor[1] = 2.0/global.field.L[1];
		for (int i=2; i<=3; i++) 
			global.field.kfactor[i] = 2*M_PI/global.field.L[i];
	}
	
	else  {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		
		if (global.program.kind == "RBC") {
			global.field.kfactor[1] =  2.0;
			global.field.kfactor[2] =  3.12;
			global.field.kfactor[3] =  3.12;
		}
		else {
			global.field.kfactor[1] = 1.0;
			global.field.kfactor[2] = 1.0;
			global.field.kfactor[3] = 1.0;
		}
		
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		global.field.L[1] = 2.0/global.field.kfactor[1];
		for (int i=2; i<=3; i++) 
			global.field.L[i] = 2*M_PI/global.field.kfactor[i];
	}
	
	
	for (int i=1; i<=3; i++) {
		global.field.Delta_x[i] =  global.field.L[i]/global.field.N[i];
		global.field.xfactor[i] = 1/global.field.kfactor[i];
	}
	
	if (master) {
        cout << "L[i] = [" << global.field.L[1] << ", "  << global.field.L[2] << ", "  << global.field.L[3] << "]" << endl;
        
        cout << "kfactor[i] = [" << global.field.kfactor[1] << ", "  << global.field.kfactor[2] << ", "  << global.field.kfactor[3] << "]" << endl;
        
        cout << "dt_fixed = " << global.time.dt_fixed << endl;
    }
    
    global.field.twobyL1 = 2.0/global.field.L[1];
    // Normalize such that L[1]=2.
    for (int i=1; i<=3; i++)
        global.field.L[i] *= global.field.twobyL1;
    
    DP kfactor_1 = global.field.kfactor[1];
    for (int i=1; i<=3; i++)
        global.field.kfactor[i] *= 1.0/kfactor_1;
    
 //   global.time.dt_fixed *= pow2(global.field.twobyL1);
 //   global.time.dt = global.time.dt_fixed;
    
    if (master) {
        cout << "Modified L[i] = [" << global.field.L[1] << ", "  << global.field.L[2] << ", "  << global.field.L[3] << "]" << endl;
		
        cout << "Modified kfactor[i] = [" << global.field.kfactor[1] << ", "  << global.field.kfactor[2] << ", "  << global.field.kfactor[3] << "]" << endl;
        
   //     cout << "Modified dt_fixed = " << global.time.dt_fixed << endl;
    }
	
	
	// Proper allocation or not?
	
	bool Nxproperdiv, Nyproperdiv, Nzproperdiv;
	
	if (global.field.Ny > 1) {
		Nxproperdiv = (global.field.Nx%global.mpi.numprocs == 0);
		Nyproperdiv = (global.field.Ny%global.mpi.numprocs == 0);
		Nzproperdiv = true;
	}
	
	else if (global.field.Ny == 1) {
		Nxproperdiv = (global.field.Nx%global.mpi.numprocs == 0);
		Nzproperdiv = ((global.field.Nz/2+1)%global.mpi.numprocs == 0);
		Nyproperdiv = true;
	}
	
	if (!(Nxproperdiv && Nyproperdiv && Nzproperdiv)) {
		cout << "N div " << global.field.Nx << " "<< global.field.Ny << " " << global.field.Nz << " " << Nxproperdiv << " " << Nyproperdiv << "  " << Nzproperdiv << endl;
		if (global.mpi.master) cerr << "ERROR in ChFF_slab.cc: Array not being divided equally.  Check dimensions" << endl;
		exit(1);
	}
	
	spectralTransform.Init("ChFF", "SLAB", Nx, Ny, Nz);
   
	global.field.maxlx = global.field.Nx-1;
	global.field.maxly = spectralTransform.local_Ny-1;
	global.field.maxlz = global.field.Nz/2;
	
	if (Ny > 1) {
		global.field.shape_complex_array = spectralTransform.local_Ny, global.field.Nz/2+1, global.field.Nx;
        global.field.shape_real_array = spectralTransform.local_Ny, global.field.Nz+2, global.field.Nx;
		
		BasicIO::shape_full_complex_array = Ny,Nz/2+1,Nx;
		BasicIO::shape_full_real_array = Ny,Nz+2,Nx;
		
		BasicIO::direction_z_complex_array = 0,1,0;
		BasicIO::direction_z_real_array = 0,1,0;
		
		BasicIO::id_complex_array = my_id,0,0;
		BasicIO::id_real_array = my_id,0,0;
		
		if (global.io.N_in_reduced.size()==3)
			BasicIO::shape_in_reduced_array = global.io.N_in_reduced[1],global.io.N_in_reduced[2]/2+1,global.io.N_in_reduced[0];
		
		if (global.io.N_out_reduced.size()==3)
			BasicIO::shape_out_reduced_array = global.io.N_out_reduced[1],global.io.N_out_reduced[2]/2+1,global.io.N_out_reduced[0];
		
		BasicIO::Fourier_directions = 1,1,0;
	}
	
	else if (Ny == 1) {
		global.field.shape_complex_array = 1,spectralTransform.local_Nz,global.field.Nx;
        global.field.shape_real_array = 1,2*spectralTransform.local_Nz,global.field.Nx;
		
		BasicIO::shape_full_complex_array = 1,Nz/2+1,Nx;
		BasicIO::shape_full_real_array = 1,Nz+2,Nx;
		
		BasicIO::direction_z_complex_array = 0,1,0;
		BasicIO::direction_z_real_array = 0,1,0;
		
		BasicIO::id_complex_array = 0,my_id,0;
		BasicIO::id_real_array = 0,my_id,0;
		
		if (global.io.N_in_reduced.size()==3)
			BasicIO::shape_in_reduced_array = 1,global.io.N_in_reduced[2]/2+1,global.io.N_in_reduced[0];
		
		if (global.io.N_out_reduced.size()==3)
			BasicIO::shape_out_reduced_array = 1,global.io.N_out_reduced[2]/2+1,global.io.N_out_reduced[0];
		
		BasicIO::Fourier_directions = 0,1,0;
	}
    
	// alias
	kfactor.resize(4);
	kfactor[0]= global.field.kfactor[0];
	kfactor[1]= global.field.kfactor[1];
	kfactor[2]= global.field.kfactor[2];
	kfactor[3]= global.field.kfactor[3];

	local_Nx=spectralTransform.local_Nx;
	local_Ny=spectralTransform.local_Ny;
	local_Nz=spectralTransform.local_Nz;
	
	local_Nx_start=spectralTransform.local_Nx_start;
	local_Ny_start=spectralTransform.local_Ny_start;
	local_Nz_start=spectralTransform.local_Nz_start;

	shape_complex_array = global.field.shape_complex_array;
    shape_real_array = global.field.shape_real_array;
	
	global.temp_array.X.resize(shape_complex_array);
    global.temp_array.X2.resize(shape_complex_array);
	
	global.temp_array.Xr.resize(shape_real_array);
    global.temp_array.Xr2.resize(shape_real_array);
	
    
	
    global.temp_array.V1_x.resize(Nx);
    
	global.temp_array.in_helm_real.resize(Nx);
	global.temp_array.out_helm_real.resize(Nx);
    global.temp_array.in_helm_complex.resize(Nx);
	global.temp_array.out_helm_complex.resize(Nx);
	
	global.temp_array.pressure_plus.resize(local_Ny, local_Nz, Nx);
	global.temp_array.pressure_minus.resize(local_Ny, local_Nz, Nx);
	global.temp_array.vx_plus.resize(local_Ny, local_Nz, Nx);
	global.temp_array.vx_minus.resize(local_Ny, local_Nz, Nx);
   
	global.temp_array.Xreal.resize(local_Ny, local_Nz, Nx);
	global.temp_array.influence_matrix.resize(local_Ny,local_Nz,2,2);

	
}