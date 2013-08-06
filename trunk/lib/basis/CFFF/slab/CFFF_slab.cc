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

/*! \file FFF_slab.cc 
 * 
 * @sa FFF_slab.cc
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "CFFF_slab.h"
#include "BasicIO.h"

CFFF_SLAB::CFFF_SLAB()
{
	
	// kfactor, xfactor, L, Delta_x
	if (global.field.kfactor.size() == 4) {
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		for (int i=1; i<=3; i++) 
			global.field.L[i] = 2*M_PI/global.field.kfactor[i];
	}
	
	else if (global.field.L.size() == 4) {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		for (int i=1; i<=3; i++) 
			global.field.kfactor[i] = 2*M_PI/global.field.L[i];
	}
    
	
	for (int i=1; i<=3; i++) {
		global.field.Delta_x[i] =  global.field.L[i]/global.field.N[i];
		global.field.xfactor[i] = 1/global.field.kfactor[i];
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
		Nzproperdiv = (global.field.Nz%global.mpi.numprocs == 0);
		Nyproperdiv = true; 
	}
	
	if (!(Nxproperdiv && Nyproperdiv && Nzproperdiv)) {
		cout << "N div " << global.field.Nx << " "<< global.field.Ny << " " << global.field.Nz << " " << Nxproperdiv << " " << Nyproperdiv << "  " << Nzproperdiv << endl;
		if (global.mpi.master) cerr << "ERROR in FFF_slab.cc: Array not being divided equally.  Check dimensions" << endl;
		exit(1);
	}
    
    spectralTransform.Init("CFFF", "SLAB", Nx, Ny, Nz);
    
    global.fft.maxlx = global.fft.slab.local_Nx-1;
	global.fft.maxly = global.field.Ny-1;
	global.fft.maxlz = global.field.Nz-1;
	
	
    if (Ny > 1) {
        global.fft.complex_arraydim_1 = spectralTransform.local_Nx;
        global.fft.complex_arraydim_2 = global.field.Ny;
        global.fft.complex_arraydim_3 = global.field.Nz;
        
        global.fft.real_arraydim_1 = spectralTransform.local_Ny;
        global.fft.real_arraydim_2 = global.field.Nx;
        global.fft.real_arraydim_3 = global.field.Nz;
    }
    
    else if (Ny == 1) { 
        global.fft.real_arraydim_1 = spectralTransform.local_Nx;
        global.fft.real_arraydim_2 = 1;
        global.fft.real_arraydim_3 = global.field.Nz;
        
        global.fft.real_arraydim_1 = spectralTransform.local_Nz;
        global.fft.real_arraydim_2 = 1;
        global.fft.real_arraydim_3 = global.field.Nx;
    }
    
    
    // alias..
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
    
    complex_arraydim_1 = global.fft.complex_arraydim_1;
    complex_arraydim_2 = global.fft.complex_arraydim_2;
    complex_arraydim_3 = global.fft.complex_arraydim_3;
    
    real_arraydim_1 = global.fft.real_arraydim_1;
    real_arraydim_2 = global.fft.real_arraydim_2;
    real_arraydim_3 = global.fft.real_arraydim_3;
	
		// temp arrays
	
    global.temp_array.X.resize(complex_arraydim_1, complex_arraydim_2, complex_arraydim_3);
	global.temp_array.X2.resize(complex_arraydim_1, complex_arraydim_2, complex_arraydim_3);
    
	global.temp_array.Xr.resize(real_arraydim_1,real_arraydim_2,real_arraydim_3);
    global.temp_array.Xr2.resize(real_arraydim_1,real_arraydim_2,real_arraydim_3);
	
    
	//*****************************************************************************************
	//set BasicIO parameters
    // FOR READIND full array
	
   // DATASPACE
	// Along X
	BasicIO::rfull_dataspace_dimension[0] = Nx;
	BasicIO::rfull_dataspace_start[0] = my_id*local_Nx;
	BasicIO::rfull_dataspace_blockcount[0] = 1;
	BasicIO::rfull_dataspace_blockstride[0] = 1;
	BasicIO::rfull_dataspace_blockdim[0] = local_Nx;
	
	BasicIO::rfull_dataspace_dimension[1] = Ny;
	BasicIO::rfull_dataspace_start[1] = 0;
	BasicIO::rfull_dataspace_blockcount[1] = 1;
	BasicIO::rfull_dataspace_blockstride[1] = 1;
	BasicIO::rfull_dataspace_blockdim[1] = Ny;
	
	BasicIO::rfull_dataspace_dimension[2] = Nz*2;
	BasicIO::rfull_dataspace_start[2] = 0;
	BasicIO::rfull_dataspace_blockcount[2] = 1;
	BasicIO::rfull_dataspace_blockstride[2] = 1;
	BasicIO::rfull_dataspace_blockdim[2] = Nz*2;
	
	// MEM SPACE (A Array)
	BasicIO::rfull_memspace_dimension[0] = local_Nx;
	BasicIO::rfull_memspace_start[0] = 0;
	BasicIO::rfull_memspace_blockcount[0] = 1;
	BasicIO::rfull_memspace_blockstride[0] = 1;
	BasicIO::rfull_memspace_blockdim[0] = local_Nx;
	
	BasicIO::rfull_memspace_dimension[1] = Ny;
	BasicIO::rfull_memspace_start[1] = 0;
	BasicIO::rfull_memspace_blockcount[1] = 1;
	BasicIO::rfull_memspace_blockstride[1] = 1;
	BasicIO::rfull_memspace_blockdim[1] = Ny;
	
	BasicIO::rfull_memspace_dimension[2] = Nz*2;
	BasicIO::rfull_memspace_start[2] = 0;
	BasicIO::rfull_memspace_blockcount[2] = 1;
	BasicIO::rfull_memspace_blockstride[2] = 1;
	BasicIO::rfull_memspace_blockdim[2] = Nz*2;

	//*****************************************************************************************
	
    // FOR WRITING FULL ARRAY
    // DATASPACE
    BasicIO::wfull_dataspace_dimension[0] = Nx;
    BasicIO::wfull_dataspace_start[0] = my_id*local_Nx;
    BasicIO::wfull_dataspace_blockstride[0] = 1;
    BasicIO::wfull_dataspace_blockcount[0] = 1;
    BasicIO::wfull_dataspace_blockdim[0] = local_Nx;
    
    BasicIO::wfull_dataspace_dimension[1] = Ny;
    BasicIO::wfull_dataspace_start[1] = 0;
    BasicIO::wfull_dataspace_blockstride[1] = 1;
    BasicIO::wfull_dataspace_blockcount[1] = 1;
    BasicIO::wfull_dataspace_blockdim[1] = Ny;
    
    BasicIO::wfull_dataspace_dimension[2] = Nz*2;
    BasicIO::wfull_dataspace_start[2] = 0;
    BasicIO::wfull_dataspace_blockstride[2] = 1;
    BasicIO::wfull_dataspace_blockcount[2] = 1;
    BasicIO::wfull_dataspace_blockdim[2] = Nz*2;
    
    // MEMSPACE
    BasicIO::wfull_memspace_dimension[0] = local_Nx;
    BasicIO::wfull_memspace_start[0] = 0;
    BasicIO::wfull_memspace_blockstride[0] = 1;
    BasicIO::wfull_memspace_blockcount[0] = 1;
    BasicIO::wfull_memspace_blockdim[0] = local_Nx;
    
    BasicIO::wfull_memspace_dimension[1] = Ny;
    BasicIO::wfull_memspace_start[1] = 0;
    BasicIO::wfull_memspace_blockstride[1] = 1;
    BasicIO::wfull_memspace_blockcount[1] = 1;
    BasicIO::wfull_memspace_blockdim[1] = Ny;
    
    BasicIO::wfull_memspace_dimension[2] = Nz*2;
    BasicIO::wfull_memspace_start[2] = 0;
    BasicIO::wfull_memspace_blockstride[2] = 1;
    BasicIO::wfull_memspace_blockcount[2] = 1;
    BasicIO::wfull_memspace_blockdim[2] = Nz*2;    
}

