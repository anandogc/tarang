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

#include "FFFW_slab.h"
#include "BasicIO.h"

FFFW_SLAB::FFFW_SLAB()
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
	
	else  {
		global.field.kfactor.resize(4);
		global.field.kfactor[0] = 1.0;
		
		if (global.program.kind == "RBC") {
			global.field.kfactor[1] = 1.0;
			global.field.kfactor[2] = 1.0;
			global.field.kfactor[3] = 1.0;
		}
		else {
			global.field.kfactor[1] = 1.0;
			global.field.kfactor[2] = 1.0;
			global.field.kfactor[3] = 1.0;
		}
		
		global.field.L.resize(4);
		global.field.L[0] = 1.0;
		for (int i=1; i<=3; i++) 
			global.field.L[i] = 2*M_PI/global.field.kfactor[i];
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
		Nzproperdiv = ((global.field.Nz/2+1)%global.mpi.numprocs == 0);
		Nyproperdiv = true;
	}
	
	if (!(Nxproperdiv && Nyproperdiv && Nzproperdiv)) {
		cout << "N div " << global.field.Nx << " "<< global.field.Ny << " " << global.field.Nz << " " << Nxproperdiv << " " << Nyproperdiv << "  " << Nzproperdiv << endl;
		if (global.mpi.master) cerr << "ERROR in FFF_slab.cc: Array not being divided equally.  Check dimensions" << endl;
		exit(1);
	}
    
    spectralTransform.Init("FFFW", "SLAB", Nx, Ny, Nz);
	
	global.fft.maxlx = global.fft.slab.local_Nx-1;
	global.fft.maxly = global.field.Ny-1;
	global.fft.maxlz = global.field.Nz/2;
    
    if (Ny > 1) {
        global.fft.complex_arraydim_1 = spectralTransform.local_Nx;
        global.fft.complex_arraydim_2 = global.field.Ny;
        global.fft.complex_arraydim_3 = global.field.Nz/2+1;
        
        global.fft.real_arraydim_1 = spectralTransform.local_Ny;
        global.fft.real_arraydim_2 = global.field.Nx;
        global.fft.real_arraydim_3 = global.field.Nz/2+1;
    }
    
    else if (Ny>1) {
        global.fft.complex_arraydim_1 = spectralTransform.local_Nx;
        global.fft.complex_arraydim_2 = 1
        global.fft.complex_arraydim_3 = global.field.Nz/2+1;
        
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
    
    // Not defined for 2D; Supposed to use original FFTW for 2D
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
	
	BasicIO::rfull_dataspace_dimension[2] = (Nz/2+1)*2;
	BasicIO::rfull_dataspace_start[2] = 0;
	BasicIO::rfull_dataspace_blockcount[2] = 1;
	BasicIO::rfull_dataspace_blockstride[2] = 1;
	BasicIO::rfull_dataspace_blockdim[2] = (Nz/2+1)*2;
	
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
	
	BasicIO::rfull_memspace_dimension[2] = (Nz/2+1)*2;
	BasicIO::rfull_memspace_start[2] = 0;
	BasicIO::rfull_memspace_blockcount[2] = 1;
	BasicIO::rfull_memspace_blockstride[2] = 1;
	BasicIO::rfull_memspace_blockdim[2] = (Nz/2+1)*2;
	
	// for kz=0 plane
	BasicIO::rkz0_full_dataspace_dimension[0] = Nx;
	BasicIO::rkz0_full_dataspace_blockdim[0] = local_Nx;
	BasicIO::rkz0_full_dataspace_start[0] = my_id*local_Nx;
	BasicIO::rkz0_full_dataspace_blockcount[0] = 1;
	BasicIO::rkz0_full_dataspace_blockstride[0] = 1;
	
	BasicIO::rkz0_full_dataspace_dimension[1] = Ny;
	BasicIO::rkz0_full_dataspace_blockdim[1] = Ny;
	BasicIO::rkz0_full_dataspace_start[1] = 0;
	BasicIO::rkz0_full_dataspace_blockcount[1] = 1;
	BasicIO::rkz0_full_dataspace_blockstride[1] = 1;
	
	BasicIO::rkz0_full_dataspace_dimension[2] = 2;
	BasicIO::rkz0_full_dataspace_blockdim[2] = 2;
	BasicIO::rkz0_full_dataspace_start[2] = 0;
	BasicIO::rkz0_full_dataspace_blockcount[2] = 1;
	BasicIO::rkz0_full_dataspace_blockstride[2] = 1;
	
	
	
	BasicIO::rkz0_full_memspace_dimension[0] = local_Nx;
	BasicIO::rkz0_full_memspace_blockdim[0] = local_Nx;
	BasicIO::rkz0_full_memspace_start[0] = 0;
	BasicIO::rkz0_full_memspace_blockcount[0] = 1;
	BasicIO::rkz0_full_memspace_blockstride[0] = 1;
	
	BasicIO::rkz0_full_memspace_dimension[1] = Ny;
	BasicIO::rkz0_full_memspace_blockdim[1] = Ny;
	BasicIO::rkz0_full_memspace_start[1] = 0;
	BasicIO::rkz0_full_memspace_blockcount[1] = 1;
	BasicIO::rkz0_full_memspace_blockstride[1] = 1;
	
	
	BasicIO::rkz0_full_memspace_dimension[2] = (Nz/2+1)*2;
	BasicIO::rkz0_full_memspace_blockdim[2] = 2;
	BasicIO::rkz0_full_memspace_start[2] = 0;
	BasicIO::rkz0_full_memspace_blockcount[2] = 1;
	BasicIO::rkz0_full_memspace_blockstride[2] = 1;
    
	//*****************************************************************************************
	
	// Reading reduced arrays
	if (global.io.N_in_reduced.size()==4) {
        // DATASPACE
        // Along X

        if (numprocs == 1) {
            BasicIO::rreduced_dataspace_dimension[0] = global.io.N_in_reduced[1];
            BasicIO::rreduced_dataspace_start[0] = 0;
            BasicIO::rreduced_dataspace_blockcount[0] = 1;
            BasicIO::rreduced_dataspace_blockstride[0] = 1;
            BasicIO::rreduced_dataspace_blockdim[0] = global.io.N_in_reduced[1];
                  
            // MEMSPACE (in Array A)
            BasicIO::rreduced_memspace_dimension[0] = Nx;
            BasicIO::rreduced_memspace_start[0] = 0;  
            BasicIO::rreduced_memspace_blockcount[0] = 2;
            BasicIO::rreduced_memspace_blockdim[0] = global.io.N_in_reduced[1]/2;
            BasicIO::rreduced_memspace_blockstride[0] = Nx-global.io.N_in_reduced[1]/2;
            // for the second block
            
			BasicIO::rkz0_reduced_dataspace_dimension[0] =  global.io.N_in_reduced[1];
			BasicIO::rkz0_reduced_dataspace_start[0] = 0;  
            BasicIO::rkz0_reduced_dataspace_blockcount[0] = 1;
            BasicIO::rkz0_reduced_dataspace_blockstride[0] = 1;
            BasicIO::rkz0_reduced_dataspace_blockdim[0] = global.io.N_in_reduced[1];
			
            BasicIO::rkz0_reduced_memspace_dimension[0] = Nx;
			BasicIO::rkz0_reduced_memspace_start[0] = 0;  
            BasicIO::rkz0_reduced_memspace_blockcount[0] = 2;
            BasicIO::rkz0_reduced_memspace_blockstride[0] = Nx-global.io.N_in_reduced[1]/2;
            BasicIO::rkz0_reduced_memspace_blockdim[0] = global.io.N_in_reduced[1]/2;
        }
        
        else {
            int my_data_size;
            int	remaining_data;
            int last_left_proc, first_right_proc;
            int first_i, xdim_filled_in_rightprocs;
			
			if (((global.io.N_in_reduced[1]/2) % local_Nx) == 0)
				last_left_proc = ((global.io.N_in_reduced[1]/2)/local_Nx) -1 ;
			else 
				last_left_proc = ((global.io.N_in_reduced[1]/2)/local_Nx);
			
			first_right_proc = numprocs-last_left_proc-1;
			
            BasicIO::rreduced_dataspace_dimension[0] = global.io.N_in_reduced[1];
            BasicIO::rkz0_reduced_dataspace_dimension[0] =  global.io.N_in_reduced[1];
			
			BasicIO::rreduced_memspace_dimension[0] = local_Nx;
			BasicIO::rkz0_reduced_memspace_dimension[0] = local_Nx;
            
            if (my_id <= last_left_proc) {
                BasicIO::rreduced_dataspace_blockcount[0] = 1;
                BasicIO::rreduced_dataspace_start[0] = my_id*local_Nx;
                
                remaining_data = global.io.N_in_reduced[1]/2 - my_id*local_Nx;
                my_data_size = minimum(local_Nx, remaining_data);
                
                BasicIO::rreduced_dataspace_blockstride[0] = 1;
                BasicIO::rreduced_dataspace_blockdim[0] = my_data_size;
                
				BasicIO::rreduced_memspace_blockcount[0] = 1;
                BasicIO::rreduced_memspace_start[0] = 0;
                BasicIO::rreduced_memspace_blockstride[0]=1;
                BasicIO::rreduced_memspace_blockdim[0] = my_data_size;
			
				
				BasicIO::rkz0_reduced_dataspace_blockcount[0] = 1;
				BasicIO::rkz0_reduced_dataspace_start[0] = my_id*local_Nx; 
				BasicIO::rkz0_reduced_dataspace_blockstride[0] = 1;
				BasicIO::rkz0_reduced_dataspace_blockdim[0] = my_data_size;
				
				BasicIO::rkz0_reduced_memspace_blockcount[0] = 1;
				BasicIO::rkz0_reduced_memspace_start[0] = 0;  
				BasicIO::rkz0_reduced_memspace_blockstride[0] = 1;
				BasicIO::rkz0_reduced_memspace_blockdim[0] = my_data_size;
            }
            
            else if (my_id >= first_right_proc) {
                
                xdim_filled_in_rightprocs = (numprocs-1-my_id)*local_Nx;
                remaining_data = (global.io.N_in_reduced[1]/2) - xdim_filled_in_rightprocs;
                
                if (remaining_data >= local_Nx) { // present proc inside neg k regime
                    BasicIO::rreduced_dataspace_blockcount[0] = 1;
                    BasicIO::rreduced_dataspace_blockdim[0] = local_Nx;
		
					BasicIO::rreduced_dataspace_start[0] = remaining_data-local_Nx+(global.io.N_in_reduced[1]/2); 
                    BasicIO::rreduced_dataspace_blockstride[0] = 1;
                    
                    
                    BasicIO::rreduced_memspace_blockcount[0] = 1;
                    BasicIO::rreduced_memspace_start[0] = 0;
                    BasicIO::rreduced_memspace_blockstride[0]=1;
                    BasicIO::rreduced_memspace_blockdim[0] = local_Nx;
					
					
					BasicIO::rkz0_reduced_dataspace_blockcount[0] = 1;
					BasicIO::rkz0_reduced_dataspace_blockdim[0] = local_Nx;
					BasicIO::rkz0_reduced_dataspace_start[0] = remaining_data-local_Nx+(global.io.N_in_reduced[1]/2); 
					BasicIO::rkz0_reduced_dataspace_blockstride[0] = 1;
					
					BasicIO::rkz0_reduced_memspace_blockcount[0] = 1;
					BasicIO::rkz0_reduced_memspace_start[0] = 0;  
					BasicIO::rkz0_reduced_memspace_blockstride[0] = 1;
					BasicIO::rkz0_reduced_memspace_blockdim[0] = local_Nx;
                }
                else { // present proc contains part of neg k regime
                    BasicIO::rreduced_dataspace_blockcount[0] = 1;
                    BasicIO::rreduced_dataspace_blockdim[0] = remaining_data;
                    BasicIO::rreduced_dataspace_start[0] = (global.io.N_in_reduced[1]/2); 
                    BasicIO::rreduced_dataspace_blockstride[0] = 1;
                    
                    BasicIO::rreduced_memspace_blockcount[0] = 1;
                    BasicIO::rreduced_memspace_start[0] = local_Nx-remaining_data;
                    BasicIO::rreduced_memspace_blockstride[0]=1;
                    BasicIO::rreduced_memspace_blockdim[0] = remaining_data;
					
					BasicIO::rkz0_reduced_dataspace_blockcount[0] = 1;
					BasicIO::rkz0_reduced_dataspace_blockdim[0] = remaining_data;
					BasicIO::rkz0_reduced_dataspace_start[0] = (global.io.N_in_reduced[1]/2);  
					BasicIO::rkz0_reduced_dataspace_blockstride[0] = 1;
					
					BasicIO::rkz0_reduced_memspace_blockcount[0] = 1;
					BasicIO::rkz0_reduced_memspace_start[0] = local_Nx-remaining_data;  
					BasicIO::rkz0_reduced_memspace_blockstride[0] = 1;
					BasicIO::rkz0_reduced_memspace_blockdim[0] = remaining_data;
                }
            }
        
            else {// Dont read
                BasicIO::rreduced_dataspace_blockcount[0] = 0;
				BasicIO::rreduced_dataspace_blockdim[0] =0;
                BasicIO::rreduced_dataspace_start[0] = 0; 
				BasicIO::rreduced_dataspace_blockstride[0] = 1;
				
				BasicIO::rkz0_reduced_dataspace_blockcount[0] = 0;
				BasicIO::rkz0_reduced_dataspace_start[0] = 0;  
				BasicIO::rkz0_reduced_dataspace_blockstride[0] = 1;
                BasicIO::rkz0_reduced_dataspace_blockdim[0] = 0;
                
                BasicIO::rreduced_memspace_blockcount[0] = 0;
                BasicIO::rreduced_memspace_start[0] = 0;
                BasicIO::rreduced_memspace_blockstride[0]=1;
                BasicIO::rreduced_memspace_blockdim[0] = 0;
                
				BasicIO::rkz0_reduced_memspace_blockcount[0] = 0;
				BasicIO::rkz0_reduced_memspace_start[0] = 0;  
				BasicIO::rkz0_reduced_memspace_blockstride[0] = 1;
                BasicIO::rkz0_reduced_memspace_blockdim[0] = 0;
            }
        }
        
        // Along Y
		if (Ny > 1) {
			BasicIO::rreduced_dataspace_dimension[1] = global.io.N_in_reduced[2];
			BasicIO::rreduced_dataspace_start[1] = 0;
			BasicIO::rreduced_dataspace_blockcount[1] = 1;
			BasicIO::rreduced_dataspace_blockstride[1] = 1;
			BasicIO::rreduced_dataspace_blockdim[1] = global.io.N_in_reduced[2];
			
			BasicIO::rreduced_memspace_dimension[1] = Ny;
			BasicIO::rreduced_memspace_start[1] = 0;  
			BasicIO::rreduced_memspace_blockcount[1] = 2;
            BasicIO::rreduced_memspace_blockstride[1] = Ny-global.io.N_in_reduced[2]/2; // for the second block
//			BasicIO::rreduced_memspace_blockstride[1] = Ny-global.io.N_in_reduced[2]/2; // for the second block
			BasicIO::rreduced_memspace_blockdim[1] = global.io.N_in_reduced[2]/2; 
			
			BasicIO::rkz0_reduced_dataspace_dimension[1] =  global.io.N_in_reduced[2];
			BasicIO::rkz0_reduced_dataspace_start[1] = 0;
			BasicIO::rkz0_reduced_dataspace_blockcount[1] = 1;
			BasicIO::rkz0_reduced_dataspace_blockstride[1] = 1;
			BasicIO::rkz0_reduced_dataspace_blockdim[1] = global.io.N_in_reduced[2];

			BasicIO::rkz0_reduced_memspace_dimension[1] = Ny;
			BasicIO::rkz0_reduced_memspace_start[1] = 0;
			BasicIO::rkz0_reduced_memspace_blockcount[1] = 2;
			BasicIO::rkz0_reduced_memspace_blockstride[1] = Ny-global.io.N_in_reduced[2]/2;
			BasicIO::rkz0_reduced_memspace_blockdim[1] = global.io.N_in_reduced[2]/2; 
		}
		else if (Ny==1) {
			BasicIO::rreduced_dataspace_dimension[1] = 1;
			BasicIO::rreduced_dataspace_start[1] = 0;
			BasicIO::rreduced_dataspace_blockcount[1] = 1;
			BasicIO::rreduced_dataspace_blockstride[1] = 1;
			BasicIO::rreduced_dataspace_blockdim[1] = 1;
			
			BasicIO::rreduced_memspace_dimension[1] = 1;
			BasicIO::rreduced_memspace_start[1] = 0;  
			BasicIO::rreduced_memspace_blockcount[1] = 1;
			BasicIO::rreduced_memspace_blockstride[1] = 1; // for the second block
			BasicIO::rreduced_memspace_blockdim[1] = 1; 
			
			BasicIO::rkz0_reduced_dataspace_dimension[1] =  global.io.N_in_reduced[2];
			BasicIO::rkz0_reduced_dataspace_start[1] = 0;
			BasicIO::rkz0_reduced_dataspace_blockcount[1] = 1;
			BasicIO::rkz0_reduced_dataspace_blockstride[1] = 1;
			BasicIO::rkz0_reduced_dataspace_blockdim[1] = global.io.N_in_reduced[2];

			BasicIO::rkz0_reduced_memspace_dimension[1] = Ny;
			BasicIO::rkz0_reduced_memspace_start[1] = 0;
			BasicIO::rkz0_reduced_memspace_blockcount[1] = 1;
			BasicIO::rkz0_reduced_memspace_blockstride[1] = 1;
			BasicIO::rkz0_reduced_memspace_blockdim[1] = 1; 
		}
        
        // Along Z
        BasicIO::rreduced_dataspace_dimension[2] = (global.io.N_in_reduced[3]/2+1)*2;
        BasicIO::rreduced_dataspace_start[2] = 0;
        BasicIO::rreduced_dataspace_blockcount[2] = 1;
        BasicIO::rreduced_dataspace_blockstride[2] = 1;
        BasicIO::rreduced_dataspace_blockdim[2] = (global.io.N_in_reduced[3]/2+1)*2;
       
        BasicIO::rreduced_memspace_dimension[2] = (Nz/2+1)*2;  // for real conversion
        BasicIO::rreduced_memspace_start[2] = 0;
        BasicIO::rreduced_memspace_blockcount[2] = 1;
        BasicIO::rreduced_memspace_blockstride[2] = 1;  
        BasicIO::rreduced_memspace_blockdim[2] = ((global.io.N_in_reduced[3]/2)+1)*2; 
        
		
		BasicIO::rkz0_reduced_dataspace_dimension[2] = 2;
		BasicIO::rkz0_reduced_dataspace_start[2] = 0;  
		BasicIO::rkz0_reduced_dataspace_blockcount[2] = 1;
		BasicIO::rkz0_reduced_dataspace_blockstride[2] = 1;
		BasicIO::rkz0_reduced_dataspace_blockdim[2] = 2;
		
		BasicIO::rkz0_reduced_memspace_dimension[2] = (Nz/2+1)*2;
		BasicIO::rkz0_reduced_memspace_start[2] = 0;  
		BasicIO::rkz0_reduced_memspace_blockcount[2] = 1;
		BasicIO::rkz0_reduced_memspace_blockstride[2] = 1;
		BasicIO::rkz0_reduced_memspace_blockdim[2] = 2; 
	}
	
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
    
    BasicIO::wfull_dataspace_dimension[2] = (Nz/2+1)*2;
    BasicIO::wfull_dataspace_start[2] = 0;
    BasicIO::wfull_dataspace_blockstride[2] = 1;
    BasicIO::wfull_dataspace_blockcount[2] = 1;
    BasicIO::wfull_dataspace_blockdim[2] = (Nz/2+1)*2;
    
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
    
    BasicIO::wfull_memspace_dimension[2] = (Nz/2+1)*2;
    BasicIO::wfull_memspace_start[2] = 0;
    BasicIO::wfull_memspace_blockstride[2] = 1;
    BasicIO::wfull_memspace_blockcount[2] = 1;
    BasicIO::wfull_memspace_blockdim[2] = (Nz/2+1)*2;
	
	// kz=0 plane
	BasicIO::wkz0_dataspace_dimension[0] = Nx;
    BasicIO::wkz0_dataspace_blockdim[0] = local_Nx;
	BasicIO::wkz0_dataspace_start[0] = my_id*local_Nx;
    BasicIO::wkz0_dataspace_blockstride[0] = 1;
    BasicIO::wkz0_dataspace_blockcount[0] = 1;
	
    BasicIO::wkz0_dataspace_dimension[1] = Ny;
    BasicIO::wkz0_dataspace_blockdim[1] = Ny;
	BasicIO::wkz0_dataspace_start[1] = 0;
    BasicIO::wkz0_dataspace_blockstride[1] = 1;
    BasicIO::wkz0_dataspace_blockcount[1] = 1;
	
    BasicIO::wkz0_dataspace_dimension[2] = 2;
    BasicIO::wkz0_dataspace_blockdim[2] = 2;
	BasicIO::wkz0_dataspace_start[2] = 0;
    BasicIO::wkz0_dataspace_blockstride[2] = 1;
    BasicIO::wkz0_dataspace_blockcount[2] = 1;
    
    BasicIO::wkz0_memspace_dimension[0] = local_Nx;
    BasicIO::wkz0_memspace_blockdim[0] = local_Nx;
	BasicIO::wkz0_memspace_start[0] = 0;
    BasicIO::wkz0_memspace_blockstride[0] = 1;
    BasicIO::wkz0_memspace_blockcount[0] = 1;
	
    BasicIO::wkz0_memspace_dimension[1] = Ny;
    BasicIO::wkz0_memspace_blockdim[1] = Ny;
	BasicIO::wkz0_memspace_start[1] = 0;
    BasicIO::wkz0_memspace_blockstride[1] = 1;
    BasicIO::wkz0_memspace_blockcount[1] = 1;
	
    BasicIO::wkz0_memspace_dimension[2] = (Nz/2+1)*2;
    BasicIO::wkz0_memspace_blockdim[2] = 2;
	BasicIO::wkz0_memspace_start[2] = 0;
    BasicIO::wkz0_memspace_blockstride[2] = 1;
    BasicIO::wkz0_memspace_blockcount[2] = 1;

	
	//*****************************************************************************************
    
	// Writing reduced Array
    // 0:Nx_red/2 (Nx_red/2+1) data (left)
    //  then Nx_red/2+1: Nx_red-1 (Nx_red/2-1) points (right)
    
	if (global.io.N_out_reduced.size()==4) {
		
		// DATASPACE
		// Along X
        if (numprocs == 1) {
            BasicIO::wreduced_dataspace_dimension[0] = global.io.N_out_reduced[1];
            BasicIO::wreduced_dataspace_start[0] = 0;
            BasicIO::wreduced_dataspace_blockcount[0] = 1;
            BasicIO::wreduced_dataspace_blockstride[0] = 1;
            BasicIO::wreduced_dataspace_blockdim[0] = global.io.N_out_reduced[1];
			
			// MEMSPACE (in Array A)
            BasicIO::wreduced_memspace_dimension[0] = Nx;
            BasicIO::wreduced_memspace_start[0] = 0;  
            BasicIO::wreduced_memspace_blockcount[0] = 2;
            BasicIO::wreduced_memspace_blockstride[0] = Nx-global.io.N_out_reduced[1]/2;
				// for the second block
            BasicIO::wreduced_memspace_blockdim[0] = global.io.N_out_reduced[1]/2;
        }
        
        else {
            int my_data_size;
            int	remaining_data;
            int last_left_proc, first_right_proc;
            int first_i, xdim_filled_in_rightprocs;
			
			if (((global.io.N_in_reduced[1]/2) % local_Nx) == 0)
				last_left_proc = ((global.io.N_out_reduced[1]/2)/local_Nx) -1 ;
			else 
				last_left_proc = ((global.io.N_out_reduced[1]/2)/local_Nx);
            
            if (((global.io.N_in_reduced[1]/2) % local_Nx) == 0)
				first_right_proc = numprocs-((global.io.N_out_reduced[1]/2)/local_Nx);
			else
				first_right_proc = numprocs-1 -((global.io.N_out_reduced[1]/2)/local_Nx);
			
            BasicIO::wreduced_dataspace_dimension[0] = global.io.N_out_reduced[1];
			BasicIO::wreduced_memspace_dimension[0] = local_Nx;
            
            if (my_id <= last_left_proc) {
                remaining_data = global.io.N_out_reduced[1]/2 - my_id*local_Nx;
                my_data_size = minimum(local_Nx, remaining_data);
                
                BasicIO::wreduced_dataspace_start[0] = my_id*local_Nx;
                BasicIO::wreduced_dataspace_blockdim[0] = my_data_size;
                BasicIO::wreduced_dataspace_blockcount[0] = 1;
                BasicIO::wreduced_dataspace_blockstride[0] = 1;
                
                BasicIO::wreduced_memspace_start[0] = 0;
                BasicIO::wreduced_memspace_blockdim[0] = my_data_size;
                BasicIO::wreduced_memspace_blockcount[0] = 1;
                BasicIO::wreduced_memspace_blockstride[0]=1;
            }
            
            else if (my_id >= first_right_proc) {
                
                xdim_filled_in_rightprocs = (numprocs-1-my_id)*local_Nx;
                remaining_data = (global.io.N_out_reduced[1]/2) - xdim_filled_in_rightprocs;
                
                if (remaining_data >= local_Nx) { // present proc inside neg k regime
					BasicIO::wreduced_dataspace_start[0] = remaining_data-local_Nx+(global.io.N_out_reduced[1]/2);
                    BasicIO::wreduced_dataspace_blockdim[0] = local_Nx;
                    BasicIO::wreduced_dataspace_blockcount[0] = 1;
                    BasicIO::wreduced_dataspace_blockstride[0] = 1;
                    
                    BasicIO::wreduced_memspace_start[0] = 0;
                    BasicIO::wreduced_memspace_blockdim[0] = local_Nx;
                    BasicIO::wreduced_memspace_blockcount[0] = 1;
                    BasicIO::wreduced_memspace_blockstride[0]=1;
                }
                else { // present proc contains part of neg k regme
                    BasicIO::wreduced_dataspace_start[0] = (global.io.N_out_reduced[1]/2);
                    BasicIO::wreduced_dataspace_blockdim[0] = remaining_data;
                    BasicIO::wreduced_dataspace_blockcount[0] = 1;
                    BasicIO::wreduced_dataspace_blockstride[0] = 1;
                    
                    BasicIO::wreduced_memspace_start[0] = local_Nx-remaining_data;
                    BasicIO::wreduced_memspace_blockdim[0] = remaining_data;
                    BasicIO::wreduced_memspace_blockcount[0] = 1;
                    BasicIO::wreduced_memspace_blockstride[0]=1;
                }
            }
			
            else {// Dont read
                BasicIO::wreduced_dataspace_blockcount[0] = 0;
				BasicIO::wreduced_dataspace_blockdim[0] =0;
                BasicIO::wreduced_dataspace_start[0] = 0; 
				BasicIO::wreduced_dataspace_blockstride[0] = 1;
                
                BasicIO::wreduced_memspace_blockcount[0] = 0;
				BasicIO::wreduced_memspace_blockdim[0] = 0;
                BasicIO::wreduced_memspace_start[0] = 0;
                BasicIO::wreduced_memspace_blockstride[0]=1;
            }
        }
        
		// Along Y
		if (Ny > 1) {
			BasicIO::wreduced_dataspace_dimension[1] = global.io.N_out_reduced[2];
			BasicIO::wreduced_dataspace_start[1] = 0;
			BasicIO::wreduced_dataspace_blockcount[1] = 1;
			BasicIO::wreduced_dataspace_blockstride[1] = 1;
			BasicIO::wreduced_dataspace_blockdim[1] = global.io.N_out_reduced[2];
			
			BasicIO::wreduced_memspace_dimension[1] = Ny;
			BasicIO::wreduced_memspace_start[1] = 0;  
			BasicIO::wreduced_memspace_blockcount[1] = 2;
			BasicIO::wreduced_memspace_blockstride[1] = Ny-global.io.N_out_reduced[2]/2; // for the second block
			BasicIO::wreduced_memspace_blockdim[1] = global.io.N_out_reduced[2]/2;
		}
		
		else if (Ny == 1) {
			BasicIO::wreduced_dataspace_dimension[1] = 1;
			BasicIO::wreduced_dataspace_start[1] = 0;
			BasicIO::wreduced_dataspace_blockcount[1] = 1;
			BasicIO::wreduced_dataspace_blockstride[1] = 1;
			BasicIO::wreduced_dataspace_blockdim[1] = 1;
			
			BasicIO::wreduced_memspace_dimension[1] = 1;
			BasicIO::wreduced_memspace_start[1] = 0;  
			BasicIO::wreduced_memspace_blockcount[1] = 1;
			BasicIO::wreduced_memspace_blockstride[1] = 1; // for the second block
			BasicIO::wreduced_memspace_blockdim[1] = 1;
		}
        
        
		// Along Z
        BasicIO::wreduced_dataspace_dimension[2] = (global.io.N_out_reduced[3]/2+1)*2;
        BasicIO::wreduced_dataspace_start[2] = 0;
        BasicIO::wreduced_dataspace_blockcount[2] = 1;
        BasicIO::wreduced_dataspace_blockstride[2] = 1;
        BasicIO::wreduced_dataspace_blockdim[2] = (global.io.N_out_reduced[3]/2+1)*2;
		
        BasicIO::wreduced_memspace_dimension[2] = (Nz/2+1)*2;  // for real conversion
        BasicIO::wreduced_memspace_start[2] = 0;
        BasicIO::wreduced_memspace_blockcount[2] = 1;
        BasicIO::wreduced_memspace_blockstride[2] = 1;  
        BasicIO::wreduced_memspace_blockdim[2] = ((global.io.N_out_reduced[3]/2)+1)*2; 
	}
	
	//*****************************************************************************************
	
	// for reading real space data
    if (Ny < 1) {
        if (!global.fft.transpose) {
            BasicIO::r_real_dataspace_dimension[0] = Nx;
            BasicIO::r_real_dataspace_blockdim[0] = local_Nx;
            BasicIO::r_real_dataspace_start[0] = my_id*local_Nx;
            BasicIO::r_real_dataspace_blockcount[0] = 1;
            BasicIO::r_real_dataspace_blockstride[0] = 1;
            
            BasicIO::r_real_dataspace_dimension[1] = Ny;
            BasicIO::r_real_dataspace_blockdim[1] = Ny;
            BasicIO::r_real_dataspace_start[1] = 0;
            BasicIO::r_real_dataspace_blockcount[1] = 1;
            BasicIO::r_real_dataspace_blockstride[1] = 1;
            
            BasicIO::r_real_memspace_dimension[0] = local_Nx;
            BasicIO::r_real_memspace_blockdim[0] = local_Nx;
            BasicIO::r_real_memspace_start[0] = 0;
            BasicIO::r_real_memspace_blockcount[0] = 1;
            BasicIO::r_real_memspace_blockstride[0] = 1;
            
            BasicIO::r_real_memspace_dimension[1] = Ny;
            BasicIO::r_real_memspace_blockdim[1] = Ny;
            BasicIO::r_real_memspace_start[1] = 0;
            BasicIO::r_real_memspace_blockcount[1] = 1;
            BasicIO::r_real_memspace_blockstride[1] = 1;
        }
        
        else {
            BasicIO::r_real_dataspace_dimension[0] = Ny;
            BasicIO::r_real_dataspace_blockdim[0] = local_Ny;
            BasicIO::r_real_dataspace_start[0] = my_id*local_Ny;
            BasicIO::r_real_dataspace_blockcount[0] = 1;
            BasicIO::r_real_dataspace_blockstride[0] = 1;
            
            BasicIO::r_real_dataspace_dimension[1] = Nx;
            BasicIO::r_real_dataspace_blockdim[1] = Nx;
            BasicIO::r_real_dataspace_start[1] = 0;
            BasicIO::r_real_dataspace_blockcount[1] = 1;
            BasicIO::r_real_dataspace_blockstride[1] = 1;
            
            BasicIO::r_real_memspace_dimension[0] = local_Ny;
            BasicIO::r_real_memspace_blockdim[0] = local_Ny;
            BasicIO::r_real_memspace_start[0] = 0;
            BasicIO::r_real_memspace_blockcount[0] = 1;
            BasicIO::r_real_memspace_blockstride[0] = 1;
            
            BasicIO::r_real_memspace_dimension[1] = Nx;
            BasicIO::r_real_memspace_blockdim[1] = Nx;
            BasicIO::r_real_memspace_start[1] = 0;
            BasicIO::r_real_memspace_blockcount[1] = 1;
            BasicIO::r_real_memspace_blockstride[1] = 1;
        }
        
        BasicIO::r_real_dataspace_dimension[2] = Nz;
        BasicIO::r_real_dataspace_blockdim[2] = Nz;
        BasicIO::r_real_dataspace_start[2] = 0;
        BasicIO::r_real_dataspace_blockcount[2] = 1;
        BasicIO::r_real_dataspace_blockstride[2] = 1;

        
        BasicIO::r_real_memspace_dimension[2] = (Nz/2+1)*2;
        BasicIO::r_real_memspace_blockdim[2] = Nz;
        BasicIO::r_real_memspace_start[2] = 0;
        BasicIO::r_real_memspace_blockcount[2] = 1;
        BasicIO::r_real_memspace_blockstride[2] = 1;
    }
    
    else if (Ny == 1) { // 2d
        if (!global.fft.transpose) {
            BasicIO::r_real_dataspace_dimension[0] = Nx;
            BasicIO::r_real_dataspace_blockdim[0] = local_Nx;
            BasicIO::r_real_dataspace_start[0] = my_id*local_Nx;
            BasicIO::r_real_dataspace_blockcount[0] = 1;
            BasicIO::r_real_dataspace_blockstride[0] = 1;
            
            BasicIO::r_real_dataspace_dimension[1] = 1;
            BasicIO::r_real_dataspace_blockdim[1] = 1;
            BasicIO::r_real_dataspace_start[1] = 0;
            BasicIO::r_real_dataspace_blockcount[1] = 1;
            BasicIO::r_real_dataspace_blockstride[1] = 1;
            
            BasicIO::r_real_dataspace_dimension[2] = Nz;
            BasicIO::r_real_dataspace_blockdim[2] = Nz;
            BasicIO::r_real_dataspace_start[2] = 0;
            BasicIO::r_real_dataspace_blockcount[2] = 1;
            BasicIO::r_real_dataspace_blockstride[2] = 1;
            
            
            BasicIO::r_real_memspace_dimension[0] = local_Nx;
            BasicIO::r_real_memspace_blockdim[0] = local_Nx;
            BasicIO::r_real_memspace_start[0] = 0;
            BasicIO::r_real_memspace_blockcount[0] = 1;
            BasicIO::r_real_memspace_blockstride[0] = 1;
            
            BasicIO::r_real_memspace_dimension[1] = 1;
            BasicIO::r_real_memspace_blockdim[1] = 1;
            BasicIO::r_real_memspace_start[1] = 0;
            BasicIO::r_real_memspace_blockcount[1] = 1;
            BasicIO::r_real_memspace_blockstride[1] = 1;
            
            BasicIO::r_real_memspace_dimension[2] = (Nz/2+1)*2;
            BasicIO::r_real_memspace_blockdim[2] = Nz;
            BasicIO::r_real_memspace_start[2] = 0;
            BasicIO::r_real_memspace_blockcount[2] = 1;
            BasicIO::r_real_memspace_blockstride[2] = 1;
        }
        
        else {
            BasicIO::r_real_dataspace_dimension[0] = Nz;
            BasicIO::r_real_dataspace_blockdim[0] = local_Nz;
            BasicIO::r_real_dataspace_start[0] = my_id*local_Nz;
            BasicIO::r_real_dataspace_blockcount[0] = 1;
            BasicIO::r_real_dataspace_blockstride[0] = 1;
            
            BasicIO::r_real_dataspace_dimension[1] = 1;
            BasicIO::r_real_dataspace_blockdim[1] = 1;
            BasicIO::r_real_dataspace_start[1] = 0;
            BasicIO::r_real_dataspace_blockcount[1] = 1;
            BasicIO::r_real_dataspace_blockstride[1] = 1;
            
            BasicIO::r_real_dataspace_dimension[2] = Nx;
            BasicIO::r_real_dataspace_blockdim[2] = Nx;
            BasicIO::r_real_dataspace_start[2] = 0;
            BasicIO::r_real_dataspace_blockcount[2] = 1;
            BasicIO::r_real_dataspace_blockstride[2] = 1;
            
            
            BasicIO::r_real_memspace_dimension[0] = local_Nz;
            BasicIO::r_real_memspace_blockdim[0] = local_Nz;
            BasicIO::r_real_memspace_start[0] = 0;
            BasicIO::r_real_memspace_blockcount[0] = 1;
            BasicIO::r_real_memspace_blockstride[0] = 1;
            
            BasicIO::r_real_memspace_dimension[1] = 1;
            BasicIO::r_real_memspace_blockdim[1] = 1;
            BasicIO::r_real_memspace_start[1] = 0;
            BasicIO::r_real_memspace_blockcount[1] = 1;
            BasicIO::r_real_memspace_blockstride[1] = 1;
            
            BasicIO::r_real_memspace_dimension[2] = Nx;
            BasicIO::r_real_memspace_blockdim[2] = Nx;
            BasicIO::r_real_memspace_start[2] = 0;
            BasicIO::r_real_memspace_blockcount[2] = 1;
            BasicIO::r_real_memspace_blockstride[2] = 1;
        }
    }
	
	//*****************************************************************************************
	
	// for writing real space
    // Note: 2D uses original FFTW so NO TRANSPOSE ORDER..
    if (Ny > 1) {
        if (!global.fft.transpose) {
            BasicIO::w_real_dataspace_dimension[0] = Nx;
            BasicIO::w_real_dataspace_blockdim[0] = local_Nx;
            BasicIO::w_real_dataspace_start[0] = my_id*local_Nx;
            BasicIO::w_real_dataspace_blockcount[0] = 1;
            BasicIO::w_real_dataspace_blockstride[0] = 1;
            
            BasicIO::w_real_dataspace_dimension[1] = Ny;
            BasicIO::w_real_dataspace_blockdim[1] = Ny;
            BasicIO::w_real_dataspace_start[1] = 0;
            BasicIO::w_real_dataspace_blockcount[1] = 1;
            BasicIO::w_real_dataspace_blockstride[1] = 1;
            
            BasicIO::w_real_memspace_dimension[0] = local_Nx;
            BasicIO::w_real_memspace_blockdim[0] = local_Nx;
            BasicIO::w_real_memspace_start[0] = 0;
            BasicIO::w_real_memspace_blockcount[0] = 1;
            BasicIO::w_real_memspace_blockstride[0] = 1;
            
            BasicIO::w_real_memspace_dimension[1] = Ny;
            BasicIO::w_real_memspace_blockdim[1] = Ny;
            BasicIO::w_real_memspace_start[1] = 0;
            BasicIO::w_real_memspace_blockcount[1] = 1;
            BasicIO::w_real_memspace_blockstride[1] = 1;
        }
        else {
            BasicIO::w_real_dataspace_dimension[0] = Ny;
            BasicIO::w_real_dataspace_blockdim[0] = local_Ny;
            BasicIO::w_real_dataspace_start[0] = my_id*local_Ny;
            BasicIO::w_real_dataspace_blockcount[0] = 1;
            BasicIO::w_real_dataspace_blockstride[0] = 1;
            
            BasicIO::w_real_dataspace_dimension[1] = Nx;
            BasicIO::w_real_dataspace_blockdim[1] = Nx;
            BasicIO::w_real_dataspace_start[1] = 0;
            BasicIO::w_real_dataspace_blockcount[1] = 1;
            BasicIO::w_real_dataspace_blockstride[1] = 1;
        
            BasicIO::w_real_memspace_dimension[0] = local_Ny;
            BasicIO::w_real_memspace_blockdim[0] = local_Ny;
            BasicIO::w_real_memspace_start[0] = 0;
            BasicIO::w_real_memspace_blockcount[0] = 1;
            BasicIO::w_real_memspace_blockstride[0] = 1;
            
            BasicIO::w_real_memspace_dimension[1] = Nx;
            BasicIO::w_real_memspace_blockdim[1] = Nx;
            BasicIO::w_real_memspace_start[1] = 0;
            BasicIO::w_real_memspace_blockcount[1] = 1;
            BasicIO::w_real_memspace_blockstride[1] = 1;
        }
        
        BasicIO::w_real_dataspace_dimension[2] = Nz;
        BasicIO::w_real_dataspace_blockdim[2] = Nz;
        BasicIO::w_real_dataspace_start[2] = 0;
        BasicIO::w_real_dataspace_blockcount[2] = 1;
        BasicIO::w_real_dataspace_blockstride[2] = 1;
        
        
        BasicIO::w_real_memspace_dimension[2] = (Nz/2+1)*2;
        BasicIO::w_real_memspace_blockdim[2] = Nz;
        BasicIO::w_real_memspace_start[2] = 0;
        BasicIO::w_real_memspace_blockcount[2] = 1;
        BasicIO::w_real_memspace_blockstride[2] = 1;
    }
    
    else if (Ny == 1) { // 2d:
        if (!global.fft.transpose) {
            BasicIO::w_real_dataspace_dimension[0] = Nx;
            BasicIO::w_real_dataspace_blockdim[0] = local_Nx;
            BasicIO::w_real_dataspace_start[0] = my_id*local_Nx;
            BasicIO::w_real_dataspace_blockcount[0] = 1;
            BasicIO::w_real_dataspace_blockstride[0] = 1;
            
            BasicIO::w_real_dataspace_dimension[1] = 1;
            BasicIO::w_real_dataspace_blockdim[1] = 1;
            BasicIO::w_real_dataspace_start[1] = 0;
            BasicIO::w_real_dataspace_blockcount[1] = 1;
            BasicIO::w_real_dataspace_blockstride[1] = 1;
            
            BasicIO::w_real_dataspace_dimension[2] = Nz;
            BasicIO::w_real_dataspace_blockdim[2] = Nz;
            BasicIO::w_real_dataspace_start[2] = 0;
            BasicIO::w_real_dataspace_blockcount[2] = 1;
            BasicIO::w_real_dataspace_blockstride[2] = 1;
            
            BasicIO::w_real_memspace_dimension[0] = local_Nx;
            BasicIO::w_real_memspace_blockdim[0] = local_Nx;
            BasicIO::w_real_memspace_start[0] = 0;
            BasicIO::w_real_memspace_blockcount[0] = 1;
            BasicIO::w_real_memspace_blockstride[0] = 1;
            
            BasicIO::w_real_memspace_dimension[1] = 1;
            BasicIO::w_real_memspace_blockdim[1] = 1;
            BasicIO::w_real_memspace_start[1] = 0;
            BasicIO::w_real_memspace_blockcount[1] = 1;
            BasicIO::w_real_memspace_blockstride[1] = 1;
            
            BasicIO::w_real_memspace_dimension[2] = (Nz/2+1)*2;
            BasicIO::w_real_memspace_blockdim[2] = Nz;
            BasicIO::w_real_memspace_start[2] = 0;
            BasicIO::w_real_memspace_blockcount[2] = 1;
            BasicIO::w_real_memspace_blockstride[2] = 1;
        }
        else {
            BasicIO::w_real_dataspace_dimension[0] = Nz;
            BasicIO::w_real_dataspace_blockdim[0] = local_Nz;
            BasicIO::w_real_dataspace_start[0] = my_id*local_Nz;
            BasicIO::w_real_dataspace_blockcount[0] = 1;
            BasicIO::w_real_dataspace_blockstride[0] = 1;
            
            BasicIO::w_real_dataspace_dimension[1] = 1;
            BasicIO::w_real_dataspace_blockdim[1] = 1;
            BasicIO::w_real_dataspace_start[1] = 0;
            BasicIO::w_real_dataspace_blockcount[1] = 1;
            BasicIO::w_real_dataspace_blockstride[1] = 1;
            
            BasicIO::w_real_dataspace_dimension[2] = Nx;
            BasicIO::w_real_dataspace_blockdim[2] = Nx;
            BasicIO::w_real_dataspace_start[2] = 0;
            BasicIO::w_real_dataspace_blockcount[2] = 1;
            BasicIO::w_real_dataspace_blockstride[2] = 1;
            
            BasicIO::w_real_memspace_dimension[0] = local_Nz;
            BasicIO::w_real_memspace_blockdim[0] = local_Nz;
            BasicIO::w_real_memspace_start[0] = 0;
            BasicIO::w_real_memspace_blockcount[0] = 1;
            BasicIO::w_real_memspace_blockstride[0] = 1;
            
            BasicIO::w_real_memspace_dimension[1] = 1;
            BasicIO::w_real_memspace_blockdim[1] = 1;
            BasicIO::w_real_memspace_start[1] = 0;
            BasicIO::w_real_memspace_blockcount[1] = 1;
            BasicIO::w_real_memspace_blockstride[1] = 1;
            
            BasicIO::w_real_memspace_dimension[2] = Nx;
            BasicIO::w_real_memspace_blockdim[2] = Nx;
            BasicIO::w_real_memspace_start[2] = 0;
            BasicIO::w_real_memspace_blockcount[2] = 1;
            BasicIO::w_real_memspace_blockstride[2] = 1;
        }
    }
    
    
}

