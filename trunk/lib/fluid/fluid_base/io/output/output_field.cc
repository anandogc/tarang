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

/*! \file  Output_field.cc
 *
 * @brief  Output_field
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "FluidIO.h"


/**********************************************************************************************
 
 Input and Output of all the components (Arrays)
 
 File name:  
 dir=time 
 
 Filenames: V1, V2, V3, same as fields of CVF, RVF, CSF etc.
 Reduced: V1reduced, V2reduced
 If plane V3kz0
 
 For frequent: USE dir=frequent.time
 
 Do not write field-frequent and regular file at the same time..
 Overwrite field frequent, but not the regular..
 May be just add frequent label with time in the dirname, i.e., frquent.time
 
 use IOdir as global variable and control...
 
 
 **********************************************************************************************/


//*********************************************************************************************

void FluidIO::Output_complex_field(FluidVF& U)
{
	U.cvf.Write_complex_field();
}


void FluidIO::Output_complex_field(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Output_complex_field_scalar(U, T);

	else if (global.program.kind == "RBC"  || global.program.kind == "STRATIFIED")
		Output_complex_field_RBC(U, T);
}

void FluidIO::Output_complex_field_scalar(FluidVF& U, FluidSF& T)
{
	// Set Output directory
	cout << "FluidIO::Output_complex_field_scalar:"
	     << "U.cvf.V1 = " << U.cvf.V1
	     << "U.cvf.V3 = " << U.cvf.V3
	     << endl;
	U.cvf.Write_complex_field();
	T.csf.Write_complex_field();
}


void FluidIO::Output_complex_field_RBC(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO") 
		Output_complex_field(U);			

	else if (global.PHYSICS.Pr_option == "PRINFTY") 
		T.csf.Write_complex_field();
		
	else
		Output_complex_field_scalar(U, T);
}


void FluidIO::Output_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	// Set Output directory
	U.cvf.Write_complex_field();
	T1.csf.Write_complex_field();
	T2.csf.Write_complex_field();
}

void FluidIO::Output_complex_field(FluidVF& U, FluidVF& W)
{
	// Set Output directory
	U.cvf.Write_complex_field();	
	W.cvf.Write_complex_field();
}


void FluidIO::Output_complex_field(FluidVF& U, FluidVF& W, FluidSF& T)
{
	// Set Output directory
	U.cvf.Write_complex_field();	
	W.cvf.Write_complex_field();
	T.csf.Write_complex_field();
}

//*********************************************************************************************

void FluidIO::Output_complex_field_frequent(FluidVF& U)
{
		// hdf5 format
		// Set Output directory
	BasicIO::Begin_frequent();
	U.cvf.Write_complex_field();
	BasicIO::End_frequent();
}


void FluidIO::Output_complex_field_frequent(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Output_complex_field_scalar_frequent(U, T);

	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_complex_field_RBC_frequent(U, T);

}

void FluidIO::Output_complex_field_scalar_frequent(FluidVF& U, FluidSF& T)
{
	BasicIO::Begin_frequent();
	U.cvf.Write_complex_field();
	T.csf.Write_complex_field();
	BasicIO::End_frequent();

}


void FluidIO::Output_complex_field_RBC_frequent(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO") 
		Output_complex_field_frequent(U);

	else if (global.PHYSICS.Pr_option == "PRINFTY"){
		BasicIO::Begin_frequent();
		T.csf.Write_complex_field();
		BasicIO::End_frequent();
	}
		
	else
		Output_complex_field_scalar_frequent(U, T);

}

void FluidIO::Output_complex_field_frequent(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	BasicIO::Begin_frequent();
	U.cvf.Write_complex_field();
	T1.csf.Write_complex_field();
	T2.csf.Write_complex_field();
	BasicIO::End_frequent();
	
}


void FluidIO::Output_complex_field_frequent(FluidVF& U, FluidVF& W)
{
	BasicIO::Begin_frequent();
	U.cvf.Write_complex_field();	
	W.cvf.Write_complex_field();
	BasicIO::End_frequent();
}


void FluidIO::Output_complex_field_frequent(FluidVF& U, FluidVF& W, FluidSF& T)
{
	BasicIO::Begin_frequent();
	U.cvf.Write_complex_field();	
	W.cvf.Write_complex_field();
	T.csf.Write_complex_field();
	BasicIO::End_frequent();
}


//*********************************************************************************************

void FluidIO::Output_reduced_complex_field(FluidVF& U)
{
	U.cvf.Write_reduced_complex_field();
}


void FluidIO::Output_reduced_complex_field(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Output_reduced_complex_field_scalar(U, T);

	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_reduced_complex_field_RBC(U, T);
}

void FluidIO::Output_reduced_complex_field_scalar(FluidVF& U, FluidSF& T)
{
	U.cvf.Write_reduced_complex_field();
	T.csf.Write_reduced_complex_field();
}


void FluidIO::Output_reduced_complex_field_RBC(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO") 
		Output_reduced_complex_field(U);			

	else if (global.PHYSICS.Pr_option == "PRINFTY") 
		T.csf.Write_reduced_complex_field();
		
	else
		Output_reduced_complex_field_scalar(U, T);
}

void FluidIO::Output_reduced_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
		// Set Output directory
	U.cvf.Write_reduced_complex_field();
	T1.csf.Write_reduced_complex_field();
	T2.csf.Write_reduced_complex_field();
}


void FluidIO::Output_reduced_complex_field(FluidVF& U, FluidVF& W)
{
	// Set Output directory
	U.cvf.Write_reduced_complex_field();	
	W.cvf.Write_reduced_complex_field();
}


void FluidIO::Output_reduced_complex_field(FluidVF& U, FluidVF& W, FluidSF& T)
{
	// Set Output directory
	U.cvf.Write_reduced_complex_field();	
	W.cvf.Write_reduced_complex_field();
	T.csf.Write_reduced_complex_field();
}


//*********************************************************************************************

void FluidIO::Output_real_field(FluidVF& U)
{
	if (global.time.now >= global.io.time.real_field_save_next) {
		U.rvf.Write_real_field();
		global.io.time.real_field_save_next += global.io.time.real_field_save_interval;
	}
}

void FluidIO::Output_real_field(FluidVF& U, FluidSF& T)
{
	if (global.time.now >= global.io.time.real_field_save_next) {
		U.rvf.Write_real_field();
		T.rsf.Write_real_field();
        
        global.io.time.real_field_save_next += global.io.time.real_field_save_interval;
	}
}

void FluidIO::Output_real_field(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	if (global.time.now >= global.io.time.real_field_save_next) {
		U.rvf.Write_real_field();
		T1.rsf.Write_real_field();
		T2.rsf.Write_real_field();
        
        global.io.time.real_field_save_next += global.io.time.real_field_save_interval;
	}
}

void FluidIO::Output_real_field(FluidVF& U, FluidVF& W)
{
	if (global.time.now >= global.io.time.real_field_save_next) {
		U.rvf.Write_real_field();	
		W.rvf.Write_real_field();
        
        global.io.time.real_field_save_next += global.io.time.real_field_save_interval;
	}
}


void FluidIO::Output_real_field(FluidVF& U, FluidVF& W, FluidSF& T)
{
	if (global.time.now >= global.io.time.real_field_save_next) {
		U.rvf.Write_real_field();	
		W.rvf.Write_real_field();
		T.rsf.Write_real_field();
        
        global.io.time.real_field_save_next += global.io.time.real_field_save_interval;
	}
}


/**********************************************************************************************

							FluidIO::Output_field_k()

***********************************************************************************************/

bool FluidIO::Spectral_space_buffer_full()
{
	if ((global.io.probes.spectral_space.buffer_index+global.io.probes.spectral_space.packet_size) >  (global.io.probes.spectral_space.buffer_size))
		return true;
    
	else
		return false;
}


//*********************************************************************************************
void FluidIO::Output_field_k(FluidVF& U)
{
	if (global.time.now < global.io.time.field_k_save_next)
		return;
		
    TinyVector<DP,3> Vk_real;
    TinyVector<complx,3> Vk_complex;
	int kx, ky, kz;
	
    bool probe_in_me_flag;
    
	for (int probe=0; probe<global.io.probes.spectral_space.number; probe++) {
		kx = global.io.probes.spectral_space.coords(probe,1);
		ky = global.io.probes.spectral_space.coords(probe,2);
		kz = global.io.probes.spectral_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me(kx, ky, kz);
        
		if (probe_in_me_flag) {
            if (Spectral_space_buffer_full()) {
				if (!field_k_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_k_out_"+To_string(my_id)+".d";
					field_k_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_k_out_file, global.io.probes.spectral_space.field_buffer, global.io.probes.spectral_space.buffer_index, global.io.probes.spectral_space.packet_size);
            }
    
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = global.time.now;
            
            for (int k=1; k<=3; k++) 
                global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = ((DP) global.io.probes.spectral_space.coords(probe,k));
            
			// Get the field now..
            if (basis_type == "SSS") {
                Vk_real = real(universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
                for (int k=1; k<=3; k++) 
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Vk_real(k-1);
            }
            else {
                Vk_complex = universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
                for (int k=1; k<=3; k++) {
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Vk_complex(k-1));
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Vk_complex(k-1));
                }
            }
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = U.Get_Tk(kx,ky,kz);            
        }
    }
	
	global.io.time.field_k_save_next += global.io.time.field_k_save_interval;
}


//*********************************************************************************************
// SCALAR
//

void FluidIO::Output_field_k(FluidVF& U, FluidSF& T)
{
	if (global.time.now < global.io.time.field_k_save_next)
		return;
	
	TinyVector<DP,3> Vk_real;
    TinyVector<complx,3> Vk_complex;
	DP Fk_real;
	complx Fk_complex;
	int kx, ky, kz;
	
    bool probe_in_me_flag;
    
	for (int probe=0; probe<global.io.probes.spectral_space.number; probe++) {
		kx = global.io.probes.spectral_space.coords(probe,1);
		ky = global.io.probes.spectral_space.coords(probe,2);
		kz = global.io.probes.spectral_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me(kx, ky, kz);
        
		if (probe_in_me_flag) {
            if (Spectral_space_buffer_full()) {
				if (!field_k_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_k_out_"+To_string(my_id)+".d";
					field_k_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_k_out_file, global.io.probes.spectral_space.field_buffer, global.io.probes.spectral_space.buffer_index, global.io.probes.spectral_space.packet_size);
            }
			
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = global.time.now;
            
            for (int k=1; k<=3; k++) 
                global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = ((DP) global.io.probes.spectral_space.coords(probe,k));
            
			// Get the field now..
            if (basis_type == "SSS") {
                Vk_real = real(universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
				Fk_real = real(universal->Get_spectral_field(kx, ky, kz, T.csf.F));
                for (int k=1; k<=3; k++) 
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Vk_real(k-1);
				 global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Fk_real;
            }
            else {
                Vk_complex = universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
				Fk_complex = universal->Get_spectral_field(kx, ky, kz, T.csf.F);
                for (int k=1; k<=3; k++) {
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Vk_complex(k-1));
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Vk_complex(k-1));
                }
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Fk_complex);
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Fk_complex);
            }
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = U.Get_Tk(kx,ky,kz);  
			global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = T.Get_Tk(kx,ky,kz);  
        }
    }
	
	global.io.time.field_k_save_next += global.io.time.field_k_save_interval;
}


//*********************************************************************************************
// SCALAR
//

void FluidIO::Output_field_k(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	if (global.time.now < global.io.time.field_k_save_next)
		return;
	
	TinyVector<DP,3> Vk_real;
    TinyVector<complx,3> Vk_complex;
	DP F1k_real, F2k_real;
	complx F1k_complex, F2k_complex;
	int kx, ky, kz;
	
    bool probe_in_me_flag;
    
	for (int probe=0; probe<global.io.probes.spectral_space.number; probe++) {
		kx = global.io.probes.spectral_space.coords(probe,1);
		ky = global.io.probes.spectral_space.coords(probe,2);
		kz = global.io.probes.spectral_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me(kx, ky, kz);
        
		if (probe_in_me_flag) {
            if (Spectral_space_buffer_full()) {
				if (!field_k_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_k_out_"+To_string(my_id)+".d";
					field_k_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_k_out_file, global.io.probes.spectral_space.field_buffer, global.io.probes.spectral_space.buffer_index, global.io.probes.spectral_space.packet_size);
            }
			
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = global.time.now;
            
            for (int k=1; k<=3; k++) 
                global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = ((DP) global.io.probes.spectral_space.coords(probe,k));
            
				// Get the field now..
            if (basis_type == "SSS") {
                Vk_real = real(universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
				F1k_real = real(universal->Get_spectral_field(kx, ky, kz, T1.csf.F));
				F2k_real = real(universal->Get_spectral_field(kx, ky, kz, T2.csf.F));
                for (int k=1; k<=3; k++) 
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Vk_real(k-1);
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = F1k_real;
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = F2k_real;
            }
			
            else {
                Vk_complex = universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
				F1k_complex = universal->Get_spectral_field(kx, ky, kz, T1.csf.F);
				F2k_complex = universal->Get_spectral_field(kx, ky, kz, T2.csf.F);
                for (int k=1; k<=3; k++) {
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Vk_complex(k-1));
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Vk_complex(k-1));
                }
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(F1k_complex);
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(F1k_complex);
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(F2k_complex);
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(F2k_complex);
            }
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = U.Get_Tk(kx,ky,kz);  
			global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = T1.Get_Tk(kx,ky,kz); 
			global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = T2.Get_Tk(kx,ky,kz); 
        }
    }
	
	global.io.time.field_k_save_next += global.io.time.field_k_save_interval;
}




//*********************************************************************************************
//		MHD
//

void FluidIO::Output_field_k(FluidVF& U, FluidVF& W)
{
	if (global.time.now < global.io.time.field_k_save_next)
		return;
	
	TinyVector<DP,3> Vk_real, Wk_real;
    TinyVector<complx,3> Vk_complex, Wk_complex;
	DP Fk_real;
	complx Fk_complex;
	int kx, ky, kz;
	
    bool probe_in_me_flag;
    
	for (int probe=0; probe<global.io.probes.spectral_space.number; probe++) {
		kx = global.io.probes.spectral_space.coords(probe,1);
		ky = global.io.probes.spectral_space.coords(probe,2);
		kz = global.io.probes.spectral_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me(kx, ky, kz);
        
		if (probe_in_me_flag) {
            if (Spectral_space_buffer_full()) {
				if (!field_k_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_k_out_"+To_string(my_id)+".d";
					field_k_out_file.open(filename.c_str());
				}
				
                global.io.Dump_buffer(field_k_out_file, global.io.probes.spectral_space.field_buffer, global.io.probes.spectral_space.buffer_index, global.io.probes.spectral_space.packet_size);
            }
			
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = global.time.now;
            
            for (int k=1; k<=3; k++) 
                global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = ((DP) global.io.probes.spectral_space.coords(probe,k));
            
				// Get the field now..
            if (basis_type == "SSS") {
                Vk_real = real(universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
				Wk_real = real(universal->Get_spectral_field(kx, ky, kz, W.cvf.V1, W.cvf.V2, W.cvf.V3));
                for (int k=1; k<=3; k++) 
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Vk_real(k-1);
				for (int k=1; k<=3; k++) 
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Wk_real(k-1);
            }
            else {
                Vk_complex = universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
				Wk_complex = universal->Get_spectral_field(kx, ky, kz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
                for (int k=1; k<=3; k++) {
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Vk_complex(k-1));
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Vk_complex(k-1));
                }
				for (int k=1; k<=3; k++) {
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Wk_complex(k-1));
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Wk_complex(k-1));
                }
            }
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = U.Get_Tk(kx,ky,kz); 
			global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = W.Get_Tk(kx,ky,kz);  
        }
    }
	
	global.io.time.field_k_save_next += global.io.time.field_k_save_interval;
}



//*********************************************************************************************
// VF + Scalar
//

void FluidIO::Output_field_k(FluidVF& U, FluidVF& W, FluidSF& T)
{
	if (global.time.now < global.io.time.field_k_save_next)
		return;
	
	TinyVector<DP,3> Vk_real, Wk_real;
    TinyVector<complx,3> Vk_complex, Wk_complex;
	DP Fk_real;
	complx Fk_complex;
	
	int kx, ky, kz;
	
    bool probe_in_me_flag;
    
	for (int probe=0; probe<global.io.probes.spectral_space.number; probe++) {
		kx = global.io.probes.spectral_space.coords(probe,1);
		ky = global.io.probes.spectral_space.coords(probe,2);
		kz = global.io.probes.spectral_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me(kx, ky, kz);
        
		if (probe_in_me_flag) {
            if (Spectral_space_buffer_full()) {
				if (!field_k_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_k_out_"+To_string(my_id)+".d";
					field_k_out_file.open(filename.c_str());
				}
				
                global.io.Dump_buffer(field_k_out_file, global.io.probes.spectral_space.field_buffer, global.io.probes.spectral_space.buffer_index, global.io.probes.spectral_space.packet_size);
            }
			
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = global.time.now;
            
            for (int k=1; k<=3; k++) 
                global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = ((DP) global.io.probes.spectral_space.coords(probe,k));
            
				// Get the field now..
            if (basis_type == "SSS") {
                Vk_real = real(universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
				Wk_real = real(universal->Get_spectral_field(kx, ky, kz, W.cvf.V1, W.cvf.V2, W.cvf.V3));
				Fk_real = real(universal->Get_spectral_field(kx, ky, kz, T.csf.F));
                for (int k=1; k<=3; k++) 
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Vk_real(k-1);
				for (int k=1; k<=3; k++) 
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Wk_real(k-1);
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = Fk_real;
            }
            else {
                Vk_complex = universal->Get_spectral_field(kx, ky, kz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
				Wk_complex = universal->Get_spectral_field(kx, ky, kz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
				Fk_complex = universal->Get_spectral_field(kx, ky, kz, T.csf.F);
                for (int k=1; k<=3; k++) {
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Vk_complex(k-1));
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Vk_complex(k-1));
                }
				for (int k=1; k<=3; k++) {
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Wk_complex(k-1));
                    global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Wk_complex(k-1));
                }
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = real(Fk_complex);
				global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = imag(Fk_complex);
            }
            global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = U.Get_Tk(kx,ky,kz); 
			global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = W.Get_Tk(kx,ky,kz); 
			global.io.probes.spectral_space.field_buffer(global.io.probes.spectral_space.buffer_index++) = T.Get_Tk(kx,ky,kz); 
        }
    }
	
	global.io.time.field_k_save_next += global.io.time.field_k_save_interval;
}


/**********************************************************************************************

							FluidIO::Output_field_r()

***********************************************************************************************/


bool FluidIO::Real_space_buffer_full()
{
	if ((global.io.probes.real_space.buffer_index+global.io.probes.real_space.packet_size) >  (global.io.probes.real_space.buffer_size))
		return true;
	
	else
		return false;
}

//*********************************************************************************************

void FluidIO::Output_field_r(FluidVF& U)
{
	if (global.time.now < global.io.time.field_r_save_next)
		return;
	
	TinyVector<DP,3> Vr;
	int rx, ry, rz;
	
	bool probe_in_me_flag;
	
	for (int probe=0; probe<global.io.probes.real_space.number; probe++) {
		rx = global.io.probes.real_space.coords(probe,1);
		ry = global.io.probes.real_space.coords(probe,2);
		rz = global.io.probes.real_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me_real_space(rx, ry, rz);
		
		if (probe_in_me_flag) {
			if (Real_space_buffer_full()) {
				if (!field_r_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_r_out_"+To_string(my_id)+".d";
					field_r_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_r_out_file, global.io.probes.real_space.field_buffer, global.io.probes.real_space.buffer_index, global.io.probes.real_space.packet_size);
			}
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = global.time.now;
            
            for (int r=1; r<=3; r++) 
                global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = ((DP) global.io.probes.real_space.coords(probe,r));
            
				// Get the field now..
                Vr = universal->Get_real_field(rx, ry, rz, U.rvf.V1r, U.rvf.V2r, U.rvf.V3r);
                for (int r=1; r<=3; r++) 
                    global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Vr(r-1);
		}
	}
	
	global.io.time.field_r_save_next += global.io.time.field_r_save_interval;
}

//*********************************************************************************************

void FluidIO::Output_field_r(FluidVF& U, FluidSF& T)
{
	if (global.time.now < global.io.time.field_r_save_next)
		return;
	
	TinyVector<DP,3> Vr;
	DP Fr;
	
	int rx, ry, rz;
	bool probe_in_me_flag;
	
	for (int probe=0; probe<global.io.probes.real_space.number; probe++) {
		rx = global.io.probes.real_space.coords(probe,1);
		ry = global.io.probes.real_space.coords(probe,2);
		rz = global.io.probes.real_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me_real_space(rx, ry, rz);
		
		if (probe_in_me_flag) {
			if (Real_space_buffer_full()) {
				if (!field_r_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_r_out_"+To_string(my_id)+".d";
					field_r_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_r_out_file, global.io.probes.real_space.field_buffer, global.io.probes.real_space.buffer_index, global.io.probes.real_space.packet_size);
			}
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = global.time.now;
            
            for (int r=1; r<=3; r++) 
                global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = ((DP) global.io.probes.real_space.coords(probe,r));
            
				// Get the field now..
			Vr = universal->Get_real_field(rx, ry, rz, U.rvf.V1r, U.rvf.V2r, U.rvf.V3r);
			Fr = universal->Get_real_field(rx, ry, rz, T.rsf.Fr);
			for (int r=1; r<=3; r++) 
				global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Vr(r-1);
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Fr;
		}
	}
	
	global.io.time.field_r_save_next += global.io.time.field_r_save_interval;
}


//*********************************************************************************************

void FluidIO::Output_field_r(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	if (global.time.now < global.io.time.field_r_save_next)
		return;
	
	TinyVector<DP,3> Vr;
	DP F1r, F2r;
	
	int rx, ry, rz;
	bool probe_in_me_flag;
	
	for (int probe=0; probe<global.io.probes.real_space.number; probe++) {
		rx = global.io.probes.real_space.coords(probe,1);
		ry = global.io.probes.real_space.coords(probe,2);
		rz = global.io.probes.real_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me_real_space(rx, ry, rz);
		
		if (probe_in_me_flag) {
			if (Real_space_buffer_full()) {
				if (!field_r_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_r_out_"+To_string(my_id)+".d";
					field_r_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_r_out_file, global.io.probes.real_space.field_buffer, global.io.probes.real_space.buffer_index, global.io.probes.real_space.packet_size);
			}
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = global.time.now;
            
            for (int r=1; r<=3; r++) 
                global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = ((DP) global.io.probes.real_space.coords(probe,r));
            
				// Get the field now..
			Vr = universal->Get_real_field(rx, ry, rz, U.rvf.V1r, U.rvf.V2r, U.rvf.V3r);
			F1r = universal->Get_real_field(rx, ry, rz, T1.rsf.Fr);
			F2r = universal->Get_real_field(rx, ry, rz, T2.rsf.Fr);
			for (int r=1; r<=3; r++) 
				global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Vr(r-1);
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = F1r;
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = F2r;
		}
	}
	
	global.io.time.field_r_save_next += global.io.time.field_r_save_interval;
}

//*********************************************************************************************


void FluidIO::Output_field_r(FluidVF& U, FluidVF& W)
{
	if (global.time.now < global.io.time.field_r_save_next)
		return;
	
	TinyVector<DP,3> Vr, Wr;
	bool probe_in_me_flag;
	int rx, ry, rz;
	
	for (int probe=0; probe<global.io.probes.real_space.number; probe++) {
		rx = global.io.probes.real_space.coords(probe,1);
		ry = global.io.probes.real_space.coords(probe,2);
		rz = global.io.probes.real_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me_real_space(rx, ry, rz);
		
		if (probe_in_me_flag) {
			if (Real_space_buffer_full()) {
				if (!field_r_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_r_out_"+To_string(my_id)+".d";
					field_r_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_r_out_file, global.io.probes.real_space.field_buffer, global.io.probes.real_space.buffer_index, global.io.probes.real_space.packet_size);
			}
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = global.time.now;
            
            for (int r=1; r<=3; r++) 
                global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = ((DP) global.io.probes.real_space.coords(probe,r));
            
				// Get the field now..
			Vr = universal->Get_real_field(rx, ry, rz, U.rvf.V1r, U.rvf.V2r, U.rvf.V3r);
			Wr = universal->Get_real_field(rx, ry, rz, W.rvf.V1r, W.rvf.V2r, W.rvf.V3r);
			for (int r=1; r<=3; r++) 
				global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Vr(r-1);
			
			for (int r=1; r<=3; r++) 
				global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Wr(r-1);
		}
	}
	
	global.io.time.field_r_save_next += global.io.time.field_r_save_interval;
}

//*********************************************************************************************


void FluidIO::Output_field_r(FluidVF& U, FluidVF& W, FluidSF& T)
{
	if (global.time.now < global.io.time.field_r_save_next)
		return;
	
	TinyVector<DP,3> Vr, Wr;
	DP Fr;
	int rx, ry, rz;
	bool probe_in_me_flag;
	
	for (int probe=0; probe<global.io.probes.real_space.number; probe++) {
		rx = global.io.probes.real_space.coords(probe,1);
		ry = global.io.probes.real_space.coords(probe,2);
		rz = global.io.probes.real_space.coords(probe,3);
		
        probe_in_me_flag = universal->Probe_in_me_real_space(rx, ry, rz);
		
		if (probe_in_me_flag) {
			if (Real_space_buffer_full()) {
				if (!field_r_out_file.is_open()){
					string filename = global.io.data_dir + "/out/field_r_out_"+To_string(my_id)+".d";
					field_r_out_file.open(filename.c_str());
				}
				
				global.io.Dump_buffer(field_r_out_file, global.io.probes.real_space.field_buffer, global.io.probes.real_space.buffer_index, global.io.probes.real_space.packet_size);
			}
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = global.time.now;
            
            for (int r=1; r<=3; r++) 
                global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = ((DP) global.io.probes.real_space.coords(probe,r));
            
				// Get the field now..
			Vr = universal->Get_real_field(rx, ry, rz, U.rvf.V1r, U.rvf.V2r, U.rvf.V3r);
			Wr = universal->Get_real_field(rx, ry, rz, W.rvf.V1r, W.rvf.V2r, W.rvf.V3r);
			Fr = universal->Get_real_field(rx, ry, rz, T.rsf.Fr);
			for (int r=1; r<=3; r++) 
				global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Vr(r-1);
			
			for (int r=1; r<=3; r++) 
				global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Wr(r-1);
			
			global.io.probes.real_space.field_buffer(global.io.probes.real_space.buffer_index++) = Fr;
		}
	}
	
	global.io.time.field_r_save_next += global.io.time.field_r_save_interval;
}


//******************************  End of output_field.cc  **************************************

