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


/*! \file  IO.h
 * 
 *	@brief  Class declaration of input/output function.
 * 
 *
 *	@author  M. K. Verma
 *	@version 2.0 MPI
 *	@date feb 2012
 *
 * @bug  No known bugs
 */
 
#include "FluidIO.h"

//*********************************************************************************************

FluidIO::FluidIO() {
	slice_file_counter = 0;
}

void FluidIO::Open_base_files()
{
	string filename;
	
	if (master) {
        
		filename = "/in/field_in.d";
		filename = global.io.data_dir+ filename;	field_in_file.open(filename.c_str());
		
		filename = "/in/force_field_in.d";
		filename = global.io.data_dir+ filename; force_field_in_file.open(filename.c_str());
		
		filename = "/out/field_out.d";
		filename = global.io.data_dir+ filename;   field_out_file.open(filename.c_str());
		
		filename = "/out/field_out_reduced.d";
		filename = global.io.data_dir+ filename;   field_out_reduced_file.open(filename.c_str());
		
		filename = "/out/realfield_out.d";
		filename = global.io.data_dir+ filename;   realfield_out_file.open(filename.c_str());
		
		filename = "/out/glob.d";
		filename = global.io.data_dir+ filename;   global_file.open(filename.c_str());
		if (!global_file.is_open()) 
			cout << "UNABLE TO OPEN FILE global_file (glob.d) " << endl;
        
        filename = "/out/glob_real.d";
		filename = global.io.data_dir+ filename;   global_real_space_file.open(filename.c_str());
		
		filename = "/out/spectrum.d";
		filename = global.io.data_dir+ filename;   spectrum_file.open(filename.c_str());
		
		filename = "/out/ring_spectrum.d";
		filename = global.io.data_dir+ filename;   ring_spectrum_file.open(filename.c_str());
		
		filename = "/out/cyl_ring_spectrum.d";
		filename = global.io.data_dir+ filename;   cylindrical_ring_spectrum_file.open(filename.c_str());
		
		filename = "/out/struct_fn.d";
		filename = global.io.data_dir+ filename;   structure_fn_file.open(filename.c_str());
		
		filename = "/out/planar_struct_fn.d";
		filename = global.io.data_dir+ filename;   planar_structure_fn_file.open(filename.c_str());
		
		filename = "/out/pressure.d";
		filename = global.io.data_dir+ filename;	  pressure_file.open(filename.c_str());
        
        filename = "/out/Tk_spectrum.d";
		filename = global.io.data_dir+ filename;	  Tk_spectrum_file.open(filename.c_str());
        
        filename = "/out/Tk_ring_spectrum.d";
		filename = global.io.data_dir+ filename;	  Tk_ring_spectrum_file.open(filename.c_str());
        
        filename = "/out/Tk_cylindrical_ring_spectrum.d";
		filename = global.io.data_dir+ filename;	  Tk_cylindrical_ring_spectrum_file.open(filename.c_str());
        
        
        filename = "/out/misc.d";
		filename = global.io.data_dir+ filename;	  misc_file.open(filename.c_str());
        
        
		
		// Set precision digit
		global_file.setf(ios::scientific);
		global_file.precision(global.io.output_precision); // sets number of decimal places
		
		field_out_file.setf(ios::scientific);
		field_out_file.precision(global.io.output_precision);
		
		pressure_file.setf(ios::scientific);
		pressure_file.precision(global.io.output_precision); 
		
		cout.setf(ios::fixed);
		cout.precision(global.io.output_precision); 
	}
		
	field_k_out_file.setf(ios::scientific);
	field_k_out_file.precision(global.io.output_precision);
	
	field_r_out_file.setf(ios::scientific);
	field_r_out_file.precision(global.io.output_precision); 
}

//*********************************************************************************************

void FluidIO::Close_base_files()
{
	
    // Output the remaining buffer from global_data and probes.
	if (master) 
		if (global.io.global_data.buffer_index > 0) 
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
	
	if (field_k_out_file.is_open() && global.io.probes.spectral_space.buffer_index > 0)
        global.io.Dump_buffer(field_k_out_file, global.io.probes.spectral_space.field_buffer, global.io.probes.spectral_space.buffer_index, global.io.probes.spectral_space.packet_size);
	
	if (field_r_out_file.is_open() && global.io.probes.real_space.buffer_index > 0)
		global.io.Dump_buffer(field_r_out_file, global.io.probes.real_space.field_buffer, global.io.probes.real_space.buffer_index, global.io.probes.real_space.packet_size);
    
	if (master)  {
		field_in_file.close();
		force_field_in_file.close();
		
		field_out_file.close();
		field_frequent_out_file.close();
		field_out_reduced_file.close();
		realfield_out_file.close();
		
		global_file.close();
        global_real_space_file.close();
		spectrum_file.close();
		ring_spectrum_file.close();
		cylindrical_ring_spectrum_file.close();
		
		structure_fn_file.close();
		planar_structure_fn_file.close();
		
		pressure_file.close();
        Tk_spectrum_file.close();
        Tk_ring_spectrum_file.close();
        Tk_cylindrical_ring_spectrum_file.close();
        
        misc_file.close();
	}
	
	if (field_k_out_file.is_open())
		field_k_out_file.close();
	
	if (field_r_out_file.is_open())
		field_r_out_file.close();
}

//************************ Class declaration of IO ends  **************************************
	

 
