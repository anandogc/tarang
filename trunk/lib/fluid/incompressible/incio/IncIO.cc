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
 
#include "IncIO.h"

extern Global global;

//*********************************************************************************************

FluidIO_incompress::~FluidIO_incompress()
{
	if (my_id == master_id)  {
		pressure_file.close();
		pressure_spectrum_file.close();
		flux_file.close();		
		shell_to_shell_file.close();
		ring_to_ring_file.close();
		cylindrical_ring_to_ring_file.close();
	}
}

//*********************************************************************************************

void FluidIO_incompress::Init_energy_transfer()
{
	if(global.energy_transfer.turnon)
		energyTr = new EnergyTr;
}


void FluidIO_incompress::Open_files()
{
	Open_base_files();
	
	if (my_id == master_id) {
		
		string filename;
		
		filename = "/out/pressure.d";
		filename = global.io.data_dir+ filename;	pressure_file.open(filename.c_str());
		
		filename = "/out/pressure_spectrum.d";
		filename = global.io.data_dir+ filename; pressure_spectrum_file.open(filename.c_str());
		
		filename = "/out/flux.d";
		filename = global.io.data_dir+ filename;   flux_file.open(filename.c_str());
		
		filename = "/out/shell_to_shell.d";
		filename = global.io.data_dir+ filename;   shell_to_shell_file.open(filename.c_str());
		
		filename = "/out/ring_to_ring.d";
		filename = global.io.data_dir+ filename;   ring_to_ring_file.open(filename.c_str());
		
		filename = "/out/cylindrical_ring_to_ring.d";
		filename = global.io.data_dir+ filename;   cylindrical_ring_to_ring_file.open(filename.c_str());
		
		// Set precision digit
		flux_file.setf(ios::scientific);
		flux_file.precision(global.io.output_precision); // sets number of decimal places
		
		shell_to_shell_file.setf(ios::scientific);
		shell_to_shell_file.precision(global.io.output_precision);
		
		ring_to_ring_file.setf(ios::scientific);
		ring_to_ring_file.precision(global.io.output_precision);
		
		cylindrical_ring_to_ring_file.setf(ios::scientific);
		cylindrical_ring_to_ring_file.precision(global.io.output_precision); 
	}
}


void FluidIO_incompress::Close_files()
{
    Close_base_files();
    
    
    pressure_file.close();
    pressure_spectrum_file.close();
    flux_file.close();
    shell_to_shell_file.close();
    ring_to_ring_file.close();
    cylindrical_ring_to_ring_file.close();
    
    
}
//************************ Class declaration of IO ends  **************************************
	

 
