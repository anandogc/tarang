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

/*! \file  Global.h
 * 
 * @brief  Class constructor of Global.
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "Global.h"



//*********************************************************************************************

void Global::Init_defaults()
{
		// program
	program.kind = "FLUID_INCOMPRESS";
	program.iter_or_diag = "ITERATION";
	program.alias_option = "DEALIAS";
	program.integration_scheme = "RK4";
	program.basis_type = "FFF";
	program.LES_switch = false;
	program.apply_strong_realitycond_alltime_switch = false;
    program.apply_weak_realitycond_alltime_switch = true;
	program.low_dimensional_switch = false;
	program.two_and_half_dimension = false;
	program.two_dimension = false;
	program.dt_option = 0;
	program.sincostr_switch = "FFF";
	
	// field
	field.incompressible = true;
	field.waveno_switch = true;
	field.anisotropy_dirn = 3;
	
	field.N[0] = 0;
	field.N[1] = 16;
	field.N[2] = 16;
	field.N[3] = 14;


	field.kfactor.resize(4);
	field.kfactor[0] = 0;
	field.kfactor[1] = 1;
	field.kfactor[2] = 1;
	field.kfactor[3] = 1;
	

	field.diss_coefficients.resize(1);
	field.diss_coefficients[0]=0.1;
	field.hyper_diss_coefficients.resize(1);
	field.hyper_diss_exponents.resize(1);
	
	// time
	time.init = 0;
	time.final = 0.010;
	time.dt_fixed = 0.001;

	// force
	force.U_switch = true;
	force.W_switch = false;
	force.T_switch = false;
	
	force.field_procedure = 0;
	
		//force.int_para;
		//force.double_para;
		//force.string_para;
	
	force.modes.number = 0;
	force.modes.number_components = no_components_table[program.kind];
	force.modes.coords.resize(force.modes.number,4);
	
	
	// IO
	io.data_dir = ".";
	io.input_field_procedure = 4;
	io.input_vx_vy_switch = false;
	io.output_vx_vy_switch = false;
	

	
	//para["io"]["int_para"] >> io.int_para;
	//para["io"]["double_para"] >> io.double_para;
	//para["io"]["string_para"] >> io.string_para;
	
	io.real_space_field_available = false;
	
	io.time.global_save_next = myconstant.INF_TIME;
	io.time.complex_field_save_next = myconstant.INF_TIME;
	io.time.field_frequent_save_next = myconstant.INF_TIME;
	io.time.field_reduced_save_next = myconstant.INF_TIME;
	io.time.real_field_save_next = myconstant.INF_TIME;
	io.time.field_k_save_next = myconstant.INF_TIME;
	io.time.field_r_save_next = myconstant.INF_TIME;
	io.time.spectrum_save_next = myconstant.INF_TIME;
	io.time.pressure_save_next = myconstant.INF_TIME;
	io.time.pressure_spectrum_save_next = myconstant.INF_TIME;
	io.time.flux_save_next = myconstant.INF_TIME;
	io.time.shell_to_shell_save_next = myconstant.INF_TIME;
	io.time.ring_spectrum_save_next = myconstant.INF_TIME;
	io.time.ring_to_ring_save_next = myconstant.INF_TIME;
	io.time.cylindrical_ring_spectrum_save_next = myconstant.INF_TIME;
	io.time.cylindrical_ring_to_ring_save_next = myconstant.INF_TIME;
	io.time.structure_fn_save_next = myconstant.INF_TIME;
	io.time.cout_save_next = myconstant.INF_TIME;
	
	io.time.global_save_interval = myconstant.INF_TIME;
	io.time.complex_field_save_interval = myconstant.INF_TIME;
	io.time.field_frequent_save_interval = myconstant.INF_TIME;
	io.time.field_reduced_save_interval = myconstant.INF_TIME;
	io.time.real_field_save_interval = myconstant.INF_TIME;
	io.time.field_k_save_interval = myconstant.INF_TIME;
	io.time.field_r_save_interval = myconstant.INF_TIME;
	io.time.spectrum_save_interval = myconstant.INF_TIME;
	io.time.pressure_save_interval = myconstant.INF_TIME;
	io.time.pressure_spectrum_save_interval = myconstant.INF_TIME;
	io.time.flux_save_interval = myconstant.INF_TIME;
	io.time.shell_to_shell_save_interval = myconstant.INF_TIME;
	io.time.ring_spectrum_save_interval = myconstant.INF_TIME;
	io.time.ring_to_ring_save_interval = myconstant.INF_TIME;
	io.time.cylindrical_ring_spectrum_save_interval = myconstant.INF_TIME;
	io.time.cylindrical_ring_to_ring_save_interval = myconstant.INF_TIME;
	io.time.structure_fn_save_interval = myconstant.INF_TIME;
	io.time.cout_save_interval = myconstant.INF_TIME;

	// SPECTRUM
	spectrum.shell.turnon = true;
	spectrum.ring.turnon = false;
	
	spectrum.cylindrical_ring.turnon = false;

	
	// EnergyTr
	energy_transfer.turnon = false;
	

			// SHELL-to-SHELL
	energy_transfer.shell_to_shell.turnon = false;
		

		
		// ring-to-ring
	energy_transfer.ring_to_ring.turnon = false;

		
		// Cylindrical ring to ring
	energy_transfer.cylindrical_ring_to_ring.turnon = false;
}









