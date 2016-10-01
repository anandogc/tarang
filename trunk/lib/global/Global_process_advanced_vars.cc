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

#include "basic_extern_vars.h"
#include "Global.h"
#include "Global_extern_vars.h"

//****************************************************************************************

void Global::Process_advanced_vars()
{
	mpi.MPI_COMM_ROW = fftk.Get_communicator("ROW");
	mpi.MPI_COMM_COL = fftk.Get_communicator("COL");

	// Time
	time.now = time.init;
	
	// force	
	force.configuration_done = false;   // default
	force.modes.read_done = false;
	
	
	//IO Output starts here
	io.global_data.buffer_index = 0;
	io.global_data.packet_size = global_data_packet_size_table[program.kind];		
	
	io.global_data.buffer_size = io.global_data.packet_size * myconstant.MAX_NO_GLOB_BUFFER_PACKETS;
	io.global_data.buffer.resize(io.global_data.buffer_size);
	
	io.global_data.buffer_size = io.global_data.packet_size * myconstant.MAX_NO_GLOB_BUFFER_PACKETS;
	io.global_data.buffer.resize(io.global_data.buffer_size);
	
	io.probes.spectral_space.buffer_index = 0;
	
	io.probes.spectral_space.packet_size = spectral_probe_packet_size_table[program.kind][program.basis_type];
	
	io.probes.spectral_space.buffer_size = io.probes.spectral_space.packet_size * myconstant.MAX_NO_PROBE_PACKETS;
	
	io.probes.spectral_space.field_buffer.resize(io.probes.spectral_space.packet_size*myconstant.MAX_NO_PROBE_PACKETS);
	
	
	io.probes.real_space.buffer_index = 0;
	io.probes.real_space.packet_size = real_probe_packet_size_table[program.kind];
	io.probes.real_space.buffer_size = io.probes.real_space.packet_size * myconstant.MAX_NO_PROBE_PACKETS;
	io.probes.real_space.field_buffer.resize(io.probes.real_space.packet_size*myconstant.MAX_NO_PROBE_PACKETS);

	// SPECTRUM
	//## commented during testing
	
	
	Real max_theta = universal->Get_max_polar_angle();
	Real dtheta;
	Real Kpll_min = 0.0;
	Real Kpll_max = universal->Anis_max_Kpll();
	Real dKpll;

	int min_radius_outside = universal->Min_radius_outside();
	int max_radius_inside =  universal->Max_radius_inside();
	int max_cylinder_radius_inside = universal->Anis_max_Krho_radius_inside();

	spectrum.shell.no_shells = min_radius_outside+1;
	
	if (spectrum.ring.turnon) {
		spectrum.ring.no_shells = max_radius_inside+1;
		
		if (spectrum.ring.sector_option == "EQUISPACED") {
			spectrum.ring.sector_angles.resize(spectrum.ring.no_sectors+1);
			dtheta = max_theta/spectrum.ring.no_sectors;
			for (int i=0; i<spectrum.ring.no_sectors; i++)
				spectrum.ring.sector_angles(i) = i*dtheta;
			
			spectrum.ring.sector_angles(spectrum.ring.no_sectors) = max_theta;
		}
		
		else if (spectrum.ring.sector_option == "EQUAL_NO_MODES") {
			;   // to implement
		}
	}
	
	if (spectrum.cylindrical_ring.turnon) {
		spectrum.cylindrical_ring.no_shells = max_cylinder_radius_inside+1;
		
		spectrum.cylindrical_ring.kpll_array.resize(spectrum.cylindrical_ring.no_slabs+1);
		
		if (spectrum.cylindrical_ring.kpll_option == "EQUISPACED") {
			spectrum.cylindrical_ring.kpll_array.resize(spectrum.cylindrical_ring.no_slabs+1);
			dKpll = (Kpll_max - Kpll_min)/spectrum.cylindrical_ring.no_slabs;
			spectrum.cylindrical_ring.kpll_array(0) = Kpll_min;
			
			for (int i=1; i<=spectrum.cylindrical_ring.no_slabs; i++) 
				spectrum.cylindrical_ring.kpll_array(i) = Kpll_min + i*dKpll;
		}
	}

	// EnergyTransfers
	
	// if the radii.size() = 0, then USER_DEFINED is off.
	if (energy_transfer.turnon) {
				
		if (energy_transfer.flux.turnon) {

			if ((energy_transfer.flux.no_spheres > 0) && (energy_transfer.flux.radii.size()==0)){
				energy_transfer.flux.radii.resize(energy_transfer.flux.no_spheres+1);
				
				energy_transfer.flux.radii(0) = 0.0;
				energy_transfer.flux.radii(1) = 2.0;
				energy_transfer.flux.radii(2) = 4.0;
				energy_transfer.flux.radii(3) = 8.0;
				
				energy_transfer.flux.radii(energy_transfer.flux.no_spheres-2) = max_radius_inside/2.0;
				energy_transfer.flux.radii(energy_transfer.flux.no_spheres-1) = max_radius_inside;
				energy_transfer.flux.radii(energy_transfer.flux.no_spheres) = INF_RADIUS;
				
				if (energy_transfer.flux.no_spheres > 6) {
					Real s = log2(max_radius_inside/16.0) / (energy_transfer.flux.no_spheres - 5);
					for (int i=4; i<=(energy_transfer.flux.no_spheres-3); i++) 
						energy_transfer.flux.radii(i) = 8*pow(2, s*(i-3));
				}
			}
		}
		
		// SHELL-to-SHELL
		
		if (energy_transfer.shell_to_shell.turnon) {
			
			if ((energy_transfer.shell_to_shell.no_shells > 0) && (energy_transfer.shell_to_shell.radii.size() ==0)) {
				energy_transfer.shell_to_shell.radii.resize(energy_transfer.shell_to_shell.no_shells+1);
				energy_transfer.shell_to_shell.radii(0) = 0.0;
				energy_transfer.shell_to_shell.radii(1) = 2.0;
				energy_transfer.shell_to_shell.radii(2) = 4.0;
				energy_transfer.shell_to_shell.radii(3) = 8.0;
				
				energy_transfer.shell_to_shell.radii(energy_transfer.shell_to_shell.no_shells-2) = max_radius_inside/2.0;
				energy_transfer.shell_to_shell.radii(energy_transfer.shell_to_shell.no_shells-1) = max_radius_inside;
				energy_transfer.shell_to_shell.radii(energy_transfer.shell_to_shell.no_shells) = INF_RADIUS;
				
				if (energy_transfer.shell_to_shell.no_shells > 6) {
					Real s = log2(max_radius_inside/16.0) / (energy_transfer.shell_to_shell.no_shells - 5);
					for (int i=4; i<=(energy_transfer.shell_to_shell.no_shells-3); i++) 
						energy_transfer.shell_to_shell.radii(i) = 8*pow(2, s*(i-3));
				}
			}
				
		}
	
		
		// RING-to-RING
		if (energy_transfer.ring_to_ring.turnon) {
			// Read ring radii
		
			if ((energy_transfer.ring_to_ring.no_shells > 0) && (energy_transfer.ring_to_ring.radii.size()==0)) {
				energy_transfer.ring_to_ring.radii.resize(energy_transfer.ring_to_ring.no_shells+1);
				energy_transfer.ring_to_ring.radii(0) = 0.0;
				energy_transfer.ring_to_ring.radii(1) = 4.0;
				energy_transfer.ring_to_ring.radii(2) = 8.0;
				
				energy_transfer.ring_to_ring.radii(energy_transfer.ring_to_ring.no_shells-2) = max_radius_inside/2.0;
				energy_transfer.ring_to_ring.radii(energy_transfer.ring_to_ring.no_shells-1) = max_radius_inside;
				energy_transfer.ring_to_ring.radii(energy_transfer.ring_to_ring.no_shells) = INF_RADIUS;
				
				if (energy_transfer.ring_to_ring.no_shells > 5) {
					Real s = log2(max_radius_inside/16.0) / (energy_transfer.ring_to_ring.no_shells - 4);
					for (int i=3; i<=(energy_transfer.ring_to_ring.no_shells-3); i++) 
						energy_transfer.ring_to_ring.radii(i) = 8*pow(2, s*(i-2));
				}
				
			}
			
			// Read sector angles
						
			if (energy_transfer.ring_to_ring.sector_option == "EQUISPACED") {
				energy_transfer.ring_to_ring.sector_angles.resize(energy_transfer.ring_to_ring.no_sectors+1);
				dtheta = max_theta/energy_transfer.ring_to_ring.no_sectors;
				for (int i=0; i<energy_transfer.ring_to_ring.no_sectors; i++)
					energy_transfer.ring_to_ring.sector_angles(i) = i*dtheta;
				
				energy_transfer.ring_to_ring.sector_angles(energy_transfer.ring_to_ring.no_sectors) = max_theta;
			}
			
			else if (energy_transfer.ring_to_ring.sector_option == "EQUAL_NO_MODES") {
				;   // to implement
			}
		}
		
		// cylindrical-ring-to-ring
		if (energy_transfer.cylindrical_ring_to_ring.turnon) {
				// Read ring radii
			
			if ((energy_transfer.cylindrical_ring_to_ring.no_shells > 0) && (energy_transfer.cylindrical_ring_to_ring.radii.size() == 0)) {
				energy_transfer.cylindrical_ring_to_ring.radii.resize(energy_transfer.cylindrical_ring_to_ring.no_shells);
				
				energy_transfer.cylindrical_ring_to_ring.radii(0) = 0.0;
				energy_transfer.cylindrical_ring_to_ring.radii(1) = 4.0;
				energy_transfer.cylindrical_ring_to_ring.radii(2) = 8.0;
				
				energy_transfer.cylindrical_ring_to_ring.radii(energy_transfer.cylindrical_ring_to_ring.no_shells-2) =max_radius_inside/2.0;
				energy_transfer.cylindrical_ring_to_ring.radii(energy_transfer.cylindrical_ring_to_ring.no_shells-1) =max_radius_inside;
				 
				energy_transfer.cylindrical_ring_to_ring.radii(energy_transfer.cylindrical_ring_to_ring.no_shells) = INF_RADIUS;
				
				if (energy_transfer.cylindrical_ring_to_ring.no_shells > 5) {
					Real s = log2(max_cylinder_radius_inside/16.0) / (energy_transfer.cylindrical_ring_to_ring.no_shells - 4);
					
					for (int i=3; i<=(energy_transfer.cylindrical_ring_to_ring.no_shells-3); i++) 
						energy_transfer.cylindrical_ring_to_ring.radii(i) = 8*pow(2, s*(i-2));
				}
				
			}
			
			// Read slabs indices
			
			if (energy_transfer.cylindrical_ring_to_ring.kpll_option == "EQUISPACED") {
				energy_transfer.cylindrical_ring_to_ring.kpll_array.resize(energy_transfer.cylindrical_ring_to_ring.no_slabs+1);
				dKpll = (Kpll_max - Kpll_min)/energy_transfer.cylindrical_ring_to_ring.no_slabs;
				
				for (int i=0; i<=energy_transfer.cylindrical_ring_to_ring.no_slabs; i++)
					energy_transfer.cylindrical_ring_to_ring.kpll_array(i) = Kpll_min + i*dKpll;
			}
			
			else if (energy_transfer.cylindrical_ring_to_ring.kpll_option == "EQUAL_NO_MODES") {
				;   // to implement
			}
		}
		
	}
	
}
