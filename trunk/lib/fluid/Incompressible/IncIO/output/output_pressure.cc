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

/*! \file  Output_energy.cc
 * 
 * @brief  Output global, spectrum (shell, rings, cylinderical rings). 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "IncIO.h" 



/**********************************************************************************************
			
						Output_pressure_spectrum
						Output written in *CS_shell_ek of CSF inherited by NLIN

***********************************************************************************************/
 
void FluidIO_incompress::Output_pressure_spectrum(Pressure& P)
{

	if (global.time.now >= global.io.time.pressure_spectrum_save_next)  {
		if (global.mpi.master)	
			pressure_file << "%% Time = " << global.time.now << endl; 	

		universal->Compute_shell_spectrum(P.F, 0, Correlation::shell_ek);
		Print_array(pressure_file, "pressure shell spectrum", Correlation::shell_ek);
		
		if (global.spectrum.ring.turnon) {
			universal->Compute_ring_spectrum(P.F, 0, Correlation::ring_ek);
			Print_array(pressure_file, "pressure ring spectrum", Correlation::ring_ek);
		}

		
		if (global.spectrum.cylindrical_ring.turnon) {		
			universal->Compute_cylindrical_ring_spectrum(P.F, 0, Correlation::cylindrical_ring_ek);
			Print_array(pressure_file, "pressure cylindrical ring spectrum", Correlation::cylindrical_ring_ek);
		}
        
        // For RBC GL scaling analysis only.  Comment else
        global.program.sincostr_switch = global.program.sincostr_switch_divergence;
        universal->Xderiv(P.F, global.temp_array.X);
        
        DP sum_rmsgradp = sqrt(sum(sqr(abs(global.temp_array.X))));
        
        universal->Yderiv(P.F, global.temp_array.X);        
        sum_rmsgradp += sqrt(sum(sqr(abs(global.temp_array.X))));
        
        universal->Zderiv(P.F, global.temp_array.X);        
        sum_rmsgradp += sqrt(sum(sqr(abs(global.temp_array.X))));
        
        misc_file << "rms(p): time, rms(|grad(p)|): " << global.time.now << " " << sum_rmsgradp << endl;
		
        // print grad(p) ended
        
		global.io.time.pressure_spectrum_save_next += global.io.time.pressure_spectrum_save_interval;
	}
}



void FluidIO_incompress::Output_pressure(Pressure& P)
{
	if (global.time.now >= global.io.time.pressure_save_next) {
		//basicIO.Write_complex_array(P.F, "Pressure");
		global.io.time.pressure_save_next += global.io.time.pressure_save_interval;
	}
}


//*******************************   End of output_energy.cc ***********************************


