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


 #include "Ifluid_main.h"
#include "shell_etc_indices.h"

//****************************************************************************************					
 
int Ifluid_main()
{
		
	// ITERATION...
	if (global.program.iter_or_diag == "ITERATION") {
        
        fluidIO_incompress.Open_files();
        fluidIO_incompress.Init_energy_transfer();
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        
		Pressure P;

		FORCE  Force;

		Time_advance_incompress  time_advance_incompress;
		
		fluidIO_incompress.Read_init_cond(U);
		
        DP total_abs_div;
        U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
        
        // true mean print nonzero div modes
        if (total_abs_div > MYEPS2) {
            cout << "abs(sum(Divergence)) of the initial field U = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
            //return (1);
        }
        		
        fluidIO_incompress.Output_all_inloop(U, P);  // for initial cond
               
		if (my_id == master_id)  
			cout << endl << "STARTING THE SIMULATION NOW" << endl;
		int  iter=0;  // iterations 

		global.time.now = global.time.init;
        
        if (basis_type=="ChFF")
            time_advance_incompress.Compute_homgeneous_soln_influence_matrix(U,P);
        
        //        cout << "influence matrix = " << global.temp_array.influence_matrix << endl;
        
		do 	{
			global.time.dt_computation_done = false;
			global.io.output_real_field_done = false;
			global.io.output_field_k_done = false;
			global.io.output_pressure_spectrum_done = false;
			global.io.output_pressure_done = false;
			
			iter++; 
			
			time_advance_incompress.Time_advance_step(U, P, Force);
			
			DP total_abs_div;
        //    U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
            // true mean print nonzero div modes
            if (total_abs_div > MYEPS2) {
                cout << "abs(sum(Divergence)) of U = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
                return (1);
            }
	
			fluidIO_incompress.Output_all_inloop(U, P);
			
            if ( (my_id == 0) && isnan(U.cvf.total_energy) )  { 
                cout << "ERROR: Numerical Overflow " << endl;  break; 
            }
		}
		while ( (global.time.now < global.time.final) && (clock() < global.time.job_time_final) );
        
        fluidIO_incompress.Output_last(U, P);
		
		fluidIO_incompress.Close_files();
	}
	
	//*******************
	// DIAGNOSTICS
	else if (global.program.iter_or_diag == "DIAGNOSTICS") {
        string filename;
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        
        Pressure P;
        
        fluidIO_incompress.Init_energy_transfer();
        fluidIO_incompress.Read_init_cond(U);
        
		
		int i=0;
		while (i <= global.io.diagnostic_procedures.size()) {
			switch (global.io.diagnostic_procedures[i])  {		
				case (0) : { 
                    filename = "/out/glob.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.global_file.open(filename.c_str());
                    if (!fluidIO_incompress.global_file.is_open()) 
                        cout << "UNABLE TO OPEN FILE global_file (glob.d) " << endl;
                    fluidIO_incompress.Output_global(U);
                    fluidIO_incompress.Close_files();
                    break;	
                }
                    
				case (1) : {
                    filename = "/out/spectrum.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_spectrum(U);
                    fluidIO_incompress.Close_files();
                    break;	
                }
                    
				case (2) : {
                    filename = "/out/ring_spectrum.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_spectrum(U);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                    
				case (3) : {
                    filename = "/out/cyl_ring_spectrum.d";
                    filename = global.io.data_dir+ filename;  
                    fluidIO_incompress.cylindrical_ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_spectrum(U);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (4) : {
                    filename = "/out/flux.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.flux_file.open(filename.c_str());
                    fluidIO_incompress.Output_flux(U, P);
                    fluidIO_incompress.Close_files();
                    break;  
                }
                    
						// force reqd for force-feed calculations
				case (5) : {
                    filename = "/out/shell_to_shell.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.shell_to_shell_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_to_shell(U, P);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                case (6) : {
                    filename = "/out/ring_to_ring.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_to_ring(U, P);
                    fluidIO_incompress.Close_files();
                    break;	
                }
                    
                    
				case (7) : {
                    filename = "/out/cylindrical_ring_to_ring.d";
                    filename = global.io.data_dir+ filename;  
                    fluidIO_incompress.cylindrical_ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_to_ring(U, P);
                    fluidIO_incompress.Close_files();
                    break;	
                }	
						//		case (7) : fluidIO_incompress.Output_structure_fn(U);  break;
						//		case (8) : fluidIO_incompress.Output_planar_structure_fn(U);  break;
				case (10) : {
                    filename = "/out/field_k_out_"+To_string(my_id)+".d";
                    filename = global.io.data_dir+ filename;	  
                    fluidIO_incompress.field_k_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_k(U);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (11) : {
                    filename = "/out/field_r_out_"+To_string(my_id)+".d";
                    filename = global.io.data_dir+ filename;	  
                    fluidIO_incompress.field_r_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_r(U);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (13) : {
                    filename = "/out/realfield_out.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.realfield_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_real_field(U); 
                    fluidIO_incompress.Close_files();
                    break;
                }
				case (14) : {
                    filename = "/out/field_out_reduced.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_out_reduced_file.open(filename.c_str());
                    fluidIO_incompress.Output_reduced_complex_field(U);
                    fluidIO_incompress.Close_files();
                    break;
                }
			}
			
			i++;
		}	

	}
  
	return(1);
  
} 


//********************************** End of Ifluid_main.cc ************************************	



