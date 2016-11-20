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


 #include "GP_main.h"


//****************************************************************************************					
 
int GP_main()
{
		
	// ITERATION...
	if (global.program.iter_or_diag == "ITERATION") {
        
        fluidIO_incompress.Open_files();
 //       fluidIO_incompress.Init_energy_transfer();
        
        FluidSF T(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.T_switch, "T");
        
        FORCE  Force;

		Time_advance_incompress  time_advance_incompress;
		
		fluidIO_incompress.Read_init_cond(T);
//		fluidIO_incompress.Output_all_inloop(T);  // for initial cond
        
        cout <<  "Initial energy " << sum(sqr(abs(T.csf.F))) << endl;
		
		if (my_id == master_id)  
			cout << endl << "STARTING THE SIMULATION NOW" << endl;

		int  iter=0;  // iterations 

		global.time.now = global.time.init;
//		fluidIO_incompress.Output_field_k(T);
        
		do 	{
			global.time.dt_computation_done = false;
			global.io.output_real_field_done = false;
			global.io.output_field_k_done = false; 
			
			iter++; 
        
            global.time.dt = global.time.dt_fixed;
            global.time.now = global.time.now + global.time.dt;
            
			time_advance_incompress.Time_advance_step(T, Force);
            
            cout << "t, total energy = " << global.time.now << "  "<< sum(sqr(abs(T.csf.F))) << endl;
	
//            fluidIO_incompress.Output_all_inloop(T);
            
 /*           if ( (my_id == 0) && isnan(U.cvf.total_energy) )  {
                cout << "ERROR: Numerical Overflow " << endl;  break; 
            } */
		} 
		while (global.time.now < global.time.final);
        
        fluidIO_incompress.Close_files();
	}
	
	
	// DIAGNOSTICS
	else if (global.program.iter_or_diag == "DIAGNOSTICS") {
/*        string filename;
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        
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
                    fluidIO_incompress.global_file.close();
                    break;	
                }
                    
				case (1) : {
                    filename = "/out/spectrum.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_spectrum(U);
                    fluidIO_incompress.spectrum_file.close();
                    break;	
                }
                    
				case (2) : {
                    filename = "/out/ring_spectrum.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_spectrum(U);
                    fluidIO_incompress.ring_spectrum_file.close();
                    break;
                }
                    
                    
				case (3) : {
                    filename = "/out/cyl_ring_spectrum.d";
                    filename = global.io.data_dir+ filename;  
                    fluidIO_incompress.cylindrical_ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_spectrum(U);
                    fluidIO_incompress.cylindrical_ring_spectrum_file.close();
                    break;
                }
                    
				case (4) : {
                    filename = "/out/flux.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.flux_file.open(filename.c_str());
                    fluidIO_incompress.Output_flux(U);
                    fluidIO_incompress.flux_file.close();
                    break;  
                }
                    
						// force reqd for force-feed calculations
				case (5) : {
                    filename = "/out/shell_to_shell.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.shell_to_shell_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_to_shell(U);
                    fluidIO_incompress.shell_to_shell_file.close();
                    break;
                }
                    
                case (6) : {
                    filename = "/out/ring_to_ring.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_to_ring(U);
                    fluidIO_incompress.ring_to_ring_file.close();
                    break;	
                }
                    
                    
				case (7) : {
                    filename = "/out/cylindrical_ring_to_ring.d";
                    filename = global.io.data_dir+ filename;  
                    fluidIO_incompress.cylindrical_ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_to_ring(U);
                    fluidIO_incompress.cylindrical_ring_to_ring_file.close();
                    break;	
                }	
						//		case (7) : fluidIO_incompress.Output_structure_fn(U);  break;
						//		case (8) : fluidIO_incompress.Output_planar_structure_fn(U);  break;
				case (10) : {
                    stringstream ss;
                    ss << my_id;
                    filename = "/out/field_k_out_"+ss.str()+".d";
                    filename = global.io.data_dir+ filename;	  
                    fluidIO_incompress.field_k_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_k(U);
                    fluidIO_incompress.field_k_out_file.close();
                    break;
                }
                    
				case (11) : {
                    stringstream ss;
                    ss << my_id;
                    filename = "/out/field_r_out_"+ss.str()+".d";
                    filename = global.io.data_dir+ filename;	  
                    fluidIO_incompress.field_r_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_r(U);
                    fluidIO_incompress.field_r_out_file.close();
                    break;
                }
                    
				case (13) : {
                    filename = "/out/realfield_out.d";
                    filename = global.io.data_dir+ filename;   
                    fluidIO_incompress.realfield_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_real_field(U); 
                    fluidIO_incompress.realfield_out_file.close();
                    break;
                }
				case (14) : {
                    filename = "/out/field_out_reduced.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_out_reduced_file.open(filename.c_str());
                    fluidIO_incompress.Output_reduced_complex_field(U);
                    fluidIO_incompress.field_out_reduced_file.close();
                    break;
                }
			}
			
			i++;
		}	*/

	}
  
	return(1);
  
} 


//********************************** End of Ifluid_main.cc ************************************	



