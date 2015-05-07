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


#include "IMHDscalar_main.h"


//****************************************************************************************					

int IMHDscalar_main()
{
	global.field.diss_coefficients[0] = 1.0/global.PHYSICS.Reynolds;
	global.field.diss_coefficients[1] = 1.0/global.PHYSICS.Reynolds_mag;
	global.field.diss_coefficients[2] = 1.0/global.PHYSICS.Peclet;
	
	global.PHYSICS.Prandtl = global.PHYSICS.Peclet/global.PHYSICS.Reynolds;
	global.PHYSICS.Prandtl_mag = global.PHYSICS.Reynolds_mag/global.PHYSICS.Reynolds;
	
    // ITERATION...
    if (global.program.iter_or_diag == "ITERATION") {
        
        fluidIO_incompress.Open_files();
        fluidIO_incompress.Init_energy_transfer();
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        
        FluidVF  B(global.field.diss_coefficients[1], global.field.hyper_diss_coefficients[1], global.field.hyper_diss_exponents[1], global.force.W_switch, "B");
		
		 FluidSF T(global.field.diss_coefficients[2], global.field.hyper_diss_coefficients[2], global.field.hyper_diss_exponents[2], global.force.T_switch, "T");
        
        Pressure P;

        FORCE  Force;

        Time_advance_incompress  time_advance_incompress;
        // EnergyTr	energytr;

        fluidIO_incompress.Read_init_cond(U, B, T);
		
		if (master) {
			B.cvf.V1(0,0,0) = Complex(1.0,0.0);
		}
		
        
        Real total_abs_div;
        U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
        // true mean print nonzero div modes
        if (total_abs_div > MYEPS2) {
            cout << "abs(sum(Divergence)) of the initial field U = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
            return (0);
        }
        
        B.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
        if (total_abs_div > MYEPS2) {
            cout << "abs(sum(Divergence)) of the initial field B = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
            return (0);
        }
		
		if (master)
			cout << "mean B = ( " << B.cvf.V1(0,0,0)  << "," << B.cvf.V1(0,0,0) << ","  << B.cvf.V1(0,0,0)  << ")" << endl;

        fluidIO_incompress.Output_all_inloop(U, B, T, P);  // for initial cond
        
        if (my_id == master_id)  
                cout << endl << "STARTING THE SIMULATION NOW" << endl;

        int  iter=0;  // iterations 

        global.time.now = global.time.init;
        
        do 	{
            global.time.dt_computation_done = false;
            global.io.output_real_field_done = false;
            global.io.output_field_k_done = false;
            global.io.output_pressure_spectrum_done = false;
            global.io.output_pressure_done = false;

            iter++;

            time_advance_incompress.Time_advance_step(U, B, T, P, Force);
			
		//	cout << "T.F: after time step = " << endl;
		//	universal->Print_large_Fourier_elements(T.csf.F);

         /*    if (global.program.kind == "KEPLERIAN")  {
                time_advance_incompress.Make_field_incompressible(U);
                time_advance_incompress.Make_field_incompressible(B);
            } */
            
            U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
            // true mean print nonzero div modes
            if (total_abs_div > MYEPS2) {
                cout << "abs(sum(Divergence)) of  U = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
                return (0);
            }
            
            B.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
            // true mean print nonzero div modes
            if (total_abs_div > MYEPS2) {
                cout << "abs(sum(Divergence)) of  B = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
                return (0);
            }

            fluidIO_incompress.Output_all_inloop(U, B, T, P);
			
		//	fluidIO_incompress.Output_cout(U,B,T);
			
		//	cout << "iter = " << iter << endl;
		//	cout  << global.time.now << " "  << U.cvf.total_energy << " " << B.cvf.total_energy << " "
		//	<< T.csf.total_energy << endl;
            
            if ( (my_id == 0) && (isnan(U.cvf.total_energy) || (isnan(B.cvf.total_energy))))  {
                cout << "ERROR: Numerical Overflow " << endl;  break;
            }
        } 
        while ( (global.time.now < global.time.final) && (clock() < global.time.job_time_final) );
        
        fluidIO_incompress.Output_last(U, B, T, P);
        
        fluidIO_incompress.Close_files();
    }


    // DIAGNOSTICS
/*    else if (global.program.iter_or_diag == "DIAGNOSTICS") {
        string filename;
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        
        FluidVF  B(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "B");
        
        Pressure P;

        fluidIO_incompress.Init_energy_transfer();
        fluidIO_incompress.Read_init_cond(U, B);
        
        int i=0;
        while (i< global.io.diagnostic_procedures.size()) {
            switch (global.io.diagnostic_procedures[i])  {
            
                case (0) : {
                    filename = "/out/glob.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.global_file.open(filename.c_str());
                    if (!fluidIO_incompress.global_file.is_open())
                        cout << "UNABLE TO OPEN FILE global_file (glob.d) " << endl;
                    fluidIO_incompress.Output_global(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (1) : {
                    filename = "/out/spectrum.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_spectrum(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (2) : {
                    filename = "/out/ring_spectrum.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_spectrum(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                    
				case (3) : {
                    filename = "/out/cyl_ring_spectrum.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.cylindrical_ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_spectrum(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (4) : {
                    filename = "/out/flux.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.flux_file.open(filename.c_str());
                    fluidIO_incompress.Output_flux(U, B, P);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                    // force reqd for force-feed calculations
				case (5) : {
                    filename = "/out/shell_to_shell.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.shell_to_shell_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_to_shell(U, B, P);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                case (6) : {
                    filename = "/out/ring_to_ring.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_to_ring(U, B, P);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                    
				case (7) : {
                    filename = "/out/cylindrical_ring_to_ring.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.cylindrical_ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_to_ring(U, B, P);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    //		case (7) : fluidIO_incompress.Output_structure_fn(U);  break;
                    //		case (8) : fluidIO_incompress.Output_planar_structure_fn(U);  break;
				case (10) : {
                    filename = "/out/field_k_out_"+To_string(my_id)+".d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_k_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_k(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (11) : {
                    filename = "/out/field_r_out_"+To_string(my_id)+".d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_r_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_r(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (13) : {
                    filename = "/out/realfield_out.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.realfield_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_real_field(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
				case (14) : {
                    filename = "/out/field_out_reduced.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_out_reduced_file.open(filename.c_str());
                    fluidIO_incompress.Output_reduced_complex_field(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
            }

                i++;
        }	

    } */

    return(1);

} 


//********************************** End of Ifluid_main.cc ************************************	



