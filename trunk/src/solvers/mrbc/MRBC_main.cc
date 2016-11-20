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


 #include "MRBC_main.h"


//****************************************************************************************					
 
int MRBC_main()
{

		//	cout << "ENTERING IFLUID MAIN my_id " << my_id << endl;
	fluidIO_incompress.Open_files();
	fluidIO_incompress.Init_energy_transfer();
	
	// Set Dissipation coefficients
	if (global.MRBC.Pr_option == "PRZERO") {
		global.field.diss_coefficients[0] = 1.0;
		global.field.diss_coefficients[1] = 0.0;
		global.field.diss_coefficients[2] = 0.0;
	}
	
	else if (global.MRBC.Pr_option == "PRLARGE") {
		if (global.MRBC.Uscaling == "USMALL") {
			global.field.diss_coefficients[0] = global.MRBC.Pr;              //  Coeff of grad^2 u
			global.field.diss_coefficients[1]  = 1.0;			// Coeff of grad^2 T
			global.field.diss_coefficients[2]  = 1.0;	
		}
		else if (global.MRBC.Uscaling == "ULARGE") {
			global.field.diss_coefficients[0] = sqrt(global.MRBC.Pr/global.MRBC.RaM);             
			global.field.diss_coefficients[1]  = 1/sqrt(global.MRBC.Pr*global.MRBC.RaM);
			global.field.diss_coefficients[2]  = 1/sqrt(global.MRBC.Pr*global.MRBC.RaM);
		}
	}
	
	else if (global.MRBC.Pr_option == "PRSMALL")  {
		if (global.MRBC.Uscaling == "USMALL")  {
			global.field.diss_coefficients[0] = 1.0;             
			global.field.diss_coefficients[1]  = 1/global.MRBC.Pr;
			global.field.diss_coefficients[2]  = 1/global.MRBC.Pr;
		}
		else if (global.MRBC.Uscaling == "ULARGE")  {
			global.field.diss_coefficients[0] = sqrt(global.MRBC.Pr/global.MRBC.RaM);             
			global.field.diss_coefficients[1]  = 1/sqrt(global.MRBC.Pr*global.MRBC.RaM);
			global.field.diss_coefficients[2]  = 1/sqrt(global.MRBC.Pr*global.MRBC.RaM);	
		}
	}
	
	else if (global.MRBC.Pr_option == "PRINFTY")  {
		if (global.MRBC.Uscaling == "USMALL")  {
			global.field.diss_coefficients[0] = global.MRBC.Pr;             
			global.field.diss_coefficients[1]  = 1.0;	
			global.field.diss_coefficients[2]  = 1.0;	
		}
		else if (global.MRBC.Uscaling == "ULARGE")  {
			global.field.diss_coefficients[0] = 1/sqrt(global.MRBC.RaM);	             
			global.field.diss_coefficients[1]  = 1/sqrt(global.MRBC.RaM);
			global.field.diss_coefficients[2]  = 1/sqrt(global.MRBC.RaM);
		}
	}

	FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
	
	FluidSF D(global.field.diss_coefficients[1], global.field.hyper_diss_coefficients[1], global.field.hyper_diss_exponents[1], global.force.T_switch, "D");
	
	FluidSF M(global.field.diss_coefficients[2], global.field.hyper_diss_coefficients[2], global.field.hyper_diss_exponents[2], global.force.T_switch, "M");

	// ITERATION...
	if (global.program.iter_or_diag == "ITERATION") {
		Pressure P;

		FORCE  Force;

		Time_advance_incompress  time_advance_incompress;
		// EnergyTr	energytr;
		
		fluidIO_incompress.Read_init_cond(U, D, M);
        
        Real total_abs_div;
        U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
        // true mean print nonzero div modes
        if (total_abs_div > MYEPS2) {
            cout << "abs(sum(Divergence)) of the initial field U = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
            return (0);
        }
        
		fluidIO_incompress.Output_all_inloop(U, D, M, P);  // for initial cond
		
		if (my_id == master_id)  
			cout << endl << "STARTING THE SIMULATION NOW" << endl;

		int  iter=0;  // iterations 

		global.time.now = global.time.init;
		fluidIO_incompress.Output_cout(U, D, M);
		fluidIO_incompress.Output_field_k(U, D, M);
		do 	{
			global.time.dt_computation_done = false;
			global.io.output_real_field_done = false;
			global.io.output_field_k_done = false;
			global.io.output_pressure_spectrum_done = false;
			global.io.output_pressure_done = false;

			iter++; 
			
			time_advance_incompress.Time_advance_step(U, D, M, P, Force);
			

            U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
            // true mean print nonzero div modes
            if (total_abs_div > MYEPS2) {
                cout << "abs(sum(Divergence)) of  U = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
                return (0);
            }
				//		fluidIO_incompress.Output_field_k(U, T);
	
				//			fluidIO_incompress.Output_all_inloop(U);
			
			
			fluidIO_incompress.Output_cout(U, D, M);
			
		} 
		while ( (global.time.now < global.time.final) && (clock() < global.time.job_time_final) );
        
        fluidIO_incompress.Output_last(U, D, M, P);
		
		fluidIO_incompress.Close_files();
	}
	
	
	// DIAGNOSTICS
	else if (global.program.iter_or_diag == "DIAGNOSTICS") {
		fluidIO_incompress.Read_init_cond(U, D, M);
        
		int i=1;
		while (global.io.diagnostic_procedures[i] <= global.io.diagnostic_procedures.size()) {
			switch (global.io.diagnostic_procedures[i])  {		
				case (0) : fluidIO_incompress.Output_global(U, D, M);		break;					
				case (1) : fluidIO_incompress.Output_shell_spectrum(U, D, M); break;			
						//				case (4) : fluidIO_incompress.Output_flux(U);		break;  
						// force reqd for force-feed calculations
				case (10) : fluidIO_incompress.Output_field_k(U, D, M);  break;
				case (11) : fluidIO_incompress.Output_field_r(U, D, M);  break;
				case (13) : fluidIO_incompress.Output_real_field(U, D, M);  break;
				case (14) : fluidIO_incompress.Output_reduced_complex_field(U, D, M); break;
			}
			
			i++;
		}	

	} 
  
	return(1);
  
} 


//********************************** End of Ifluid_main.cc ************************************	



