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

/*! \file  init_cond_main.cc
 * 
 * @brief  Initial condition main file.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug No known bugs
 */

#include "../FluidIO.h"


/**********************************************************************************************

		Sets up initial conditions of the Vector and Scalar field

***********************************************************************************************/


void FluidIO::Read_init_cond(FluidVF& U)
{
    switch (global.io.input_field_procedure)	{
		
		case (1) : Init_cond_complex_field(U);	break;					
		// read from field_in_file
		
		case (2) : Init_cond_reduced_complex_field(U);	break;	
		// read from field_in_file with Nreduced D
		
		case (3): Init_cond_real_field(U); break; 
		
		case (4) : Init_cond_modes(U);	break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
		case (5) : Init_cond_energy_helicity_spectrum(U); break;	
		// given energy and hel spectrum
		
		case (6) : Init_cond_Taylor_Green(U); break;
		
		case (7) : Init_cond_ABC(U); break;
        
        case (400): Init_cond_vortex(U); break;
            
        case (420) : Init_cond_channel_flow(U); break;
            
            
	}
	
	if (master)
		cout  << "Reading of field configurations ended successfully" << endl;

    // preprocess the data
    int input_proc = global.io.input_field_procedure;

    
/*    if ((input_proc==1) || (input_proc==2) || (input_proc==3) || (input_proc==4) || (input_proc==5)) {
      U.Satisfy_strong_reality_condition_field();

        if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS")) 
            universal->Zero_modes(U.cvf.V1, U.cvf.V2, U.cvf.V3);
    }*/
}

//
//

void FluidIO::Read_init_cond(FluidVF& U, FluidSF& T)
{
	
	switch (global.io.input_field_procedure) {
		
		case (1) : Init_cond_complex_field(U, T);	break;			
		// read from field_in_file
		
		case (2) : Init_cond_reduced_complex_field(U, T);	break;		
		// read from field_in_file with Nreduced D
		
		case (3):Init_cond_real_field(U, T); break;	
		
		case (4) : Init_cond_modes(U, T); break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
		case (5) : Init_cond_energy_helicity_spectrum(U, T); break;	
		// given energy and hel spectrum
		
		case (6) : Init_cond_Taylor_Green(U, T); break;	
		// initialize only V field
		
		case (7) : Init_cond_ABC(U, T); break;			
		// initialize only V field
			
		case (100) : Init_cond_Rayleigh_Taylor(U, T); break;			
		// initialize For Rayleigh Taylor instability.
            
        case (400): Init_cond_vortex(U, T); break;
	}
	
	if (master)
		cout  << "Reading of field configurations ended successfully" << endl;
    
    // preprocess the data
    int input_proc = global.io.input_field_procedure;
    
    if ((input_proc==1) || (input_proc==2) || (input_proc==3) || (input_proc==4) || (input_proc==5)) {
        if (global.program.apply_strong_realitycond_alltime_switch){
	        U.Satisfy_strong_reality_condition_field();
	        T.Satisfy_strong_reality_condition_field();
    	}	
        
        if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS")) {
            universal->Zero_modes(U.cvf.V1,U.cvf.V2, U.cvf.V3);
            universal->Zero_modes(T.csf.F);
        } 
    }
}

void FluidIO::Read_init_cond(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	
	switch (global.io.input_field_procedure) {
			
		case (1) : Init_cond_complex_field(U, T1, T2);	break;			
				// read from field_in_file
			
		case (2) : Init_cond_reduced_complex_field(U, T1, T2);	break;		
				// read from field_in_file with Nreduced D
			
		case (3):Init_cond_real_field(U, T1, T2); break;	
			
		case (4) : Init_cond_modes(U, T1, T2); break;		
				// Modes - ki, Vx, (Vy:3D),Theta
			
	}
	
	if (global.mpi.master)
		cout  << "Reading of field configurations ended successfully" << endl;
    
    // preprocess the data
    int input_proc = global.io.input_field_procedure;
    
    if ((input_proc==1) || (input_proc==2) || (input_proc==3) || (input_proc==4)) {
        U.Satisfy_strong_reality_condition_field();
        T1.Satisfy_strong_reality_condition_field();
        T2.Satisfy_strong_reality_condition_field();
        
        if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS")) {
            universal->Zero_modes(U.cvf.V1, U.cvf.V2, U.cvf.V3);
            universal->Zero_modes(T1.csf.F);
            universal->Zero_modes(T2.csf.F);
        }
    }
}

//
//

void FluidIO::Read_init_cond(FluidVF& U, FluidVF& W)
{
	switch (global.io.input_field_procedure) {
		
		case (1) : Init_cond_complex_field(U, W);	break;			
		// read from field_in_file
		
		case (2) : Init_cond_reduced_complex_field(U, W); break;	
		// read from field_in_file with Nreduced D
		
		case (3): Init_cond_real_field(U, W); break;	
		
		case (4) : Init_cond_modes(U, W);	break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
		case (5) : Init_cond_energy_helicity_spectrum(U, W); break;	
		// given energy and hel spectrum
		
		case (6) : Init_cond_Taylor_Green(U, W); break;	
		// initialize V, W field
		
		case (7) : Init_cond_ABC(U, W); break;			
		// initialize V,W  field
		
		case (101) : Init_cond_DYNAMO_SIX_MODE(U, W); break;
		// For Verma et al.'s six mode model
			
		case (102) : Init_cond_dynamo_full_velocity_field(U, W); break;	
		// Full velocity field for a given Re, then introduce random B field in a band
            
        case (400): Init_cond_vortex(U, W); break;
	}
	
	if (global.mpi.master)
		cout  << "Reading of field configurations ended successfully" << endl;
    
    // preprocess the data
    int input_proc = global.io.input_field_procedure;
    
    if ((input_proc==1) || (input_proc==2) || (input_proc==3) || (input_proc==4) || (input_proc==5) || (input_proc==102)) {
        
        U.Satisfy_strong_reality_condition_field();
        W.Satisfy_strong_reality_condition_field();
        
        if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS")) {
            universal->Zero_modes(U.cvf.V1, U.cvf.V2, U.cvf.V3);
            universal->Zero_modes(W.cvf.V1, W.cvf.V2, W.cvf.V3);
        }
    }
}


//
//

void FluidIO::Read_init_cond(FluidVF& U, FluidVF& W, FluidSF& T)
{
	switch (global.io.input_field_procedure) {
		
	/*	case (1) : Init_cond_complex_field(U, W, T);	break;
		// read from field_in_file
		
		case (2) : Init_cond_reduced_complex_field(U, W, T);	break;	
		// read from field_in_file with Nreduced D
		
		case (3):Init_cond_real_field(U, W, T); break;	*/
		
		case (4) : Init_cond_modes(U, W, T);	break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
	/*	case (5) : Init_cond_energy_helicity_spectrum(U, W, T); break;
		// given energy and hel spectrum
		
		case (6) : Init_cond_Taylor_Green(U, W, T); break;	
		// initialize V, W field
		
		case (7) : Init_cond_ABC(U, W, T); break;	*/		
		// initialize V, W field
	}
	
	if (global.mpi.master)
		cout  << "Reading of field configurations ended successfully" << endl;
    
    // preprocess the data
    int input_proc = global.io.input_field_procedure;
    
    if ((input_proc==1) || (input_proc==2) || (input_proc==3) || (input_proc==4) || (input_proc==5)) {
        
        U.Satisfy_strong_reality_condition_field();
        W.Satisfy_strong_reality_condition_field();
        T.Satisfy_strong_reality_condition_field();
        
        if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS")) {
            universal->Zero_modes(U.cvf.V1, U.cvf.V2, U.cvf.V3);
            universal->Zero_modes(W.cvf.V1, W.cvf.V2, W.cvf.V3);
            universal->Zero_modes(T.csf.F);
        }
    }
}



void FluidIO::Read_init_cond(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	switch (global.io.input_field_procedure) {
			
	/*	case (1) : Init_cond_complex_field(U, W, T);	break;
			// read from field_in_file
			
		case (2) : Init_cond_reduced_complex_field(U, W, T);	break;
			// read from field_in_file with Nreduced D
			
		case (3):Init_cond_real_field(U, W, T); break; */
			
		case (4) : Init_cond_modes(U, W, T, C);	break;
			// Modes - ki, Vx, (Vy:3D),Theta
			
	/*	case (5) : Init_cond_energy_helicity_spectrum(U, W, T); break;
			// given energy and hel spectrum
			
		case (6) : Init_cond_Taylor_Green(U, W, T); break;
			// initialize V, W field
			
		case (7) : Init_cond_ABC(U, W, T); break; */
			// initialize V, W field
	}
	
	if (global.mpi.master)
		cout  << "Reading of field configurations ended successfully" << endl;
    
    // preprocess the data
 /*   int input_proc = global.io.input_field_procedure;
    
    if ((input_proc==1) || (input_proc==2) || (input_proc==3) || (input_proc==5)) {
        
        U.Satisfy_strong_reality_condition_field();
        W.Satisfy_strong_reality_condition_field();
        T.Satisfy_strong_reality_condition_field();
        
        if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS")) {
            universal->Zero_modes(U.cvf.V1, U.cvf.V2, U.cvf.V3);
            universal->Zero_modes(W.cvf.V1, W.cvf.V2, W.cvf.V3);
            universal->Zero_modes(T.csf.F);
        }
    } */
}


//******************* For GP

void FluidIO::Read_init_cond(FluidSF& T)
{
    
    switch (global.io.input_field_procedure)	{
            
        case (400): Init_cond_vortex(T); break;
            
	}
	
	if (my_id==master_id)
		cout  << "Reading of field configurations ended successfully" << endl;
    
    // preprocess the data
    int input_proc = global.io.input_field_procedure;
    
/*    if ((input_proc==1) || (input_proc==2) || (input_proc==3) || (input_proc==5)) {
        U.Satisfy_strong_reality_condition_field();
        
    if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
            universal->Zero_modes(U.cvf.V1, U.cvf.V2, U.cvf.V3);
    } */
}


//******************************** End of Init_cond_main.cc ***********************************




  
