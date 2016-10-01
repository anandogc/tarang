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

//*********************************************************************************************

void Global::Process_basic_vars()
{
		// Alias....
	basis_type = program.basis_type;
	
	N[0]=field.N[0];
	N[1]=field.N[1];
	N[2]=field.N[2];
	N[3]=field.N[3];
	
	field.Nx=field.N[1];
	field.Ny=field.N[2];
	field.Nz=field.N[3];
	
	Nx=field.Nx;
	Ny=field.Ny;
	Nz=field.Nz;
	
	//job_time
	if (time.job_time.length() > 0){
		int hh, mm, ss;
		sscanf(time.job_time.c_str(), "%d:%d:%d", &hh, &mm, &ss);
		time.job_time_final = (hh*60*60 + mm*60 + ss)*CLOCKS_PER_SEC;
	}
	else
		time.job_time_final =  numeric_limits<clock_t>::max();


	//mpi
	my_id=mpi.my_id;
	numprocs=mpi.numprocs;
	
	master_id=mpi.master_id;
	master=mpi.master;

	num_p_rows = mpi.num_p_rows;
	num_p_cols = mpi.num_p_cols;

	I = myconstant.I;		
	
		/// minusI = -sqrt(-1).			
	minusI = myconstant.minusI;
	
		/// minus2I = -2*sqrt(-1).			
	minus2I = myconstant.minus2I;
	
	MYEPS = myconstant.MYEPS;
	MYEPS2 = myconstant.MYEPS2;
	
	MY_MAX_INT = myconstant.MY_MAX_INT;
		// cut off while reading diagnostic_procedure() array from input file and similar ops
	
		/// Infinite radius.. All the modes outside -- for flux and shelltr calc
	INF_RADIUS = myconstant.INF_RADIUS; 
	INF_TIME = myconstant.INF_TIME;

		
	// program
	program.T_exists = ((program.kind == "SCALAR_INCOMPRESS") || (program.kind == "RBC") || (program.kind == "STRATIFIED") || (program.kind == "MRBC")) || (program.kind == "MHD_SCALAR_INCOMPRESS") || (program.kind == "MHD_SCALAR_INCOMPRESS_ALL") || (program.kind == "MHD_ASTRO_INCOMPRESS") || (program.kind == "MHD_ASTRO_INCOMPRESS_ALL");
	program.W_exists = ((program.kind == "MHD_INCOMPRESS") || (program.kind == "KEPLERIAN")) || (program.kind == "MHD_SCALAR_INCOMPRESS") || (program.kind == "MHD_SCALAR_INCOMPRESS_ALL") || (program.kind == "MHD_ASTRO_INCOMPRESS") || (program.kind == "MHD_ASTRO_INCOMPRESS_ALL");

	program.sincostr_switch_Vx = basis_table[program.basis_type][program.sincostr_switch]["Vx"];
	program.sincostr_switch_Vy = basis_table[program.basis_type][program.sincostr_switch]["Vy"];
	program.sincostr_switch_Vz = basis_table[program.basis_type][program.sincostr_switch]["Vz"];
	program.sincostr_switch_F = basis_table[program.basis_type][program.sincostr_switch]["Vx"];
	program.sincostr_switch_Visqr = basis_table[program.basis_type][program.sincostr_switch]["Vi2"];
	program.sincostr_switch_VxVy = basis_table[program.basis_type][program.sincostr_switch]["VxVy"];
	program.sincostr_switch_VxVz = basis_table[program.basis_type][program.sincostr_switch]["VxVz"];
	program.sincostr_switch_VyVz = basis_table[program.basis_type][program.sincostr_switch]["VyVz"];
	
	program.sincostr_switch_FVx = basis_table[program.basis_type][program.sincostr_switch]["Vi2"];
	program.sincostr_switch_FVy = basis_table[program.basis_type][program.sincostr_switch]["VxVy"];
	program.sincostr_switch_FVz = basis_table[program.basis_type][program.sincostr_switch]["VxVz"];
	program.sincostr_switch_divergence = basis_table[program.basis_type][program.sincostr_switch]["div"];
	
		// ALIAS
	sincostr_switch=program.sincostr_switch;     // SINX, COSX, SCC, CSC, CCS, SSC, CSS, SCS
	sincostr_switch_Vx=program.sincostr_switch_Vx;
	sincostr_switch_Vy=program.sincostr_switch_Vy;
	sincostr_switch_Vz=program.sincostr_switch_Vz;
	sincostr_switch_F=program.sincostr_switch_F;
	sincostr_switch_Visqr=program.sincostr_switch_Visqr;
	sincostr_switch_VxVy=program.sincostr_switch_VxVy;
	sincostr_switch_VxVz=program.sincostr_switch_VxVz;
	sincostr_switch_VyVz=program.sincostr_switch_VyVz;
	sincostr_switch_FVx=program.sincostr_switch_FVx;
	sincostr_switch_FVy=program.sincostr_switch_FVy;
	sincostr_switch_FVz=program.sincostr_switch_FVz;
	sincostr_switch_divergence=program.sincostr_switch_divergence;
		// ALIAS END..

	if ((field.N[2] != 1) && (program.two_dimension || program.two_and_half_dimension)) {
		cerr << "N[2]=1 for 2D flows" << endl;
		exit(1);
	}
	
	int no_diss_coeff;
	if (program.kind == "FLUID_INCOMPRESS") {
		field.howmany = 3;
		no_diss_coeff = 1;
	}
	
	else if ((program.kind == "SCALAR_INCOMPRESS") || (program.kind =="RBC") || (program.kind == "STRATIFIED")) {
		field.howmany = 4;
		no_diss_coeff = 2;
	}
	
	else if ((program.kind == "MHD_INCOMPRESS") || (program.kind == "KEPLERIAN"))     {
		field.howmany = 6;
		no_diss_coeff = 2;
	}
	
	else if ((program.kind == "MHD_SCALAR_INCOMPRESS") || (program.kind == "MHD_SCALAR_INCOMPRESS_ALL"))     {
		field.howmany = 7;
		no_diss_coeff = 3;
	}
	
	else if ((program.kind == "MHD_ASTRO_INCOMPRESS") || (program.kind == "MHD_ASTRO_INCOMPRESS_ALL"))     {
		field.howmany = 8;
		no_diss_coeff = 4;
	}
	
	else if (program.kind == "MRBC") { // moist RBC
		field.howmany = 5;
		no_diss_coeff = 3;
	}
    
    else if (program.kind == "GP") { // Gross-Pateveskii
		field.howmany = 1;
		no_diss_coeff = 1;
	}

	if ((program.kind.length()>0) && (field.diss_coefficients.size() != no_diss_coeff) ) 
		Show_error("Number of Dissipation coefficeints must be equal to "+To_string(no_diss_coeff));

	
}
