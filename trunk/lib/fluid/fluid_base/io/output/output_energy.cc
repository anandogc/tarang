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

#include "FluidIO.h" 
#include <algorithm>
//#include "../../scft/scft_energy.h"  
//#include "../Output.h"


//*********************************************************************************************
//class IncFluid;

/* class Output{
	public;
	IncFluid *incFluid;
	void global();
	
	connect(IncFluid *incFluid){
		this->incFluid=incFluid;
	}
}; */


bool FluidIO::Global_data_buffer_full()
{
	if ((global.io.global_data.buffer_index+global.io.global_data.packet_size) >  (global.io.global_data.buffer_size))
		return true;
    
	
	return false;
}



void FluidIO::Output_global(FluidVF& U)
{	
	U.cvf.Compute_total_energy();
	U.cvf.Compute_total_k2energy();
	U.cvf.Compute_entropy(); 
    
    if (global.program.helicity_switch)
        U.cvf.Compute_total_helicity();
    
    DP total_dissipation;
    // energy dissipation due to hyperviscosity.
    if (U.hyper_dissipation_switch) {
        DP kn_energy;
        U.cvf.Compute_total_kn_energy(U.hyper_dissipation_exponent, kn_energy);
        total_dissipation = U.hyper_dissipation_coefficient*kn_energy;
    }
	
	if (global.mpi.master)  {
		DP kmax = universal->Max_radius_inside();
		
		total_dissipation += U.dissipation_coefficient*U.cvf.total_k2energy;
		
		DP kolm_scale_u = sqrt(sqrt((my_pow(U.dissipation_coefficient, 3) /total_dissipation)));
		
		DP kmax_eta = kmax*kolm_scale_u;
		
		DP Rlambda = 2*U.cvf.total_energy* sqrt(15/total_dissipation); 			
		
		if (Global_data_buffer_full()) {
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
		}
		
			// Total 15 fields, so packetsize = 15
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.now;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Rlambda;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.dt;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E3;
	}
}

//*********************************************************************************************


void FluidIO::Output_global(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Output_global_scalar(U, T);
	
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_global_RBC(U, T);
}


void FluidIO::Output_global_scalar(FluidVF &U, FluidSF& T)
{

	U.cvf.Compute_total_energy();
	U.cvf.Compute_total_k2energy();

	T.csf.Compute_total_energy();
	T.csf.Compute_total_k2energy();
    
	if (global.program.helicity_switch)
        U.cvf.Compute_total_helicity();
    
	U.cvf.Compute_entropy(); 
	T.csf.Compute_entropy();
    
    DP total_dissipation, Ttotal_dissipation;
    
    // energy dissipation due to hyperviscosity.
    if (U.hyper_dissipation_switch) {
        DP kn_energy;
        U.cvf.Compute_total_kn_energy(U.hyper_dissipation_exponent, kn_energy);
        total_dissipation = U.hyper_dissipation_coefficient*kn_energy;
    }
    
    if (T.hyper_diffusion_switch) {
        DP kn_energy;
        T.csf.Compute_total_kn_energy(T.hyper_diffusion_exponent, kn_energy);
        Ttotal_dissipation = T.hyper_diffusion_coefficient*kn_energy;
    }
	
	if (global.mpi.master) {
		DP kmax = universal->Max_radius_inside();
		
		total_dissipation += U.dissipation_coefficient*U.cvf.total_k2energy;
		
		DP kolm_scale_u = sqrt(sqrt((my_pow(U.dissipation_coefficient,3) / total_dissipation)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = U.dissipation_coefficient/T.diffusion_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
	
		DP Rlambda = 2*U.cvf.total_energy* sqrt(15/total_dissipation);
		
        Ttotal_dissipation += T.diffusion_coefficient*T.csf.total_k2energy;
		
		if (Global_data_buffer_full()) {
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
		}
		
			// Total 20 fields
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.now;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Ttotal_dissipation;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.diffusion_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Rlambda;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.dt;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E3;
	}	
}

//*********************************************************************************************

// RB Convection //

void FluidIO::Output_global_RBC(FluidVF& U, FluidSF& T)
{
    
	static DP nusselt_no;
	
	U.cvf.Compute_total_energy();
	U.cvf.Compute_total_k2energy();
	
/*	DP residual_T_energy;    // Exclude |T(kx,0,0)|^2
	DP peclet_nu;			 
	DP nusselt_nu2;
	
	residual_T_energy = Get_total_energy_residual_SCFT(T.F);
	peclet_nu = globalvar_Ra*sqrt(2*Get_total_Sn_anis_SCFT(T.F, 6));
	// For P=infty, Pe= R*sum_K(theta(k)^2 Kperp^2/K^6)
		

	nusselt_nu2 = 1+global.PHYSICS.Ra*2*Get_total_Sn_anis_SCFT(T.F, 4);
	// For P=infty, Nu = 1+ R*sum_K(theta(k)^2 Kperp^2/K^4) */
	
	T.csf.Compute_total_energy();
	T.csf.Compute_total_k2energy();
	
	U.cvf.Compute_entropy(); 
	T.csf.Compute_entropy(); 
    
	if (global.program.helicity_switch)
        U.cvf.Compute_total_helicity();
	
	nusselt_no = Correlation::Get_Nusselt_no(U, T);
    
    DP total_dissipation, Ttotal_dissipation;
    
    // Energy dissipation due to hyperviscosity.
    if (U.hyper_dissipation_switch) {
        DP kn_energy;
        U.cvf.Compute_total_kn_energy(U.hyper_dissipation_exponent, kn_energy);
        total_dissipation = U.hyper_dissipation_coefficient*kn_energy;
    }
    
    if (T.hyper_diffusion_switch) {
        DP kn_energy;
        T.csf.Compute_total_kn_energy(T.hyper_diffusion_exponent, kn_energy);
        Ttotal_dissipation = T.hyper_diffusion_coefficient*kn_energy;
    }
	
	if (global.mpi.master) {
		DP kmax = universal->Max_radius_inside();
		
        total_dissipation += U.dissipation_coefficient*U.cvf.total_k2energy;
		
		DP kolm_scale_u = sqrt(sqrt((my_pow(U.dissipation_coefficient,3) / total_dissipation)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = U.dissipation_coefficient/T.diffusion_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*U.cvf.total_energy* sqrt(15/total_dissipation); 
		
        Ttotal_dissipation += T.diffusion_coefficient*T.csf.total_k2energy;
		
		if (Global_data_buffer_full()) {
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
		}
		
			// Total 21 fields
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.now;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Ttotal_dissipation;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = nusselt_no;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.diffusion_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Rlambda;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.dt;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E3;
	}
}

//*********************************************************************************************

void FluidIO::Output_global(FluidVF &U, FluidSF& T1, FluidSF& T2)
{	
	U.cvf.Compute_total_energy();
	U.cvf.Compute_total_k2energy();
	
	T1.csf.Compute_total_energy();
	T1.csf.Compute_total_k2energy();
	
	T2.csf.Compute_total_energy();
	T2.csf.Compute_total_k2energy();
	
    if (global.program.helicity_switch)
        U.cvf.Compute_total_helicity();
	
	U.cvf.Compute_entropy(); 
	T1.csf.Compute_entropy(); 
	T1.csf.Compute_entropy();
    
    DP total_dissipation, T1total_dissipation, T2total_dissipation;
    
    // Add energy dissipation due to hyperviscosity.
    if (U.hyper_dissipation_switch) {
        DP kn_energy;
        U.cvf.Compute_total_kn_energy(U.hyper_dissipation_exponent, kn_energy);
        total_dissipation = U.hyper_dissipation_coefficient*kn_energy;
    }
    
    if (T1.hyper_diffusion_switch) {
        DP kn_energy;
        T1.csf.Compute_total_kn_energy(T1.hyper_diffusion_exponent, kn_energy);
        T1total_dissipation = T1.hyper_diffusion_coefficient*kn_energy;
    }
    
    if (T2.hyper_diffusion_switch) {
        DP kn_energy;
        T2.csf.Compute_total_kn_energy(T2.hyper_diffusion_exponent, kn_energy);
        T2total_dissipation = T2.hyper_diffusion_coefficient*kn_energy;
    }
	
	if (global.mpi.master) {
		DP kmax = universal->Max_radius_inside();
		
		total_dissipation += U.dissipation_coefficient*U.cvf.total_k2energy;
		
		DP kolm_scale_u = sqrt(sqrt((my_pow(U.dissipation_coefficient,3) / total_dissipation)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = U.dissipation_coefficient/T1.diffusion_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*U.cvf.total_energy* sqrt(15/total_dissipation);
		
		T1total_dissipation += T1.diffusion_coefficient*T1.csf.total_k2energy;
		
		T2total_dissipation += T2.diffusion_coefficient*T2.csf.total_k2energy;
		
		
		if (Global_data_buffer_full()) {
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
		}
		
			// Total 24 fields
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.now;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T1.csf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T2.csf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T1total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T2total_dissipation;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T1.csf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T2.csf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.dissipation_coefficient;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T1.diffusion_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T2.diffusion_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Rlambda;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.dt;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E3;
	}	
}

//*********************************************************************************************

void FluidIO::Output_global(FluidVF& U, FluidVF& W)
{
	static DP Hc;

	U.cvf.Compute_total_energy();
	U.cvf.Compute_total_k2energy();

	
	W.cvf.Compute_total_energy();
	W.cvf.Compute_total_k2energy();

	Hc = Correlation::Get_cross_helicity(U, W);
	U.cvf.Compute_total_k2Hc(W.cvf.V1, W.cvf.V2, W.cvf.V3);

    if (global.program.helicity_switch) {
        U.cvf.Compute_total_helicity();
        W.cvf.Compute_total_helicity();
    }

	U.cvf.Compute_entropy();
	W.cvf.Compute_entropy();

    DP total_dissipation, Wtotal_dissipation;
    
    // Energy dissipation due to hyperviscosity.
    if (U.hyper_dissipation_switch) {
        DP kn_energy;
        U.cvf.Compute_total_kn_energy(U.hyper_dissipation_exponent, kn_energy);
        total_dissipation += U.hyper_dissipation_coefficient*kn_energy;
    }
    
    if (W.hyper_dissipation_switch) {
        DP kn_energy;
        W.cvf.Compute_total_kn_energy(W.hyper_dissipation_exponent, kn_energy);
        Wtotal_dissipation += W.hyper_dissipation_coefficient*kn_energy;
    }
    
		
	if (global.mpi.master) {
		DP kmax = universal->Max_radius_inside();
		
		total_dissipation += U.dissipation_coefficient*U.cvf.total_k2energy;
        
		DP total_dissipation_Hc = ((U.dissipation_coefficient+W.dissipation_coefficient)/2)*U.cvf.total_k2Hc;
		
		DP kolm_scale_u = sqrt(sqrt((my_pow(U.dissipation_coefficient,3) / total_dissipation)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = U.dissipation_coefficient/W.dissipation_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*U.cvf.total_energy* sqrt(15/total_dissipation);
		
		Wtotal_dissipation += W.dissipation_coefficient*W.cvf.total_k2energy;
		
		
		if (Global_data_buffer_full()) {
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
		}
		
			// Total 29 fields
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.now;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Wtotal_dissipation;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Hc;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation_Hc;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_k2H2;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Rlambda;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.dt;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E3;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E3; 
        
	}
}


//**************************************************************************************

void FluidIO::Output_global(FluidVF& U, FluidVF& W, FluidSF& T)
{
	static DP Hc;
	
	U.cvf.Compute_total_energy();
	U.cvf.Compute_total_k2energy();
	
	W.cvf.Compute_total_energy();
	W.cvf.Compute_total_k2energy();
	
	T.csf.Compute_total_energy();
	T.csf.Compute_total_k2energy();

	
	Hc = Correlation::Get_cross_helicity(U, W);
	U.cvf.Compute_total_k2Hc(W.cvf.V1, W.cvf.V2, W.cvf.V3);
	
    if (global.program.helicity_switch) {
        U.cvf.Compute_total_helicity();
        W.cvf.Compute_total_helicity();
    }
	
	DP nusselt_no = Correlation::Get_Nusselt_no(U, T);
	
	U.cvf.Compute_entropy();
	W.cvf.Compute_entropy();
	
    DP total_dissipation, Wtotal_dissipation, Ttotal_dissipation;
    
    // Energy dissipation due to hyperviscosity.
    if (U.hyper_dissipation_switch) {
        DP kn_energy;
        U.cvf.Compute_total_kn_energy(U.hyper_dissipation_exponent, kn_energy);
        total_dissipation += U.hyper_dissipation_coefficient*kn_energy;
    }
    
    if (W.hyper_dissipation_switch) {
        DP kn_energy;
        W.cvf.Compute_total_kn_energy(W.hyper_dissipation_exponent, kn_energy);
        Wtotal_dissipation += W.hyper_dissipation_coefficient*kn_energy;
    }
    
	if (T.hyper_diffusion_switch) {
        DP kn_energy;
        T.csf.Compute_total_kn_energy(T.hyper_diffusion_exponent, kn_energy);
        Ttotal_dissipation = T.hyper_diffusion_coefficient*kn_energy;
    }
	
	if (global.mpi.master) {
		DP kmax = universal->Max_radius_inside();
		
		total_dissipation += U.dissipation_coefficient*U.cvf.total_k2energy;
        
		DP total_dissipation_Hc = ((U.dissipation_coefficient+W.dissipation_coefficient)/2)*U.cvf.total_k2Hc;
		
		DP kolm_scale_u = sqrt(sqrt((my_pow(U.dissipation_coefficient,3) / total_dissipation)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = U.dissipation_coefficient/W.dissipation_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*U.cvf.total_energy* sqrt(15/total_dissipation);
		
		Wtotal_dissipation += W.dissipation_coefficient*W.cvf.total_k2energy;
		
		
		if (Global_data_buffer_full()) {
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
		}
		
		// Total 34 fields
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.now;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Wtotal_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Ttotal_dissipation;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Hc;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = nusselt_no;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation_Hc;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_k2H2;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.entropy;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.diffusion_coefficient;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Rlambda;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.dt;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E3;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E3;
        
	}
}


//**************************************************************************************

void FluidIO::Output_global(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	static DP Hc;
	
	U.cvf.Compute_total_energy();
	U.cvf.Compute_total_k2energy();
	
	W.cvf.Compute_total_energy();
	W.cvf.Compute_total_k2energy();
	
	T.csf.Compute_total_energy();
	T.csf.Compute_total_k2energy();
	
	C.csf.Compute_total_energy();
	C.csf.Compute_total_k2energy();
	
	Hc = Correlation::Get_cross_helicity(U, W);
	U.cvf.Compute_total_k2Hc(W.cvf.V1, W.cvf.V2, W.cvf.V3);
	
    if (global.program.helicity_switch) {
        U.cvf.Compute_total_helicity();
        W.cvf.Compute_total_helicity();
    }
	
	DP nusselt_no = Correlation::Get_Nusselt_no(U, T);
	
	U.cvf.Compute_entropy();
	W.cvf.Compute_entropy();
	
    DP total_dissipation, Wtotal_dissipation, Ttotal_dissipation, Ctotal_dissipation;
    
    // Energy dissipation due to hyperviscosity.
    if (U.hyper_dissipation_switch) {
        DP kn_energy;
        U.cvf.Compute_total_kn_energy(U.hyper_dissipation_exponent, kn_energy);
        total_dissipation += U.hyper_dissipation_coefficient*kn_energy;
    }
    
    if (W.hyper_dissipation_switch) {
        DP kn_energy;
        W.cvf.Compute_total_kn_energy(W.hyper_dissipation_exponent, kn_energy);
        Wtotal_dissipation += W.hyper_dissipation_coefficient*kn_energy;
    }
    
	if (T.hyper_diffusion_switch) {
        DP kn_energy;
        T.csf.Compute_total_kn_energy(T.hyper_diffusion_exponent, kn_energy);
        Ttotal_dissipation = T.hyper_diffusion_coefficient*kn_energy;
    }
	
	if (C.hyper_diffusion_switch) {
        DP kn_energy;
        C.csf.Compute_total_kn_energy(C.hyper_diffusion_exponent, kn_energy);
        Ctotal_dissipation = C.hyper_diffusion_coefficient*kn_energy;
    }
	
	
	if (global.mpi.master) {
		DP kmax = universal->Max_radius_inside();
		
		total_dissipation += U.dissipation_coefficient*U.cvf.total_k2energy;
        
		DP total_dissipation_Hc = ((U.dissipation_coefficient+W.dissipation_coefficient)/2)*U.cvf.total_k2Hc;
		
		DP kolm_scale_u = sqrt(sqrt((my_pow(U.dissipation_coefficient,3) / total_dissipation)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = U.dissipation_coefficient/W.dissipation_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*U.cvf.total_energy* sqrt(15/total_dissipation);
		
		Wtotal_dissipation += W.dissipation_coefficient*W.cvf.total_k2energy;
		
		
		if (Global_data_buffer_full()) {
			global.io.Dump_buffer(global_file, global.io.global_data.buffer, global.io.global_data.buffer_index, global.io.global_data.packet_size);
		}
		
		// Total 38 fields
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.now;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = C.csf.total_energy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Wtotal_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Ttotal_dissipation;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Ctotal_dissipation;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Hc;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = nusselt_no;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = total_dissipation_Hc;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_helicity1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_helicity2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_k2H2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_k2H1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_k2H2;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.csf.entropy;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = C.csf.entropy;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = kmax_eta2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.dissipation_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = T.diffusion_coefficient;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = C.diffusion_coefficient;
		
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = Rlambda;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = global.time.dt;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = U.cvf.total_E3;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E1;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E2;
		global.io.global_data.buffer(global.io.global_data.buffer_index++) = W.cvf.total_E3;
        
	}
}




//*********************************************************************************************

// Outputs global.time.now and total energy to cout.
void FluidIO::Output_cout(FluidVF &U)
{	
	U.cvf.Compute_total_energy();

	if (my_id == master_id)
		cout  << global.time.now << " "  << U.cvf.total_energy << endl; 
}


void FluidIO::Output_cout(FluidVF &U, FluidSF& T)
{
	U.cvf.Compute_total_energy();
	T.csf.Compute_total_energy();
	
	if (my_id == master_id)
		cout  << global.time.now << " "  << U.cvf.total_energy << " " << T.csf.total_energy  << endl;
}

void FluidIO::Output_cout(FluidVF &U, FluidSF& T1, FluidSF& T2)
{
	U.cvf.Compute_total_energy();
	T1.csf.Compute_total_energy();
	T2.csf.Compute_total_energy();
	
	if (my_id == master_id)
		cout  << global.time.now << " "  << U.cvf.total_energy << " " << T1.csf.total_energy << " " << T2.csf.total_energy  << endl;
}

void FluidIO::Output_cout(FluidVF &U, FluidVF &W)
{
	U.cvf.Compute_total_energy();
	W.cvf.Compute_total_energy();
	
	if (my_id == master_id)
		cout  << global.time.now << " "  << U.cvf.total_energy << " " << W.cvf.total_energy << endl;
}

void FluidIO::Output_cout(FluidVF &U, FluidVF &W, FluidSF & T)		
{
	U.cvf.Compute_total_energy();
	W.cvf.Compute_total_energy();
	T.csf.Compute_total_energy();
	
	if (my_id == master_id)
		cout  << global.time.now << " "  << U.cvf.total_energy << " " << W.cvf.total_energy << " " 
		<< T.csf.total_energy << endl;
}


void FluidIO::Output_cout(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	U.cvf.Compute_total_energy();
	W.cvf.Compute_total_energy();
	T.csf.Compute_total_energy();
	C.csf.Compute_total_energy();
	
	if (my_id == master_id)
		cout  << global.time.now << " "  << U.cvf.total_energy << " " << W.cvf.total_energy << " "
		<< T.csf.total_energy << " " << C.csf.total_energy << endl;
}

//*********************************************************************************************
//*********************************************************************************************
 
void FluidIO::Output_shell_spectrum(FluidVF& U)
{
	if (global.spectrum.shell.turnon) {
		if (global.mpi.master) 
			spectrum_file << "%% Time = " << global.time.now << "\n \n"; 	

		
		Correlation::Compute_shell_spectrum(U);
		Print_array(spectrum_file, "Uek", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		Print_array(spectrum_file, "UDk", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
		
		if (U.force_switch) {
			Correlation::Compute_force_shell_spectrum(U);
			Print_array(spectrum_file, "(Fv.v)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_shell_spectrum_helicity(U);
			Print_array(spectrum_file, "hk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (master)
			spectrum_file.flush();

	}
}  

//*********************************************************************************************  
// scalar
void FluidIO::Output_shell_spectrum(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Output_shell_spectrum_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_shell_spectrum_RBC(U, T);
}


  
void FluidIO::Output_shell_spectrum_scalar(FluidVF& U, FluidSF& T)
{
	
	if (global.spectrum.shell.turnon) {
		if (global.mpi.master)
			spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_shell_spectrum(U);
		Print_array(spectrum_file, "Uek", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		Print_array(spectrum_file, "UDk", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
		
		Correlation::Compute_shell_spectrum(T);
		Print_array(spectrum_file, "Tek", Correlation::shell_ek);
		
		Correlation::Compute_shell_spectrum(U,T);
		Print_array(spectrum_file, "(U.T)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		if (U.force_switch) {
			Correlation::Compute_force_shell_spectrum(U);
			Print_array(spectrum_file, "(Fv.v)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (T.force_switch) {
			Correlation::Compute_force_shell_spectrum(T);
			Print_array(spectrum_file, "(FT.T)(k)", Correlation::shell_ek);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_shell_spectrum_helicity(U);
			Print_array(spectrum_file, "hk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (master)
			spectrum_file.flush();

	}
}

//  RB Convection	


void FluidIO::Output_shell_spectrum_RBC(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO")
		Output_shell_spectrum(U);
	
	else
		Output_shell_spectrum_scalar(U, T);
}


void FluidIO::Output_shell_spectrum(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	
	if (global.spectrum.shell.turnon) {
		if (global.mpi.master)
			spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_shell_spectrum(U);
		Print_array(spectrum_file, "Uek", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		Print_array(spectrum_file, "UDk", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
		
		Correlation::Compute_shell_spectrum(T1);
		Print_array(spectrum_file, "T1ek", Correlation::shell_ek);
		
		Correlation::Compute_shell_spectrum(T2);
		Print_array(spectrum_file, "T2ek", Correlation::shell_ek);
		
		Correlation::Compute_shell_spectrum(U,T1);
		Print_array(spectrum_file, "(U.T1)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		Correlation::Compute_shell_spectrum(U,T2);
		Print_array(spectrum_file, "(U.T2)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		if (U.force_switch) {
			Correlation::Compute_force_shell_spectrum(U);
			Print_array(spectrum_file, "(Fv.v)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (T1.force_switch) {
			Correlation::Compute_force_shell_spectrum(T1);
			Print_array(spectrum_file, "(FT1.T1)(k)", Correlation::shell_ek);
		}
		
		if (T2.force_switch) {
			Correlation::Compute_force_shell_spectrum(T2);
			Print_array(spectrum_file, "(FT2.T2)(k)", Correlation::shell_ek);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_shell_spectrum_helicity(U);
			Print_array(spectrum_file, "hk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (master)
			spectrum_file.flush();

	}
}
 
//********************************************************************************************* 
// Vector
  
void FluidIO::Output_shell_spectrum(FluidVF& U, FluidVF& W)
{
	if (global.spectrum.shell.turnon) {
		if (global.mpi.master)
			spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_shell_spectrum(U);
		Print_array(spectrum_file, "Uek", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		Print_array(spectrum_file, "UDk", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
		
		Correlation::Compute_shell_spectrum(W);
		Print_array(spectrum_file, "Wek", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		Print_array(spectrum_file, "WDk", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
		
		Correlation::Compute_shell_spectrum(U,W);
		Print_array(spectrum_file, "(U.W)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		if (U.force_switch) {
			Correlation::Compute_force_shell_spectrum(U);
			Print_array(spectrum_file, "(Fv.v)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (W.force_switch) {
			Correlation::Compute_force_shell_spectrum(W);
			Print_array(spectrum_file, "(Fw.w)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_shell_spectrum_helicity(U);
			Print_array(spectrum_file, "hk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
			
			Correlation::Compute_shell_spectrum_helicity(W);
			Print_array(spectrum_file, "Whk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (master)
			spectrum_file.flush();
	}
}


//*********************************************************************************************
// Magneto+scalar
  
void FluidIO::Output_shell_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	if (global.spectrum.shell.turnon) {
		if (global.mpi.master)
			spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_shell_spectrum(U);
		Print_array(spectrum_file, "Uek", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		Print_array(spectrum_file, "UDk", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
		
		Correlation::Compute_shell_spectrum(W);
		Print_array(spectrum_file, "Wek", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		Print_array(spectrum_file, "WDk", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
		
		Correlation::Compute_shell_spectrum(T);
		Print_array(spectrum_file, "Tek", Correlation::shell_ek);
		
		Correlation::Compute_shell_spectrum(U,W);
		Print_array(spectrum_file, "(U.W)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		Correlation::Compute_shell_spectrum(U,T);
		Print_array(spectrum_file, "(U.T)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		if (U.force_switch) {
			Correlation::Compute_force_shell_spectrum(U);
			Print_array(spectrum_file, "(Fv.v)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (W.force_switch) {
			Correlation::Compute_force_shell_spectrum(W);
			Print_array(spectrum_file, "(Fw.W)(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (T.force_switch) {
			Correlation::Compute_force_shell_spectrum(T);
			Print_array(spectrum_file, "(FT.T)(k)", Correlation::shell_ek);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_shell_spectrum_helicity(U);
			Print_array(spectrum_file, "hk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
			
			Correlation::Compute_shell_spectrum_helicity(W);
			Print_array(spectrum_file, "Whk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		}
		
		if (master)
			spectrum_file.flush();
	}
}



 
/**********************************************************************************************
	
						Output_ring_spectrum()

***********************************************************************************************/
 
 
void FluidIO::Output_ring_spectrum(FluidVF& U)
{
	if (global.spectrum.ring.turnon) {
	
		if (global.mpi.master)
			ring_spectrum_file << "%% Time = " << global.time.now << endl; 	

		Correlation::Compute_ring_spectrum(U);
		Print_array(ring_spectrum_file, "Uek", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		Print_array(ring_spectrum_file, "UDk", Correlation::ring_dissk1, Correlation::ring_dissk2, Correlation::ring_dissk3);
		
		if (U.force_switch) {
			Correlation::Compute_force_ring_spectrum(U);
			Print_array(ring_spectrum_file, "(Fv.v)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_ring_spectrum_helicity(U);
			Print_array(ring_spectrum_file, "hk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (master)
			spectrum_file.flush();
	}
}  
 

//*********************************************************************************************  
// scalar
void FluidIO::Output_ring_spectrum(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR") 
		Output_ring_spectrum_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_ring_spectrum_RBC(U, T);
}
  
void FluidIO::Output_ring_spectrum_scalar(FluidVF& U, FluidSF& T)
{

	if (global.spectrum.ring.turnon) {
		
		if (global.mpi.master)
			ring_spectrum_file << "%% Time = " << global.time.now << endl; 	

		Correlation::Compute_ring_spectrum(U);
		Print_array(ring_spectrum_file, "Uek", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		Print_array(ring_spectrum_file, "UDk", Correlation::ring_dissk1, Correlation::ring_dissk2, Correlation::ring_dissk3);
		
		Correlation::Compute_ring_spectrum(T);
		Print_array(ring_spectrum_file, "Tek", Correlation::ring_ek);
		
		Correlation::Compute_ring_spectrum(U,T);
		Print_array(ring_spectrum_file, "(U.T)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		if (U.force_switch) {
			Correlation::Compute_force_ring_spectrum(U);
			Print_array(ring_spectrum_file, "(Fv.v)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (T.force_switch) {
			Correlation::Compute_force_ring_spectrum(T);
			Print_array(ring_spectrum_file, "(FT.T)(k)", Correlation::ring_ek);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_ring_spectrum_helicity(U);
			Print_array(ring_spectrum_file, "hk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (master)
			ring_spectrum_file.flush();
	}
}

//  RB Convection	//

void FluidIO::Output_ring_spectrum_RBC(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO")
		Output_ring_spectrum(U);
	
	else
		Output_ring_spectrum_scalar(U, T);
}
 
//********************************************************************************************* 
// Vector
  
void FluidIO::Output_ring_spectrum(FluidVF& U, FluidVF& W)
{
	if (global.spectrum.ring.turnon) {
		
		if (global.mpi.master)
			ring_spectrum_file << "%% Time = " << global.time.now << endl;
		
		Correlation::Compute_ring_spectrum(U);
		Print_array(ring_spectrum_file, "Uek", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		Print_array(ring_spectrum_file, "UDk", Correlation::ring_dissk1, Correlation::ring_dissk2, Correlation::ring_dissk3);
		
		Correlation::Compute_ring_spectrum(W);
		Print_array(ring_spectrum_file, "Wek", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		Print_array(ring_spectrum_file, "WDk", Correlation::ring_dissk1, Correlation::ring_dissk2, Correlation::ring_dissk3);
		
		Correlation::Compute_ring_spectrum(U,W);
		Print_array(ring_spectrum_file, "(U.W)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		if (U.force_switch) {
			Correlation::Compute_force_ring_spectrum(U);
			Print_array(ring_spectrum_file, "(Fv.v)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (W.force_switch) {
			Correlation::Compute_force_ring_spectrum(W);
			Print_array(ring_spectrum_file, "(Fw.W)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}

		if (global.program.helicity_switch) {
			Correlation::Compute_ring_spectrum_helicity(U);
			Print_array(ring_spectrum_file, "hk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
			
			Correlation::Compute_ring_spectrum_helicity(W);
			Print_array(ring_spectrum_file, "Whk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (master)
			ring_spectrum_file.flush();
	}

}


//*********************************************************************************************
// Magneto+scalar
  
void FluidIO::Output_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{
	if (global.spectrum.ring.turnon) {
		
		if (global.mpi.master)
			ring_spectrum_file << "%% Time = " << global.time.now << endl; 	

		Correlation::Compute_ring_spectrum(U);
		Print_array(ring_spectrum_file, "Uek", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		Print_array(ring_spectrum_file, "UDk", Correlation::ring_dissk1, Correlation::ring_dissk2, Correlation::ring_dissk3);
		
		Correlation::Compute_ring_spectrum(W);
		Print_array(ring_spectrum_file, "Wek", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		Print_array(ring_spectrum_file, "WDk", Correlation::ring_dissk1, Correlation::ring_dissk2, Correlation::ring_dissk3);
		
		Correlation::Compute_ring_spectrum(T);
		Print_array(ring_spectrum_file, "Tek", Correlation::ring_ek);
		
		Correlation::Compute_ring_spectrum(U,W);
		Print_array(ring_spectrum_file, "(U.W)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_ring_spectrum(U,T);
		Print_array(ring_spectrum_file, "(U.T)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		if (U.force_switch) {
			Correlation::Compute_force_ring_spectrum(U);
			Print_array(ring_spectrum_file, "(Fv.v)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (W.force_switch) {
			Correlation::Compute_force_ring_spectrum(W);
			Print_array(ring_spectrum_file, "(Fw.W)(k)", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (T.force_switch) {
			Correlation::Compute_force_ring_spectrum(T);
			Print_array(ring_spectrum_file, "(FT.T)(k)", Correlation::ring_ek);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_ring_spectrum_helicity(U);
			Print_array(ring_spectrum_file, "hk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
			
			Correlation::Compute_ring_spectrum_helicity(W);
			Print_array(ring_spectrum_file, "Whk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		}
		
		if (master)
			ring_spectrum_file.flush();
	}

}


/**********************************************************************************************
	
						Output_cylindrical_ring_spectrum()

***********************************************************************************************/
 
 
void FluidIO::Output_cylindrical_ring_spectrum(FluidVF& U)
{
	
	if (global.spectrum.cylindrical_ring.turnon) {
		
		if (global.mpi.master)
			cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << endl;
		
		Correlation::Compute_cylindrical_ring_spectrum(U);
		Print_array(cylindrical_ring_spectrum_file, "Uek", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		Print_array(cylindrical_ring_spectrum_file, "UDk", Correlation::cylindrical_ring_dissk1, Correlation::cylindrical_ring_dissk2);
		
		if (U.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(U);
			Print_array(cylindrical_ring_spectrum_file, "(Fv.v)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_cylindrical_ring_spectrum_helicity(U);
			Print_array(cylindrical_ring_spectrum_file, "hk", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (master)
			cylindrical_ring_spectrum_file.flush();
	}								

}  
  
//
//

//
//  scalar
void FluidIO::Output_cylindrical_ring_spectrum(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR") 
		Output_cylindrical_ring_spectrum_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_cylindrical_ring_spectrum_RBC(U, T);
}


void FluidIO::Output_cylindrical_ring_spectrum_scalar(FluidVF& U, FluidSF& T)
{
	
	if (global.spectrum.cylindrical_ring.turnon){
		if (global.mpi.master)
			cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << endl;
		
		Correlation::Compute_cylindrical_ring_spectrum(U);
		Print_array(cylindrical_ring_spectrum_file, "Uek", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		Print_array(cylindrical_ring_spectrum_file, "UDk", Correlation::cylindrical_ring_dissk1, Correlation::cylindrical_ring_dissk2);
		
		Correlation::Compute_cylindrical_ring_spectrum(T);
		Print_array(cylindrical_ring_spectrum_file, "Tek", Correlation::cylindrical_ring_ek);
		
		Correlation::Compute_cylindrical_ring_spectrum(U,T);
		Print_array(cylindrical_ring_spectrum_file, "(U.T)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		
		if (U.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(U);
			Print_array(cylindrical_ring_spectrum_file, "(Fv.v)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (T.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(T);
			Print_array(cylindrical_ring_spectrum_file, "(FT.T)(k)", Correlation::cylindrical_ring_ek);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_cylindrical_ring_spectrum_helicity(U);
			Print_array(cylindrical_ring_spectrum_file, "hk", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (master)
			cylindrical_ring_spectrum_file.flush();
	}															
								
} 

//  RB Convection	//

void FluidIO::Output_cylindrical_ring_spectrum_RBC(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO")
		Output_cylindrical_ring_spectrum(U);
	
	else
		Output_cylindrical_ring_spectrum_scalar(U, T);
}

//
//

void FluidIO::Output_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W)
{

	if (global.spectrum.cylindrical_ring.turnon) {

		if (global.mpi.master)
			cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << endl;
		
		
		Correlation::Compute_cylindrical_ring_spectrum(U);
		Print_array(cylindrical_ring_spectrum_file, "Uek", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		Print_array(cylindrical_ring_spectrum_file, "UDk", Correlation::cylindrical_ring_dissk1, Correlation::cylindrical_ring_dissk2);
		
		Correlation::Compute_cylindrical_ring_spectrum(W);
		Print_array(cylindrical_ring_spectrum_file, "Wek", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		Print_array(cylindrical_ring_spectrum_file, "WDk", Correlation::cylindrical_ring_dissk1, Correlation::cylindrical_ring_dissk2);
		
		Correlation::Compute_cylindrical_ring_spectrum(U,W);
		Print_array(cylindrical_ring_spectrum_file, "(U.W)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
	
		
		if (U.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(U);
			Print_array(cylindrical_ring_spectrum_file, "(Fv.v)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (W.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(W);
			Print_array(cylindrical_ring_spectrum_file, "(Fw.W)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}

		
		if (global.program.helicity_switch) {
			Correlation::Compute_cylindrical_ring_spectrum_helicity(U);
			Print_array(cylindrical_ring_spectrum_file, "hk", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
			
			Correlation::Compute_cylindrical_ring_spectrum_helicity(W);
			Print_array(cylindrical_ring_spectrum_file, "Whk", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (master)
			cylindrical_ring_spectrum_file.flush();
				
	}																
									
}  
  
//
//

void FluidIO::Output_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{	
	
	if (global.spectrum.cylindrical_ring.turnon) { 	

		if (global.mpi.master)
			cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << endl;
		
		Correlation::Compute_cylindrical_ring_spectrum(U);
		Print_array(cylindrical_ring_spectrum_file, "Uek", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		Print_array(cylindrical_ring_spectrum_file, "UDk", Correlation::cylindrical_ring_dissk1, Correlation::cylindrical_ring_dissk2);
		
		Correlation::Compute_cylindrical_ring_spectrum(W);
		Print_array(cylindrical_ring_spectrum_file, "Wek", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		Print_array(cylindrical_ring_spectrum_file, "WDk", Correlation::cylindrical_ring_dissk1, Correlation::cylindrical_ring_dissk2);
		
		Correlation::Compute_cylindrical_ring_spectrum(T);
		Print_array(cylindrical_ring_spectrum_file, "Tek", Correlation::cylindrical_ring_ek);
		
		Correlation::Compute_cylindrical_ring_spectrum(U,W);
		Print_array(cylindrical_ring_spectrum_file, "(U.W)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		
		Correlation::Compute_cylindrical_ring_spectrum(U,T);
		Print_array(cylindrical_ring_spectrum_file, "(U.T)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		
		if (U.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(U);
			Print_array(cylindrical_ring_spectrum_file, "(Fv.v)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (W.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(W);
			Print_array(cylindrical_ring_spectrum_file, "(Fw.W)(k)", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (T.force_switch) {
			Correlation::Compute_force_cylindrical_ring_spectrum(T);
			Print_array(cylindrical_ring_spectrum_file, "(FT.T)(k)", Correlation::cylindrical_ring_ek);
		}
		
		if (global.program.helicity_switch) {
			Correlation::Compute_cylindrical_ring_spectrum_helicity(U);
			Print_array(cylindrical_ring_spectrum_file, "hk", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
			
			Correlation::Compute_cylindrical_ring_spectrum_helicity(W);
			Print_array(cylindrical_ring_spectrum_file, "Whk", Correlation::cylindrical_ring_ek1, Correlation::cylindrical_ring_ek2);
		}
		
		if (master)
			cylindrical_ring_spectrum_file.flush();
										
	}
}


//*********************************************************************************************
//*********************************************************************************************

void FluidIO::Output_Tk_shell_spectrum(FluidVF& U)
{
	if (global.time.now >= global.io.time.Tk_shell_spectrum_save_next) {
       
		if (global.mpi.master)
			Tk_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		
		Correlation::Compute_Tk_shell_spectrum(U);
		Print_array(Tk_spectrum_file, "UTk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		if (master)
			Tk_spectrum_file.flush();
    
		global.io.time.Tk_shell_spectrum_save_next += global.io.time.Tk_shell_spectrum_save_interval;
	}
		
}

void FluidIO::Output_Tk_shell_spectrum(FluidVF& U, FluidSF& T)
{
	
	if (global.time.now >= global.io.time.Tk_shell_spectrum_save_next) {
		
		if (global.mpi.master)
			Tk_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_shell_spectrum(U);
		Print_array(Tk_spectrum_file, "UTk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		Correlation::Compute_Tk_shell_spectrum(T);
		Print_array(Tk_spectrum_file, "TTk", Correlation::shell_ek);
		
		if (master)
			Tk_spectrum_file.flush();
		
		global.io.time.Tk_shell_spectrum_save_next += global.io.time.Tk_shell_spectrum_save_interval;
	}
}

void FluidIO::Output_Tk_shell_spectrum(FluidVF& U, FluidVF& W)
{
	if (global.time.now >= global.io.time.Tk_shell_spectrum_save_next) {
		
		if (global.mpi.master)
			Tk_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_shell_spectrum(U);
		Print_array(Tk_spectrum_file, "UTk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		Correlation::Compute_Tk_shell_spectrum(U);
		Print_array(Tk_spectrum_file, "WTk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		if (master)
			Tk_spectrum_file.flush();
		
		global.io.time.Tk_shell_spectrum_save_next += global.io.time.Tk_shell_spectrum_save_interval;
	}
}

void FluidIO::Output_Tk_shell_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{
	if (global.time.now >= global.io.time.Tk_shell_spectrum_save_next) {
		
		if (global.mpi.master)
			Tk_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_shell_spectrum(U);
		Print_array(Tk_spectrum_file, "UTk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		Correlation::Compute_Tk_shell_spectrum(W);
		Print_array(Tk_spectrum_file, "WTk", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
		
		Correlation::Compute_Tk_shell_spectrum(T);
		Print_array(Tk_spectrum_file, "TTk", Correlation::shell_ek);
		
		if (master)
			Tk_spectrum_file.flush();
		
		global.io.time.Tk_shell_spectrum_save_next += global.io.time.Tk_shell_spectrum_save_interval;
	}
}

//************

void FluidIO::Output_Tk_ring_spectrum(FluidVF& U)
{
	if ((global.spectrum.ring.turnon) && (global.time.now >= global.io.time.Tk_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_ring_spectrum(U);
		Print_array(Tk_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		if (master)
			Tk_ring_spectrum_file.flush();
		
		global.io.time.Tk_ring_spectrum_save_next += global.io.time.Tk_ring_spectrum_save_interval;
		
	}
}

void FluidIO::Output_Tk_ring_spectrum(FluidVF& U, FluidSF& T)
{
	
	if ((global.spectrum.ring.turnon) && (global.time.now >= global.io.time.Tk_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_ring_spectrum(U);
		Print_array(Tk_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_ring_spectrum(T);
		Print_array(Tk_ring_spectrum_file, "TTk", Correlation::ring_ek);
		
		if (master)
			Tk_ring_spectrum_file.flush();
		
		global.io.time.Tk_ring_spectrum_save_next += global.io.time.Tk_ring_spectrum_save_interval;
		
	}
}

void FluidIO::Output_Tk_ring_spectrum(FluidVF& U, FluidVF& W)
{
	
	if ((global.spectrum.ring.turnon) && (global.time.now >= global.io.time.Tk_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_ring_spectrum(U);
		Print_array(Tk_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_ring_spectrum(W);
		Print_array(Tk_ring_spectrum_file, "WTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		if (master)
			Tk_ring_spectrum_file.flush();
		
		global.io.time.Tk_ring_spectrum_save_next += global.io.time.Tk_ring_spectrum_save_interval;
		
	}
}

void FluidIO::Output_Tk_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	if ((global.spectrum.ring.turnon) && (global.time.now >= global.io.time.Tk_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_ring_spectrum(U);
		Print_array(Tk_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_ring_spectrum(W);
		Print_array(Tk_ring_spectrum_file, "WTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_ring_spectrum(T);
		Print_array(Tk_ring_spectrum_file, "TTk", Correlation::ring_ek);
		
		if (master)
			Tk_ring_spectrum_file.flush();
		
		global.io.time.Tk_ring_spectrum_save_next += global.io.time.Tk_ring_spectrum_save_interval;
		
	}
}



//************

void FluidIO::Output_Tk_cylindrical_ring_spectrum(FluidVF& U)
{
    
	if ((global.spectrum.cylindrical_ring.turnon) && (global.time.now >= global.io.time.Tk_cylindrical_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(U);
		Print_array(Tk_cylindrical_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		if (master)
			Tk_cylindrical_ring_spectrum_file.flush();
		
		global.io.time.Tk_cylindrical_ring_spectrum_save_next += global.io.time.Tk_cylindrical_ring_spectrum_save_interval;		
	}
}

void FluidIO::Output_Tk_cylindrical_ring_spectrum(FluidVF& U, FluidSF& T)
{
    
	if ((global.spectrum.cylindrical_ring.turnon) && (global.time.now >= global.io.time.Tk_cylindrical_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(U);
		Print_array(Tk_cylindrical_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(T);
		Print_array(Tk_cylindrical_ring_spectrum_file, "TTk", Correlation::ring_ek);
		
		if (master)
			Tk_cylindrical_ring_spectrum_file.flush();
		
		global.io.time.Tk_cylindrical_ring_spectrum_save_next += global.io.time.Tk_cylindrical_ring_spectrum_save_interval;
	}
    
}

void FluidIO::Output_Tk_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W)
{
    
	if ((global.spectrum.cylindrical_ring.turnon) && (global.time.now >= global.io.time.Tk_cylindrical_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(U);
		Print_array(Tk_cylindrical_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(W);
		Print_array(Tk_cylindrical_ring_spectrum_file, "WTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		if (master)
			Tk_cylindrical_ring_spectrum_file.flush();
		
		global.io.time.Tk_cylindrical_ring_spectrum_save_next += global.io.time.Tk_cylindrical_ring_spectrum_save_interval;
	}
    
}

void FluidIO::Output_Tk_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T)
{
    
	if ((global.spectrum.cylindrical_ring.turnon) && (global.time.now >= global.io.time.Tk_cylindrical_ring_spectrum_save_next)) {
        
		if (global.mpi.master)
			Tk_cylindrical_ring_spectrum_file << "%% Time = " << global.time.now << "\n \n";
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(U);
		Print_array(Tk_cylindrical_ring_spectrum_file, "UTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(W);
		Print_array(Tk_cylindrical_ring_spectrum_file, "WTk", Correlation::ring_ek1, Correlation::ring_ek2, Correlation::ring_ek3);
		
		Correlation::Compute_Tk_cylindrical_ring_spectrum(T);
		Print_array(Tk_cylindrical_ring_spectrum_file, "TTk", Correlation::ring_ek);
		
		if (master)
			Tk_cylindrical_ring_spectrum_file.flush();
		
		global.io.time.Tk_cylindrical_ring_spectrum_save_next += global.io.time.Tk_cylindrical_ring_spectrum_save_interval;
	}
    
}


//********
