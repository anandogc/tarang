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


/*! \file  Output_ET.cc
 * 
 * @brief  Output flux, shell-to-shell, ring-to-ring (spherical & cylinderical), and
 *			Skpq.
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "IncIO.h"   


//*********************************************************************************************

void FluidIO_incompress::Output_flux(FluidVF& U, Pressure& P)
{
	
	if (master)
		flux_file << "%% Time = " << global.time.now << '\n'; 	

	energyTr->Compute_flux(U);
	Print_array(flux_file, "flux: U2U ", energyTr->flux_self);
	
	energyTr->Power_supply_within_sphere(U); 
	Print_array(flux_file, "sum(Fv.v)", energyTr->sphere_force_x_field);
	
/*	if (global.energy_transfer.helicity_flux_switch) {
		if (Ny > 1) // 3D
			energyTr->Compute_kinetic_helicity_flux(U);
		else if (Ny==1)	// 2D
			energyTr->Compute_enstrophy_flux(U);
		
		Print_array(flux_file, "flux: U2U_hk", energyTr->flux_hk);
	} */
    
    
    // Vpll to Vpll flux
    if (global.energy_transfer.Vpll_switch) {
        
        energyTr->Compute_flux_Vpll(U, P);
		
		Print_array(flux_file, "flux: Upll to Upll ", energyTr->flux_self);
        
		Print_array(flux_file, "sum(grad(p).v)", energyTr->sphere_force_x_field);
    }
	
	if (master)
		flux_file.flush();
} 


//*********************************************************************************************
// scalar

void FluidIO_incompress::Output_flux(FluidVF& U, FluidSF& T, Pressure& P)
{
	
	if (master)
		flux_file << "%% Time = " << global.time.now << '\n';
	
	energyTr->Compute_flux(U,T);
	Print_array(flux_file, "flux: U2U ", energyTr->flux_self);
	Print_array(flux_file, "flux: T2T ", energyTr->flux_SF);
	
	energyTr->Power_supply_within_sphere(U);
	Print_array(flux_file, "sum(Fv.v)", energyTr->sphere_force_x_field);
	
	energyTr->Power_supply_within_sphere(T);
	Print_array(flux_file, "sum(FT.T)", energyTr->sphere_force_x_field);
	
/*	if (global.energy_transfer.helicity_flux_switch) {
		if (Ny > 1) // 3D
			energyTr->Compute_kinetic_helicity_flux(U);
		else if (Ny==1)	// 2D
			energyTr->Compute_enstrophy_flux(U);
		
		Print_array(flux_file, "flux: U2U_hk", energyTr->flux_hk);
	} */
    
    
    // Vpll to Vpll flux
    if (global.energy_transfer.Vpll_switch) {
        
        energyTr->Compute_flux_Vpll(U, P);
		
		Print_array(flux_file, "flux: Upll to Upll ", energyTr->flux_self);
        
		Print_array(flux_file, "sum(grad(p).v)", energyTr->sphere_force_x_field);
    }
	
	if (master)
		flux_file.flush();
} 


//*********************************************************************************************
// Vector
  
void FluidIO_incompress::Output_flux(FluidVF& U, FluidVF& W, Pressure& P)
{
	
	if (master)
		flux_file << "%% Time = " << global.time.now << '\n';
	
	energyTr->Compute_flux(U, W);
	Print_array(flux_file, "flux: U2U ", energyTr->flux_self);
	Print_array(flux_file, "flux: flux_VF_Uin_Wout ", energyTr->flux_VF_Uin_Wout);
	Print_array(flux_file, "flux: flux_VF_Uin_Win ", energyTr->flux_VF_Uin_Win);
	Print_array(flux_file, "flux: flux_VF_Win_Wout ", energyTr->flux_VF_Win_Wout);
	Print_array(flux_file, "flux: flux_VF_Win_Uout ", energyTr->flux_VF_Win_Uout);
	Print_array(flux_file, "flux: flux_VF_Uout_Wout ", energyTr->flux_VF_Uout_Wout);
	Print_array(flux_file, "flux: flux_Elsasser_plus ", energyTr->flux_Elsasser_plus);
	Print_array(flux_file, "flux: flux_Elsasser_minus ", energyTr->flux_Elsasser_minus);

	
	energyTr->Power_supply_within_sphere(U);
	Print_array(flux_file, "sum(Fv.v)", energyTr->sphere_force_x_field);
	
	energyTr->Power_supply_within_sphere(W);
	Print_array(flux_file, "sum(Fw.w)", energyTr->sphere_force_x_field);
	
/*	if (global.energy_transfer.helicity_flux_switch) {
		if (Ny > 1) {// 3D
			energyTr->Compute_kinetic_helicity_flux(U);
		}
		
		else if (Ny==1)	// 2D
			energyTr->Compute_enstrophy_flux(U);
		
		Print_array(flux_file, "flux: U2U_hk", energyTr->flux_hk);
	}
    */
    
    // Vpll to Vpll flux
    if (global.energy_transfer.Vpll_switch) {
        
        energyTr->Compute_flux_Vpll(U, W, P);

		Print_array(flux_file, "pll2pll: U2U ", energyTr->flux_self);
		Print_array(flux_file, "pll2pll: flux_VF_Uin_Wout ", energyTr->flux_VF_Uin_Wout);
		Print_array(flux_file, "pll2pll: flux_VF_Uin_Win ", energyTr->flux_VF_Uin_Win);
		Print_array(flux_file, "pll2pll: flux_VF_Win_Wout ", energyTr->flux_VF_Win_Wout);
		Print_array(flux_file, "pll2pll: flux_VF_Win_Uout ", energyTr->flux_VF_Win_Uout);
		Print_array(flux_file, "pll2pll: flux_VF_Uout_Wout ", energyTr->flux_VF_Uout_Wout);
		Print_array(flux_file, "pll2pll: flux_Elsasser_plus ", energyTr->flux_Elsasser_plus);
		Print_array(flux_file, "pll2pll: flux_Elsasser_minus ", energyTr->flux_Elsasser_minus);
        
		Print_array(flux_file, "sum(grad(p).v)", energyTr->sphere_force_x_field);
    }
	
	if (master)
		flux_file.flush();
}

//*********************************************************************************************
// Magnetoconvection
  
void FluidIO_incompress::Output_flux(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P)
{
	
	if (master)
		flux_file << "%% Time = " << global.time.now << '\n';
	
	energyTr->Compute_flux(U, W, T);
	Print_array(flux_file, "flux: U2U ", energyTr->flux_self);
	Print_array(flux_file, "flux: flux_VF_Uin_Wout ", energyTr->flux_VF_Uin_Wout);
	Print_array(flux_file, "flux: flux_VF_Uin_Win ", energyTr->flux_VF_Uin_Win);
	Print_array(flux_file, "flux: flux_VF_Win_Wout ", energyTr->flux_VF_Win_Wout);
	Print_array(flux_file, "flux: flux_VF_Win_Uout ", energyTr->flux_VF_Win_Uout);
	Print_array(flux_file, "flux: flux_VF_Uout_Wout ", energyTr->flux_VF_Uout_Wout);
	Print_array(flux_file, "flux: flux_Elsasser_plus ", energyTr->flux_Elsasser_plus);
	Print_array(flux_file, "flux: flux_Elsasser_minus ", energyTr->flux_Elsasser_minus);
	Print_array(flux_file, "flux: T2T ", energyTr->flux_SF);
	
	energyTr->Power_supply_within_sphere(U);
	Print_array(flux_file, "sum(Fv.v)", energyTr->sphere_force_x_field);
	
	energyTr->Power_supply_within_sphere(W);
	Print_array(flux_file, "sum(Fw.w)", energyTr->sphere_force_x_field);
	
	energyTr->Power_supply_within_sphere(T);
	Print_array(flux_file, "sum(FT.T)", energyTr->sphere_force_x_field);
	
    
    // Vpll to Vpll flux
    if (global.energy_transfer.Vpll_switch) {
        
        cerr << "Vpll to Vpll flux for Output_flux(U,W,T) Not implemented presently " << endl;
    }
	
	if (master)
		flux_file.flush();
}


//*********************************************************************************************
//*********************************************************************************************
// Shell-to-shell

void FluidIO_incompress::Output_shell_to_shell(FluidVF& U, Pressure& P)
{
	static Range ra1(1,global.energy_transfer.shell_to_shell.no_shells-1);
	static Range ra2(1,global.energy_transfer.shell_to_shell.no_shells);
	
	if (master)
		shell_to_shell_file << "%% Time = " << global.time.now << '\n';
	
	energyTr->Compute_shell_tr(U);
	Print_array(shell_to_shell_file, "shell_to_shell: U2U ", energyTr->shelltoshell_self(ra1,ra2));
	
    
    // Vpll to Vpll shell-to-shell
    if (global.energy_transfer.Vpll_switch) {
        energyTr->Compute_shell_tr_Vpll(U, P);
		Print_array(shell_to_shell_file, "pll2pll shell_to_shell: U2U ", energyTr->shelltoshell_self(ra1,ra2));
    }
	
	if (master)
		shell_to_shell_file.flush();
}
 
//*********************************************************************************************   
// scalar


void FluidIO_incompress::Output_shell_to_shell(FluidVF& U, FluidSF& T, Pressure& P)
{
	
	static Range ra1(1,global.energy_transfer.shell_to_shell.no_shells-1);
	static Range ra2(1,global.energy_transfer.shell_to_shell.no_shells);
	
	if (master)
		shell_to_shell_file << "%% Time = " << global.time.now << '\n';
	
	energyTr->Compute_shell_tr(U, T);
	Print_array(shell_to_shell_file, "shell_to_shell: U2U ", energyTr->shelltoshell_self(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: T2T ", energyTr->shelltoshell_SF(ra1,ra2));
	
    
    // Vpll to Vpll shell-to-shell
    if (global.energy_transfer.Vpll_switch) {
        energyTr->Compute_shell_tr_Vpll(U, P);
		Print_array(shell_to_shell_file, "pll2pll shell_to_shell: U2U ", energyTr->shelltoshell_self(ra1,ra2));
    }
	
	if (master)
		shell_to_shell_file.flush();
} 


//*********************************************************************************************
// MHD

void FluidIO_incompress::Output_shell_to_shell(FluidVF& U, FluidVF& W, Pressure& P)
{

	
	static Range ra1(1,global.energy_transfer.shell_to_shell.no_shells-1);
	static Range ra2(1,global.energy_transfer.shell_to_shell.no_shells);
	
	if (master)
		shell_to_shell_file << "%% Time = " << global.time.now << '\n';
	
	energyTr->Compute_shell_tr(U, W);
	Print_array(shell_to_shell_file, "shell_to_shell: U2U ", energyTr->shelltoshell_self(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_VF_WtoW ", energyTr->shelltoshell_VF_WtoW(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_VF_UtoW ", energyTr->shelltoshell_VF_UtoW(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_Elsasser_plus ", energyTr->shelltoshell_Elsasser_plus(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_Elsasser_minus ", energyTr->shelltoshell_Elsasser_minus(ra1,ra2));
	
    
    // Vpll to Vpll shell-to-shell
    if (global.energy_transfer.Vpll_switch) {
        energyTr->Compute_shell_tr_Vpll(U, W, P);
		Print_array(shell_to_shell_file, "pll2pll shell_to_shell: U2U ", energyTr->shelltoshell_self(ra1,ra2));
		Print_array(shell_to_shell_file, "pll2pll shell_to_shell: shelltoshell_VF_WtoW ", energyTr->shelltoshell_VF_WtoW(ra1,ra2));
		Print_array(shell_to_shell_file, "pll2pll shell_to_shell: shelltoshell_VF_UtoW ", energyTr->shelltoshell_VF_UtoW(ra1,ra2));
		Print_array(shell_to_shell_file, "pll2pll shell_to_shell: shelltoshell_Elsasser_plus ", energyTr->shelltoshell_Elsasser_plus(ra1,ra2));
		Print_array(shell_to_shell_file, "pll2pll shell_to_shell: shelltoshell_Elsasser_minus ", energyTr->shelltoshell_Elsasser_minus(ra1,ra2));
    }
	
	
	DP B0mag = W.Get_mag_V0();
	
	if (B0mag > MYEPS) {
		energyTr->Compute_shell_ET_B0(U, W);
		Print_array(shell_to_shell_file, "b to u due to B0 ", energyTr->energy_tr_shell_B0(ra2));
	}
    
	if (master)
		shell_to_shell_file.flush();
				
} 


//*********************************************************************************************
//
//  MHD+Scalar

void FluidIO_incompress::Output_shell_to_shell(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P)
{

	static Range ra1(1,global.energy_transfer.shell_to_shell.no_shells-1);
	static Range ra2(1,global.energy_transfer.shell_to_shell.no_shells);
	
	if (master)
		shell_to_shell_file << "%% Time = " << global.time.now << '\n';
	
	energyTr->Compute_shell_tr(U, W, T);
	Print_array(shell_to_shell_file, "shell_to_shell: U2U ", energyTr->shelltoshell_self(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_VF_WtoW ", energyTr->shelltoshell_VF_WtoW(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_VF_UtoW ", energyTr->shelltoshell_VF_UtoW(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_Elsasser_plus ", energyTr->shelltoshell_Elsasser_plus(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: shelltoshell_Elsasser_minus ", energyTr->shelltoshell_Elsasser_minus(ra1,ra2));
	Print_array(shell_to_shell_file, "shell_to_shell: T2T ", energyTr->shelltoshell_SF(ra1,ra2));
	
	DP B0mag = W.Get_mag_V0();
	
	if (B0mag > MYEPS) {
		energyTr->Compute_shell_ET_B0(U, W);
		Print_array(shell_to_shell_file, "b to u due to B0 ", energyTr->energy_tr_shell_B0(ra2));
	}
	
    
    // Vpll to Vpll shell-to-shell
    if (global.energy_transfer.Vpll_switch) {
        cerr << "Vpll to Vpll flux for Output_shell2shell(U,W,T) Not implemented presently " << endl;
    }
    
	if (master)
		shell_to_shell_file.flush();
} 


//********************************************************************************************* 

// ring-to-ring (spherical)

void FluidIO_incompress::Output_ring_to_ring(FluidVF& U, Pressure& P)
{
	if (global.energy_transfer.ring_to_ring.turnon) {
		
		static Range ra1(1, global.energy_transfer.ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.ring_to_ring.no_sectors);
		
		if (master)
			ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_ring_tr(U);
		Print_array(ring_to_ring_file, "ring_to_ring: U2U ", energyTr->ring_to_ring_self(ra1,ra2,ra1,ra2));
		
		energyTr->Power_supply_ring(U);
		Print_array(ring_to_ring_file, "sum(Fv.v)", energyTr->ring_force_x_field(ra1,ra2));

		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			energyTr->Compute_ring_tr_Vpll(U, P);
			Print_array(ring_to_ring_file, "pll2pll ring_to_ring: U2U ", energyTr->ring_to_ring_self(ra1,ra2,ra1,ra2));
			
			Print_array(ring_to_ring_file, "grad(p)*Vpll", energyTr->ring_force_x_field(ra1,ra2));
		}
		
		if (master)
			ring_to_ring_file.flush();
        
	}
}  

//*********************************************************************************************
// scalar


void FluidIO_incompress::Output_ring_to_ring(FluidVF& U, FluidSF& T, Pressure& P)
{
	if (global.program.kind == "INC_SCALAR")
		Output_ring_to_ring_scalar(U, T, P);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_ring_to_ring_RBC(U, T, P);
}	

 
void FluidIO_incompress::Output_ring_to_ring_scalar(FluidVF& U, FluidSF& T, Pressure& P)
{
	if (global.energy_transfer.ring_to_ring.turnon  == 1) {
		static Range ra1(1, global.energy_transfer.ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.ring_to_ring.no_sectors);
		
		if (master)
			ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_ring_tr(U, T);
		Print_array(ring_to_ring_file, "ring_to_ring: U2U ", energyTr->ring_to_ring_self(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: T2T ", energyTr->ring_to_ring_SF(ra1,ra2,ra1,ra2));

		
		energyTr->Power_supply_ring(U);
		Print_array(ring_to_ring_file, "sum(Fv.v)", energyTr->ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_ring(T);
		Print_array(ring_to_ring_file, "sum(FT.T)", energyTr->ring_force_x_field(ra1,ra2));
		
		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			energyTr->Compute_ring_tr_Vpll(U, P);
			Print_array(ring_to_ring_file, "pll2pll ring_to_ring: U2U ", energyTr->ring_to_ring_self(ra1,ra2,ra1,ra2));
			
			Print_array(ring_to_ring_file, "grad(p)*Vpll", energyTr->ring_force_x_field(ra1,ra2));
		}
		
		if (master)
			ring_to_ring_file.flush();
	}	
}

//	RB Convection	//

void FluidIO_incompress::Output_ring_to_ring_RBC(FluidVF& U, FluidSF& T, Pressure& P)
{
	if (global.PHYSICS.Pr_option == "PRZERO")
		Output_ring_to_ring(U, P);
	
	else
		Output_ring_to_ring_scalar(U, T, P);
}

//*********************************************************************************************
void FluidIO_incompress::Output_ring_to_ring(FluidVF& U, FluidVF& W, Pressure& P)
{
	if (global.energy_transfer.ring_to_ring.turnon) {
		
		static Range ra1(1, global.energy_transfer.ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.ring_to_ring.no_sectors);
		
		if (master)
			ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_ring_tr(U, W);
		Print_array(ring_to_ring_file, "ring_to_ring: U2U ", energyTr->ring_to_ring_self(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: W2W ", energyTr->ring_to_ring_VF_WtoW(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: U2W ", energyTr->ring_to_ring_VF_UtoW(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: Z+ ", energyTr->ring_to_ring_Elsasser_plus(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: Z- ", energyTr->ring_to_ring_Elsasser_minus(ra1,ra2,ra1,ra2));

		
		
		energyTr->Power_supply_ring(U);
		Print_array(ring_to_ring_file, "sum(Fv.v)", energyTr->ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_ring(W);
		Print_array(ring_to_ring_file, "sum(Fw.w)", energyTr->ring_force_x_field(ra1,ra2));
		
		
		DP B0mag = W.Get_mag_V0();
		
		if (B0mag > MYEPS) {
			energyTr->Compute_ring_ET_B0(U, W);
			Print_array(ring_to_ring_file, "b to u due to B0 ", energyTr->energy_tr_ring_B0(ra1,ra2));
		}
		
		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			energyTr->Compute_ring_tr_Vpll(U, P);
			Print_array(ring_to_ring_file, "pll2pll ring_to_ring: U2U ", energyTr->ring_to_ring_self(ra1,ra2,ra1,ra2));
			Print_array(ring_to_ring_file, "pll2pll ring_to_ring: W2W ", energyTr->ring_to_ring_VF_WtoW(ra1,ra2,ra1,ra2));
			Print_array(ring_to_ring_file, "pll2pll ring_to_ring: U2W ", energyTr->ring_to_ring_VF_UtoW(ra1,ra2,ra1,ra2));
			Print_array(ring_to_ring_file, "pll2pll ring_to_ring: Z+ ", energyTr->ring_to_ring_Elsasser_plus(ra1,ra2,ra1,ra2));
			Print_array(ring_to_ring_file, "pll2pll ring_to_ring: Z- ", energyTr->ring_to_ring_Elsasser_minus(ra1,ra2,ra1,ra2));
			
			Print_array(ring_to_ring_file, "grad(p)*Vpll", energyTr->ring_force_x_field(ra1,ra2));
		}
		
		if (master)
			ring_to_ring_file.flush();
	}	
} 


//*********************************************************************************************
//

void FluidIO_incompress::Output_ring_to_ring(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P)
{

	if (global.energy_transfer.ring_to_ring.turnon  == 1) {
		static Range ra1(1, global.energy_transfer.ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.ring_to_ring.no_sectors);
		
		if (master)
			ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_ring_tr(U, W, T);
		Print_array(ring_to_ring_file, "ring_to_ring: U2U ", energyTr->ring_to_ring_self(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: W2W ", energyTr->ring_to_ring_VF_WtoW(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: U2W ", energyTr->ring_to_ring_VF_UtoW(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: Z+ ", energyTr->ring_to_ring_Elsasser_plus(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: Z- ", energyTr->ring_to_ring_Elsasser_minus(ra1,ra2,ra1,ra2));
		Print_array(ring_to_ring_file, "ring_to_ring: T2T ", energyTr->ring_to_ring_SF(ra1,ra2,ra1,ra2));
		
		
		energyTr->Power_supply_ring(U);
		Print_array(ring_to_ring_file, "sum(Fv.v)", energyTr->ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_ring(W);
		Print_array(ring_to_ring_file, "sum(Fw.w)", energyTr->ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_ring(T);
		Print_array(ring_to_ring_file, "sum(FT.T)", energyTr->ring_force_x_field(ra1,ra2));
		
		
		DP B0mag = W.Get_mag_V0();
		
		if (B0mag > MYEPS) {
			energyTr->Compute_ring_ET_B0(U, W);
			Print_array(ring_to_ring_file, "b to u due to B0 ", energyTr->energy_tr_ring_B0(ra1,ra2));
		}
		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			cerr << "Vpll to Vpll for Output_ring2ring(U,W,T) Not implemented presently " << '\n';
		}
		
		if (master)
			ring_to_ring_file.flush();
	
	}										
} 


//*********************************************************************************************
//*********************************************************************************************
// ring-to-ring (cylinderical)

void FluidIO_incompress::Output_cylindrical_ring_to_ring(FluidVF& U, Pressure& P)
{

	if (global.energy_transfer.cylindrical_ring_to_ring.turnon == 1) {
		
		static Range ra1(1, global.energy_transfer.cylindrical_ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.cylindrical_ring_to_ring.no_slabs);
		
		if (master)
			cylindrical_ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_cylindrical_ring_tr(U);
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: U2U ", energyTr->cylindrical_ring_to_ring_self(ra1,ra2,ra1,ra2));
		
		energyTr->Power_supply_cylindrical_ring(U);
		Print_array(cylindrical_ring_to_ring_file, "sum(Fv.v)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			energyTr->Compute_cylindrical_ring_tr_Vpll(U, P);
			Print_array(cylindrical_ring_to_ring_file, "pll2pll cyl_ring_to_ring: U2U ", energyTr->cylindrical_ring_to_ring_self(ra1,ra2,ra1,ra2));
			
			Print_array(cylindrical_ring_to_ring_file, "grad(p)*Vpll", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		}
		
		if (master)
			cylindrical_ring_to_ring_file.flush();
	}
	
}

//*********************************************************************************************
// scalar

void FluidIO_incompress::Output_cylindrical_ring_to_ring(FluidVF& U, FluidSF& T, Pressure& P)
{
	if (global.program.kind == "INC_SCALAR")
		Output_cylindrical_ring_to_ring_scalar(U, T, P);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Output_cylindrical_ring_to_ring_RBC(U, T, P);
}


 
void FluidIO_incompress::Output_cylindrical_ring_to_ring_scalar(FluidVF& U, FluidSF& T, Pressure& P)
{
	if (global.energy_transfer.cylindrical_ring_to_ring.turnon == 1) {
		
		static Range ra1(1, global.energy_transfer.cylindrical_ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.cylindrical_ring_to_ring.no_slabs);
		
		if (master)
			cylindrical_ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_cylindrical_ring_tr(U, T);
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: U2U ", energyTr->cylindrical_ring_to_ring_self(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: T2T ", energyTr->cylindrical_ring_to_ring_SF(ra1,ra2,ra1,ra2));
		
		energyTr->Power_supply_cylindrical_ring(U);
		Print_array(cylindrical_ring_to_ring_file, "sum(Fv.v)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_cylindrical_ring(T);
		Print_array(cylindrical_ring_to_ring_file, "sum(FT.T)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			energyTr->Compute_cylindrical_ring_tr_Vpll(U, P);
			Print_array(cylindrical_ring_to_ring_file, "pll2pll cyl_ring_to_ring: U2U ", energyTr->cylindrical_ring_to_ring_self(ra1,ra2,ra1,ra2));
			
			Print_array(cylindrical_ring_to_ring_file, "grad(p)*Vpll", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		}
		
		if (master)
			cylindrical_ring_to_ring_file.flush();
	}	
} 

//	RB Convection	//

void FluidIO_incompress::Output_cylindrical_ring_to_ring_RBC(FluidVF& U, FluidSF& T, Pressure& P)
{
	if (global.PHYSICS.Pr_option == "PRZERO")
		Output_cylindrical_ring_to_ring(U, P);
	
	else
		Output_cylindrical_ring_to_ring_scalar(U, T, P);
}

//*********************************************************************************************
void FluidIO_incompress::Output_cylindrical_ring_to_ring(FluidVF& U, FluidVF& W, Pressure& P)
{
	if (global.energy_transfer.cylindrical_ring_to_ring.turnon == 1) {
		
		static Range ra1(1, global.energy_transfer.cylindrical_ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.cylindrical_ring_to_ring.no_slabs);
		
		if (master)
			cylindrical_ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_cylindrical_ring_tr(U, W);
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: U2U ", energyTr->cylindrical_ring_to_ring_self(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: W2W", energyTr->cylindrical_ring_to_ring_VF_WtoW(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: U2W", energyTr->cylindrical_ring_to_ring_VF_UtoW(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: Z+", energyTr->cylindrical_ring_to_ring_Elsasser_plus(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: Z-", energyTr->cylindrical_ring_to_ring_Elsasser_minus(ra1,ra2,ra1,ra2));

		
		
		energyTr->Power_supply_cylindrical_ring(U);
		Print_array(cylindrical_ring_to_ring_file, "sum(Fv.v)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_cylindrical_ring(W);
		Print_array(cylindrical_ring_to_ring_file, "sum(Fw.w)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			energyTr->Compute_cylindrical_ring_tr_Vpll(U, W, P);
			Print_array(cylindrical_ring_to_ring_file, "pll2pll cyl_ring_to_ring: U2U ", energyTr->cylindrical_ring_to_ring_self(ra1,ra2,ra1,ra2));
			
			Print_array(cylindrical_ring_to_ring_file, "grad(p)*Vpll", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		}
		
		if (master)
			cylindrical_ring_to_ring_file.flush();
	}			
} 

//*********************************************************************************************
//

void FluidIO_incompress::Output_cylindrical_ring_to_ring(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P)
{

	if (global.energy_transfer.cylindrical_ring_to_ring.turnon == 1) {
		
		static Range ra1(1, global.energy_transfer.cylindrical_ring_to_ring.no_shells-1);
		static Range ra2(1, global.energy_transfer.cylindrical_ring_to_ring.no_slabs);
		
		if (master)
			cylindrical_ring_to_ring_file << "%% Time = " << global.time.now << '\n';
		
		energyTr->Compute_cylindrical_ring_tr(U, W, T);
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: U2U ", energyTr->cylindrical_ring_to_ring_self(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: W2W", energyTr->cylindrical_ring_to_ring_VF_WtoW(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: U2W", energyTr->cylindrical_ring_to_ring_VF_UtoW(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: Z+", energyTr->cylindrical_ring_to_ring_Elsasser_plus(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: Z-", energyTr->cylindrical_ring_to_ring_Elsasser_minus(ra1,ra2,ra1,ra2));
		Print_array(cylindrical_ring_to_ring_file, "cyl_ring_to_ring: T2T ", energyTr->cylindrical_ring_to_ring_SF(ra1,ra2,ra1,ra2));
		
		energyTr->Power_supply_cylindrical_ring(U);
		Print_array(cylindrical_ring_to_ring_file, "sum(Fv.v)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_cylindrical_ring(W);
		Print_array(cylindrical_ring_to_ring_file, "sum(Fw.w)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		energyTr->Power_supply_cylindrical_ring(T);
		Print_array(cylindrical_ring_to_ring_file, "sum(FT.T)", energyTr->cylindrical_ring_force_x_field(ra1,ra2));
		
		
		// Vpll to Vpll shell-to-shell
		if (global.energy_transfer.Vpll_switch) {
			cerr << "Vpll to Vpll for Output_cyl_ring2ring(U,W,T) Not implemented presently " << '\n';
		}
		
		if (master)
			cylindrical_ring_to_ring_file.flush();
		
	}						
}



//*******************************  End of Output_ET.cc  ***************************************




