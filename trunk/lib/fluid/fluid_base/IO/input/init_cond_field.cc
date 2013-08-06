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

/*! \file  init_cond_field.cc
 *
 * @brief   Read field vars from file field_in_file.
 * Design issue:  These functions are kept here so that we can use them for compressible as 
 *  well as incompressible fields.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 */


#include "FluidIO.h"


/**********************************************************************************************

						Input from a file: field_in_file

***********************************************************************************************/


void  FluidIO::Init_cond_complex_field(FluidVF& U)
{

	U.cvf.Read_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
}

//*********************************************************************************************

void  FluidIO::Init_cond_complex_field(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Init_cond_scalar_complex_field(U, T);

	else if (global.program.kind == "RBC")
		Init_cond_RBC_complex_field(U, T);
}


void  FluidIO::Init_cond_scalar_complex_field(FluidVF& U, FluidSF& T)
{
	U.cvf.Read_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	T.csf.Read_complex_field();
}


void  FluidIO::Init_cond_RBC_complex_field(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO") {
		U.cvf.Read_complex_field();
        
        if (global.io.input_vx_vy_switch && global.field.incompressible)
            universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
            // Vz(lx,ly,lz=0) already read; construct for lz>=1.
        
		U.Zero_Prandtl_number_compute_temperature(T);
	}
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
		T.csf.Read_complex_field();
		
		U.Infinite_Prandtl_number_compute_velocity(T);
	}
	else
		Init_cond_scalar_complex_field(U, T);
}

void  FluidIO::Init_cond_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	U.cvf.Read_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	T1.csf.Read_complex_field();
	T2.csf.Read_complex_field();
}


//
//
//*********************************************************************************************
void  FluidIO::Init_cond_complex_field(FluidVF& U, FluidVF& W)
{
	U.cvf.Read_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	W.cvf.Read_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(W.cvf.V1, W.cvf.V2, W.cvf.V3);
        // W.Vz(lx,ly,lz=0) already read; construct for lz>=1.
}

//
//*********************************************************************************************
void  FluidIO::Init_cond_complex_field(FluidVF& U, FluidVF& W, FluidSF& T)
{
	U.cvf.Read_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	W.cvf.Read_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(W.cvf.V1, W.cvf.V2, W.cvf.V3);
        // W.Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	T.csf.Read_complex_field();
}


/**********************************************************************************************

						Input from a file: field_in_file(N_in_reduced[])

***********************************************************************************************/


// Fluid
void  FluidIO::Init_cond_reduced_complex_field(FluidVF& U)
{
    U.cvf.Read_reduced_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
}

//*********************************************************************************************
// Passive scalar + RB convection


void  FluidIO::Init_cond_reduced_complex_field(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Init_cond_reduced_complex_field_scalar(U, T);

	else if (global.program.kind == "RBC")
		Init_cond_reduced_complex_field_RBC(U, T);
}


void  FluidIO::Init_cond_reduced_complex_field_scalar(FluidVF& U, FluidSF& T)
{
	U.cvf.Read_reduced_complex_field();	
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	T.csf.Read_reduced_complex_field();	
}

void  FluidIO::Init_cond_reduced_complex_field_RBC(FluidVF& U, FluidSF& T)
{

	if (global.PHYSICS.Pr_option == "PRZERO") {
		U.cvf.Read_reduced_complex_field();
        
        if (global.io.input_vx_vy_switch && global.field.incompressible)
            universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
            // Vz(lx,ly,lz=0) already read; construct for lz>=1.
        
		U.Zero_Prandtl_number_compute_temperature(T);
	}

	if (global.PHYSICS.Pr_option == "PRINFTY") {
		T.csf.Read_reduced_complex_field();
		U.Infinite_Prandtl_number_compute_velocity(T);
	}

	else
		Init_cond_reduced_complex_field_scalar(U, T);

		//	U.Zero_modes_RB_slip(T);
}

void  FluidIO::Init_cond_reduced_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	U.cvf.Read_reduced_complex_field();	
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	T1.csf.Read_reduced_complex_field();
	T2.csf.Read_reduced_complex_field();
}

//*********************************************************************************************
// MHD

void  FluidIO::Init_cond_reduced_complex_field(FluidVF& U, FluidVF& W)
{
	U.cvf.Read_reduced_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	W.cvf.Read_reduced_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(W.cvf.V1, W.cvf.V2, W.cvf.V3);
        // W.Vz(lx,ly,lz=0) already read; construct for lz>=1.
}


//*********************************************************************************************
// Magnetoconvection

void  FluidIO::Init_cond_reduced_complex_field(FluidVF& U, FluidVF& W, FluidSF& T)
{
	U.cvf.Read_reduced_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(U.cvf.V1, U.cvf.V2, U.cvf.V3);
        // Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	W.cvf.Read_reduced_complex_field();
    
    if (global.io.input_vx_vy_switch && global.field.incompressible)
        universal->Fill_Vz(W.cvf.V1, W.cvf.V2, W.cvf.V3);
        // W.Vz(lx,ly,lz=0) already read; construct for lz>=1.
    
	T.csf.Read_reduced_complex_field();
}


/**********************************************************************************************
 
 Input real-field from a file: field_in_file
 
 ***********************************************************************************************/


void  FluidIO::Init_cond_real_field(FluidVF& U)
{
	U.rvf.Read_real_field();
	U.Forward_transform();
}



void  FluidIO::Init_cond_real_field(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Init_cond_real_field_scalar(U, T);
	
	else if (global.program.kind == "RBC")
		Init_cond_real_field_RBC(U, T);
}


void  FluidIO::Init_cond_real_field_scalar(FluidVF& U, FluidSF& T)
{

	U.rvf.Read_real_field();
	T.rsf.Read_real_field();
	
	U.Forward_transform();
	T.Forward_transform();
}

void  FluidIO::Init_cond_real_field_RBC(FluidVF& U, FluidSF& T)
{
	if (global.PHYSICS.Pr_option == "PRZERO") {
		Init_cond_real_field(U);
		U.Zero_Prandtl_number_compute_temperature(T);
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
		T.rsf.Read_real_field();
		T.Forward_transform();
		U.Infinite_Prandtl_number_compute_velocity(T);
	}
	
	else
		Init_cond_real_field_scalar(U, T);
}


void  FluidIO::Init_cond_real_field(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	
	U.rvf.Read_real_field();
	T1.rsf.Read_real_field();
	T2.rsf.Read_real_field();
	
	U.Forward_transform();
	T1.Forward_transform();
	T2.Forward_transform();
}

//
//
//*********************************************************************************************
void  FluidIO::Init_cond_real_field(FluidVF& U, FluidVF& W)
{
	U.rvf.Read_real_field();
	W.rvf.Read_real_field();
	
	U.Forward_transform();
	W.Forward_transform();
}

//
//*********************************************************************************************
void  FluidIO::Init_cond_real_field(FluidVF& U, FluidVF& W, FluidSF& T)
{
	U.rvf.Read_real_field();
	W.rvf.Read_real_field();
	T.rsf.Read_real_field();
	
	U.Forward_transform();
	W.Forward_transform();
	T.Forward_transform();
}



//******************************** End of Init_cond_field.cc **********************************

