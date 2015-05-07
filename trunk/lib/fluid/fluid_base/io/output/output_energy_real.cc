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
// For Chebyshev basis fns

void FluidIO::Output_global_real_space(FluidVF& U)
{
	if ( (global.time.now >= global.io.time.global_save_next) && (basis_type.find("Ch") != string::npos) ){
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		if (master)
			global_file << global.time.now << " " << U.cvf.total_energy << "   " << U.cvf.total_E1 << " "<< U.cvf.total_E2 << " " << U.cvf.total_E3 << endl;
		
		global.io.time.global_save_next +=  global.io.time.global_save_interval;
	}
}

void FluidIO::Output_global_real_space(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Output_global_real_space_scalar(U, T);
	
	
	else if (global.program.kind == "RBC")
		Output_global_real_space_RBC(U, T);
}

void FluidIO::Output_global_real_space_scalar(FluidVF& U, FluidSF& T)
{
	if ( (global.time.now >= global.io.time.global_save_next) && (basis_type.find("Ch") != string::npos) ) {
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			T.Inverse_transform();	// Fr = Inv_transform(F)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		T.csf.total_energy = universal->Get_total_energy_real_space(T.rsf.Fr);
		
		if (master)
			global_file << global.time.now << " " << U.cvf.total_energy << " " << T.csf.total_energy <<  "   " << U.cvf.total_E1 << " "<< U.cvf.total_E2 << " " << U.cvf.total_E3 << endl;
		
		global.io.time.global_save_next +=  global.io.time.global_save_interval;
	}
}

void FluidIO::Output_global_real_space_RBC(FluidVF& U, FluidSF& T)
{
	if ( (global.time.now >= global.io.time.global_save_next) && (basis_type.find("Ch") != string::npos) ) {
	
		static Real nusselt_no;
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			T.Inverse_transform();	// Fr = Inv_transform(F)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		T.csf.total_energy = universal->Get_total_energy_real_space(T.rsf.Fr);
		
		nusselt_no = Correlation::Get_Nusselt_no(U, T);
		
		// total dissipation..
		
		
		if (master)
			global_file << global.time.now << " " << U.cvf.total_energy << " " << T.csf.total_energy <<  "   " << nusselt_no << "    "<< U.cvf.total_E1 << " "<< U.cvf.total_E2 << " " << U.cvf.total_E3 << endl;
		
		global.io.time.global_save_next +=  global.io.time.global_save_interval;
	}
}



void FluidIO::Output_global_real_space(FluidVF& U, FluidVF& W)
{
	if ( (global.time.now >= global.io.time.global_save_next) && (basis_type.find("Ch") != string::npos) ) {
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			W.Inverse_transform();	// Wir = Inv_transform(W)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		W.cvf.total_E1 = universal->Get_total_energy_real_space(W.rvf.V1r);
		W.cvf.total_E2 = universal->Get_total_energy_real_space(W.rvf.V2r);
		W.cvf.total_E3 = universal->Get_total_energy_real_space(W.rvf.V3r);
		
		W.cvf.total_energy = W.cvf.total_E1 + W.cvf.total_E2 + W.cvf.total_E3;
		
		if (master)
			global_file << global.time.now << " " << U.cvf.total_energy << " " << W.cvf.total_energy << "   " << U.cvf.total_E1 << " "<< U.cvf.total_E2 << " " << U.cvf.total_E3 << "   " << W.cvf.total_E1 << " "<< W.cvf.total_E2 << " " << W.cvf.total_E3 << endl;
		
		global.io.time.global_save_next +=  global.io.time.global_save_interval;
	}
}

void FluidIO::Output_global_real_space(FluidVF& U, FluidVF& W, FluidSF& T)
{
	if ( (global.time.now >= global.io.time.global_save_next) && (basis_type.find("Ch") != string::npos) ) {
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			W.Inverse_transform();	// Wir = Inv_transform(W)
			T.Inverse_transform();	// Fr = Inv_transform(F)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		W.cvf.total_E1 = universal->Get_total_energy_real_space(W.rvf.V1r);
		W.cvf.total_E2 = universal->Get_total_energy_real_space(W.rvf.V2r);
		W.cvf.total_E3 = universal->Get_total_energy_real_space(W.rvf.V3r);
		
		W.cvf.total_energy = W.cvf.total_E1 + W.cvf.total_E2 + W.cvf.total_E3;
		
		T.csf.total_energy = universal->Get_total_energy_real_space(T.rsf.Fr);
		
		if (master)
			global_file << global.time.now << " " << U.cvf.total_energy << " " << W.cvf.total_energy  << " " << T.csf.total_energy << "   " << U.cvf.total_E1 << " "<< U.cvf.total_E2 << " " << U.cvf.total_E3 << "   " << W.cvf.total_E1 << " "<< W.cvf.total_E2 << " " << W.cvf.total_E3 << endl;
		
		global.io.time.global_save_next +=  global.io.time.global_save_interval;
	}
}


//**********************

void FluidIO::Output_cout_real_space(FluidVF& U)
{
	if ( (global.time.now >= global.io.time.cout_save_next) && (basis_type.find("Ch") != string::npos) ) {
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		if (master)
			cout  << global.time.now << " "  << U.cvf.total_energy << endl;
		
		global.io.time.cout_save_next +=  global.io.time.cout_save_interval;
	}
}

void FluidIO::Output_cout_real_space(FluidVF& U, FluidSF& T)
{
	if ( (global.time.now >= global.io.time.cout_save_next) && (basis_type.find("Ch") != string::npos) ) {
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			T.Inverse_transform();	// Fr = Inv_transform(F)
		}
	
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		T.csf.total_energy = universal->Get_total_energy_real_space(T.rsf.Fr);
		
		if (master)
			cout << global.time.now << " " << U.cvf.total_energy << " " << T.csf.total_energy << endl;
		
		global.io.time.cout_save_next +=  global.io.time.cout_save_interval;
	}
}

void FluidIO::Output_cout_real_space(FluidVF& U, FluidVF& W)
{
	if ( (global.time.now >= global.io.time.cout_save_next) && (basis_type.find("Ch") != string::npos) ) {
	
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			W.Inverse_transform();	// Wir = Inv_transform(W)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		W.cvf.total_E1 = universal->Get_total_energy_real_space(W.rvf.V1r);
		W.cvf.total_E2 = universal->Get_total_energy_real_space(W.rvf.V2r);
		W.cvf.total_E3 = universal->Get_total_energy_real_space(W.rvf.V3r);
		
		W.cvf.total_energy = W.cvf.total_E1 + W.cvf.total_E2 + W.cvf.total_E3;
		
		if (master)
			cout << global.time.now << " " << U.cvf.total_energy << " " << W.cvf.total_energy << endl;
		
		global.io.time.cout_save_next +=  global.io.time.cout_save_interval;
	}
}

void FluidIO::Output_cout_real_space(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	if ( (global.time.now >= global.io.time.cout_save_next) && (basis_type.find("Ch") != string::npos) ) {
		
		if (!global.io.real_space_field_available) {
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			W.Inverse_transform();	// Wir = Inv_transform(W)
			T.Inverse_transform();	// Fr = Inv_transform(F)
		}
		
		U.cvf.total_E1 = universal->Get_total_energy_real_space(U.rvf.V1r);
		U.cvf.total_E2 = universal->Get_total_energy_real_space(U.rvf.V2r);
		U.cvf.total_E3 = universal->Get_total_energy_real_space(U.rvf.V3r);
		
		U.cvf.total_energy = U.cvf.total_E1 + U.cvf.total_E2 + U.cvf.total_E3;
		
		W.cvf.total_E1 = universal->Get_total_energy_real_space(W.rvf.V1r);
		W.cvf.total_E2 = universal->Get_total_energy_real_space(W.rvf.V2r);
		W.cvf.total_E3 = universal->Get_total_energy_real_space(W.rvf.V3r);
		
		W.cvf.total_energy = W.cvf.total_E1 + W.cvf.total_E2 + W.cvf.total_E3;
		
		T.csf.total_energy = universal->Get_total_energy_real_space(T.rsf.Fr);
		
		if (master)
			cout << global.time.now << " " << U.cvf.total_energy << " " << W.cvf.total_energy  << " " << T.csf.total_energy << endl;
		
		global.io.time.cout_save_next +=  global.io.time.cout_save_interval;
	}
}





 
//*******************************   End of output_energy.cc ***********************************


