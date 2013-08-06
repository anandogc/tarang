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

/*! \file  compute_nlin_SF.cc
 * 
 * @brief  Compute nonlinear \f$ N_i^V = \mathcal{F} (D_j V_j V_i) \f$ and 
 *		\f$ N_i^\psi = \mathcal{F} (D_j V_j \psi) \f$.
 *
 *	Steps <BR>
 * # Inverse transform of V & T: \f$ \vec{Vr} = \mathcal{F}(\vec{V}) \f$  & 
 *			\f$ Fr = \mathcal{F}(F) \f$  <BR>
 * # T.nlin \f$ N^\psi = \mathcal{F} (D_j V_j \psi) \f$. <BR>
 * # nlin \f$ N_i  = (Vr_i)^2 \f$  <BR>
 * # nlin \f$ N_i  = \mathcal{F}(D_i (Vr_i)^2) \f$ <BR>
 * # compute off-diagonal terms and put them in \f$ \vec{Vr} \f$. <BR>
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j V_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * # Dealias \f$ \vec{N} \f$ if DEALIASE switch is on. <BR>
 *  
 * # RB Convection: If Pr==0: nlin(), else nlin(T).
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   Transpose order needs to tested
 */

  
#include "Nlin.h"
#include "IncIO.h"
#include "extern_var_incompress.h"

//**************************************************************************************************

void Nlin_incompress::Compute_nlin(FluidVF& U, FluidSF& T) 
{
	if (global.program.kind == "INC_SCALAR")
		Compute_nlin_scalar(U, T);
	
	else if (global.program.kind == "RBC") 
		Compute_nlin_RBC(U, T);
}

//*********************************************************************************************									

void Nlin_incompress::Compute_nlin_scalar(FluidVF& U, FluidSF& T) 
{
	
	// Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS") {
        U.cvf.Dealias();
        T.csf.Dealias();
    }
    
	
	U.Inverse_transform();	// Vir = Inv_transform(Vi)
	T.Inverse_transform();	// Fr = Inv_transform(F)

	global.io.real_space_field_available = true;
	
    // Output real field
    if (!global.io.output_real_field_done) {
		fluidIO_incompress.Output_real_field(U, T);
		fluidIO_incompress.Output_field_r(U, T);
		fluidIO_incompress.Output_global_real_space(U, T);  // Only for Cheby basis
		fluidIO_incompress.Output_cout_real_space(U, T);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	}
	
	if (!global.time.dt_computation_done) {
		global.time.dt = U.Get_dt(T);
		global.time.now = global.time.now + global.time.dt;
		global.time.dt_computation_done = true;
	}
	global.io.real_space_field_available = false;
	
	Compute_nlin_VxT(U, T);  // Computes only T.nlin = FT(Dj(ujT))
    
	Compute_nlin_diag(U);
	
	Compute_nlin_offdiag(U);
	
	if (!global.io.output_field_k_done) {
		fluidIO_incompress.Output_field_k(U, T);
        fluidIO_incompress.Output_Tk_shell_spectrum(U, T);
        fluidIO_incompress.Output_Tk_ring_spectrum(U, T);
        fluidIO_incompress.Output_Tk_cylindrical_ring_spectrum(U, T);
		global.io.output_field_k_done = true;
	} 
    
	

/*	global.io.real_space_field_available = false;
	
    // Output real field
    if (!global.io.output_real_field_done) {
	//	fluidIO_incompress.Output_real_field(U, T);
	//	fluidIO_incompress.Output_field_r(U, T);
		fluidIO_incompress.Output_global_real_space(U, T);  // Only for Cheby basis
		fluidIO_incompress.Output_cout_real_space(U, T);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	}
	
	global.time.now = global.time.now + global.time.dt;
	
    U.nlin1 = 0.0;
    U.nlin2 = 0.0;
    U.nlin3 = 0.0;
    T.nlin = 0.0; 
 */
}


//*********************************************************************************************
// For RB convection
//

void Nlin_incompress::Compute_nlin_RBC(FluidVF& U, FluidSF& T) 
{
	// For chebyshev
	if (global.program.basis_type.find("Ch") != string::npos)
	{
		if (global.PHYSICS.Pr_option == "PRZERO") {
			Compute_nlin(U);
			T.nlin = 0;
		}
		// if real-field is to be written, only V(r)
		
		else if (global.PHYSICS.Pr_option == "PRINFTY")	{
			// nlin1 = 0.0; nlin2 = 0.0; nlin3 = 0.0;
			
			// Dealiasing reqd before the real-space product.
			if (global.program.alias_option == "DEALIAS")
				T.csf.Dealias();
			
			T.Inverse_transform();	// Fr = Inv_transform(F)
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			
			if (!global.time.dt_computation_done) {
				global.time.dt = U.Get_dt();
				global.time.now = global.time.now + global.time.dt;
				global.time.dt_computation_done = true;
			}
			
			if (!global.io.output_real_field_done) {
				fluidIO_incompress.Output_real_field(U, T);
				fluidIO_incompress.Output_field_r(U, T);
				fluidIO_incompress.Output_global_real_space(U);  // Only for Cheby basis
				fluidIO_incompress.Output_cout_real_space(U);	 // Only for Cheby basis
				global.io.output_real_field_done = true;
			}
			
			Compute_nlin_VxT(U, T);		// For T.nlin only
			
			U.nlin1 = 0.0;
			U.nlin2 = 0.0;
			U.nlin3 = 0.0;
		}
		
		else
			Compute_nlin_scalar(U, T);
		
#ifdef GROSSMANN_LOHSE	
		if (!global.io.output_nlin_magnitude_done) {
			fluidIO_incompress.Output_field_k(U, T);
			
			global.io.output_field_k_done = true;

			// For RBC GL scaling analysis only.  Comment else
			fluidIO_incompress.misc_file << "rms(nlin): time, U.nlin1, U.nlin2, U.nlin3, T.nlin: " << global.time.now << " " << sqrt(sum(sqr(abs(U.nlin1)))) << " " << sqrt(sum(sqr(abs(U.nlin2)))) << " " << sqrt(sum(sqr(abs(U.nlin3)))) << " " << sqrt(sum(sqr(abs(T.nlin)))) << endl;
			
			global.io.output_nlin_magnitude_done = true;
		}
#endif
	}
	
	//**********
	// Other than ChFF basis
	else {
		if (global.PHYSICS.Pr_option == "PRZERO")
			Compute_nlin(U);   
			// if real-field is to be written, only V(r)
		
		else if (global.PHYSICS.Pr_option == "PRINFTY")	{
			// nlin1 = 0.0; nlin2 = 0.0; nlin3 = 0.0;
			
			// Dealiasing reqd before the real-space product.
			if (global.program.alias_option == "DEALIAS") 
				T.csf.Dealias();
		
			T.Inverse_transform();	// Fr = Inv_transform(F)
			U.Infinite_Prandtl_number_compute_velocity(T); // find U
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			
			if (!global.time.dt_computation_done) {
				global.time.dt = U.Get_dt();
				global.time.now = global.time.now + global.time.dt;
				global.time.dt_computation_done = true;
			}
			
			if (!global.io.output_real_field_done) {
				fluidIO_incompress.Output_real_field(U, T);
				fluidIO_incompress.Output_field_r(U, T);
				global.io.output_real_field_done = true;
			}
			
			Compute_nlin_VxT(U, T);		// For T.nlin only
		}
		
		else
			Compute_nlin_scalar(U, T);

#ifdef GROSSMANN_LOHSE
		if (!global.io.output_nlin_magnitude_done) {
			fluidIO_incompress.Output_field_k(U, T);
			
			global.io.output_field_k_done = true;
			
			// For RBC GL scaling analysis only.  Comment else
			fluidIO_incompress.misc_file << "rms(nlin): time, U.nlin1, U.nlin2, U.nlin3, T.nlin: " << global.time.now << " " << sqrt(sum(sqr(abs(U.nlin1)))) << " " << sqrt(sum(sqr(abs(U.nlin2)))) << " " << sqrt(sum(sqr(abs(U.nlin3)))) << " " << sqrt(sum(sqr(abs(T.nlin)))) << endl;
			
			global.io.output_nlin_magnitude_done = true;
		}
#endif
	}
	
}		

//***************************** fn compute_nlin_SF ends ***************************************


