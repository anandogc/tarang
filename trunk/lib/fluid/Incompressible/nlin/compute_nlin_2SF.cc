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
#include "extern_var_incompress.h"


//*********************************************************************************************									

void Nlin_incompress::Compute_nlin(FluidVF& U, FluidSF& T1, FluidSF& T2) 
{
    // Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS") {
        U.cvf.Dealias();
        T1.csf.Dealias();
        T2.csf.Dealias();
    }
    
	U.Inverse_transform();	// Vir = Inv_transform(Vi)
	T1.Inverse_transform();	// Fr = Inv_transform(F)
	T2.Inverse_transform();
	
	if (!global.time.dt_computation_done) {
		global.time.dt = U.Get_dt(T1);
		global.time.now = global.time.now + global.time.dt;
		global.time.dt_computation_done = true;
	}	
	
    // Output real field
/*    if (!global.io.output_real_field_done) {
		fluidIO_incompress.Output_real_field(U, T);
		fluidIO_incompress.Output_field_r(U, T);
//	fluidIO_incompress.Output_global_real_space(U, T1, T2);  // Only for Cheby basis
 //	fluidIO_incompress.Output_cout_real_space(U, T1, T2);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	} */
	
		Compute_nlin_VxT(U, T1);  // Computes only T.nlin = FT(Dj(ujT))
	
		Compute_nlin_VxT(U, T2); 
	
		Compute_nlin_diag(U);
	
		Compute_nlin_offdiag(U);
	
/*	if (!global.io.output_field_k_done) {
		fluidIO_incompress.Output_field_k(U, T);
		global.io.output_field_k_done = true;
	} */
}


//***************************** fn compute_nlin_SF ends ***************************************


