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


/*! \file  compute_nlin.cc
 * 
 * @brief  Compute nonlinear \f$ N_i = \mathcal{F} (D_j V_j V_i) \f$.
 *
 *	Steps <BR>
 * # Inverse transform of  V: \f$ \vec{Vr} = \mathcal{F}(\vec{V}) \f$  <BR>
 * # nlin \f$ N_i  = (Vr_i)^2 \f$  <BR>
 * # nlin \f$ N_i  = D_i \mathcal{F}((Vr_i)^2) \f$ <BR>
 * # compute off-diagonal terms and put them in \f$ \vec{Vr} \f$. <BR>
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j V_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * # Dealias \f$ \vec{N} \f$ if DEALIASE switch is on. <BR>
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   TRANSPOSE order needs to be tested.
 */

#include "Nlin.h"
#include "extern_var_incompress.h"

//*********************************************************************************************
	    
void Nlin_incompress::Compute_nlin(FluidVF& U) 
{
	// Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS")
        U.cvf.Dealias();
    	

	U.Inverse_transform();  // Vir = Inv_transform(Vi)

    global.io.real_space_field_available = true;
  
    if (!global.time.dt_computation_done) {
		global.time.dt = U.Get_dt();
		global.time.now = global.time.now + global.time.dt;
		global.time.dt_computation_done = true;
	}
    
	if (!global.io.output_real_field_done) {
		fluidIO_incompress.Output_real_field(U);
		fluidIO_incompress.Output_real_field_slice(U);
		fluidIO_incompress.Output_field_r(U);
		fluidIO_incompress.Output_global_real_space(U);  // Only for Cheby basis
		fluidIO_incompress.Output_cout_real_space(U);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	} 
	
    global.io.real_space_field_available = false;
    
	Compute_nlin_diag(U);
	
	Compute_nlin_offdiag(U);

	
	/*if (master) cout << "nlin1 " << sum(abs(U.nlin1)) << endl;
	universal->Print_large_Fourier_elements(U.nlin1);
	
	if (master) cout << "nlin2 " << sum(abs(U.nlin2));
	universal->Print_large_Fourier_elements(U.nlin2);
	
	if (master) cout << "nlin3 " << sum(abs(U.nlin3)) ;
	universal->Print_large_Fourier_elements(U.nlin3); */
 
	
	
	
	if (!global.io.output_field_k_done) {
		fluidIO_incompress.Output_field_k(U);
        fluidIO_incompress.Output_Tk_shell_spectrum(U);
        fluidIO_incompress.Output_Tk_ring_spectrum(U);
        fluidIO_incompress.Output_Tk_cylindrical_ring_spectrum(U);
		global.io.output_field_k_done = true;
	}
}



//*********************************  End of compute_nlin.cc  **********************************




