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


/*! \file  compute_nlin_VF_MHD.cc
 * 
 * @brief  Compute nonlinear \f$ N_i^V = \mathcal{F} (D_j (V_j V_i - B_j B_i)) \f$ and 
 *		\f$ N_i^B = \mathcal{F} (D_j (V_j B_i - B_j V_i)) \f$.
 *
 *  We use \f$ \vec{Z}^\pm \f$ variables that decreases the no of required FFTs.
 *
 * @note \f$ N_i^V = \mathcal{F} (D_j (V_j V_i - B_j B_i)) 
 *				= \mathcal{F} (D_j (Z_j^- Z_i^+ + Z^+_j Z^-_i))\f$
 *
 * @note \f$ N_i^B = \mathcal{F} (D_j (V_j B_i - B_j V_i))
 *				= \mathcal{F} (D_j (Z_j^- Z_i^+ - Z^+_j Z^-_i))\f$
 *
 *	Steps <BR>
 *  # (V,B) <- \f$ \vec{Z}^+,  \vec{Z}^- \f$. <BR>
 *  # Inverse transform of \f$ \vec{Z}^\pm \f$ and store in real arrays. <BR>
 * # nlin \f$ N_i  = (Vr_i * Wr_i) \f$  <BR>
 * # nlin \f$ N_i  = D_i \mathcal{F}(Vr_i Wr_i) \f$ <BR>
 * # compute off-diagonal terms \f$ (Vr_i Wr_j) \f$ and put them in 
 *			\f$ \vec{Vr},  \vec{Wr}\f$. <BR>
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j W_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * # Go back (V,B) vars from \f$ (Z^+, Z^-) \f$.
 * # Dealias \f$ \vec{N} \f$ if DEALIASE switch is on. <BR>
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Sept 2008
 *
 * @bug   Transpose order needs to be tested.
 */


#include "Nlin.h"
#include "extern_var_incompress.h"

//*********************************************************************************************	

void Nlin_incompress::Compute_nlin(FluidVF& U, FluidVF& W) 
{
    // Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS") {
        U.cvf.Dealias();
        W.cvf.Dealias();
    }
    
	U.Inverse_transform();		// Vir = Inv_transform(Vi)
	W.Inverse_transform();		// Wir = Inv_transform(Wi)
								// Vi, Wi are unaffected in this operation.

    
    // Output real field
   if (!global.io.output_real_field_done) {
		fluidIO_incompress.Output_real_field(U, W);
		fluidIO_incompress.Output_real_field_slice(U, W);
		fluidIO_incompress.Output_field_r(U, W);
		fluidIO_incompress.Output_global_real_space(U, W);  // Only for Cheby basis
		fluidIO_incompress.Output_cout_real_space(U, W);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	}
    
    if (!global.time.dt_computation_done) {
		global.time.dt = U.Get_dt(W);
		global.time.now = global.time.now + global.time.dt;
		global.time.dt_computation_done = true;
	}
    
	// U=Zp=(U+B); B=Zm=(U-B);
	MHD::UB_to_Elsasser_real_field(U, W);	
	
	Compute_nlin_diag(U, W);
	
	Compute_nlin_offdiag(U, W);
	
	// U.nlin=(Zp.nlin+Zm.nlin)/2; B.nlin=(Zp.nlin-Zm.nlin)/2;
	MHD::Elsasser_to_UB_nlin(U, W);	
	
	if (!global.io.output_field_k_done) {
		fluidIO_incompress.Output_field_k(U, W);
        fluidIO_incompress.Output_Tk_shell_spectrum(U, W);
        fluidIO_incompress.Output_Tk_ring_spectrum(U, W);
        fluidIO_incompress.Output_Tk_cylindrical_ring_spectrum(U, W);
		global.io.output_field_k_done = true;
	}
    
  /*  U.nlin1 = 0;
    U.nlin2 = 0;
    U.nlin3 = 0;
    W.nlin1 = 0;W.nlin2 = 0;W.nlin3 = 0; */
 
}

//***************************** fn compute_nlin_VF_MHD ends ***********************************



