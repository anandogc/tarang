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


//*********************************************************************************************	
// U.nlin_i = FT(Dj (Vj Wi)
void Nlin_incompress::Compute_nlin(FluidVF& U, PlainFluidVF& W) 
{
    // Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS") {
        U.cvf.Dealias();
        W.cvf.Dealias();
    }
    
	U.Inverse_transform();		// Vir = Inv_transform(Vi)
	W.Inverse_transform();		// Wir = Inv_transform(Wi)
								// Vi, Wi are unaffected in this operation.

	
	Compute_nlin_diag(U, W);
	
	Compute_nlin_offdiag(U, W);	// U.nlin_i = FT(Dj (Vj Wi)
}

// U.nlin_1 = FT(Dj (Vj W1)
void Nlin_incompress::Compute_nlin_first_component(FluidVF& U, PlainFluidVF& W) 
{
	
    // Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS") {
        U.cvf.Dealias();
        universal->Dealias(W.cvf.V1);
    }
    
	U.Inverse_transform();		// Vir = Inv_transform(Vi)
	W.Inverse_transform_first_component();		// W1r = Inv_transform(W1)
								// Vi, Wi are unaffected in this operation.
	
	Compute_nlin_first_component_product(U, W); // U.nlin_1 = FT(Dj (Vj W1)
}

//***************************** fn compute_nlin_VF_MHD ends ***********************************



