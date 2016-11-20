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


/*! \file  compute_nlin_VF_SF.cc
 * 
 * @brief  Compute nonlinear \f$ N_i^V = \mathcal{F} (D_j (V_j V_i - B_j B_i)) \f$ and 
 *		\f$ N_i^B = \mathcal{F} (D_j (V_j B_i - B_j V_i)) \f$,
 *		\f$ N_i^\psi = \mathcal{F} (D_j V_j \psi) \f$.
 *
 *  We use compute_nlin(W) to compute MHD nlin terms, and use the method of 
 *	compute_nlin(T) to compute the temperature nlin term.
 *
 *	@sa compute_nlin_VF_MHD.cc
 *  @sa compute_nlin_SF.cc
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @note  The basis functions of V and B are the same in the present implementation.
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Sept 2008
 *
 * @bug  Transpose order
 */

#include "Nlin.h"
#include "extern_var_incompress.h"


//*********************************************************************************************
// C is the normalized density of gas
void Nlin_incompress::Compute_nlin(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	if (global.program.kind == "MHD_ASTRO_INCOMPRESS")
		Compute_nlin_MHD_ASTRO_INCOMPRESS(U, W,T, C);
	
	// Rayleigh-Benard magnetoconvection 
	else if (global.program.kind == "MHD_ASTRO_INCOMPRESS_ALL")
		Compute_nlin_MHD_ASTRO_INCOMPRESS_ALL(U, W,T,C);
}


void Nlin_incompress::Compute_nlin_MHD_ASTRO_INCOMPRESS(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C) 
{
    // Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS") {
        U.cvf.Dealias();
        W.cvf.Dealias();
        T.csf.Dealias();
		C.csf.Dealias();
    }
    
	T.Inverse_transform();		// Fr = Inv_transform(F)
	C.Inverse_transform();
	U.Inverse_transform();		// Vir = Inv_transform(Vi)
	W.Inverse_transform();		// Wir = Inv_transform(Wi)
								// F, Vi, Wi are unaffected in this operation.
    
    // Output real field
    if (!global.io.output_real_field_done) {
	//	fluidIO_incompress.Output_real_field(U, W, T, C);
	//	fluidIO_incompress.Output_field_r(U, W, T, C);
		//	fluidIO_incompress.Output_global_real_space(U, W, T, C);  // Only for Cheby basis
		//	fluidIO_incompress.Output_cout_real_space(U, W, T, C);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	}
    
    if (!global.time.dt_computation_done) {
		global.time.dt = U.Get_dt(W);
		global.time.now = global.time.now + global.time.dt;
		global.time.dt_computation_done = true;
	}
								
	Compute_nlin_VxT(U, T);
	Compute_nlin_VxT(U, C);
	
	// U=Zp=(U+B); B=Zm=(U-B);
	MHD::UB_to_Elsasser_real_field(U, W);	
	
	Compute_nlin_diag(U, W);
	
	Compute_nlin_offdiag(U, W);
	
	// U.nlin=(Zp.nlin+Zm.nlin)/2; B.nlin=(Zp.nlin-Zm.nlin)/2;
	MHD::Elsasser_to_UB_nlin(U, W);
	
	if (!global.io.output_field_k_done) {
	//	fluidIO_incompress.Output_field_k(U, W, T, C);
    //    fluidIO_incompress.Output_Tk_shell_spectrum(U, W, T, C);
	//  fluidIO_incompress.Output_Tk_ring_spectrum(U, W, T, C);
    //    fluidIO_incompress.Output_Tk_cylindrical_ring_spectrum(U, W, T, C);
		global.io.output_field_k_done = true;
	}
}

//*********************************************************************************************

void Nlin_incompress::Compute_nlin_MHD_ASTRO_INCOMPRESS_ALL(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C) 
{
	
	Compute_nlin_MHD_ASTRO_INCOMPRESS(U, W, T, C);
}

//***************************** fn compute_nlin_VF_SF ends ************************************


