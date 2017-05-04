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

/*! \file  IncVF.h
 * 
 * @brief  Class declaration of IncVF, Incompressible Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   No known bugs
 */
 
//*********************************************************************************************

#ifndef _H_NLIN
#define _H_NLIN

#include "fluid_base.h"

//! @brief Incompressible vector field IncVF 
/*!
 *  Inherits CVF that contains the complex vector field. <BR>
 *  RVF that contains real vector field, typically Inverse tranform of CVF. <BR>
 *  NLIN that contains the nonlinear term \f$ N_i = \mathcal{F} (D_j V_j V_i) \f$.<BR>
 *  EnergyTr that contains energy transfer vars like flux, shell-to-shell transfers. <BR>
 * 
 *	Compute nonlinear terms <BR>
 *  Compute energy transfer functions: <BR>
 *	Isotropic: flux, shell-to-shell <BR>
 *  Anisotropic: ring-to-ring in spherical shell and in cylinderical shells.
 *
 *	@sa IncSF.h
 *	@sa Nlin.h
 *  @sa EnergyTr.h
 */
 
//*********************************************************************************************	

class Nlin_incompress
{
	/// For computing FT[nabla.(V V)]

public:
	static void Compute_nlin_diag(FluidVF& U);						///  Vr[i] -> Vr[i]^2 stored in nlin[i]
	static void Compute_nlin_offdiag(FluidVF& U);					/// Computes Vr[i]*Vr[j] 
	
	static void Compute_nlin_diag(FluidVF& U, FluidVF& W);			///  Vr[i] -> Vr[i]^2 stored in nlin[i]
	static void Compute_nlin_offdiag(FluidVF& U, FluidVF& W);		/// Computes Vr[i]*Vr[j]  
	
	static void Compute_nlin_VxT(FluidVF& U, FluidSF& T);
	
	static void Compute_nlin_vector_potential(FluidVF& U, PlainFluidVF& W);
	
	static void Compute_nlin(FluidVF& U);							/// nlin[i] = FT[v.grad(v)][i]
	static void Compute_nlin(FluidVF& U, FluidSF& T);					///  nlin[i] = FT[v.grad(v)][i]; T.nlin= FT[v.grad(T)]
	
	static void Compute_nlin_scalar(FluidVF& U, FluidSF& T);				// scalar convection
	static void Compute_nlin_RBC(FluidVF& U, FluidSF& T);					// RB convection
	
	static void Compute_nlin(FluidVF& U, FluidSF& T1, FluidSF& T2);

	static void Compute_nlin(FluidVF& U, FluidVF& W);
  
	static void Compute_nlin(FluidVF& U, FluidVF& W, FluidSF& T);
	static void Compute_nlin_MHD_SCALAR_INCOMPRESS(FluidVF& U, FluidVF& W, FluidSF& T);
	static void Compute_nlin_MHD_SCALAR_INCOMPRESS_ALL(FluidVF& U, FluidVF& W, FluidSF& T);

	static void Compute_nlin(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	static void Compute_nlin_MHD_ASTRO_INCOMPRESS(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	static void Compute_nlin_MHD_ASTRO_INCOMPRESS_ALL(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);

	
    
    //GP
    static void Compute_nlin(FluidSF& T);
    static void Compute_nlin_Tsqr(FluidSF& T);
	
		// FOR Energy Transfers
	
	static void Compute_nlin_diag(FluidVF& U, PlainFluidVF& W);
	static void Compute_nlin_offdiag(FluidVF& U, PlainFluidVF& W);
	static void Compute_nlin_first_component_product(FluidVF& U, PlainFluidVF& W);
	
	static void Compute_nlin(FluidVF& U, PlainFluidVF& W);
	static void Compute_nlin_first_component(FluidVF& U, PlainFluidVF& W);
    static void Compute_nlin_helical(FluidVF& U, PlainFluidVF& W);
};

#endif

//************************************  End of IncVF.h  ***************************************


