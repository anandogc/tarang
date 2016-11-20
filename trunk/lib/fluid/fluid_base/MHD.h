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
 
#ifndef _H_MHD
#define _H_MHD

#include "FluidVF.h"
#include "FluidSF.h"

//*********************************************************************************************


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


class MHD 
{
public:
	static void UB_to_Elsasser_field(FluidVF& U, FluidVF& W);
	static void Elsasser_to_UB_field(FluidVF& U, FluidVF& W);
	
	static void UB_to_Elsasser_real_field(FluidVF& U, FluidVF& W);
	static void Elsasser_to_UB_real_field(FluidVF& U, FluidVF& W);
	
	static void UB_to_Elsasser_force(FluidVF& U, FluidVF& W);
	static void Elsasser_to_UB_force(FluidVF& U, FluidVF& W);
	
	static void UB_to_Elsasser_nlin(FluidVF& U, FluidVF& W);
	static void Elsasser_to_UB_nlin(FluidVF& U, FluidVF& W);
};


#endif

//************************************  End of IncVF.h  ***************************************


