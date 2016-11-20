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

/*! \file  IncSF.h
 * 
 * @brief  Class declaration of IncSF, Incompressible scalar Field 
 *		example: passive scalar, temperature in RB convection
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug  No known bugs
 */

//*********************************************************************************************

#ifndef _Pressure_H
#define _Pressure_H

#include "fluid_base.h"

//! @brief Incompressible scalar field IncSF 
/*!
 *  Inherits CSF that contains the complex scalar field. <BR>
 *  RSF that contains real vector field, typically Inverse tranform of CSF. <BR>
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

class Pressure: public CSF 
{ 
	
 public:
	
	void Compute_pressure(FluidVF &U);
	
};

#endif

//************************** End of IncSF.h   *************************************************	



