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

/*! \file  Cvf.h
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008 
 * @bugs  No known bugs
 */

#ifndef _H_PLAIN_FluidVF
#define _H_PLAIN_FluidVF

#include "fields.h" 


//*********************************************************************************************

//! Complex vector field
/*!
 * 3 dimensional complex vector field has 3 components with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2,N_3/2+1] \f$.   <BR>
 * These arrays can also store real values, and their dimension are \f$[N_1, N_2, N_3] \f$. 
 *
 * Inverse transform: Complex -> Real.  <BR>
 * Forward transform: Real -> Complex. <BR>
 * The implemented basis functions: <BR>
 *  FOUR: Fourier along all directions. <BR>
 *  SCFT: Sin/Cos along x, and Fourier along all perpendicular directions. <BR>
 * 
 *  The modal energy \f$ C(\vec{k}) \f$ is defined in each basis function.  
 *  The dissipation rate for the mode is \f$ 2 K^2 C(\vec{k}) \f$.
 *  
 *  Helicity1 = \f$ H1(K) = \vec{K} . [\vec{Vr} \times \vec{Vi}] \f$. <BR>
 *  Helicity2 = \f$ H2(K) = \vec{K} . [(\vec{Vr} \times \vec{Vi})] /K^2 \f$. <BR>
 * 
 *  Entropy = \f$ \sum p(k) log(1/p(k)) \f$ where probability  
 *				\f$ p(k) =  E(k)/E_{total) \f$ with \f$ E(k) = |A(k)|^2/2 \f$.
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */

//*********************************************************************************************

class PlainFluidVF
{
	
public:
	CVF cvf;
	RVF rvf;
	PlainFluidVF(string field_name="");
	
	void Inverse_transform();
	void Forward_transform();
	void Inverse_transform_first_component();
	void Forward_transform_first_component();
};
#endif

//******************************** End of CVF.h ***********************************************	






