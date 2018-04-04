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

#ifndef _H_Correlation
#define _H_Correlation


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
 
#include "FluidVF.h"
#include "FluidSF.h"

//*********************************************************************************************	

class Correlation
{ 
public:
	
	//!  Energy of Vx in shell k
	static  Array<Real,1>		shell_ek1;
	
	//!  Energy of Vy in shell k
	static Array<Real,1>		shell_ek2;
	
	//!  Energy of Vz in shell k
	static Array<Real,1>		shell_ek3;
	//shubhadeep
	static  Array<Real,1>		shell_ek1_force;
	static Array<Real,1>		shell_ek2_force;
	static Array<Real,1>		shell_ek3_force;
	//shubhadeep
	
	//!  Energy Dissipation rate in shell k (without \f$ \nu \f$).
	static Array<Real,1>		shell_dissk1;
	static Array<Real,1>		shell_dissk2;
	static Array<Real,1>		shell_dissk3;
	
	//! Sum \f$ \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K})) \f$ in shell k
	//! components along 1,2,3 directions.
	
	//!  Energy in ring(m,n) along e1 (toroidal direction).
	static Array<Real,2>		ring_ek1;
	
	//!  Energy in ring(m,n) along e2 (poloidal direction).
	static Array<Real,2>		ring_ek2;
	
	static Array<Real,2>		ring_ek3;

	static Array<Real,2>		ring_hk1;
	static Array<Real,2>		ring_hk2;
	static Array<Real,2>		ring_hk3;
	
	//!  Energy Dissipation rate in ring(m,n) (without \f$ \nu \f$).
	static Array<Real,2>		ring_dissk1, ring_dissk2, ring_dissk3;
	
	
	
	//!  Energy spectrum along the anisotropy direction.
	static Array<Real,2>		cylindrical_ring_ek1;
	
	//!  Energy spectrum perpendicular to the anisotropy direction.
	static Array<Real,2>		cylindrical_ring_ek2;
	
	//!  Energy Dissipation rate in ring(m,n) (without \f$ \nu \f$).
	static Array<Real,2>		cylindrical_ring_dissk1, cylindrical_ring_dissk2;
	
	//! \f$ \sum \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K})) \f$
	/// in ring(m,n)
	
	//!  For scalars
	static Array<Real,1>		shell_ek;
	static Array<Real,1>		shell_dissk;
	
	static Array<Real,2>		ring_ek;
	static Array<Real,2>		ring_dissk;
	
	static Array<Real,2>		cylindrical_ring_ek;
	static Array<Real,2>		cylindrical_ring_dissk;
	
	
	//! for forcing and init-cond
	static Array<Real,1>		shell_ek_temp1;
	static Array<Real,1>		shell_ek_temp2;

	
	static Array<Real,2> ring_spectrum;
	static Array<Real,2> cylindrical_ring_spectrum;

	static void Initialize();
	static Real Get_Nusselt_no(FluidVF& U, FluidSF& T);
	static Real Get_cross_helicity(FluidVF &U, FluidVF& W);
	
	static void Compute_shell_spectrum(FluidVF& U);
	static void Compute_shell_spectrum(FluidVF& U, FluidVF& W);
	static void Compute_shell_spectrum_dissipation(FluidVF& U, FluidVF& W); //shubhadeep
	static void Compute_shell_spectrum(FluidVF& U, FluidSF& T);
	static void Compute_force_shell_spectrum_helical(FluidVF& U,FluidVF& helicalU);//shubhadeep
	
	static void Compute_ring_spectrum(FluidVF& U);
	static void Compute_ring_spectrum(FluidVF& U, FluidVF& W);
	static void Compute_ring_spectrum(FluidVF&U, FluidSF& T);
	static void Compute_helical_ring_spectrum(FluidVF& U, FluidVF& helicalU);
	
	static void Compute_cylindrical_ring_spectrum(FluidVF& U);
	static void Compute_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W);
	static void Compute_cylindrical_ring_spectrum(FluidVF& U, FluidSF& T);
	
	
	static void Compute_shell_spectrum_helicity(FluidVF& U);
    static void Compute_shell_spectrum_helicity2(FluidVF& U);
	static void Compute_ring_spectrum_helicity(FluidVF& U);
	static void Compute_cylindrical_ring_spectrum_helicity(FluidVF& U);
	
	static void Compute_shell_spectrum(FluidSF& T);
	static void Compute_shell_spectrum(FluidSF& T1, FluidSF& T2);
	
	static void Compute_ring_spectrum(FluidSF& T);
	static void Compute_ring_spectrum(FluidSF& T1, FluidSF& T2);
	
	static void Compute_cylindrical_ring_spectrum(FluidSF& T);
	static void Compute_cylindrical_ring_spectrum(FluidSF& T1, FluidSF& T2);

	static void Compute_force_shell_spectrum(FluidVF& U);
	static void Compute_force_ring_spectrum(FluidVF& U);
	static void Compute_force_cylindrical_ring_spectrum(FluidVF& U);
	
	static void Compute_force_shell_spectrum(FluidSF& T);
	static void Compute_force_ring_spectrum(FluidSF& T);
	static void Compute_force_cylindrical_ring_spectrum(FluidSF& T);
	
	
	static void Compute_Tk_shell_spectrum(FluidVF& U);
	static void Compute_Tk_ring_spectrum(FluidVF& U);
	static void Compute_Tk_cylindrical_ring_spectrum(FluidVF& U);
	
	static void Compute_Tk_shell_spectrum(FluidSF& T);
	static void Compute_Tk_ring_spectrum(FluidSF& T);
	static void Compute_Tk_cylindrical_ring_spectrum(FluidSF& T);
};
	

#endif

//************************************  End of IncVF.h  ***************************************


