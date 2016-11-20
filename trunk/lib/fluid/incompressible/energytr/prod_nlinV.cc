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
 * Tarang-2 is distriEnergyTr::buted in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-2; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */
/*! \file  prod_nlinV.cc
 * 
 * @brief  Returns \f$\sum \Re(N_i^*(k) V_i(k)) \f$  or \f$ \sum \Re(N_i^*(k) W_i(k)) \f$
 *			for the specified region  (inside sphere or outside sphere).  
 *			\f$ \vec{N} \f$ is the nlin term. 
 * 
 * We also compute ans_real = \f$ \sum \Re(N_i(k)) * \Re(V_i(k)) \f$ and 
 *				   ans_imag = \f$ \sun \Im(N_i(k)) * \Im(V_i(k)) \f$ 
 *					(Could be \fR W_i \f$ as well).
 *
 * Sphere_index is passed to the function as a parameter.  The function
 *		picks the radius and computes real( A[k]*conj(B[k]) ) using shell-mult functions. <BR>
 *		For inside sphere, the surface modes are included.
 *		Fou outside sphere, the surface modes are excluded.
 *
 * @sa void universal->Shell_mult_single(...)
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "EnergyTr.h"

//********************************************************************************************* 

/** @brief For fluid: compute \f$\sum \Re(N_i^*(k) V_i(k)) \f$ outside sphere.  
 *	
 *	@param	sphere_index 1:no_sphere, the last sphere contains all the modes.
 *
 *	@return  \f$\sum \Re(U.N_i^*(k) W.cvf.V_i(k)) \f$ outside sphere.
 */
Real EnergyTr::Prod_out_sphere_nlinV(int sphere_index, FluidVF& U, FluidVF& W)
{
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = global.myconstant.INF_RADIUS;
	
	return  ( universal->Shell_mult_single(U.nlin1, W.cvf.V1, inner_radius, outer_radius)
			+ universal->Shell_mult_single(U.nlin2, W.cvf.V2, inner_radius, outer_radius)
			+ universal->Shell_mult_single(U.nlin3, W.cvf.V3, inner_radius, outer_radius) ); 
}

// FOR PASSIVE SCALAR --- OUTSIDE SPHERE

Real EnergyTr::Prod_out_sphere_nlinT(int sphere_index, FluidVF& U, FluidSF& T)
{
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = global.myconstant.INF_RADIUS;			
	
	return ( universal->Shell_mult_single(U.nlin1, T.csf.F, inner_radius, outer_radius) );	
}

// Vpll to Vpll ET
Real EnergyTr::Prod_out_sphere_nlin_Vpll(int sphere_index, FluidVF& U, FluidVF& W)
{
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = global.myconstant.INF_RADIUS;
	
    if (global.field.anisotropy_dirn == 1)
        return ( universal->Shell_mult_single(U.nlin1, W.cvf.V1, inner_radius, outer_radius) );
    
    else if (global.field.anisotropy_dirn == 2)
        return ( universal->Shell_mult_single(U.nlin1, W.cvf.V2, inner_radius, outer_radius) );
    
    else if (global.field.anisotropy_dirn == 3)
        return ( universal->Shell_mult_single(U.nlin1, W.cvf.V3, inner_radius, outer_radius) );

    return 0;
}

//*********************************************************************************************

/** @brief For fluid: compute \f$\sum \Re(N_i^*(k) V_i(k)) \f$ inside sphere.  
 *	
 *	@param	sphere_index 1:no_sphere, the last sphere contains all the modes.
 *
 *	@return  \f$\sum \Re(N_i^*(k) V_i(k)) \f$ inside sphere.
 */
Real EnergyTr::Prod_in_sphere_nlinV(int sphere_index, FluidVF& U, FluidVF& W)
{
	Real inner_radius = 0;
	Real outer_radius = global.energy_transfer.flux.radii(sphere_index);

	return  ( universal->Shell_mult_single(U.nlin1, W.cvf.V1, inner_radius, outer_radius)
			 + universal->Shell_mult_single(U.nlin2, W.cvf.V2, inner_radius, outer_radius)
			 + universal->Shell_mult_single(U.nlin3, W.cvf.V3, inner_radius, outer_radius) ); 
}

// FOR PASSIVE SCALAR --- INSIDE SPHERE

Real EnergyTr::Prod_in_sphere_nlinT(int sphere_index, FluidVF& U, FluidSF& T)
{
	Real inner_radius = 0;
	Real outer_radius = global.energy_transfer.flux.radii(sphere_index);
	
	return ( universal->Shell_mult_single(U.nlin1, T.csf.F, inner_radius, outer_radius) );		
}

Real EnergyTr::Prod_in_sphere_nlin_Vpll(int sphere_index, FluidVF& U, FluidVF& W)
{
	Real inner_radius = 0;
	Real outer_radius = global.energy_transfer.flux.radii(sphere_index);
	
    if (global.field.anisotropy_dirn == 1)
        return ( universal->Shell_mult_single(U.nlin1, W.cvf.V1, inner_radius, outer_radius) );
    
    else if (global.field.anisotropy_dirn == 2)
        return ( universal->Shell_mult_single(U.nlin1, W.cvf.V2, inner_radius, outer_radius) );
    
    else if (global.field.anisotropy_dirn == 3)
        return ( universal->Shell_mult_single(U.nlin1, W.cvf.V3, inner_radius, outer_radius) );

    return 0;
}

//*********************************************************************************************
// Helicity flux
Real EnergyTr::Prod_out_sphere_nlin_vorticity(int sphere_index, FluidVF& U, FluidVF& W)
{
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = global.myconstant.INF_RADIUS;	
	
	return universal->Shell_mult_vorticity(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, inner_radius, outer_radius);
	
}


// HM Helicity flux
Real EnergyTr::Prod_out_sphere_nlin_vector_potential(int sphere_index, FluidVF& U, FluidVF& W)
{
	Real inner_radius = global.energy_transfer.flux.radii(sphere_index);
	Real outer_radius = global.myconstant.INF_RADIUS;	
	
	return universal->Shell_mult_vector_potential(U.nlin1, U.nlin2, U.nlin3, W.cvf.V1, W.cvf.V2, W.cvf.V3, inner_radius, outer_radius);	
}

//******************************  End of Prod_nlinV.cc ****************************************






